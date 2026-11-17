# shellcheck shell=bash
# =============================================================================
# haps-guard.sh -- sourced by the 3b arms (RFMix and FLARE). One check, used
# on every reference .haps/.sample pair extract_haps caches and on the
# simulated cohort's pair before it is read.
#
# WHY. In August, chr13's NAT_PEL_chr13.sample in the RFMix work folder was 62
# bytes -- the two header lines and no individual -- beside an intact .haps.
# The cache test was "the .haps exists", so every later task reused the pair,
# and the NAT_PEL panel ran on chr13 with a reference sample list that did not
# match its haplotypes. Nothing failed.
#
# The likely cause is in the callers: `[[ -s list ]] && extract_haps ... || true`
# turns off `set -e` inside extract_haps, so a plink2 that died part way still
# had its output stamped as cached and its .sample rewritten into a header.
# The callers now use `if`, plink2's exit status is checked, and this check
# runs on every cached pair and again after every extraction.
#
# A SHAPEIT .haps has 5 leading columns and 2 per individual, so a pair is
# consistent only if every line has 5 + 2 x (individuals in the .sample)
# columns. The first and the LAST line are checked: a write cut off part way
# leaves a short last line, and reading the last line costs nothing.
#
#   haps_pair_check <file.haps> <file.sample>   -> 0, or prints why and returns 1
# =============================================================================
haps_pair_check () {
    local haps="$1" samp="$2" n c1 cn
    if [[ ! -s "$haps" ]]; then echo "${haps} is missing or empty"; return 1; fi
    if [[ ! -s "$samp" ]]; then echo "${samp} is missing or empty"; return 1; fi
    n=$(awk 'NR > 2 && NF' "$samp" | wc -l)
    if [[ "$n" -lt 1 ]]; then
        echo "${samp} lists no individual ($(wc -c < "$samp") bytes: the header only)"
        return 1
    fi
    c1=$(awk 'NR == 1 {print NF; exit}' "$haps")
    cn=$(tail -n 1 "$haps" | awk '{print NF}')
    if [[ "${c1:-0}" -ne $((5 + 2 * n)) || "${cn:-0}" -ne $((5 + 2 * n)) ]]; then
        echo "${haps} has ${c1:-0} columns on its first line and ${cn:-0} on its last, but ${samp}" \
             "lists ${n} individuals, so both should be $((5 + 2 * n))"
        return 1
    fi
    return 0
}
