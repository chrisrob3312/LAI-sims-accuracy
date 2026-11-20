# Christina's preferences for this repo

## Copy-paste block hygiene (IMPORTANT)

When giving shell commands where a **later step depends on verifying an
earlier step's output**, put them in **separate copy-paste blocks**. Never
combine "diagnose" and "act" (rm / sbatch / git push / anything that
mutates state) into the same block she might paste all at once.

Concretely:

- **Block A** = read-only diagnostics (`ls`, `cat`, `wc`, `sacct`,
  `squeue`, `find`, `git log`, `awk` inspections, etc.).
- **Explain what to look for** in Block A's output.
- **Block B** = the destructive / submission commands (`rm`, `sbatch`,
  `scancel`, `git commit`, `git push`, `scontrol update`).
- Preface Block B with an explicit "**only run if <verification>**" line.

Rule of thumb: if a wrong answer to the diagnostic would make the second
block harmful (delete good data, submit a bad job, push wrong commit),
they MUST be in separate blocks.

## Diagnose from the logs first (IMPORTANT)

When a job fails or a result looks wrong, **read the job's own `.out` and
`.err` in full (`cat`, not a filtered `grep`) before forming a theory**.
A diagnostic written around an assumption can confirm the wrong
assumption -- the first lines of a 4a `.err` held the real bug
(`line 264: sp: unbound variable`) while a targeted grep looked elsewhere.

- Never suppress stderr in a diagnostic (`2>/dev/null`, `|| true` on the
  thing being checked). A check that fails silently looks like a check
  that passed.
- **Look paths up from the scripts** (`grep -n` the variable that writes
  the file) instead of guessing them. Say where each path came from.
- Run cluster tools inside the env the pipeline uses:
  `module load anaconda3/2024.06`, `source "$(conda info --base)/etc/profile.d/conda.sh"`,
  `conda activate shapeit5` -- plink2 and bcftools are not on `(base)`.
- Show raw output (headers, first rows, counts) alongside any summary, so
  a parsing error in the check itself is visible.

## Other running preferences

- Long RFMix / accuracy pilots run on the mhgcp cluster; slides + tables
  live in this repo and get regen'd on the main-loop side, not the cluster
  (shapeit5 env lacks matplotlib).
- Speaker notes in the committee pptx carry the depth — slides stay at
  headline + one dominant figure + tiny callouts (Spring TAC style).
- Deck palette = Spring TAC theme3 colors (navy #0E2841, teal #156082 for
  Mexican-like cohort, orange #E97132 for Brazilian-like cohort, forest
  #196B24 EUR, cyan #0F9ED5 SAS, plum #A02B93 AMR).
- HOMOG panel = HGDP + 1KG + MXB joint-called and Q≥0.95 filtered by
  supervised ADMIXTURE. It **includes MXB samples** — do not describe it
  as MXB-free / donor-safe.
