#!/usr/bin/env python3
"""
failure-report.py -- which pipeline tasks failed, how, and their own logs.
READ-ONLY on the pipeline. Run by scripts/failure-report.sh, which supplies:
  --sacct   `sacct -n -P -o JobID,JobIDRaw,JobName,State,ExitCode,Elapsed,
             Timelimit,MaxRSS,ReqMem,NodeList,Submit` text
  --squeue  `squeue -h -o '%i|%j|%T|%r'` text
  --logdir  where every stage writes <prefix>_<jobid>.out/.err
Writes the whole report (no log cut) to --report and the same to stdout with
logs longer than --cap lines cut to their first and last lines (--cap 0 = never).
"""
import argparse
import glob
import os
import re
import sys
from collections import Counter, OrderedDict, defaultdict
from datetime import datetime, timedelta

# Every job name the repo's stage scripts set (#SBATCH --job-name), in run
# order. Any other name in sacct (an interactive `srun --pty bash` is "bash")
# is reported as "not a pipeline job".
ORDER = ["sgdp_fetch", "1b_sgdp_raw", "1b_sgdp_liftover", "strip_almarri", "1a_prep_mxb_liftover",
         "1d_almarri", "1c_phase", "1b_prep", "1b_chunk", "1b_ligate", "1e_king_dedup", "4a_homog",
         "4a_homog_kN", "2b_simu", "3b_rfmix", "3b_homog", "3b_chip", "3b_chip_homog", "3b_flare",
         "finalize_pilot", "5a_lift19", "6b_simset", "6c_real", "6d_cells", "6e_kinship", "site_report"]
NOT_PIPELINE = "(not a pipeline job)"
OK = {"COMPLETED", "RUNNING", "PENDING", "REQUEUED", "RESIZING", "SUSPENDED"}
# Submissions of one name made within this window of its latest one are the
# same run (the chunk stage is two arrays, main + heavy, submitted together).
SAME_RUN = timedelta(minutes=10)


def parse_time(s):
    try:
        return datetime.strptime(s, "%Y-%m-%dT%H:%M:%S")
    except (ValueError, TypeError):
        return None


def task_key(jobid):
    """'3754713_5' -> (3754713, '5'); '3754282' -> (3754282, ''); pending '3754713_[16-22%15]' too."""
    base, _, idx = jobid.partition("_")
    return int(base), idx


def idx_sort(i):
    m = re.match(r"\d+", (i or "").lstrip("["))
    return (int(m.group(0)) if m else -1, i)


def n_tasks(idx):
    """Tasks one sacct row stands for. A pending array range is ONE row,
    '3754839_[212-594%40]' -> 383; a single task or a plain job -> 1."""
    m = re.match(r"^\[(.*)\]$", idx or "")
    if not m:
        return 1
    n = 0
    for part in m.group(1).split("%")[0].split(","):
        rng, _, step = part.partition(":")
        lo, _, hi = rng.partition("-")
        try:
            lo, hi, step = int(lo), int(hi or lo), int(step or 1)
        except ValueError:
            return 1
        if hi >= lo and step > 0:
            n += (hi - lo) // step + 1
    return max(n, 1)


FLAG = re.compile(r"warn|error|fail|fatal|abort|killed|segfault|core dump|exception|traceback|not found|no such"
                  r"|cannot|can't|unable|invalid|unbound|stale|mismatch|refus|denied|out of memory|oom|truncat"
                  r"|corrupt|skipp|\[[WE]::", re.I)   # [W::...] / [E::...]: htslib's own warnings and errors


def norm(line):
    """A log line with its numbers made generic, so 22 chromosomes' copies group."""
    s = re.sub(r"\x1b\[[0-9;]*m", "", line).strip()
    s = re.sub(r"0x[0-9a-fA-F]+", "0x#", s)
    s = re.sub(r"\d+(?:\.\d+)?(?:[eE][-+]?\d+)?", "#", s)
    return re.sub(r"\s+", " ", s)[:240]


def sweep(a, names, by_name, current_by_name, emit, out_full, out_screen):
    """Section 6: every pipeline task of the current run that has started, its .err
    and .out read in full, lines grouped across tasks by their generic form."""
    emit("\n== 6. every log of the current run: lines with a warning/error word, then lines only a"
         " few tasks print ==")
    emit("  (numbers are shown as # to group the chromosomes' copies; each group gives how many tasks"
         " print it, how often, and one raw example)")

    def show(lines):
        for i, s in enumerate(lines):
            if i < a.sweep_cap:
                emit(s)
            else:
                out_full.append(s)
        if len(lines) > a.sweep_cap:
            out_screen.append(f"      [... {len(lines) - a.sweep_cap} more groups in the report file ...]")

    for n in names:
        if n not in ORDER:
            continue
        tasks = [r for b in current_by_name.get(n, []) for r in by_name[n][b] if r["state"] != "PENDING"]
        if not tasks:
            continue
        groups, n_logs, n_err, n_err_lines, missing = {}, 0, 0, 0, 0
        for r in tasks:
            for ext in ("err", "out"):
                hits = sorted(glob.glob(os.path.join(a.logdir, f"*_{r['raw']}.{ext}")))
                if not hits:
                    missing += 1
                    continue
                p = hits[0]
                n_logs += 1
                try:
                    with open(p, errors="replace") as fh:
                        lines = fh.read().splitlines()
                except OSError as e:
                    emit(f"  {n}: {p}: {e}")
                    continue
                if ext == "err" and any(l.strip() for l in lines):
                    n_err += 1
                    n_err_lines += sum(1 for l in lines if l.strip())
                for l in lines:
                    if not l.strip():
                        continue
                    k = (ext, norm(l))
                    g = groups.get(k)
                    if g is None:
                        g = groups[k] = [set(), 0, f"{os.path.basename(p)}: {l.strip()[:200]}"]
                    g[0].add(r["jobid"])
                    g[1] += 1
        nt = len(tasks)
        emit(f"\n  -- {n}: {nt} task(s) started; {n_logs} log files read ({missing} not found);"
             f" .err has text in {n_err} ({n_err_lines} lines)")
        flagged = sorted(((k, g) for k, g in groups.items() if FLAG.search(k[1])),
                         key=lambda kg: (-len(kg[1][0]), -kg[1][1], kg[0]))
        if flagged:
            emit(f"     lines with a warning/error word: {len(flagged)} group(s)")
            show([f"      [{k[0]}] {len(g[0])}/{nt} tasks x{g[1]}  {k[1]}\n"
                  f"             e.g. {g[2]}" for k, g in flagged])
        else:
            emit("     lines with a warning/error word: none")
        # Tools write their warnings to stderr, worded however they like: every
        # other .err line too, however many tasks print it.
        other_err = sorted(((k, g) for k, g in groups.items() if k[0] == "err" and not FLAG.search(k[1])),
                           key=lambda kg: (-len(kg[1][0]), -kg[1][1], kg[0]))
        if other_err:
            emit(f"     other .err lines: {len(other_err)} group(s)")
            show([f"      [err] {len(g[0])}/{nt} tasks x{g[1]}  {k[1]}\n"
                  f"             e.g. {g[2]}" for k, g in other_err])
        if nt >= 5:
            few = max(1, nt // 10)
            rare = sorted(((k, g) for k, g in groups.items()
                           if k[0] == "out" and not FLAG.search(k[1]) and len(g[0]) <= few),
                          key=lambda kg: (len(kg[1][0]), kg[0]))
            if rare:
                emit(f"     lines only {few} or fewer of the {nt} tasks print (no warning word): {len(rare)} group(s)")
                show([f"      [{k[0]}] {len(g[0])}/{nt} tasks x{g[1]}  {k[1]}\n"
                      f"             e.g. {g[2]}" for k, g in rare])


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--sacct", required=True)
    ap.add_argument("--squeue", required=True)
    ap.add_argument("--logdir", required=True)
    ap.add_argument("--since", default="")
    ap.add_argument("--report", required=True)
    ap.add_argument("--cap", type=int, default=150)
    ap.add_argument("--sweep", action="store_true",
                    help="section 6: every log of the current run, completed tasks too, grouped")
    ap.add_argument("--sweep-cap", type=int, default=15, help="groups per stage and kind shown on screen")
    a = ap.parse_args()

    alloc, rss = {}, {}
    with open(a.sacct) as fh:
        for line in fh:
            f = line.rstrip("\n").split("|")
            if len(f) < 11 or not f[0]:
                continue
            jobid, raw, name, state, ec, el, tl, maxrss, reqmem, node, submit = f[:11]
            if "." in jobid:                       # a step: keep the batch step's MaxRSS
                if jobid.endswith(".batch") and maxrss:
                    rss[jobid.rsplit(".", 1)[0]] = maxrss
                continue
            if not re.match(r"^\d+", jobid):
                continue
            alloc[jobid] = dict(jobid=jobid, raw=raw, name=name, state=state.split()[0] if state else "",
                                state_full=state, ec=ec, el=el, tl=tl, reqmem=reqmem, node=node,
                                submit=parse_time(submit))
    for j, r in alloc.items():
        r["maxrss"] = rss.get(j, "")

    by_name = defaultdict(lambda: defaultdict(list))   # name -> base id -> tasks
    for r in alloc.values():
        base, idx = task_key(r["jobid"])
        r["base"], r["idx"] = base, idx
        by_name[r["name"]][base].append(r)

    out_full, out_screen = [], []

    def emit(s="", screen=None):
        out_full.append(s)
        out_screen.append(s if screen is None else screen)

    emit(f"failure-report  {datetime.now():%Y-%m-%d %H:%M}  sacct since {a.since or '?'}  logs: {a.logdir}")
    names = [n for n in ORDER if n in by_name] + sorted(n for n in by_name if n not in ORDER)
    failures = []
    current_by_name = {}
    emit("\n== 1. task states per stage (current run = latest submission of each name; earlier ones listed as superseded) ==")
    emit("  (a pending array range is one sacct row; it is counted here as the tasks it holds)")
    for n in names:
        subs = by_name[n]
        latest = max(subs)
        t_latest = min((r["submit"] for r in subs[latest] if r["submit"]), default=None)
        current = [b for b in subs if b == latest or (t_latest and any(
            r["submit"] and abs(r["submit"] - t_latest) <= SAME_RUN for r in subs[b]))]
        old = sorted(b for b in subs if b not in current)
        current_by_name[n] = current
        cnt = Counter()
        for b in current:
            for r in subs[b]:
                cnt[r["state"]] += n_tasks(r["idx"])
        states = "  ".join(f"{s} {c}" for s, c in sorted(cnt.items()))
        emit(f"  {n:16s} {'+'.join(str(b) for b in sorted(current)):24s} {states}"
             + (f"   (superseded: {' '.join(str(b) for b in old)})" if old else "")
             + ("" if n in ORDER else f"   {NOT_PIPELINE}"))
        for b in sorted(current):
            for r in sorted(subs[b], key=lambda r: idx_sort(r["idx"])):
                if r["state"] not in OK:
                    failures.append(r)

    n_other = sum(1 for r in failures if r["name"] not in ORDER)
    emit(f"\n== 2. tasks that did not complete in the current run: {len(failures)}"
         f" ({len(failures) - n_other} pipeline, {n_other} {NOT_PIPELINE.strip('()')}) ==")
    if not failures:
        emit("  none")
    logs = {}
    for r in failures:
        found = {}
        for ext in ("err", "out"):
            hits = sorted(glob.glob(os.path.join(a.logdir, f"*_{r['raw']}.{ext}")))
            found[ext] = hits[0] if hits else ""
        logs[r["jobid"]] = found
        emit(f"  {r['jobid']:18s} {r['name']:14s} {r['state_full']:20s} exit {r['ec']:6s} "
             f"ran {r['el']} of {r['tl']}  RSS {r['maxrss'] or '-'} of {r['reqmem']}  node {r['node']}"
             + ("" if r["name"] in ORDER else f"  {NOT_PIPELINE}"))
        emit(f"  {'':18s} .err {found['err'] or '(no *_' + r['raw'] + '.err in ' + a.logdir + ')'}")
        emit(f"  {'':18s} .out {found['out'] or '(no *_' + r['raw'] + '.out)'}")

    groups = OrderedDict()
    for r in failures:
        groups.setdefault((r["name"], r["state"], r["ec"]), []).append(r)
    if failures:
        emit("\n== 3. failures grouped: stage, state, exit code -> how many, on which nodes ==")
        for (n, s, ec), rs in groups.items():
            nodes = Counter(r["node"] for r in rs)
            idxs = ",".join(r["idx"] or r["jobid"] for r in rs)
            emit(f"  {n:16s} {s:14s} exit {ec:6s} x{len(rs):<4d} nodes: "
                 + ", ".join(f"{k} {v}" for k, v in nodes.most_common()) + f"   tasks: {idxs}")

        emit("\n== 4. the logs of the FIRST failure of each kind, in full ==")
        for (n, s, ec), rs in groups.items():
            r = rs[0]
            emit(f"\n---------- {n} task {r['jobid']} ({s}, exit {ec}, node {r['node']}) ----------")
            for ext in ("err", "out"):
                p = logs[r["jobid"]][ext]
                if not p:
                    emit(f"[.{ext}: not found]")
                    continue
                try:
                    with open(p, errors="replace") as fh:
                        lines = fh.read().splitlines()
                except OSError as e:
                    emit(f"[.{ext}: {p}: {e}]")
                    continue
                emit(f"[.{ext}: {p} -- {len(lines)} lines]")
                if a.cap and len(lines) > a.cap:
                    head, tail = a.cap // 5, a.cap - a.cap // 5
                    for i, l in enumerate(lines):
                        if i < head or i >= len(lines) - tail:
                            emit(l)
                        else:
                            out_full.append(l)
                    out_screen.insert(len(out_screen) - tail,
                                      f"[... {len(lines) - head - tail} lines not shown on screen; all of them are in {a.report} ...]")
                else:
                    for l in lines:
                        emit(l)

    emit("\n== 5. still queued (squeue): stage, state, reason -> how many ==")
    q = Counter()
    with open(a.squeue) as fh:
        for line in fh:
            f = line.rstrip("\n").split("|")
            if len(f) >= 4:
                q[(f[1], f[2], f[3])] += 1
    if not q:
        emit("  nothing queued or running")
    for (n, s, rsn), c in sorted(q.items(), key=lambda kv: (ORDER.index(kv[0][0]) if kv[0][0] in ORDER else 99, kv[0])):
        flag = "   <-- will never start: its dependency failed" if "DependencyNeverSatisfied" in rsn else ""
        emit(f"  {n:16s} {s:10s} {rsn:28s} x{c}{flag}")

    if a.sweep:
        sweep(a, names, by_name, current_by_name, emit, out_full, out_screen)

    with open(a.report, "w") as g:
        g.write("\n".join(out_full) + "\n")
    print("\n".join(out_screen))
    print(f"\nfull report (no log cut): {a.report}")


if __name__ == "__main__":
    main()
