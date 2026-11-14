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

ORDER = ["sgdp_fetch", "1b_sgdp_raw", "1a_prep_mxb_liftover", "1d_almarri", "1c_phase", "1b_chunk",
         "1b_ligate", "1e_king_dedup", "4a_homog_kN", "2b_simu", "3b_rfmix", "3b_homog", "3b_chip",
         "3b_chip_homog", "3b_flare", "finalize_pilot", "site_report"]
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
    m = re.match(r"\d+", i or "")
    return (int(m.group(0)) if m else -1, i)


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--sacct", required=True)
    ap.add_argument("--squeue", required=True)
    ap.add_argument("--logdir", required=True)
    ap.add_argument("--since", default="")
    ap.add_argument("--report", required=True)
    ap.add_argument("--cap", type=int, default=150)
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
    emit("\n== 1. task states per stage (current run = latest submission of each name; earlier ones listed as superseded) ==")
    for n in names:
        subs = by_name[n]
        latest = max(subs)
        t_latest = min((r["submit"] for r in subs[latest] if r["submit"]), default=None)
        current = [b for b in subs if b == latest or (t_latest and any(
            r["submit"] and abs(r["submit"] - t_latest) <= SAME_RUN for r in subs[b]))]
        old = sorted(b for b in subs if b not in current)
        cnt = Counter(r["state"] for b in current for r in subs[b])
        states = "  ".join(f"{s} {c}" for s, c in sorted(cnt.items()))
        emit(f"  {n:16s} {'+'.join(str(b) for b in sorted(current)):24s} {states}"
             + (f"   (superseded: {' '.join(str(b) for b in old)})" if old else ""))
        for b in sorted(current):
            for r in sorted(subs[b], key=lambda r: idx_sort(r["idx"])):
                if r["state"] not in OK:
                    failures.append(r)

    emit(f"\n== 2. tasks that did not complete in the current run: {len(failures)} ==")
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
             f"ran {r['el']} of {r['tl']}  RSS {r['maxrss'] or '-'} of {r['reqmem']}  node {r['node']}")
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

    with open(a.report, "w") as g:
        g.write("\n".join(out_full) + "\n")
    print("\n".join(out_screen))
    print(f"\nfull report (no log cut): {a.report}")


if __name__ == "__main__":
    main()
