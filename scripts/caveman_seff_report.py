#!/usr/bin/env python3
"""Read Nextflow trace TSVs and collect seff metrics. Python 3 standard library only.

Run in caveman-full-scatter-test:
  python3 caveman_seff_report.py
  python3 caveman_seff_report.py --trace trace-22346168.tsv trace-NEW.tsv

Default: read all trace-*.tsv in --project-dir (default current directory).
Outputs: jobs.tsv, summary.txt, raw_seff/*.txt. Does not modify pipeline outputs.
Each Slurm job is a group of shards, NOT one original splitList region.
"""
import argparse
import csv
import math
import os
from pathlib import Path
import re
import shutil
import statistics
import subprocess
import sys


def seconds(value):
    match = re.fullmatch(r"(?:(\d+)-)?(\d+):(\d+):(\d+(?:\.\d+)?)", value.strip())
    if not match:
        return None
    d, h, m, s = match.groups()
    return int(d or 0)*86400 + int(h)*3600 + int(m)*60 + float(s)


def gib(value):
    match = re.match(r"([\d.]+)\s*([KMGTPE]?)(?:i?B)?\b", value.strip(), re.I)
    if not match:
        return None
    size, unit = match.groups()
    return float(size) * 1024 ** (" KMGTPE".index(unit.upper()) - 3)


def parse_seff(raw):
    fields = {}
    for line in raw.splitlines():
        if ":" in line:
            k, v = line.split(":", 1)
            fields[k.strip()] = v.strip()
    def pct(key):
        match = re.match(r"([\d.]+)%", fields.get(key, ""))
        return float(match[1]) if match else None
    def integer(key):
        match = re.match(r"\d+", fields.get(key, ""))
        return int(match[0]) if match else None
    nodes = integer("Nodes")
    cpus = integer("Cores")
    if cpus is None and nodes is not None and integer("Cores per node") is not None:
        cpus = nodes * integer("Cores per node")
    wall = seconds(fields.get("Job Wall-clock time", ""))
    cpu = seconds(fields.get("CPU Utilized", ""))
    mem = fields.get("Memory Efficiency", "")
    match = re.search(r"% of\s+([\d.]+\s*[KMGTPE]?i?B)\b", mem, re.I)
    result = dict(slurm_state=fields.get("State", "UNKNOWN"), allocated_cpus=cpus,
                  wall_hours=wall/3600 if wall is not None else None,
                  cpu_hours=cpu/3600 if cpu is not None else None,
                  average_used_cores=cpu/wall if cpu is not None and wall else None,
                  cpu_efficiency_pct=pct("CPU Efficiency"),
                  peak_memory_gib=gib(fields.get("Memory Utilized", "")),
                  allocated_memory_gib=gib(match[1]) if match else None,
                  memory_efficiency_pct=pct("Memory Efficiency"))
    issues = []
    if any(result[k] is None for k in ("wall_hours", "cpu_hours", "peak_memory_gib", "allocated_cpus")):
        issues.append("missing_accounting_fields")
    if cpu == 0 or result["peak_memory_gib"] == 0:
        issues.append("zero_accounting_value_verify_before_optimizing")
    if result["cpu_efficiency_pct"] is not None and result["cpu_efficiency_pct"] > 105:
        issues.append("cpu_efficiency_above_105_verify_accounting")
    if result["memory_efficiency_pct"] is not None and result["memory_efficiency_pct"] > 100:
        issues.append("reported_memory_exceeds_allocation_verify")
    result["notes"] = ";".join(issues)
    return result


def load_jobs(paths):
    jobs = {}
    skipped = 0
    for path in paths:
        with path.open(encoding="utf-8-sig", newline="") as handle:
            reader = csv.DictReader(handle, delimiter="\t")
            if not {"native_id", "name", "status"} <= set(reader.fieldnames or []):
                raise ValueError("Trace needs native_id/name/status columns: " + str(path))
            for row in reader:
                name = row.get("name", "")
                # Allow workflow-qualified names such as workflow:Mstep (group_1).
                match = re.search(r"(?:^|:)(?:CaVEMan_)?(Mstep|Estep|ModelMerge|Gather)(?:\s|$)", name)
                if not match:
                    continue
                job = row.get("native_id", "").strip()
                if not re.fullmatch(r"\d+(?:_\d+)?", job):
                    skipped += 1
                    continue
                if job not in jobs:
                    jobs[job] = dict(job_id=job, stage=match[1], task=name,
                                     trace_statuses=set(), trace_files=set(),
                                     nextflow_hash=row.get("hash", ""))
                jobs[job]["trace_statuses"].add(row["status"])
                jobs[job]["trace_files"].add(str(path))
    for row in jobs.values():
        row["trace_statuses"] = ",".join(sorted(row["trace_statuses"]))
        row["trace_files"] = ";".join(sorted(row["trace_files"]))
    return list(jobs.values()), skipped


def p95(values):
    return sorted(values)[max(0, math.ceil(len(values)*0.95)-1)]


def summarize(rows):
    lines = ["CaVEMan Slurm/seff 资源报告", "",
             "统计单位：一个 Slurm 分组任务；不是单个原始 genomic shard。",
             "同一作业号跨 trace 去重；旧作业和失败重跑分别保留。",
             "缓存任务查询原 Slurm 作业的实际资源，不把缓存命中时间当作运行时间。",
             "分组汇总仅使用 Slurm COMPLETED 且记账完整、非零、无异常标记的任务。",
             "seff 内存是 Slurm 采样估计，可能漏掉瞬时峰值。", ""]
    for stage in ("Mstep", "ModelMerge", "Estep", "Gather"):
        all_rows = [r for r in rows if r["stage"] == stage]
        good = [r for r in all_rows if r.get("slurm_state", "").startswith("COMPLETED")
                and not r.get("notes")]
        lines.append("{}: 作业记录 {}，用于汇总 {}".format(stage, len(all_rows), len(good)))
        if not good:
            lines.append("  无可靠数据；查看 jobs.tsv 的 notes 和 raw_seff。\n")
            continue
        for key, label, unit in (("wall_hours", "实际耗时", "h"),
                                 ("average_used_cores", "平均实际使用核心", ""),
                                 ("cpu_efficiency_pct", "CPU效率", "%"),
                                 ("peak_memory_gib", "峰值内存", "GiB"),
                                 ("memory_efficiency_pct", "内存效率", "%")):
            values = [r[key] for r in good if r.get(key) is not None]
            if values:
                lines.append("  {}: 中位数 {:.2f}{}，P95 {:.2f}{}，最大 {:.2f}{}".format(
                    label, statistics.median(values), unit, p95(values), unit, max(values), unit))
        allocated = sorted(set(r["allocated_cpus"] for r in good))
        lines.append("  申请核心数: " + ", ".join(str(x) for x in allocated))
        allocated_hours = sum(r["wall_hours"]*r["allocated_cpus"] for r in good)
        used_hours = sum(r["cpu_hours"] for r in good)
        lines.append("  合计申请核时 {:.2f} h，实际CPU核时 {:.2f} h，加权CPU效率 {:.2f}%".format(
            allocated_hours, used_hours, 100*used_hours/allocated_hours if allocated_hours else 0))
        peak = max(r["peak_memory_gib"] for r in good)
        lines.append("  内存试跑参考下限: ceil(本批最大采样峰值 × 1.3) = {} GiB".format(math.ceil(peak*1.3)))
        if stage in ("Mstep", "Estep"):
            lines.append("  请结合 result/pool.json 的 workers 解读平均核心数；并行池不能按单线程组直接减少内存。")
        lines.append("")
    lines += ["优化解释：",
              "1. CPU效率低不一定是算法慢：可能申请核心过多，也可能等待CRAM读取/文件系统。",
              "2. 当前一个组内含多个区域；无法从seff获取每个区域单独的耗时和内存。",
              "3. 内存参考下限只针对这批样本；按集群每核内存约束同时调整CPU/内存。",
              "4. ModelMerge/Gather 各只有一个任务，建议先小规模验证再改资源。",
              "5. 多轮 trace 中同一组的多个成功作业会各自保留；限定 --trace 可只分析选定轮次。",
              "6. 缺失/零记账与失败作业不能用于推断安全内存；未自动回填为零。"]
    return "\n".join(lines) + "\n"


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--project-dir", type=Path, default=Path.cwd())
    parser.add_argument("--trace", nargs="+", type=Path, help="Explicit trace files; default all trace-*.tsv")
    parser.add_argument("--outdir", type=Path, help="Default PROJECT/seff_report")
    parser.add_argument("--timeout", type=int, default=30, help="Seconds per seff query")
    args = parser.parse_args()
    if not shutil.which("seff"):
        parser.error("seff 未找到，请在有 Slurm/seff 的登录节点运行。")
    paths = args.trace or sorted(args.project_dir.glob("trace-*.tsv"), key=lambda p: p.stat().st_mtime)
    if not paths:
        parser.error("找不到 trace-*.tsv；用 --trace 指定 trace 文件。")
    jobs, skipped = load_jobs(paths)
    if not jobs:
        parser.error("trace 中没有可查询的 CaVEMan Slurm 作业号。")
    out = args.outdir or args.project_dir / "seff_report"
    rawdir = out / "raw_seff"
    rawdir.mkdir(parents=True, exist_ok=True)
    print("读取 {} 个trace，{} 个唯一作业；缺失作业号记录 {}".format(len(paths), len(jobs), skipped), flush=True)
    rows = []
    for i, row in enumerate(jobs, 1):
        try:
            env = dict(os.environ, LC_ALL="C", LANG="C")
            result = subprocess.run(["seff", row["job_id"]], capture_output=True, text=True,
                                    timeout=args.timeout, env=env)
            raw = result.stdout + ("\nSTDERR:\n"+result.stderr if result.stderr else "")
            metrics = parse_seff(result.stdout)
            if result.returncode:
                metrics["notes"] = "seff_exit_{};".format(result.returncode) + metrics["notes"]
        except subprocess.TimeoutExpired:
            raw = "seff timed out after {} seconds\n".format(args.timeout)
            metrics = dict(slurm_state="UNKNOWN", notes="seff_timeout")
        (rawdir / (row["job_id"] + ".txt")).write_text(raw, encoding="utf-8")
        rows.append(dict(row, **metrics))
        if i % 20 == 0 or i == len(jobs):
            print("seff: {}/{}".format(i, len(jobs)), flush=True)
    columns = ["stage", "task", "job_id", "trace_statuses", "slurm_state", "allocated_cpus",
               "wall_hours", "cpu_hours", "average_used_cores", "cpu_efficiency_pct",
               "peak_memory_gib", "allocated_memory_gib", "memory_efficiency_pct",
               "notes", "nextflow_hash", "trace_files"]
    rows.sort(key=lambda r: (r["stage"], r["task"], r["job_id"]))
    with (out / "jobs.tsv").open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=columns, delimiter="\t")
        writer.writeheader()
        for row in rows:
            writer.writerow({k: ("{:.4f}".format(v) if isinstance(v, float) else v)
                             for k, v in row.items()})
    summary = summarize(rows)
    (out / "summary.txt").write_text(summary, encoding="utf-8")
    print(summary)
    print("报告目录:", out.resolve())


if __name__ == "__main__":
    try:
        main()
    except (OSError, ValueError) as exc:
        sys.exit("错误: " + str(exc))
