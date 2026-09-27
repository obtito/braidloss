"""Run one native COMSOL job, retaining every input, log and resource sample.

Only the COMSOL batch command saves the MPH. No COMSOL security settings change.
Heavy jobs are serial; unrelated user processes are never altered.
"""
from pathlib import Path
import hashlib
import json
import subprocess
import sys
import time
import psutil

BIN = Path(r"D:/tools/comosol/COMSOL62/Multiphysics/bin/win64")


def run(java: Path, threads=4):
    java = java.resolve()
    job = java.parent
    record = {"java": java.name, "sha256": hashlib.sha256(java.read_bytes()).hexdigest(),
              "threads": threads, "started_unix": time.time(), "status": "compiling"}
    compilation = subprocess.run([str(BIN / "comsolcompile.exe"), str(java)], capture_output=True)
    (job / "compile.log").write_bytes(compilation.stdout + compilation.stderr)
    if compilation.returncode or not java.with_suffix(".class").exists():
        record.update(status="compile_failed", exit_code=compilation.returncode)
        (job / "run.json").write_text(json.dumps(record, indent=2), encoding="utf-8")
        raise RuntimeError(f"Compilation failed: {job / 'compile.log'}")
    t0 = time.monotonic()
    samples = []
    with (job / "stdout.log").open("wb") as stream:
        p = subprocess.Popen([str(BIN / "comsolbatch.exe"), "-np", str(threads),
                              "-inputfile", str(java.with_suffix('.class')),
                              "-outputfile", str(job / "model.mph"),
                              "-batchlog", str(job / "solver.log")],
                             stdout=stream, stderr=subprocess.STDOUT)
        parent = psutil.Process(p.pid)
        low_memory_since = None
        while p.poll() is None:
            rss = private = 0
            try:
                family = [parent] + parent.children(recursive=True)
                for proc in family:
                    try:
                        mi = proc.memory_info()
                        rss += mi.rss
                        private += getattr(mi, "private", mi.vms)
                    except psutil.Error:
                        pass
            except psutil.Error:
                family = []
            free = psutil.virtual_memory().available
            samples.append([time.monotonic() - t0, rss, private, free])
            if free < 400 * 1024**2:
                low_memory_since = low_memory_since or time.monotonic()
                if time.monotonic() - low_memory_since > 20:
                    for proc in reversed(family):
                        try:
                            proc.terminate()
                        except psutil.Error:
                            pass
                    record["resource_stop"] = "Available RAM below 400 MiB for 20 s; stopped this job only"
                    p.wait()
                    break
            else:
                low_memory_since = None
            time.sleep(1)
    log = (job / "stdout.log").read_text(encoding="utf-8", errors="replace")
    record.update(elapsed_seconds=time.monotonic() - t0, exit_code=p.returncode,
                  peak_family_rss_GiB=max((s[1] for s in samples), default=0) / 1024**3,
                  peak_family_private_GiB=max((s[2] for s in samples), default=0) / 1024**3,
                  minimum_available_GiB=min((s[3] for s in samples), default=0) / 1024**3,
                  status="solved" if "SOLVE_COMPLETE" in log else "incomplete",
                  mph=[f.name for f in job.glob("*.mph")])
    (job / "resources.json").write_text(json.dumps({"columns": ["seconds", "rss_bytes", "private_bytes", "available_bytes"], "samples": samples}), encoding="utf-8")
    (job / "run.json").write_text(json.dumps(record, indent=2), encoding="utf-8")
    print(json.dumps(record, indent=2), flush=True)
    print("\n".join(line for line in log.splitlines() if line.startswith("METRIC|") or
        (line.startswith("STAGE|") and not line.startswith(("STAGE|curve|","STAGE|solid|","STAGE|selection|")))), flush=True)
    return record


if __name__ == "__main__":
    run(Path(sys.argv[1]))
