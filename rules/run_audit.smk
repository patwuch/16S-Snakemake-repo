# Writes how this run was launched into the Snakemake log itself, so a run can
# be audited from its log alone: launcher, cgroup memory limits, code state,
# and the container images with their staleness check.
import subprocess
from datetime import datetime
from snakemake.logging import logger


def _audit_cmd(*cmd):
    p = subprocess.run(cmd, cwd=WORKFLOW_DIR, capture_output=True, text=True)
    return (p.stdout + p.stderr).rstrip(), p.returncode


def _audit_cgroup():
    rel = Path("/proc/self/cgroup").read_text().strip().split("::", 1)[-1]
    base = Path("/sys/fs/cgroup") / rel.lstrip("/")
    limits = {f: (base / f).read_text().strip() if (base / f).exists() else "n/a"
              for f in ("memory.max", "memory.swap.max")}
    return rel, limits


onstart:
    lines = ["RUN AUDIT", "========="]

    wrapper = os.environ.get("RUN_SNAKEMAKE_WRAPPER")
    lines.append(f"Launcher: {wrapper if wrapper else 'snakemake called directly (not scripts/run_snakemake.sh)'}")
    lines.append(f"Working directory: {os.getcwd()}")
    cgroup, limits = _audit_cgroup()
    lines.append(f"Cgroup: {cgroup}")
    lines.append(f"  memory.max={limits['memory.max']}  memory.swap.max={limits['memory.swap.max']}")
    for k in sorted(os.environ):
        if k.startswith(("SNAKEMAKE_", "APPTAINER", "SINGULARITY")) or k == "SKIP_CONTAINER_CHECK":
            lines.append(f"  env {k}={os.environ[k]}")

    head, _ = _audit_cmd("git", "rev-parse", "HEAD")
    branch, _ = _audit_cmd("git", "rev-parse", "--abbrev-ref", "HEAD")
    dirty, _ = _audit_cmd("git", "status", "--porcelain")
    lines.append(f"Git: {branch} @ {head}")
    lines.append("  uncommitted: " + ("none" if not dirty else ""))
    lines += [f"    {l}" for l in dirty.splitlines()]

    lines.append("Containers:")
    for img in sorted({r.container_img for r in workflow.rules if r.container_img}):
        p = Path(img)
        if p.exists():
            st = p.stat()
            built = datetime.fromtimestamp(st.st_mtime).isoformat(timespec="seconds")
            lines.append(f"  {img}  (built {built}, {st.st_size} bytes)")
        else:
            lines.append(f"  {img}  (MISSING)")
    check, rc = _audit_cmd("scripts/check_containers.sh")
    lines.append("  staleness check: " + ("all match their sources" if rc == 0 else ""))
    lines += [f"    {l}" for l in check.splitlines() if rc != 0]

    logger.info("\n".join(lines))
