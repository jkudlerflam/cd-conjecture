#!/usr/bin/env python3
"""Stage working-tree files and run one calculation per IAS Slurm job."""
import argparse
import hashlib
import json
from pathlib import Path
import re
import shlex
import subprocess
import sys
import tarfile
import tempfile


HERE = Path(__file__).resolve().parent
ROOT = HERE.parent.parent
LOCAL = ROOT / "runs/typhon"
REMOTE = "cd-conjecture/runs"  # Relative to the IAS user's home.
HOST = "ias-cluster"
SSH = ["ssh", "-o", "BatchMode=yes", "-o", "ForwardAgent=no", "-o", "ConnectTimeout=15"]
DEFAULTS = {"cpus-per-task": 1, "mem": "2G", "time": "00:05:00", "qos": "short"}
ACCOUNTING = ('sacct -X -P -j "$job" --format='
              'JobIDRaw,JobName,State,ExitCode,Elapsed,AllocCPUS,ReqMem,Timelimit,NodeList,Partition,Account,QOS,Submit,Start,End')


def ssh(command, **options):
    return subprocess.run([*SSH, HOST, command], check=True, **options)


def project_file(name):
    """Only individually named regular files from the checkout can enter a snapshot."""
    if not isinstance(name, str):
        raise ValueError("Each file name must be a string")
    rel = Path(name)
    path = ROOT / rel
    if rel.is_absolute() or ".." in rel.parts or path.resolve() != path or not path.is_file():
        raise ValueError(f"Expected a repository-relative regular file without symlinks: {name}")
    excluded = {".git", ".ssh", ".aws", ".venv", "venv", "site-packages", "__pycache__", ".DS_Store"}
    binary = {".so", ".dylib", ".dll", ".pyd", ".a", ".o", ".pyc", ".pyo", ".whl", ".egg"}
    if (LOCAL in path.parents or excluded.intersection(rel.parts) or binary.intersection(path.suffixes)
            or any((p / "pyvenv.cfg").is_file() or (p / "conda-meta").is_dir()
                   for p in (ROOT / parent for parent in rel.parents))):
        raise ValueError(f"Environments, compiled packages and local run files cannot be staged: {name}")
    return rel


def digest(path):
    with path.open("rb") as f:
        return hashlib.file_digest(f, "sha256").hexdigest()


def submit(name, recipe):
    recipe = project_file(recipe)
    recipe_bytes = (ROOT / recipe).read_bytes()
    spec = json.loads(recipe_bytes)
    if not isinstance(spec, dict):
        raise ValueError("Expected a JSON object for the job recipe")
    extra = set(spec) - {"files", "command", "resources", "seed"}
    if extra:
        raise ValueError(f"Unknown recipe fields: {sorted(extra)}")
    command = spec.get("command")
    if not isinstance(command, list) or not command:
        raise ValueError("command must be a nonempty argument list")
    if not all(isinstance(x, str) for x in command) or not command[0]:
        raise ValueError("command must be a nonempty argument list")
    if not isinstance(spec.get("files"), list):
        raise ValueError("files must list the source and inputs relative to the repository root")
    files = sorted({recipe, *(project_file(path) for path in spec["files"])})
    overrides = spec.get("resources", {})
    if not isinstance(overrides, dict) or set(overrides) - DEFAULTS.keys():
        raise ValueError(f"Configurable Slurm resources are {', '.join(DEFAULTS)}")
    resources = {"account": "physics", "partition": "typhon", "nodes": 1, "ntasks": 1, **DEFAULTS, **overrides}
    local = LOCAL / name
    if local.exists():
        raise ValueError(f"{local} already exists. Use status or fetch to revisit that run.")

    helpers = [HERE / "batch.sh", HERE / "environment.sh"]
    selected = [*files, *(p.relative_to(ROOT) for p in helpers)]
    git = ["git", "-C", str(ROOT)]
    manifest = {
        "run": name, "host": HOST, "remote_dir": f"~/{REMOTE}/{name}",
        "git_head": subprocess.check_output([*git, "rev-parse", "HEAD"], text=True).strip(),
        "git_status_selected": subprocess.check_output(
            [*git, "status", "--short", "--untracked-files=all", "--", *selected], text=True).splitlines(),
        "recipe": str(recipe), "files": [str(p) for p in files],
        "command": command, "cwd": "src", "seed": spec.get("seed"), "resources": resources,
    }
    sbatch = ["sbatch", "--parsable", f"--job-name=cd-{name}",
              "--output=logs/%j.out", "--error=logs/%j.err",
              *(f"--{key}={value}" for key, value in resources.items()), "batch.sh"]
    with tempfile.TemporaryDirectory() as tmp:
        # Archive the same recipe bytes from which the command and resources were read.
        recipe_copy = Path(tmp) / "recipe"
        recipe_copy.write_bytes(recipe_bytes)
        manifest_path = Path(tmp) / "manifest.json"
        manifest_text = json.dumps(manifest, indent=2) + "\n"
        manifest_path.write_text(manifest_text)
        archive = Path(tmp) / "snapshot.tar.gz"
        with tarfile.open(archive, "w:gz") as tar:
            for path in files:
                source = recipe_copy if path == recipe else ROOT / path
                tar.add(source, arcname=str(Path("src") / path), recursive=False)
            for path in [manifest_path, *helpers]:
                tar.add(path, arcname=path.name, recursive=False)
        checksum = f"{digest(archive)}  snapshot.tar.gz"
        # Exclusive mkdir prevents a retry or a second checkout from submitting the same run twice.
        stage = f"""set -e
umask 077
mkdir -p {REMOTE}
mkdir {REMOTE}/{name}
cd {REMOTE}/{name}
cat > snapshot.tar.gz
echo '{checksum}' > snapshot.sha256
sha256sum -c snapshot.sha256 >&2
tar -xzf snapshot.tar.gz
mkdir logs results
{shlex.join(sbatch)} > job-id
cat job-id"""
        with archive.open("rb") as stream:
            receipt = ssh(stage, stdin=stream, stdout=subprocess.PIPE, text=True).stdout.strip()
    print(f"{name}  Slurm {receipt}  {HOST}:~/{REMOTE}/{name}", flush=True)
    local.mkdir(parents=True, exist_ok=True)
    (local / "job-id").write_text(receipt + "\n")
    (local / "manifest.json").write_text(manifest_text)
    (local / "snapshot.sha256").write_text(checksum + "\n")


def job_context(name):
    return f"""set -e
cd {REMOTE}/{name}
job=$(cut -d';' -f1 job-id)
case "$job" in ''|*[!0-9]*) echo 'No Slurm job ID recorded in this run' >&2; exit 1;; esac
"""


def fetch(name):
    # Query first. If a job finishes during a live copy, its receipt still labels it as live.
    accounting = ssh(job_context(name) + ACCOUNTING, stdout=subprocess.PIPE, text=True).stdout
    local = LOCAL / name
    local.mkdir(parents=True, exist_ok=True)
    includes = ("job-id", "manifest.json", "snapshot.tar.gz", "snapshot.sha256", "logs/", "logs/**", "results/", "results/**")
    subprocess.run(["rsync", "-rtc", "-e", shlex.join(SSH), *(f"--include=/{p}" for p in includes),
                    "--exclude=*", f"{HOST}:{REMOTE}/{name}/", str(local) + "/"], check=True)
    if digest(local / "snapshot.tar.gz") != (local / "snapshot.sha256").read_text().split()[0]:
        raise ValueError("Downloaded snapshot differs from its recorded SHA-256")
    (local / "accounting.psv").write_text(accounting)
    print(accounting.rstrip(), local, sep="\n")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    actions = parser.add_subparsers(dest="action", required=True)
    actions.add_parser("setup", help="initialize the Linux environment once on IAS")
    for action in ("submit", "status", "logs", "fetch"):
        sub = actions.add_parser(action)
        sub.add_argument("name", help="unique run directory name")
        if action == "submit":
            sub.add_argument("recipe", help="repository-relative JSON job specification")
    args = parser.parse_args()
    if args.action == "setup":
        with (HERE / "environment.sh").open("rb") as script:
            ssh("bash -se -- setup", stdin=script)
        print("Typhon environment ready.")
        return
    if not re.fullmatch(r"[A-Za-z0-9][A-Za-z0-9_.-]{0,63}", args.name):
        parser.error("use a run name of 1–64 letters, digits, underscores, dots or dashes, starting with a letter or digit")
    if args.action == "submit":
        submit(args.name, args.recipe)
    elif args.action == "fetch":
        fetch(args.name)
    elif args.action == "status":
        ssh(job_context(args.name) + f'echo "Slurm $job"\nsqueue --me --name=cd-{args.name} -o "%.18i %.12T %R"\n' + ACCOUNTING)
    else:
        ssh(job_context(args.name) + """if [ -f "logs/$job.out" ] || [ -f "logs/$job.err" ]; then
    for log in "logs/$job.out" "logs/$job.err"; do
        if [ -f "$log" ]; then tail -n 60 -v "$log"; fi
    done
else
    echo "Slurm $job has no output logs yet."
fi""")


if __name__ == "__main__":
    try:
        main()
    except (OSError, ValueError, subprocess.CalledProcessError) as error:
        sys.exit(str(error))
