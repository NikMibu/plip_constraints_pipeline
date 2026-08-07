#!/usr/bin/env python3
"""Setup validation.

Checks that everything the pipeline actually calls is in place, before a run
starts. Reads config.yaml so the environment names and the DiffDock path are
the ones the pipeline will use, rather than assumptions baked in here.

Checks are grouped:
  core      - needed for every run
  crystal   - needed for --steps crystal (fetch, clean, plip, generate)
  diffdock  - needed for --steps diffdock_full

A missing diffdock requirement is reported as a warning, not a failure: the
crystal workflow is useful on its own and needs no GPU.

Exit codes: 0 all required checks passed, 1 a required check failed.
"""
import argparse
import os
import shutil
import subprocess
import sys

OK = "[ ok ]"
FAIL = "[fail]"
WARN = "[warn]"


class Report:
    def __init__(self):
        self.required = []
        self.optional = []

    def add(self, name, passed, detail="", optional=False, hint=""):
        """detail is shown always, hint only when the check did not pass."""
        mark = OK if passed else (WARN if optional else FAIL)
        suffix = detail if passed else (hint or detail)
        print(f"  {mark} {name}" + (f" - {suffix}" if suffix else ""))
        (self.optional if optional else self.required).append((name, passed))
        return passed


def load_config(path):
    try:
        import yaml
    except ImportError:
        return None, "PyYAML is not installed, cannot read the config"
    if not os.path.exists(path):
        return None, f"{path} not found"
    try:
        with open(path) as f:
            return yaml.safe_load(f), ""
    except Exception as e:
        return None, str(e)


def micromamba_env_exists(exe, env_name):
    """True if `exe` knows an environment called `env_name`."""
    try:
        result = subprocess.run([exe, "env", "list"],
                                capture_output=True, timeout=30, text=True)
        if result.returncode != 0:
            return False
        for line in result.stdout.splitlines():
            fields = line.split()
            if fields and fields[0] == env_name:
                return True
        return False
    except Exception:
        return False


def tool_in_env(exe, env_name, tool):
    """True if `tool` runs inside the micromamba environment."""
    try:
        result = subprocess.run([exe, "run", "-n", env_name, tool, "--help"],
                                capture_output=True, timeout=120)
        return result.returncode in (0, 1, 2)  # --help conventions vary
    except Exception:
        return False


def check_core(rep, config_path):
    print("\nCore")
    v = sys.version_info
    rep.add("Python >= 3.8", v >= (3, 8), f"{v.major}.{v.minor}.{v.micro}")

    for label, module in [("pandas", "pandas"), ("PyYAML", "yaml"),
                          ("requests", "requests"), ("BioPython", "Bio")]:
        try:
            __import__(module)
            rep.add(label, True)
        except ImportError:
            rep.add(label, False, hint="pip install -r requirements.txt")

    cfg, err = load_config(config_path)
    rep.add(f"{os.path.basename(config_path)}", cfg is not None, hint=err)
    return cfg


def check_crystal(rep, cfg):
    print("\nCrystal workflow")
    mm = (cfg or {}).get("micromamba", {})
    exe = mm.get("executable", "micromamba")
    env = mm.get("plip_env", "plip")

    # Mirrors plip_pipeline.plip.resolve_plip_command: micromamba first,
    # then a plain plip on PATH.
    if shutil.which(exe):
        rep.add(f"micromamba ({exe})", True, detail=shutil.which(exe))
        if rep.add(f"environment '{env}'", micromamba_env_exists(exe, env),
                   hint=f"{exe} create -n {env} -c conda-forge python=3.9 openbabel"):
            rep.add(f"plip inside '{env}'", tool_in_env(exe, env, "plip"),
                    hint=f"{exe} run -n {env} pip install plip==3.0.0")
    else:
        rep.add(f"micromamba ({exe})", False, detail="not on PATH", optional=True)
        rep.add("plip on PATH (fallback)", shutil.which("plip") is not None,
                detail=shutil.which("plip") or "",
                hint="pip install plip==3.0.0, or install micromamba")

    for host in ["https://files.rcsb.org", "https://rest.uniprot.org"]:
        reachable, detail = check_url(host)
        rep.add(host, reachable, detail)


def check_diffdock(rep, cfg):
    print("\nDiffDock workflow (optional)")
    mm = (cfg or {}).get("micromamba", {})
    exe = mm.get("executable", "micromamba")
    env = mm.get("diffdock_env", "diffdock")

    repo = os.environ.get("DIFFDOCK_HOME") or os.path.expandvars(
        (cfg or {}).get("diffdock", {}).get("repo_path", "") or ""
    )
    # expandvars leaves unknown variables untouched, so "${DIFFDOCK_HOME}"
    # surviving means the variable is simply not set.
    if "${" in repo or repo.startswith("$"):
        rep.add("DiffDock repository", False,
                hint="DIFFDOCK_HOME is not set - export DIFFDOCK_HOME=/path/to/DiffDock",
                optional=True)
    elif repo:
        exists = os.path.isdir(repo)
        rep.add("DiffDock repository", exists, detail=repo,
                hint=f"{repo} does not exist - set diffdock.repo_path "
                     f"(or $DIFFDOCK_HOME)", optional=True)
        if exists:
            script = os.path.join(repo, "inference.py")
            rep.add("inference.py", os.path.exists(script), detail=script, optional=True)
    else:
        rep.add("DiffDock repository", False, hint="diffdock.repo_path is not set",
                optional=True)

    if shutil.which(exe):
        rep.add(f"environment '{env}'", micromamba_env_exists(exe, env),
                hint="see MICROMAMBA_SETUP.md", optional=True)

    try:
        import rdkit  # noqa: F401
        rep.add("RDKit", True, optional=True)
    except ImportError:
        rep.add("RDKit", False, hint="needed to build 3D ligands", optional=True)


def check_url(url):
    try:
        import requests
        r = requests.get(url, timeout=10)
        return r.status_code < 400, f"HTTP {r.status_code}"
    except Exception as e:
        return False, type(e).__name__


def main():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--config", default="config.yaml", help="pipeline config to validate against")
    p.add_argument("--skip-network", action="store_true",
                   help="do not contact RCSB and UniProt")
    args = p.parse_args()

    print("=" * 66)
    print("  PLIP Constraints Pipeline - setup validation")
    print("=" * 66)

    rep = Report()
    cfg = check_core(rep, args.config)

    if args.skip_network:
        global check_url
        check_url = lambda url: (True, "skipped")

    check_crystal(rep, cfg)
    check_diffdock(rep, cfg)

    passed = sum(1 for _, ok in rep.required if ok)
    total = len(rep.required)
    opt_failed = [n for n, ok in rep.optional if not ok]

    print("\n" + "=" * 66)
    if passed == total:
        print(f"{OK} {passed}/{total} required checks passed")
        if opt_failed:
            print(f"{WARN} optional, not available: {', '.join(opt_failed)}")
            print("      The crystal workflow does not need these:")
            print("        python pipeline.py --steps crystal")
        else:
            print("\nReady:  python pipeline.py --config " + args.config)
        sys.exit(0)
    else:
        print(f"{FAIL} {passed}/{total} required checks passed")
        print("\nFailed:")
        for name, ok in rep.required:
            if not ok:
                print(f"  - {name}")
        sys.exit(1)


if __name__ == "__main__":
    main()
