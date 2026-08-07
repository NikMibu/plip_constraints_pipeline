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

try:
    # Share the pipeline's own resolver so the two cannot drift apart.
    from plip_pipeline.utils import resolve_micromamba
except Exception:  # missing deps are reported by the core checks below
    resolve_micromamba = None

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
    configured = mm.get("executable", "micromamba")
    env = mm.get("plip_env", "plip")

    # Uses the pipeline's own resolver rather than a second implementation.
    # The previous version of this script checked Docker while the pipeline
    # had long since moved to micromamba; sharing the code is what stops that
    # from happening again.
    exe = resolve_micromamba(configured) if resolve_micromamba else shutil.which(configured)

    if exe:
        rep.add("micromamba", True, detail=exe)
        if rep.add(f"environment '{env}'", micromamba_env_exists(exe, env),
                   hint=f"micromamba create -n {env} -c conda-forge python=3.9 openbabel=3.1.1"):
            rep.add(f"plip inside '{env}'", tool_in_env(exe, env, "plip"),
                    hint=f"micromamba run -n {env} pip install plip==3.0.0")
    else:
        rep.add("micromamba", False, optional=True,
                hint="not found. A standard install is a shell function, invisible to "
                     "Python - export $MAMBA_EXE or set micromamba.executable")
        found = shutil.which("plip")
        rep.add("plip on PATH (fallback)", found is not None,
                detail=f"{found} - WARNING: not the pinned 3.0.0 environment" if found else "",
                hint="no micromamba and no plip - see SETUP.md section 2")

    # Probe the endpoints the pipeline actually calls, not the bare hosts.
    # rest.uniprot.org answers 403 at its root while the FASTA endpoint works,
    # so checking the host would report a failure that does not exist.
    uniprot_id = (cfg or {}).get("target", {}).get("uniprot_id", "P00918")
    for label, url in [
        ("RCSB entry data", "https://data.rcsb.org/rest/v1/core/entry/1CRN"),
        ("RCSB file download", "https://files.rcsb.org/download/1CRN.pdb"),
        (f"UniProt {uniprot_id}", f"https://rest.uniprot.org/uniprotkb/{uniprot_id}.fasta"),
    ]:
        reachable, detail = check_url(url)
        rep.add(label, reachable, detail=detail, hint=f"{url} - {detail}")

    ok, detail = check_rcsb_search(uniprot_id)
    rep.add("RCSB search", ok, detail=detail,
            hint=f"search query failed or returned nothing: {detail}")


def boltz_version(exe, env_name, explicit=None):
    """Return the reported Boltz version and how it was found, or (None, None).

    Boltz does not have to live in a micromamba environment - a plain venv is
    the route that was actually tested - so an explicit executable path comes
    first. Order: $BOLTZ_EXE, boltz.executable, the micromamba environment,
    PATH.
    """
    attempts = []

    for candidate in (os.environ.get("BOLTZ_EXE"), explicit):
        if candidate:
            path = os.path.expanduser(os.path.expandvars(candidate))
            if os.path.exists(path):
                attempts.append([path, "--version"])

    if exe and env_name:
        attempts.append([exe, "run", "-n", env_name, "boltz", "--version"])
    if shutil.which("boltz"):
        attempts.append(["boltz", "--version"])

    for cmd in attempts:
        try:
            r = subprocess.run(cmd, capture_output=True, timeout=180, text=True)
            text = (r.stdout + r.stderr).strip()
            if r.returncode == 0 and text:
                return text.splitlines()[-1].strip(), cmd[0]
        except Exception:
            continue
    return None, None


def check_boltz(rep, cfg, exe):
    """Boltz-2 is not called by the pipeline, but it consumes its output.

    Reported as optional: generating YAMLs is useful without it.
    """
    print("\nBoltz-2 prediction (optional)")
    env = (cfg or {}).get("micromamba", {}).get("boltz_env", "boltz")
    explicit = (cfg or {}).get("boltz", {}).get("executable")

    version, via = boltz_version(exe, env, explicit)
    if version:
        expected = "2.2.1"
        rep.add("boltz", True, detail=f"{version} (via {via})", optional=True)
        if expected not in version:
            rep.add(f"boltz == {expected}", False, optional=True,
                    hint=f"thesis used {expected}; output format differs between versions")
    else:
        rep.add("boltz", False, optional=True,
                hint="not found. Set boltz.executable to the binary (a venv works), "
                     "or micromamba.boltz_env to an environment that has it. "
                     "Needed only to run the generated YAMLs.")


def check_diffdock(rep, cfg):
    print("\nDiffDock workflow (optional)")
    mm = (cfg or {}).get("micromamba", {})
    configured = mm.get("executable", "micromamba")
    env = mm.get("diffdock_env", "diffdock")
    exe = resolve_micromamba(configured) if resolve_micromamba else shutil.which(configured)

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

    if exe:
        rep.add(f"environment '{env}'", micromamba_env_exists(exe, env),
                hint="see SETUP.md section 3", optional=True)

    try:
        import rdkit  # noqa: F401
        rep.add("RDKit", True, optional=True)
    except ImportError:
        rep.add("RDKit", False, hint="needed to build 3D ligands", optional=True)


def check_url(url):
    """GET probe for one endpoint."""
    try:
        import requests
        r = requests.get(url, timeout=15)
        return r.status_code < 400, f"HTTP {r.status_code}"
    except Exception as e:
        return False, type(e).__name__


def check_rcsb_search(uniprot_id):
    """Run a real search query, the way PDBFetcher._search_pdb_ids does.

    The endpoint only accepts POST and answers a bare GET with 400, so probing
    it like a normal URL says nothing. Sending an actual query - narrowed to a
    single hit - confirms the service works and the query shape is still valid.
    """
    try:
        import requests
        query = {
            "query": {
                "type": "terminal",
                "service": "text",
                "parameters": {
                    "attribute": "rcsb_polymer_entity_container_identifiers"
                                 ".reference_sequence_identifiers.database_accession",
                    "operator": "exact_match",
                    "value": uniprot_id,
                },
            },
            "request_options": {"paginate": {"start": 0, "rows": 1}},
            "return_type": "entry",
        }
        r = requests.post("https://search.rcsb.org/rcsbsearch/v2/query",
                          json=query, timeout=30)
        if r.status_code >= 400:
            return False, f"HTTP {r.status_code}"
        hits = r.json().get("total_count", 0)
        return hits > 0, f"{hits} structures for {uniprot_id}"
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
        global check_url, check_rcsb_search
        check_url = lambda url: (True, "skipped")
        check_rcsb_search = lambda uniprot_id: (True, "skipped")

    check_crystal(rep, cfg)

    configured = (cfg or {}).get("micromamba", {}).get("executable", "micromamba")
    mamba = resolve_micromamba(configured) if resolve_micromamba else shutil.which(configured)
    check_boltz(rep, cfg, mamba)
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
