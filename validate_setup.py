#!/usr/bin/env python3
"""
Setup Validation Script

Checks if all requirements are met before running the pipeline.
"""
import subprocess
import sys
import os


def check_python_version():
    """Check Python version."""
    print("Checking Python version...", end=" ")
    version = sys.version_info
    if version.major >= 3 and version.minor >= 8:
        print(f"✓ {version.major}.{version.minor}.{version.micro}")
        return True
    else:
        print(f"✗ {version.major}.{version.minor}.{version.micro} (need 3.8+)")
        return False


def check_module(module_name, import_name=None):
    """Check if Python module is installed."""
    import_name = import_name or module_name
    print(f"Checking {module_name}...", end=" ")
    try:
        __import__(import_name)
        print("✓")
        return True
    except ImportError:
        print(f"✗ (run: pip install {module_name})")
        return False


def check_docker():
    """Check if Docker is available."""
    print("Checking Docker...", end=" ")
    try:
        result = subprocess.run(
            ["docker", "--version"],
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            timeout=5
        )
        if result.returncode == 0:
            version = result.stdout.decode().strip()
            print(f"✓ {version}")
            return True
        else:
            print("✗ Docker not found")
            return False
    except FileNotFoundError:
        print("✗ Docker not installed")
        return False
    except subprocess.TimeoutExpired:
        print("✗ Docker timeout")
        return False


def check_plip_image():
    """Check if PLIP Docker image exists."""
    print("Checking PLIP Docker image...", end=" ")
    try:
        result = subprocess.run(
            ["docker", "images", "-q", "pharmai/plip"],
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            timeout=10
        )
        if result.returncode == 0 and result.stdout.strip():
            print("✓")
            return True
        else:
            print("✗ (run: docker pull pharmai/plip)")
            return False
    except Exception:
        print("✗ Failed to check")
        return False


def check_config_file():
    """Check if config file exists."""
    print("Checking config.yaml...", end=" ")
    if os.path.exists("config.yaml"):
        print("✓")
        return True
    else:
        print("✗ (create config.yaml)")
        return False


def check_internet():
    """Check internet connectivity."""
    print("Checking internet connection...", end=" ")
    try:
        import requests
        response = requests.get("https://www.rcsb.org", timeout=5)
        if response.status_code == 200:
            print("✓")
            return True
        else:
            print(f"✗ RCSB PDB unreachable (HTTP {response.status_code})")
            return False
    except Exception:
        print("✗ No internet connection")
        return False


def main():
    """Run all checks."""
    print("=" * 60)
    print("  PLIP Constraints Pipeline - Setup Validation")
    print("=" * 60)
    print()
    
    checks = [
        ("Python Version", check_python_version),
        ("BioPython", lambda: check_module("biopython", "Bio")),
        ("Pandas", lambda: check_module("pandas")),
        ("Requests", lambda: check_module("requests")),
        ("PyYAML", lambda: check_module("pyyaml", "yaml")),
        ("Docker", check_docker),
        ("PLIP Image", check_plip_image),
        ("Config File", check_config_file),
        ("Internet", check_internet),
    ]
    
    results = []
    
    for name, check_func in checks:
        try:
            result = check_func()
            results.append((name, result))
        except Exception as e:
            print(f"✗ Error: {e}")
            results.append((name, False))
    
    print()
    print("=" * 60)
    
    passed = sum(1 for _, result in results if result)
    total = len(results)
    
    if passed == total:
        print(f"✓ All checks passed ({passed}/{total})")
        print()
        print("You're ready to run the pipeline!")
        print("  python pipeline.py --config config.yaml")
        sys.exit(0)
    else:
        print(f"✗ Some checks failed ({passed}/{total} passed)")
        print()
        print("Failed checks:")
        for name, result in results:
            if not result:
                print(f"  - {name}")
        print()
        print("Fix the issues above before running the pipeline.")
        sys.exit(1)


if __name__ == "__main__":
    main()

