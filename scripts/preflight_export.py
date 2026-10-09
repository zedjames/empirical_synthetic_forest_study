#!/usr/bin/env python3
"""Preflight a candidate public export. Does not copy, publish, or grant licenses.

Example:
  python scripts/preflight_export.py dist/HFD06 --inventory review/export_inventory.json
This is a conservative file-integrity and disclosure *aid*, not a substitute
for an owner reviewing the code and all third-party licenses.
"""
import argparse
import hashlib
import json
from pathlib import Path
import re
import sys

BLOCKED_PARTS = {".git", ".ssh", ".aws", ".cache", ".inputs", "__pycache__", ".venv", "venv"}
BLOCKED_SUFFIXES = {".pem", ".p12", ".pfx", ".key", ".env", ".sqlite"}
TEXT_SUFFIXES = {".py", ".md", ".txt", ".json", ".toml", ".yaml", ".yml", ".tex", ".sh", ".cff", ".csv"}
SUSPICIOUS = [
    re.compile(r"gh[pousr]_[A-Za-z0-9]{25,}"),
    re.compile(r"github_pat_[A-Za-z0-9_]{25,}"),
    re.compile(r"AKIA[0-9A-Z]{16}"),
    re.compile(r"-----BEGIN (?:OPENSSH|RSA |EC |DSA |)PRIVATE KEY-----"),
    re.compile(r"(?i)(?:/Users/|/home/)[^\s,\"']{4,}"),
]
RAW_INPUT_BASENAMES = {
    "hf253-05-stems-2014.csv", "hf253-06-stems-2019.csv",
    "hf253-04-trees-2014.csv", "hf253-02-species-codes.csv",
}

def sha256_file(path):
    h = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()

def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("export_dir", type=Path)
    parser.add_argument("--inventory", type=Path, help="Local review inventory with files mapping")
    parser.add_argument("--report", type=Path, help="Write JSON audit report")
    a = parser.parse_args()
    root = a.export_dir.resolve(strict=True)
    if not root.is_dir():
        parser.error("Export root must be a directory")
    errors, warnings, observed = [], [], {}
    manifest = {}
    if a.inventory:
        obj = json.loads(a.inventory.read_text(encoding="utf8"))
        manifest = obj.get("files", {})
        if not isinstance(manifest, dict):
            parser.error("Inventory files must be a filename-to-metadata map")
    for path in sorted(root.rglob("*")):
        rel = path.relative_to(root).as_posix()
        if path.is_symlink():
            errors.append(f"SYMLINK: {rel}")
            continue
        if path.is_dir():
            continue
        parts = Path(rel).parts
        if any(part in BLOCKED_PARTS for part in parts):
            errors.append(f"PRIVATE/CACHE PATH: {rel}")
        if path.name.lower() in RAW_INPUT_BASENAMES:
            errors.append(f"EXTERNAL RAW DATA: {rel}; retrieve from publisher instead")
        if path.suffix.lower() in BLOCKED_SUFFIXES or path.name.lower().startswith(".env"):
            errors.append(f"SENSITIVE SUFFIX: {rel}")
        if path.stat().st_size > 100_000_000:
            warnings.append(f"VERY LARGE FILE: {rel} ({path.stat().st_size} bytes)")
        digest = sha256_file(path)
        observed[rel] = {"sha256": digest, "bytes": path.stat().st_size}
        expected = manifest.get(rel)
        if expected is None and manifest:
            errors.append(f"NOT ALLOWLISTED: {rel}")
        if isinstance(expected, dict):
            reference = expected.get("sha256") or expected.get("export_sha256")
            if reference and reference.lower() != digest.lower():
                errors.append(f"SHA256 MISMATCH: {rel}")
            status = expected.get("redistribution", "")
            if "NOT_AUTHORIZED" in status or "pending" in status.lower() or "Review only" in status:
                errors.append(f"LICENSE/DISCLOSURE UNAPPROVED: {rel}: {status}")
        if path.suffix.lower() in TEXT_SUFFIXES and path.stat().st_size <= 5_000_000:
            content = path.read_text(encoding="utf-8", errors="replace")
            for rule in SUSPICIOUS:
                if rule.search(content):
                    errors.append(f"SENSITIVE TEXT PATTERN: {rel}: {rule.pattern[:35]}")
    for key in manifest:
        if key not in observed:
            errors.append(f"ALLOWLIST FILE MISSING: {key}")
    report = {
        "status": "FAILED" if errors else "PREFLIGHT_PASS_PENDING_HUMAN_APPROVAL",
        "files": len(observed), "errors": errors, "warnings": warnings,
        "file_sha256": observed,
        "note": "Passing a static preflight does not authorize publication or validate science.",
    }
    if a.report:
        a.report.parent.mkdir(parents=True, exist_ok=True)
        a.report.write_text(json.dumps(report, indent=2) + "\n", encoding="utf8")
    print(json.dumps({k:v for k,v in report.items() if k != "file_sha256"}, indent=2))
    return 1 if errors else 0

if __name__ == "__main__":
    sys.exit(main())
