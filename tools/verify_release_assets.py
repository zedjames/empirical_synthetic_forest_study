#!/usr/bin/env python3
"""Independent, read-only hash and manifest check for published Paper V ZIP assets.

Usage: python scripts/verify_release_assets.py /path/to/downloaded-release
This verifies byte identity and file inventory. It does NOT execute the science.
"""
import hashlib
import json
from pathlib import Path
import sys
import zipfile

EXPECTED = {"HFD05": ("HFD05_Artifact_Manifest.json", 492),
            "HFD06": ("HFD06_Artifact_Manifest.json", 67)}


def sha(stream):
    h = hashlib.sha256()
    while True:
        chunk = stream.read(1024*1024)
        if not chunk:
            return h.hexdigest()
        h.update(chunk)


def verify(root):
    sums = root / "SHA256SUMS"
    lines = sums.read_text().splitlines()
    if len(lines) != 2:
        raise RuntimeError("Expected exactly two frozen scientific ZIP identities")
    reports = {}
    for line in lines:
        digest, delim, filename = line.partition("  ")
        if not delim or len(digest) != 64:
            raise RuntimeError("Invalid SHA256SUMS entry")
        if Path(filename).name != filename or not filename.endswith(".zip"):
            raise RuntimeError("Unsafe asset name")
        kind = "HFD05" if "HFD05" in filename else "HFD06" if "HFD06" in filename else None
        if not kind or kind in reports:
            raise RuntimeError("Ambiguous package identity")
        with (root/filename).open("rb") as handle:
            if sha(handle) != digest:
                raise RuntimeError("Asset SHA mismatch: "+filename)
        manifest_name, expected_count = EXPECTED[kind]
        with zipfile.ZipFile(root/filename) as archive:
            names = archive.namelist()
            hits = [n for n in names if n == manifest_name or n.endswith("/"+manifest_name)]
            if len(hits) != 1:
                raise RuntimeError("Missing package manifest: "+filename)
            prefix = hits[0][:-len(manifest_name)]
            manifest = json.loads(archive.read(hits[0]))
            files = manifest["files"]
            if len(files) != expected_count:
                raise RuntimeError("Unexpected file count: "+filename)
            for path, entry in files.items():
                with archive.open(prefix+path) as stream:
                    actual = sha(stream)
                expected = entry.get("sha256") or entry.get("export_sha256")
                if actual != expected:
                    raise RuntimeError("Member SHA mismatch: "+path)
        reports[kind] = {"asset":filename,"files":len(files),"asset_sha256":digest}
    return reports


if __name__ == "__main__":
    folder = Path(sys.argv[1]) if len(sys.argv)==2 else None
    if not folder or not folder.is_dir():
        raise SystemExit("Usage: verify_release_assets.py DOWNLOAD_DIRECTORY")
    print(json.dumps({"status":"PACKAGE_BYTES_VERIFIED_NOT_SCIENTIFIC_REPLAY",
                      "packages":verify(folder)},indent=2))
