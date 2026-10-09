#!/usr/bin/env python3
"""Read-only SHA-256 and archived-member verification of Paper V release assets.
Passing this does not prove scientific replay, model adequacy or release approval.
"""
import hashlib
import json
from pathlib import Path
import sys
import zipfile

KINDS = {"HFD05": ("HFD05_Artifact_Manifest.json", 492),
         "HFD06": ("HFD06_Artifact_Manifest.json", 67)}
def digest(handle):
    h = hashlib.sha256()
    for block in iter(lambda: handle.read(1024 * 1024), b""):
        h.update(block)
    return h.hexdigest()

def verify(root):
    lines = (root / "SHA256SUMS").read_text().splitlines()
    if len(lines) != 2:
        raise ValueError("Exactly two release package checksums required")
    completed = {}
    for row in lines:
        wanted, sep, filename = row.partition("  ")
        if not sep or len(wanted) != 64 or Path(filename).name != filename or not filename.endswith(".zip"):
            raise ValueError("Invalid checksum record")
        kind = "HFD05" if "HFD05" in filename else "HFD06" if "HFD06" in filename else None
        if not kind or kind in completed:
            raise ValueError("Missing or repeated package identity")
        with (root / filename).open("rb") as stream:
            if digest(stream) != wanted:
                raise ValueError("Archive hash mismatch: " + filename)
        marker, expected_size = KINDS[kind]
        with zipfile.ZipFile(root / filename) as archive:
            hits = [p for p in archive.namelist() if p == marker or p.endswith("/"+marker)]
            if len(hits) != 1:
                raise ValueError("Manifest missing or ambiguous: "+filename)
            prefix = hits[0][:-len(marker)]
            files = json.loads(archive.read(hits[0]))["files"]
            if len(files) != expected_size:
                raise ValueError("Frozen inventory size changed")
            for name, meta in files.items():
                with archive.open(prefix+name) as stream:
                    if digest(stream) != (meta.get("sha256") or meta.get("export_sha256")):
                        raise ValueError("Scientific member differs: "+name)
        completed[kind] = {"filename":filename,"files":len(files),"sha256":wanted}
    return completed

if __name__ == "__main__":
    if len(sys.argv) != 2:
        raise SystemExit("Usage: verify_release_assets.py DOWNLOAD_DIR")
    print(json.dumps({"status":"BYTE_INTEGRITY_ONLY","packages":verify(Path(sys.argv[1]))},indent=2))
