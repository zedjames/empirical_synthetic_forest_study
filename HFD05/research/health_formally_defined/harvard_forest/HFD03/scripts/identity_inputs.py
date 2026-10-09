"""Small, separately provenanced public inputs. Never writes HFD01/HFD02S."""
import hashlib
import json
from pathlib import Path
import urllib.request

ROOT = Path(__file__).resolve().parents[1]
INPUTS = [
    ("hf253-04-trees-2014.csv", 9699549, "53bd19e300a99c0677a265a5ac58f9a6"),
    ("hf253-02-species-codes.csv", 5740, "8ab6408dd633a821508e26e9beae1337"),
]


def main():
    folder = ROOT / ".inputs"
    folder.mkdir(exist_ok=True)
    records = []
    for name, size, md5 in INPUTS:
        path = folder / name
        url = "https://harvardforest.fas.harvard.edu/data/p25/hf253/" + name
        if not path.exists():
            request = urllib.request.Request(url, headers={"User-Agent": "HFD03-reproducible-research/1"})
            with urllib.request.urlopen(request, timeout=60) as response:
                raw = response.read(size + 1)
            if len(raw) != size or hashlib.md5(raw).hexdigest() != md5:
                raise ValueError("Publisher length/MD5 mismatch: " + name)
            path.write_bytes(raw)
        raw = path.read_bytes()
        if len(raw) != size or hashlib.md5(raw).hexdigest() != md5:
            raise ValueError("Input identity changed")
        records.append(dict(name=name, url=url, bytes=size, publisher_md5=md5,
                            sha256=hashlib.sha256(raw).hexdigest(), rows=len(raw.splitlines())-1,
                            provenance="OBSERVED", retrieval_date="2026-10-07",
                            metadata_anchor="HFD01 HF253-EML 4469f4b3f1950612",
                            permission="Public Harvard Forest research archive; attribution required"))
    (ROOT / "provenance").mkdir(exist_ok=True)
    (ROOT / "provenance/identity_inputs.json").write_text(json.dumps(dict(
        sources=records, unavailable="No E1 tree table listed in HF253 version 6; E1 tree grouping is DERIVED",
        licensing_metadata="Read-only frozen HF253 EML; no new whole-plant death inference"), indent=2)+"\n")
    print(json.dumps(records, indent=2))


if __name__ == "__main__":
    main()
