"""Hash-checked external CSV parser with lossless raw byte preservation."""
import csv,io,json
import h05_common as h
def load(dataset,prefix):
    h.guard("external_protocol");h.guard("external_encoding")
    sources=json.loads((h.ROOT/"external_sources/source_inventory.json").read_text())["sources"]
    selected=[r for r in sources if r["dataset"]==dataset and r["name"].startswith(prefix)]
    if len(selected)!=1:raise ValueError("Ambiguous external source")
    record=selected[0];path=h.ROOT/record["local_input"]
    if h.sha(path)!=record["sha256"]:raise ValueError("External input changed")
    raw=path.read_bytes();encoding="latin-1" if dataset=="HF453" else "ascii";text=raw.decode(encoding)
    if text.encode(encoding)!=raw:raise ValueError("Raw byte parser not reversible")
    rows=list(csv.DictReader(io.StringIO(text,newline="")))
    expected=set(record["schema"])
    if not rows or set(rows[0])!=expected:raise ValueError("Undeclared external schema")
    for row in rows:
        if None in row or any(v is None for v in row.values()):raise ValueError("Malformed external CSV row")
        fields=["StemTag","year","status"] if dataset=="HF453" else [k for k in ["seedlingID","yearOfOb","status","flag","sampled","alive","date","taxonCode","coordinate"] if k in row]
        for key in fields:row[key].encode("ascii")
    return rows,record
def inspect():
    from collections import Counter
    outputs=[]
    for dataset,prefix in [("HF453","hf453-01"),("HF355","hf355-01"),("HF355","hf355-02"),("HF355","hf355-05")]:
        rows,record=load(dataset,prefix)
        value=dict(dataset=dataset,name=record["name"],source_sha256=record["sha256"],rows=len(rows),schema=list(rows[0]),
            field_support={k:dict(Counter(r[k] for r in rows)) for k in ["year","status","flag","sampled","alive","trueRecrt"] if k in rows[0]},
            first_outcome_inspection_after_protocol_commit="35c4c384ad",parser_amendment="Lossless explicit byte mapping; no model/target changes",model_fitting=False,biological_scores_computed=False)
        outputs.append(value);print(json.dumps(value,sort_keys=True),flush=True)
    h.write_json("external_sources/record_support_inventory.json",outputs)
if __name__=="__main__":inspect()
