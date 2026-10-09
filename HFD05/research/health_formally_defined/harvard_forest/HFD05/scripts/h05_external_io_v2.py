"""Second explicit source parser: narrow documented aliases, no raw rewrites."""
import csv,io,json
import h05_common as h
def load(dataset,prefix):
    for phase in ["external_protocol","external_encoding","external_schema"]:h.guard(phase)
    sources=json.loads((h.ROOT/"external_sources/source_inventory.json").read_text())["sources"]
    selected=[r for r in sources if r["dataset"]==dataset and r["name"].startswith(prefix)]
    if len(selected)!=1:raise ValueError("Ambiguous external source")
    record=selected[0];path=h.ROOT/record["local_input"]
    if h.sha(path)!=record["sha256"]:raise ValueError("External input changed")
    raw=path.read_bytes();encoding="latin-1" if dataset=="HF453" else "ascii";text=raw.decode(encoding)
    if text.encode(encoding)!=raw:raise ValueError("Raw byte parser not reversible")
    parser=csv.DictReader(io.StringIO(text,newline=""));header=parser.fieldnames;rows=list(parser)
    aliases=h.config("external_schema_protocol.json")["aliases"].get(prefix,{})
    renamed=[aliases.get(k,k) for k in header]
    if not rows or len(renamed)!=len(set(renamed)) or set(renamed)!=set(record["schema"]):raise ValueError("Undeclared external schema")
    for row in rows:
        if None in row or any(v is None for v in row.values()):raise ValueError("Malformed external CSV row")
    rows=[{aliases.get(k,k):v for k,v in row.items()} for row in rows]
    for row in rows:
        fields=["StemTag","year","status"] if dataset=="HF453" else [k for k in ["seedlingID","yearOfOb","status","flag","sampled","alive","date","taxonCode","coordinate"] if k in row]
        for key in fields:row[key].encode("ascii")
    return rows,dict(record,actual_header=header,explicit_aliases=aliases)
def inspect():
    from collections import Counter
    outputs=[]
    for dataset,prefix in [("HF453","hf453-01"),("HF355","hf355-01"),("HF355","hf355-02"),("HF355","hf355-05")]:
        rows,record=load(dataset,prefix)
        value=dict(dataset=dataset,name=record["name"],source_sha256=record["sha256"],rows=len(rows),schema=list(rows[0]),actual_header=record["actual_header"],explicit_aliases=record["explicit_aliases"],
            field_support={k:dict(Counter(r[k] for r in rows)) for k in ["year","status","flag","sampled","alive","trueRecrt"] if k in rows[0]},
            first_outcome_inspection_after_protocol_commit="35c4c384ad",parser_amendments=["Lossless byte mapping","Exact metadata header aliases"],model_fitting=False,biological_scores_computed=False)
        outputs.append(value);print(json.dumps(value,sort_keys=True),flush=True)
    h.write_json("external_sources/record_support_inventory.json",outputs)
if __name__=="__main__":inspect()
