"""Narrow canonical external metadata/data retrieval; no model selection."""
import datetime,json,subprocess,xml.etree.ElementTree as ET
from pathlib import Path
import h05_common as h

def fetch_metadata():
    c=h.config("external_protocol.json");records=[]
    for identity,item in c["sources"].items():
        payload=subprocess.check_output(["curl","--fail","--silent","--show-error","--location","--max-time","90",item["metadata"]])
        path=h.ROOT/"external_sources"/identity/"metadata.xml";path.parent.mkdir(parents=True,exist_ok=True)
        if path.exists() and path.read_bytes()!=payload:raise ValueError("Metadata changed; retain original before versioning")
        path.write_bytes(payload)
        root=ET.fromstring(payload)
        if root.get("packageId")!="knb-lter-hfr."+identity[2:]+"."+str(item["version"]):raise ValueError("Dataset version differs from declared source")
        def local(element):return element.tag.split("}")[-1]
        datasets=[e for e in root.iter() if local(e)=="dataTable"]
        for table in datasets:
            name=next((e.text for e in table if local(e)=="entityName"),"")
            urls=[e.text for e in table.iter() if local(e)=="url"]
            fields=[e.text for e in table.iter() if local(e)=="attributeName"]
            if not any(name.startswith(prefix) for prefix in item["tables"]):continue
            if len(urls)!=1 or not urls[0].startswith("https://"):raise ValueError("Unresolved canonical public data URL "+name)
            records.append(dict(dataset=identity,version=item["version"],DOI=item["DOI"],license=item["license"],name=name,url=urls[0],
                schema=fields,metadata_url=item["metadata"],metadata_sha256=h.sha(path),access_date="2026-10-08",
                raw_observations_opened=False))
    h.write_json("external_sources/metadata_inventory.json",dict(sources=records,observation_outcomes_opened=False))
    print(json.dumps([dict(dataset=r["dataset"],name=r["name"],url=r["url"],fields=r["schema"]) for r in records],indent=2),flush=True)

def fetch_data():
    h.guard("external_protocol")
    inventory=json.loads((h.ROOT/"external_sources/metadata_inventory.json").read_text())
    records=[]
    for row in inventory["sources"]:
        path=h.ROOT/".inputs"/row["dataset"]/row["name"];path.parent.mkdir(parents=True,exist_ok=True)
        if path.exists():raise ValueError("External raw input already present; do not overwrite")
        payload=subprocess.check_output(["curl","--fail","--silent","--show-error","--location","--max-time","120",row["url"]]);path.write_bytes(payload)
        records.append(dict(row,sha256=h.sha(path),bytes=len(payload),local_input=str(path.relative_to(h.ROOT)),raw_observations_opened=False))
    h.write_json("external_sources/source_inventory.json",dict(sources=records,protocol_sha256=h.sha(h.ROOT/"config/external_protocol.json"),outcome_selection=False))
    print("Exact external inputs retrieved",sum(r["bytes"] for r in records),"bytes; no observations parsed",flush=True)

if __name__=="__main__":
    import sys
    fetch_metadata() if sys.argv[1]=="metadata" else fetch_data()
