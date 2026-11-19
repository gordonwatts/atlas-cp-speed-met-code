import sys, json, pathlib, hashlib
sys.path.insert(0,"/home/gwatts/atlas-cp-speed/muon-20261101/muon/install/python313")
sys.path.append("/home/gwatts/atlas-cp-speed/muon-20261101/venv/lib/python3.13/site-packages")
import correctionlib
from correctionlib import _core
p=pathlib.Path(sys.argv[1]); out=pathlib.Path(sys.argv[2]); out.mkdir(parents=True,exist_ok=False)
c=correctionlib.CorrectionSet.from_file(str(p/'tables.json')); probes=json.loads((p/'probes.json').read_text()); bad=[]
for i,t in enumerate(probes):
    actual=c[t['correction']].evaluate(*t['inputs'])
    if actual!=t['expected']: bad.append(dict(t,index=i,actual=actual))
core=pathlib.Path(_core.__file__)
r={"checks":len(probes),"mismatches":len(bad),"failures":bad,"binding":str(core),"binding_sha256":hashlib.sha256(core.read_bytes()).hexdigest(),"payload_sha256":hashlib.sha256((p/'tables.json').read_bytes()).hexdigest(),"scope":"ROOT table contents and axis boundary/flow probes only; not MET recipe acceptance"}
(out/'comparison.json').write_text(json.dumps(r,indent=2));print(json.dumps(r,indent=2));sys.exit(bool(bad))
