import ROOT, json, pathlib, hashlib, os
out=pathlib.Path('/home/gwatts/atlas/correctionlib_tests/met/results/discovery-20261104-01')
ROOT.xAOD.Init().ignore()
tool=ROOT.met.METSystematicsTool('METDiscovery')
assert tool.initialize().isSuccess()
result={'release':'25.2.117','recommended':[s.name() for s in tool.recommendedSystematics()], 'affecting':[s.name() for s in tool.affectingSystematics()]}
f=ROOT.TFile.Open(pathlib.Path('/home/gwatts/atlas/correctionlib_tests/met/config/input.txt').read_text().strip())
t=f.Get('CollectionTree')
result['input']={'uuid':str(f.GetUUID().AsString()), 'events':int(t.GetEntries()), 'bytes':int(f.GetSize())}
result['met_branches']=[b.GetName() for b in t.GetListOfBranches() if 'MET' in b.GetName()]
(out/'discovery.json').write_text(json.dumps(result,indent=2))
print(json.dumps(result,indent=2))
