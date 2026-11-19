import ROOT, json, pathlib, hashlib, math, itertools, sys
ROOT.gInterpreter.Declare('#define XAOD_STANDALONE 1\n#include "PathResolver/PathResolver.h"')
out=pathlib.Path(sys.argv[1]); out.mkdir(parents=True,exist_ok=False)
config=pathlib.Path(str(ROOT.PathResolverFindCalibFile("METUtilities/R22_PreRecs/TrackSoftTerms-pflow_Dec24.config")))
env=ROOT.TEnv(str(config)); source=pathlib.Path(str(ROOT.PathResolverFindCalibFile(env.GetValue("Conf.InputFile",""))))
f=ROOT.TFile.Open(str(source)); suffix=env.GetValue("PtHardSyst.suffix",""); directory=env.GetValue("PtHardSyst.sourcedir","")
corrections=[]; inventory=[]; probes=[]
for name in ("shiftpara", "resopara", "resoperp"):
    key=directory+"/"+name+"_"+suffix
    h=f.Get(key)
    if not h or h.ClassName()!="TH3D": raise RuntimeError(key)
    axes=[h.GetXaxis(),h.GetYaxis(),h.GetZaxis()]
    edges=[[float(a.GetBinLowEdge(i)) for i in range(1,a.GetNbins()+2)] for a in axes]
    # Include ROOT flow bins. Only pt-hard overflow is clamped by the original tool.
    content=[float(h.GetBinContent(min(i,h.GetNbinsX()),j,k)) for i in range(h.GetNbinsX()+2) for j in range(h.GetNbinsY()+2) for k in range(h.GetNbinsZ()+2)]
    corrections.append({"name":name,"version":1,"inputs":[{"name":n,"type":"real"} for n in ("pthard","njet","mu")],"output":{"name":"value","type":"real"},"data":{"nodetype":"multibinning","inputs":["pthard","njet","mu"],"edges":[["-inf"]+e+["inf"] for e in edges],"content":content,"flow":"error"}})
    inventory.append({"correction":name,"object":key,"class":h.ClassName(),"axis_bins":[a.GetNbins() for a in axes],"edges":edges,"errors_exported":False,"errors_consumed":False,"flow":"ROOT under/overflow retained; pt-hard overflow clamps to last ordinary bin"})
    defaults=[float(a.GetBinCenter(1)) for a in axes]
    for dim,a in enumerate(axes):
        points=[float(a.GetBinCenter(i)) for i in range(0,a.GetNbins()+2)]
        points += [p for edge in edges[dim] for p in (math.nextafter(edge,-math.inf),edge,math.nextafter(edge,math.inf))]
        for p in points:
            v=defaults.copy(); v[dim]=p
            bins=[axes[d].FindFixBin(v[d]) for d in range(3)]; bins[0]=min(bins[0],h.GetNbinsX())
            probes.append({"correction":name,"inputs":v,"expected":float(h.GetBinContent(*bins)),"bins":bins,"dimension":dim})
    for bins in itertools.product(*[range(a.GetNbins()+2) for a in axes]):
        v=[float(a.GetBinCenter(i)) for a,i in zip(axes,bins)]
        b=list(bins); b[0]=min(b[0],h.GetNbinsX())
        probes.append({"correction":name,"inputs":v,"expected":float(h.GetBinContent(*b)),"bins":b,"dimension":"all-bin-centers"})
payload={"schema_version":2,"corrections":corrections}
(out/'tables.json').write_text(json.dumps(payload,allow_nan=False,separators=(',',':')))
(out/'probes.json').write_text(json.dumps(probes,allow_nan=False,separators=(',',':')))
provenance={"config":str(config),"config_sha256":hashlib.sha256(config.read_bytes()).hexdigest(),"source":str(source),"source_sha256":hashlib.sha256(source.read_bytes()).hexdigest(),"inventory":inventory,"probe_count":len(probes),"status":"candidate table conversion only; no complete MET recipe"}
(out/'provenance.json').write_text(json.dumps(provenance,indent=2))
print(json.dumps({"probe_count":len(probes),"objects":inventory,"source":str(source)},indent=2))
