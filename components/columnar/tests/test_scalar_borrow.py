"""Test native scalar borrowing against the generalized linked-column path."""
from unittest import TestCase, mock

import awkward as ak
import numpy as np

from atlas_met_prepare import Data, Resolved, original_scalar_view
from test_empty_projection import branch, dataset, link

CLUSTER = "egammaClusters"
GSF = "GSFTrackParticles"
MU = "CombinedMuonTrackParticles"
TAU = "TauTracks"
PDG = "TruthTaus"


def scalar(counts, values, dtype="f4"):
    offsets = np.array([0,*np.cumsum(counts)],dtype="i8")
    content = np.array(values,dtype=dtype)
    return ak.Array(ak.contents.ListOffsetArray(ak.index.Index64(offsets),ak.contents.NumpyArray(content)))


def same_columns(test, resolved, field, *args, **kwargs):
    fast = resolved.column(field,*args,**kwargs)
    with mock.patch.object(Data,"scalar_view",return_value=None):
        ordinary = resolved.column(field,*args,**kwargs)
    for observed,expected in zip(fast,ordinary,strict=True):
        test.assertEqual((observed.shape,observed.dtype),(expected.shape,expected.dtype))
        test.assertEqual(observed.tobytes(),expected.tobytes())
    return fast


class ScalarBorrow(TestCase):
    def test_routes_used_by_electron_photon_muon_and_tau(self):
        counts = [2,1]
        raw = {
            branch(CLUSTER,"calEta"):scalar(counts,[1.,2.,3.]),
            branch(CLUSTER,"ETACALOFRAME"):scalar(counts,[np.nan,-0.,0.75]),
            branch(GSF,"qOverP"):scalar(counts,[-1.,3.,2.]),
            branch(GSF,"phi"):scalar(counts,[0.125,-0.25,0.5]),
            branch(MU,"qOverP"):scalar(counts,[0.125,-0.25,0.5]),
            branch(MU,"numberOfPixelHits"):scalar(counts,[6,7,9],"u1"),
            branch(TAU,"flagSet"):scalar(counts,[0,32,0xffffffff],"u4"),
            branch(PDG,"pdgId"):scalar(counts,[15,-15,22],"i4"),
            branch(PDG,"IsHadronicTau"):scalar(counts,[False,True,False],"?"),
        }
        data = dataset(raw)
        data.keys = {target:i+101 for i,target in enumerate((CLUSTER,GSF,MU,TAU,PDG))}
        events = np.array([0,0,1,1],"i8")
        roles = (
            (CLUSTER,"ETACALOFRAME","f4",0,[0,1,0,3]),  # electron first cluster
            (CLUSTER,"ETACALOFRAME","f4",0,[1,0,0,3]),  # photon first cluster
            (GSF,"phi","f4",0,[0,1,0,3]),             # electron GSF track
            (MU,"numberOfPixelHits","i8",-1,[0,1,0,3]),
            (TAU,"flagSet","i8",0,[0,1,0,3]),
            (PDG,"IsHadronicTau","i8",-1,[0,1,0,3]),
        )
        for target,field,dtype,default,indices in roles:
            resolved = Resolved(data,link([data.keys[target]]*4,indices),events,(target,))
            with self.subTest(target=target,field=field):
                _,present = same_columns(self,resolved,field,dtype,default)
                np.testing.assert_array_equal(present,[True,True,True,False])
                name = branch(target,field)
                owned = ak.to_layout(raw[name]).content.data
                self.assertIs(data.scalar_views[name],owned)
        self.assertEqual(len(data.scalar_views),5)  # shared clusters use only one borrow

    def test_options_nested_empty_shapes_and_bad_rows_keep_generic_behavior(self):
        raw = {
            branch(GSF,"qOverP"):scalar([2,1],[2.,3.,4.]),
            branch(GSF,"option"):ak.Array([[None,7.],[9.]]),
            branch(GSF,"nested"):ak.Array([[[3.,4.],[]],[[8.,None]]]),
            branch(GSF,"unaligned"):scalar([1,2],[2.,3.,4.]),
            branch(GSF,"short"):scalar([1,1],[5.,8.]),
        }
        data = dataset(raw)
        resolved = Resolved(data,link([12,12,12],[0,1,0]),np.array([0,0,1]),(GSF,))
        _,present = same_columns(self,resolved,"option","f4",-1.)
        np.testing.assert_array_equal(present,[False,True,True])
        same_columns(self,resolved,"nested","f4",-1.,slot=1)
        same_columns(self,resolved,"nested","i8",-1,lengths=True)
        same_columns(self,resolved,"unaligned","f4",0.)
        self.assertIsNone(data.scalar_views[branch(GSF,"option")])
        self.assertIsNone(data.scalar_views[branch(GSF,"unaligned")])
        self.assertNotIn(branch(GSF,"nested"),data.scalar_views)
        with self.assertRaises(IndexError):
            resolved.column("short")

        empty = Resolved(data,link([12],[0xffffffff]),np.array([0]),(GSF,))
        _,present = same_columns(self,empty,"nested","f4",-1.,slot=1)
        self.assertFalse(present.any())
        self.assertIn((branch(GSF,"nested"),"f4",float,-1.,1,False),data.validated_empty_columns)

        # An invalid resolved row must take the original branch and raise.
        force = Resolved(data,link([12],[0]),np.array([0]),(GSF,))
        force.locations = [(GSF,np.array([0],"i8"),np.array([123],"i8"))]
        with self.assertRaises(IndexError):
            force.column("qOverP")
        self.assertIsNotNone(data.scalar_views[branch(GSF,"qOverP")])

    def test_slices_tagged_and_noncontiguous_arrays_do_not_borrow(self):
        leaf = scalar([1,2,1],[1.,2.,3.,4.])
        self.assertIsNone(original_scalar_view(leaf[1:],np.array([0,2,3],"i8")))
        padded_empty = scalar([0,1],[3.])
        self.assertIsNone(original_scalar_view(padded_empty[1:],np.array([0,1],"i8")))
        self.assertIsNone(original_scalar_view(ak.Array([[],[]]),np.array([0,0,0],"i8")))
        tagged = ak.with_parameter(leaf,"test","unrecognized vector")
        self.assertIsNone(original_scalar_view(tagged,np.array([0,1,3,4],"i8")))
        source = np.arange(8.,dtype="f4")[::2]
        stride = ak.Array(ak.contents.ListOffsetArray(ak.index.Index64(np.array([0,1,3,4],"i8")),
                                                       ak.contents.NumpyArray(source)))
        self.assertIsNone(original_scalar_view(stride,np.array([0,1,3,4],"i8")))
