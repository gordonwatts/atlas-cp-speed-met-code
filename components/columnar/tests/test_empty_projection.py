"""Exercise empty linked-column validation and fallback without a physics input."""
from unittest import TestCase, mock

import awkward as ak
import numpy as np
from numpy.exceptions import AxisError

from atlas_met_columnar import OWNERS
from atlas_met_prepare import Data, Link, Resolved

ID = "InDetTrackParticles"
GSF = "GSFTrackParticles"
OWNERS_BY_KEY = {ID: 11, GSF: 12}


def branch(target, field):
    return target + "AuxDyn." + field


def link(keys, indices):
    return Link(ak.zip({"m_persKey": np.asarray(keys, dtype="u4"),
                        "m_persIndex": np.asarray(indices, dtype="u4")}))


def dataset(raw, events=2):
    output = {name+".o": np.zeros(events+1, dtype="i8") for name in OWNERS}
    return Data(raw,output,OWNERS_BY_KEY,None,False)


class EmptyProjection(TestCase):
    def test_empty_target_repeat_retains_default_and_false_presence(self):
        scalar = ak.unflatten(np.empty(0,"f4"),[0,0])
        vector = ak.unflatten(ak.unflatten(np.empty(0,"f4"),np.empty(0,"i8")),[0,0])
        data = dataset({branch(ID,"qOverP"):scalar, branch(ID,"cov"):vector})
        event_rows = np.array([0,1])
        elements = link([11,11],[0,0xffffffff])
        first = Resolved(data,elements,event_rows,(ID,))
        again = Resolved(data,elements,event_rows,(ID,))
        self.assertFalse(first.valid.any())
        with mock.patch("atlas_met_prepare.ak.flatten", wraps=ak.flatten) as flatten:
            for resolved in (first,again):
                got, present = resolved.column("cov","f4",-1.,slot=1)
                np.testing.assert_array_equal(got,np.array([-1.,-1.],"f4"))
                np.testing.assert_array_equal(present,[False,False])
            self.assertEqual(flatten.call_count,1)

    def test_unreachable_first_target_keeps_alternate_defaults(self):
        data = dataset({
            branch(ID,"qOverP"):ak.Array([[3.],[]]),
            branch(ID,"cov"):ak.Array([[[2.,7.]],[]]),
            branch(GSF,"qOverP"):ak.Array([[],[8.,9.]]),
            branch(GSF,"cov"):ak.Array([[],[[12.,None],[18.,19.]]]),
        })
        events = np.array([1,1,1,0])
        elements = link([11,12,12,12],[0,0,1,0])
        resolved = Resolved(data,elements,events,(ID,GSF))
        np.testing.assert_array_equal(resolved.valid,[False,True,True,False])
        for _ in range(2):
            got, present = resolved.column("cov","f4",-1.,slot=1)
            np.testing.assert_array_equal(got,[-1.,-1.,19.,-1.])
            np.testing.assert_array_equal(present,[False,False,True,False])
            sizes, present = resolved.column("cov","i8",-1,lengths=True)
            np.testing.assert_array_equal(sizes,[-1,2,2,-1])
            np.testing.assert_array_equal(present,[False,True,True,False])
        missing, present = resolved.column("absent","i8",-7)
        np.testing.assert_array_equal(missing,[-7,-7,-7,-7])
        self.assertFalse(present.any())

    def test_empty_projection_still_rejects_malformed_axis_and_scatter(self):
        data = dataset({
            branch(ID,"qOverP"):ak.Array([[3.],[]]),
            branch(ID,"theta"):ak.Array([[0.5],[]]),
            branch(ID,"cov"):ak.unflatten(ak.Array(np.array([[3.,4.]])),[1,0]),
        })
        resolved = Resolved(data,link([11],[0]),np.array([1]),(ID,))
        ordinary, present = resolved.column("theta","f4",-1.)
        np.testing.assert_array_equal(ordinary,[-1.])
        self.assertFalse(present.any())
        # Successful scalar or slot validation must not certify a different view.
        resolved.column("cov","f4",-1.,slot=0)
        for _ in range(2):
            with self.assertRaises(AxisError):
                resolved.column("theta","f4",-1.,slot=0)
            with self.assertRaises(ValueError):
                resolved.column("cov","f4",-1.)
