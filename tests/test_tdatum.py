"""Regression tests for mathematical validation and mutation-loop operations."""

import unittest

from sage.all import QQ, ZZ, LaurentPolynomialRing, Permutation
from sage.all import diagonal_matrix, identity_matrix, matrix, zero_matrix

from tdatum import MutationLoop, TDatum
from tdatum.examples import LengthOne, RSG, Rank2, SG, UntwistedAffine


class TDatumTests(unittest.TestCase):
    def setUp(self):
        self.ring = LaurentPolynomialRing(QQ, "z")
        self.z = self.ring.gen()
        self.ap = diagonal_matrix(self.ring, [1+self.z**2]*2)
        self.am = matrix(self.ring, 2, 2, [1+self.z**2, -self.z, -2*self.z, 1+self.z**2])
        self.d = diagonal_matrix([1, 2])

    def test_reject_invalid_data(self):
        z, ring = self.z, self.ring
        scalar = matrix(ring, 1, 1, [1+z**2])
        cases = [
            (scalar, scalar, zero_matrix(ZZ, 1)),
            (scalar, scalar, matrix(ZZ, 1, 1, [-1])),
            (self.ap, self.ap, matrix(ZZ, 2, 2, [1, 1, 0, 1])),
            (scalar, matrix(ring, 1, 1, [1+z**2-z/2]), 'identity'),
            (scalar, 2*scalar, 'identity'),
            (matrix(ring, 1, 1, [2]), matrix(ring, 1, 1, [2]), 'identity'),
            (matrix(ring, 1, 1, [1+1/z]), matrix(ring, 1, 1, [1+1/z]), 'identity'),
            (self.ap, self.am, 'identity'),
            (matrix(ring, 1, 1, [1+z**2-z]), matrix(ring, 1, 1, [1+z**2-z]), 'identity'),
            (matrix(ring, 1, 1, [1+z**2-z**3]), scalar, 'identity'),
            (zero_matrix(ring, 0), zero_matrix(ring, 0), 'identity'),
        ]
        for ap, am, d in cases:
            with self.subTest(ap=ap, am=am, d=d), self.assertRaises(ValueError):
                TDatum(ap, am, d)

    def test_symmetrizer_compatibility(self):
        z = self.z
        swap = matrix(self.ring, 2, 2, [1, z, z, 1])
        with self.assertRaisesRegex(ValueError, "commute"):
            TDatum(swap, swap, self.d)
        # This violates dual integrality before the symplectic identity is tested.
        am = matrix(self.ring, 2, 2, [1+z**2, 0, -z, 1+z**2])
        with self.assertRaisesRegex(ValueError, "integer coefficients"):
            TDatum(self.ap, am, self.d)

    def test_sign_and_langlands_duals(self):
        td = TDatum(self.ap, self.am, self.d)
        dual = td.sign_dual()
        self.assertEqual(dual.pair(), (td.pair()[1], td.pair()[0]))
        self.assertEqual(dual.symmetrizer(), (1, 2))
        self.assertEqual(dual.sign_dual().pair(), td.pair())
        self.assertEqual(td.langlands_dual().langlands_dual().pair(), td.pair())

    def test_nonsymmetric_round_trip(self):
        td = TDatum(self.ap, self.am, self.d)
        loop = td.mutation_loop()
        self.assertTrue(loop.is_complete())
        self.assertEqual(loop.symmetrizer(), (1, 1, 2, 2))
        recovered = loop.t_datum()
        self.assertEqual(recovered.pair(), td.pair())
        self.assertEqual(recovered.symmetrizer(), td.symmetrizer())
        self.assertEqual(loop.inverse().symmetrizer(), loop.symmetrizer())
        self.assertEqual(loop.inverse().inverse(), loop)

    def test_variable_name(self):
        loop = TDatum(self.ap, self.am, self.d).mutation_loop()
        self.assertEqual(str(loop.t_datum("q").variable()), "q")

    def test_input_matrices_cannot_change_datum(self):
        td = TDatum(self.ap, self.am, self.d)
        self.ap[0, 0] = 0
        self.d[0, 0] = 5
        self.assertEqual(td.pair()[0][0, 0], 1+self.z**2)
        self.assertEqual(td.symmetrizer(), (1, 2))
        with self.assertRaises(ValueError):
            td.pair()[0][0, 0] = 0

    def test_rsg_cases(self):
        cases = [([2,1], (4,)), ([3,1], (2,6)), ([4,1], (2,2,8)),
                 ([3,2,2], (2,6,6,14,14)), ([4,2,3], (2,2,8,8,18,18,18)),
                 ([3,2,2,2], (2,6,6,14,14,34,34))]
        for params, degrees in cases:
            with self.subTest(params=params):
                td = RSG(params).t_datum()
                self.assertEqual(td.degrees(), degrees)
                ap, am = td.pair()
                z = td.variable()
                self.assertEqual(ap*am.transpose().subs({z: 1/z}), am*ap.transpose().subs({z: 1/z}))

    def test_rank_two_catalogue(self):
        for label in range(1, 7):
            with self.subTest(label=label):
                td = Rank2(label).t_datum()
                self.assertEqual(td.size(), 2)
                loop = td.mutation_loop()
                self.assertTrue(loop.is_complete())
                self.assertEqual(loop.t_datum().pair(), td.pair())

    def test_additional_constructors(self):
        self.assertEqual(LengthOne([1, -2, 1]).t_datum().degrees(), (4,))
        self.assertGreater(SG([3, 2]).t_datum().size(), 0)
        self.assertEqual(UntwistedAffine('A', 2, 2).t_datum().degrees(), (2, 2))


class MutationLoopTests(unittest.TestCase):
    def test_equality_preserves_time_steps(self):
        b, nu = zero_matrix(ZZ, 2), Permutation([1, 2])
        self.assertNotEqual(MutationLoop(b, [[0], [0]], nu), MutationLoop(b, [[1], [1]], nu))
        self.assertNotEqual(MutationLoop(b, [[0,1]], nu), MutationLoop(b, [[0], [1]], nu))
        self.assertEqual(MutationLoop(b, [[0,1]], nu), MutationLoop(b, [[1,0]], nu))
        self.assertNotEqual(MutationLoop(b, [[0]], nu), object())

    def test_invalid_loops(self):
        b = matrix(ZZ, 2, 2, [0, 1, -1, 0])
        nu = Permutation([1, 2])
        with self.assertRaisesRegex(ValueError, "commute"):
            MutationLoop(b, [[0,1], [1,0]], nu)
        with self.assertRaisesRegex(ValueError, "distinct vertices"):
            MutationLoop(b, [[0,0]], nu)
        with self.assertRaisesRegex(ValueError, "skew-symmetric"):
            MutationLoop(b, [[0], [0]], nu, [1,2])
        with self.assertRaisesRegex(ValueError, "complete"):
            MutationLoop(zero_matrix(ZZ, 2), [[0]], nu).t_datum()

    def test_loop_input_is_copied(self):
        b, blocks = zero_matrix(ZZ, 2), [[0], [1]]
        loop = MutationLoop(b, blocks, Permutation([1, 2]))
        blocks[0][0] = 1
        loop.vertices()[0][0] = 1
        loop.indices(0)[0] = 1
        self.assertEqual(loop.vertices(), [[0], [1]])


if __name__ == '__main__':
    unittest.main()
