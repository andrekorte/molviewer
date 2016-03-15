import unittest

from atom import Atom
from molecule import Molecule
import numpy as np


class AtomTest(unittest.TestCase):
    def test_missingarguments(self):
        with self.assertRaises(TypeError):
            x = Atom()

    def test_Z_as_symbol(self):
        li = Atom('Li', [0, 0, 0])
        self.assertEqual(li.symbol, 'Li')

    def test_wrong_Z(self):
        with self.assertRaises(KeyError):
            x = Atom('Z', [0, 0, 0])
        with self.assertRaises(KeyError):
            ne = Atom('NE', [0, 0, 0])

    def test_wrong_r_length(self):
        with self.assertRaises(AssertionError):
            x = Atom('Z', [0, 0])
        with self.assertRaises(AssertionError):
            x = Atom('Z', [0, 0, 0, 0])

    def test_isinstance(self):
        li = Atom(3, [0.0, 0.0, 1.0])
        self.assertIsInstance(li, Atom)

    def test_symbol(self):
        li = Atom(3, [0.0, 0.0, 1.0])
        self.assertEqual(li.symbol, 'Li')

    def test_charge(self):
        li = Atom(3, [0.0, 0.0, 1.0])
        self.assertEqual(li.charge, 0.0)

    def test_Z(self):
        li = Atom(3, [0.0, 0.0, 1.0])
        self.assertEqual(li.Z, 3)

    def test_mass(self):
        li = Atom(3, [0.0, 0.0, 1.0])
        self.assertEqual(li.mass, 6.941)

    def test_name(self):
        li = Atom(3, [0.0, 0.0, 1.0])
        self.assertEqual(li.name, 'lithium')

    def test_color(self):
        c = Atom(6, [1.0, 1.0, 1.0])
        self.assertTupleEqual(c.color(), (0.4392156862745098, 0.5019607843137255, 0.5647058823529412))

    def test_distance(self):
        h = Atom(1, [1.0, 1.0, 1.0])
        he = Atom(2, [0.0, 0.0, 0.0])
        self.assertEqual(he.distance(h), 1.7320508075688772)

    def test_distance_wrong_position_vector_type(self):
        h = Atom(1, [1.0, 1.0, 1.0])
        with self.assertRaises(AssertionError):
            h.distance([3, 0, 1])

    def test_radius(self):
        he = Atom(2, [0.0, 0.0, 0.0])
        self.assertEqual(he.radius(), 1.4)

    def test_xyz(self):
        he = Atom(2, [0.0, 0.0, 0.0])
        self.assertEqual(he.xyz(), 'He      0.000000    0.000000    0.000000')

    def test_translate(self):
        n = Atom(4, [0.0, 0.0, 0.0])
        n.translate(np.asarray([1.0, 1.0, 1.0]))
        self.assertTrue(np.allclose(n.r, [1.0, 1.0, 1.0], rtol=1e-05, atol=1e-08))

    def test_translate_wrong_vector_type(self):
        n = Atom(4, [0.0, 0.0, 0.0])
        with self.assertRaises(AssertionError):
            n.translate([0, 0, 1])
        with self.assertRaises(AssertionError):
            n.translate((0, 0, 1))

    def test_translate_wrong_vector_length(self):
        n = Atom(4, [0.0, 0.0, 0.0])
        too_short = np.asarray([0, 2])
        too_long = np.asarray([0, 1, 2, 3])
        with self.assertRaises(ValueError):
            n.translate(too_short)
        with self.assertRaises(ValueError):
            n.translate(too_long)

if __name__ == '__main__':
    unittest.main()
