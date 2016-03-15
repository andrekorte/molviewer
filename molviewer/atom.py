# -*- coding: utf-8 -*-
'''
This module contains the Atom class

'''
import numpy as np
from .moldata import radius, name, sym2no, symbol, mass, floatcolor


class Atom(object):
    '''Class that provides an Atom object

    The atomic charge can be given as a float. The mass is computed
    based on the nuclear charge. The nuclear charge cann be given as the
    elements symbol (see example below).

    Example:
        >>> h = Atom(1, [1.0, 1.0, 1.0])
        >>> he = Atom('He', [0.0, 0.0, 0.0])
        >>> print(he)
        He   4.0026    0.000000    0.000000    0.000000
        >>> he.distance(h)
        1.7320508075688772

    :param r: Array of length 3, position of atom
    :type r: numpy.ndarray
    :param Z: Nuclear charge, defines element (can be given as element symbol)
    :type Z: int
    :param charge: Atomic charge (default: 0.0, optional)
    :type charge: float

    :raises: KeyError, TypeError
    '''

    def __init__(self, Z, r, **kwargs):
        assert len(r) == 3
        if isinstance(Z, int):
            self.Z = Z
        elif sym2no[Z]:
            # raises KeyError if Z not in sym2no dict
            self.Z = sym2no[Z]
        self.symbol = symbol[self.Z]
        self.name = name[self.Z]
        self.r = np.asarray(r)
        self.mass = mass[self.Z]
        self.charge = float(kwargs.get('charge', 0.0))

    def __repr__(self):
        string = '{:4s}{:7.4f}{:12.6f}{:12.6f}{:12.6f}'.format(
            self.symbol, self.mass, self.r[0], self.r[1], self.r[2])
        return string

    def color(self):
        '''Returns the atom color as RGB values.

        The color is used in the molecular viewer.

        :returns: tuple -- A tuple of RGB values.

        '''
        return floatcolor[self.Z]

    def distance(self, atom):
        '''Returns the euclidian distance between two atoms.

        :param other: The second atom.
        :type other: atom.Atom
        :returns: float -- The euclidian distance in angstrom.

        '''
        assert isinstance(atom, Atom)
        return np.linalg.norm(self.r - atom.r)

    def radius(self):
        '''Returns the atomic (vdW) radius in angstrom.'''
        return radius[self.Z]

    def xyz(self):
        '''Returns a string for atom to be written in xyz file format.'''
        string = '{:4s}{:12.6f}{:12.6f}{:12.6f}'.format(
            self.symbol, self.r[0], self.r[1], self.r[2])
        return string

    def translate(self, delta_r):
        '''Translate the atom by vector delta_r.

        Values of delta_r are given in angstroms.

        :param delta_r: The x,y,z values.
        :type delta_r: numpy.ndarray
        :raises: AssertionError, ValueError

        '''
        assert isinstance(delta_r, np.ndarray)
        self.r += delta_r

if __name__ == '__main__':
    hydrogen = Atom(1, [1.0, 1.0, 1.0])
    helium = Atom('He', [0.0, 0.0, 0.0])
    print('Distance between atoms is {0:f}'.format(helium.distance(hydrogen)))
    print()
    print(helium)
