# -*- coding: utf-8 -*-
'''
This module contains the Molecule class.

'''
import numpy as np
from bond import Bond
from moldata import symbol, bohr2angstrom
from atom import Atom


class Molecule(object):
    '''Class that provides a molecule object.

    :param atoms: The list of atoms which comprise the molecule.
    :type atoms: list
    :param name: The name of the molecule (optional).
    :type name: str
    :param charge: The molecular charge (optional).
    :type charge: int

    '''

    def __init__(self, atoms=[], **kwargs):
        self.filename = kwargs.get('filename', None)
        self.atoms = atoms
        self.name = kwargs.get('name', 'molecule')
        self.charge = int(kwargs.get('charge', 0))

    def __repr__(self):
        return self.name

    def __str__(self):
        lines = []
        first_line = 'Name: {}'.format(self.name)
        second_line = 'Mass: {}'.format(self.mass())
        lines.append(first_line)
        lines.append(second_line)
        for atom in self.atoms:
            lines.append(atom.__repr__())
        record = '\n'.join(lines)
        return record

    def bonds(self, scale=1.1):
        '''Calculates the bonds based on the vdW-radii.

        :param scale: The scaling factor (default: 1.1).
        :type scale: float

        :returns: list -- A list of Bond objects.

        '''
        bonds = []
        for i, ati in enumerate(self.atoms):
            for j in range(i):
                atj = self.atoms[j]
                if ati.distance(atj) < bohr2angstrom * scale * (ati.radius() + atj.radius()):
                    bonds.append(Bond(i, j, ati, atj))
        return bonds

    def add_atom(self, atom):
        '''Add an atom to the molecule.

        :param atom: The atom to add.
        :type atom: atom.Atom

        :returns: None

        '''
        try:
            assert isinstance(atom, Atom)
            self.atoms.append(atom)
        except AssertionError:
            print("Can't add non-atom type object to atomlist")

    def remove_atom(self, atom_index):
        '''Remove atom with atom_index from the molecule.

        :param atom_index: The index of the atom to remove.
        :type atom_index: int

        '''
        self.atoms.pop(atom_index)

    def mass(self):
        '''Returns the mass of the molecule.

        :returns: float -- The mass in amu.

        '''
        return sum([atom.mass for atom in self.atoms])

    def nel(self):
        '''Returns number of electrons of the molecule.

        :returns: int -- The number of electrons.

        '''
        return sum(atom.Z for atom in self.atoms) - self.charge

    def stoichiometry(self):
        '''Computes the stoichiometry of the molecule.

        :returns: str -- The stoichiometry of the molecule.

        '''
        from collections import Counter
        cnt = Counter()
        for atom in self.atoms:
            cnt[atom.Z] += 1
        keys = sorted(cnt.keys())
        s = []
        for key in keys:
            if cnt[key] == 1:
                s.append(symbol[key])
            else:
                s.append('{}{}'.format(symbol[key], cnt[key]))
        return ''.join(s)

    def translate(self, dr):
        '''Translate molecule by dr.

        Input is given in angstrom.

        :param dr: A 1x3 array which translates the molecule.
        :type dr: numpy.ndarray

        '''
        for atom in self.atoms:
            atom.translate(dr)

    def center_of_mass(self):
        '''Returns the center of mass.

        :returns: numpy.array -- The center of mass of the molecule. A 1x3 array.

        '''
        X = sum([atom.mass * atom.r[0] for atom in self.atoms]) / self.mass()
        Z = sum([atom.mass * atom.r[2] for atom in self.atoms]) / self.mass()
        Y = sum([atom.mass * atom.r[1] for atom in self.atoms]) / self.mass()
        return np.asarray([X, Y, Z])

    def xyz(self):
        '''Returns the molecule in xyz format.

        :returns: str -- The molecule in xyz format.

        '''
        lines = []
        first_line = '{}'.format(len(self.atoms))
        second_line = '{}'.format(self.name)
        lines.append(first_line)
        lines.append(second_line)
        for atom in self.atoms:
            lines.append(atom.xyz())
        record = '\n'.join(lines)
        return record

if __name__ == '__main__':
    h1 = Atom(1, [0.0, 0.0, 0.36])
    h2 = Atom(1, [0.0, 0.0, -0.36])
    h3 = Atom(1, [0.0, 0.0, 0.0])
    mol = Molecule(atoms=[h1, h2])
    mol.add_atom(123)
    mol.add_atom(h3)
    print(mol)
    print(mol.mass())
