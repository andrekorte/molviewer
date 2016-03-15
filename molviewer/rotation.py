# -*- coding: utf-8 -*-
'''
This module contains the RotationConstant class.

'''
import numpy as np
from moldata import amu
from molutils import isclose


class RotationConstant(object):
    '''Class to calculate the rotation constants of a molecule.

    :param molecule: The molecule for which to calculate the constants.
    :type molecule: molecule.Molecule

    '''
    def __init__(self):
        pass

    def inertia_tensor(self, mol):
        '''Returns the inertia tensor.

        Input is in amu and angstrom; output in kg*m**2.

        :param molecule: The molecule for which to calculate the constants.
        :type molecule: molecule.Molecule

        :returns: numpy.ndarray -- The inertia tensor. A 3x3 matrix.

        '''
        # calculate center of mass
        cm = mol.center_of_mass       # angstrom
        # shift molecule to new reference frame
        mol.translate(-cm)
        tensor = np.zeros([3, 3])
        # Diagonal elements
        I_xx = np.sum([atom.mass * amu * ((1e-10 * atom.r[1])**2 + (1e-10 * atom.r[2])**2) for atom in mol.atoms])
        I_yy = np.sum([atom.mass * amu * ((1e-10 * atom.r[0])**2 + (1e-10 * atom.r[2])**2) for atom in mol.atoms])
        I_zz = np.sum([atom.mass * amu * ((1e-10 * atom.r[0])**2 + (1e-10 * atom.r[1])**2) for atom in mol.atoms])

        tensor[0, 0] = I_xx
        tensor[1, 1] = I_yy
        tensor[2, 2] = I_zz

        tensor[0, 1] = np.sum([atom.mass * amu * 1e-10 * atom.r[0] * 1e-10 * atom.r[1] for atom in mol.atoms])
        tensor[0, 2] = np.sum([atom.mass * amu * 1e-10 * atom.r[0] * 1e-10 * atom.r[2] for atom in mol.atoms])
        tensor[1, 2] = np.sum([atom.mass * amu * 1e-10 * atom.r[1] * 1e-10 * atom.r[2] for atom in mol.atoms])
        tensor[1, 0] = tensor[0, 1]
        tensor[2, 0] = tensor[0, 2]
        tensor[2, 1] = tensor[1, 2]

        return tensor

    def inertia_tensor_amu(self, mol):
        '''Returns the inertia tensor.

        Input is in amu and angstrom; output in amu and angstrom.

        :param molecule: The molecule for which to calculate the constants.
        :type molecule: molecule.Molecule

        :returns: numpy.ndarray -- The inertia tensor. A 3x3 matrix.

        '''
        # calculate center of mass
        cm = mol.center_of_mass  # angstrom
        # shift molecule to new reference frame
        mol.translate(-cm)
        tensor = np.zeros([3, 3])
        # Diagonal elements
        I_xx = np.sum([atom.mass * (atom.r[1]**2 + atom.r[2]**2) for atom in mol.atoms])
        I_yy = np.sum([atom.mass * (atom.r[0]**2 + atom.r[2]**2) for atom in mol.atoms])
        I_zz = np.sum([atom.mass * (atom.r[0]**2 + atom.r[1]**2) for atom in mol.atoms])

        tensor[0, 0] = I_xx
        tensor[1, 1] = I_yy
        tensor[2, 2] = I_zz

        tensor[0, 1] = np.sum([atom.mass * atom.r[0] * atom.r[1] for atom in mol.atoms])
        tensor[0, 2] = np.sum([atom.mass * atom.r[0] * atom.r[2] for atom in mol.atoms])
        tensor[1, 2] = np.sum([atom.mass * atom.r[1] * atom.r[2] for atom in mol.atoms])
        tensor[1, 0] = tensor[0, 1]
        tensor[2, 0] = tensor[0, 2]
        tensor[2, 1] = tensor[1, 2]

        return tensor

    def symmetry(self, mol):
        '''Returns the symmetry of the molecule based on inertia tensor.

        :param molecule: The molecule for which to calculate the constants.
        :type molecule: molecule.Molecule

        :returns: str -- The rotational symmetry.

        '''
        tensor = self.inertia_tensor_amu(mol)
        if tensor.shape == (3, 3):
            # eigenvalues are principal moments of inertia
            vals, vecs = np.linalg.eig(tensor)
        elif tensor.shape == (3,):
            vals = tensor

        if isclose(vals[0], 0) and isclose(vals[1], vals[2]):
            return 'linear'
        elif isclose(vals[0], vals[1]) and isclose(vals[0], vals[2]):
            return 'spherical rotor'
        elif isclose(vals[0], vals[1]) and not isclose(vals[0], vals[2]):
            return 'symmetric rotor'
        else:
            return 'asymmetric rotor'

    def rotational_constants(self, mol, unit='wavenumber'):
        '''Returns the rotational constants given the inertia tensor.

        :param molecule: The molecule for which to calculate the constants.
        :type molecule: molecule.Molecule
        :param unit: The unit in which the constants are returned. Valid
                    choices are 'wavenumber' (default) and 'GHz'.
        :type str:

        :returns: numpy.array -- The rotational constants. A 1x3 matrix.

        '''
        tensor = self.inertia_tensor(mol)
        if tensor.shape == (3, 3):
            # eigenvalues are principal moments of inertia
            vals, vecs = np.linalg.eig(tensor)
        elif tensor.shape == (3,):
            vals = tensor

        # We want A > B > C
        vals = sorted(vals)

        # rotational constant = hbar / (4*pi*c*I)
        hbar = 1.054571800e-34  # Js
        c = 299792458  # m/s

        A = hbar / (4 * np.pi * c * vals[0])  # 1/m
        B = hbar / (4 * np.pi * c * vals[1])  # 1/m
        C = hbar / (4 * np.pi * c * vals[2])  # 1/m

        A = A / 100.  # 1/cm
        B = B / 100.  # 1/cm
        C = C / 100.  # 1/cm

        if unit == 'GHz':
            from moldata import wavenumber2ghz
            A = A * wavenumber2ghz
            B = B * wavenumber2ghz
            C = C * wavenumber2ghz
            return np.asarray([A, B, C])

        return np.asarray([A, B, C])


if __name__ == '__main__':
    from molio import read_xyz

    mol = read_xyz('molecules/xyz/methane.xyz')
    rotconst = RotationConstant()
    print(rotconst.symmetry(mol))
