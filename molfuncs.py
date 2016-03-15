# -*- coding: utf-8 -*-
'''
This module contains all functions needed for:
    * geometrical calculations and transformations
    * retrieval of molecular information

'''
import numpy as np


def moleculeinfo(mol):
    '''Returns basic information about the molecule

    This function is used in the viewer class to create the info tree view.

    Args:
        mol (molecule.molecule): The molecule instance

    Returns:
        data (dict): The information about the molecule
    '''
    bonddict = {}
    for i, b in enumerate(mol.bonds()):
        bonddict[i] = b.asdict()
    data = {
            'General': {
            'File': mol.filename,
            'Name': mol.name,
            'Mass': mol.mass(),
            'Electrons': mol.nel(),
            'Charge': mol.charge,
            'Formula': mol.stoichiometry()
            },
            'Bond table': [{'ishbond': b.ishbond, 'Atoms': b.atoms} for b in mol.bonds()]
           }
    return data

#
# geometrical functions
#
def distance_matrix(mol):
    '''Calculate the distance matrix for a molecule

    Returned values are in angstrom.

    Args:
        mol (molecule.Molecule): The molecule instance.

    Returns:
        mat (numpy.array): A N*N matrix with the atomic distances, where N is
        the number of atoms in the molecule.
    '''
    atoms = mol.atoms
    mat = np.zeros([len(atoms), len(atoms)])
    for i, ati in enumerate(atoms):
        for j, atj in enumerate(atoms):
            mat[i, j] = ati.distance(atj)
    return mat

def unit_vector(v):
    '''Returns the unit vector of vector v

    Args:
        vector (numpy.array): The vector

    Returns:
        unit_vector (numpy.array): The unit vector.
    '''
    return v / np.linalg.norm(v)


def angle_between_vectors(v1, v2):
    '''Returns the angle in radians between vector v1 and v2

    Args:
        vector1 (numpy.array): The first vector.
        vector2 (numpy.array): The second vector.

    Returns:
        angle (float): The angle in radians.
    '''
    v1u = unit_vector(v1)
    v2u = unit_vector(v2)
    angle = np.arccos(np.dot(v1u, v2u))

    if np.isnan(angle):
        if (v1u == v2u).all():
            return 0.0
        else:
            return np.pi
    return angle


def angle(i, j, k):
    '''Returns the angle in degrees between the points i, j, k

    Point i is the vertex.

    Args:
        i (numpy.array): The first point.
        j (numpy.array): The second point.
        k (numpy.array): The third point.

    Returns:
        angle (float): The angle in degrees.
    '''
    v1 = j - i
    v2 = k - i
    angle = angle_between_vectors(v1, v2)
    return angle * 180 / np.pi


def oop_angle(i, j, k, l):
    '''Returns the out-of-plane angle between points i,j,k,l

    Args:
        i (numpy.array): The first point.
        j (numpy.array): The second point.
        k (numpy.array): The third point.
        l (numpy.array): The fourth point.

    Returns:
        angle (float): The angle in degrees.
    '''
    elj = unit_vector(l - j)
    elk = unit_vector(l - k)
    eli = unit_vector(l - i)

    perp = np.cross(elj, elk)
    projection = np.dot(perp, eli)
    square_dist = np.dot(elj, elk)

    sine = projection / np.sqrt(square_dist)
    sine = min(1.0, max(-1.0, sine))
    angle = np.arcsin(sine)
    return angle


def torsion_angle(i, j, k, l):
    '''Returns the torsion angle between four points

    Args:
        i (numpy.array): The first point.
        j (numpy.array): The second point.
        k (numpy.array): The third point.
        l (numpy.array): The fourth point.

    Returns:
        angle (float): The angle in degrees.
    '''
    eij = unit_vector(i - j)
    ejk = unit_vector(j - k)
    ekl = unit_vector(k - l)

    perp_v1 = np.cross(eij, ejk)
    perp_v2 = np.cross(ejk, ekl)
    projection = np.dot(perp_v1, perp_v2)
    square_dist1 = np.dot(perp_v1, perp_v1)
    square_dist2 = np.dot(perp_v2, perp_v2)

    cosine = projection / np.sqrt(square_dist1 * square_dist2)
    cosine = min(1.0, max(-1.0, cosine))
    angle = np.arccos(cosine)
    return angle
