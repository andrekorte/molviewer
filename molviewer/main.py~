import sys
from tabulate import tabulate
from molio import read_xyz, read_mol
from rotation import inertia_tensor, inertia_tensor_amu, find_symmetry, rotational_constants
from molfuncs import distance_matrix, angle, oop_angle, torsion_angle
from molio import view
from moldata import wavenumber2ghz


# Functions for pretty printing
#
#
def print_seperator(s):
    '''Prints a seperator

    Args:
        s (str): A string to insert in the seperator.

    The seperator string looks like this:

    ************************************************************
    *                         string                           *
    ************************************************************
    '''
    print('')
    print(60*'*')
    string = str(s)
    print('*' + string.center(58, ' ') + '*')
    print(60*'*')


def print_center(s, n=60):
    '''Prints a centerd string

    Args:
        s (str): The string to print.
        n (int): Line length (default:60).
    '''
    print(s.center(n, ' '))


def print_end_of_output():
    '''Prints the last line of the output'''
    print('')
    string = ' END OF OUTPUT '
    print(string.center(60, '*'))


def print_distance_matrix(matrix):
    '''Prints a nicely formatted distance matrix table

    Printing is skipped when there are more than twelve atoms in the molecule.

    Args:
        matrix (numpy.array): The distance matrix array.
    '''
    if matrix.shape > (12, 12):
        print('Skipping printing of distance matrix')
        return
    headers = [i for i, _ in enumerate(matrix)]
    print(tabulate(matrix, headers, showindex=True))


def run_all():
    # General Information
    print_seperator('General information')
    print('Name: ', m.name)
    print('Formula: {0:s}\nMass: {1:12.6f}'.format(m.stoichiometry, m.mass))

    # Coordinates
    print_seperator('Cartesian coordinates in xyz format (angstrom)')
    print(m.xyz)

    # Distance matrix
    print_seperator('Distance matrix (angstrom)')
    dist_mat = distance_matrix(m)
    print_distance_matrix(dist_mat)

    # Rotational constants
    I_amu = inertia_tensor_amu(m)
    I = inertia_tensor(m)
    A_cm = rotational_constants(I)[0]
    B_cm = rotational_constants(I)[1]
    C_cm = rotational_constants(I)[2]
    A_GHz = rotational_constants(I)[0] * wavenumber2ghz
    B_GHz = rotational_constants(I)[1] * wavenumber2ghz
    C_GHz = rotational_constants(I)[2] * wavenumber2ghz

    print_seperator('Rotational constants')
    print('Molecule is {:s}'.format(find_symmetry(I_amu)))
    print('')
    print_center('Rotational constants in cm-1')
    print_center('A {0:6.2f}     B {1:6.2f}     C {2:6.2f}'.format(A_cm, B_cm, C_cm))
    print('')
    print_center('Rotational constants in GHz')
    print_center('A {0:6.2f}     B {1:6.2f}     C {2:6.2f}'.format(A_GHz, B_GHz, C_GHz))

    print_end_of_output()


if __name__ == '__main__':
    # File Information
    print_seperator('I/O Information')
    infile = sys.argv[1]
    ext = infile.split('.')[-1]
    if ext == 'xyz':
        m = read_xyz(infile)
    elif ext == 'mol':
        m = read_mol(infile)
    elif ext == 'sdf':
        m = read_mol(infile)
    else:
        sys.exit('Unknown file format')

    if isinstance(m, list):
        for i, m in enumerate(m):
            print('! Begin frame {:d}'.format(i+1))
            run_all()
            print('! End frame {:d}'.format(i+1))
    else:
        run_all()
        #view(m, 'avogadro')
