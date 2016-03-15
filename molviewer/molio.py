'''
I/O functions

'''
import sys
import os
from os.path import basename, splitext
import requests
from atom import Atom
from molecule import Molecule
from moldata import sym2no
from config import data_dir, valid_extensions


####################
# Read files
####################
def count_frames(infile):
    '''Returns the number of frames in a xyz files

    All frames must be the same molecule

    Args:
        infile (str): filename
    '''
    with open(infile, 'r') as inf:
        num_atoms = int(inf.readline().strip())
        inf.seek(0)
        for i, l in enumerate(inf):
            pass
    num_frames = int((i + 1) / (num_atoms + 2))
    return num_frames


def read_xyz(infile):
    '''Read frames from a xyz format file

    Returns Molecule if only one frame in file, else retruns a list of Molecules.

    Args:
        infile (str): filename

    Returns:
        mol (molecule.Molecule): A molecule instance
    '''
    num_frames = count_frames(infile)
    print('Found {:d} frame(s) in file {:s}'.format(num_frames, os.path.realpath(infile)))

    mols = []
    with open(infile, 'r') as inf:
        filename = basename(infile)
        print('Reading file: {:s}'.format(infile))
        for i in range(num_frames):
            print('Reading frame {:d}'.format(i+1))
            num_atoms = int(next(inf).strip())
            name = next(inf).strip()
            atoms = []
            for i in range(num_atoms):
                line = next(inf)
                # 4 fields --> symbol x y z
                # 5 fields --> symbol x y z charge
                # 7 fields --> symbol x y z vec_x vec_y vec_z
                # 8 fields --> symbol x y z charge vec_x vec_y vec_z
                num_fields = len(line.split())
                if num_fields == 4:
                    symbol, x, y, z = line.split()
                    x = float(x)
                    y = float(y)
                    z = float(z)
                    atom = Atom(sym2no[symbol], [x, y, z])
                elif num_fields == 5:
                    symbol, x, y, z, q = line.split()
                    x = float(x)
                    y = float(y)
                    z = float(z)
                    q = float(q)
                    atom = Atom(sym2no[symbol], [x, y, z], charge=q)
                elif num_fields == 7:
                    symbol, x, y, z, vec_y, vec_y, vec_z = line.split()
                    x = float(x)
                    y = float(y)
                    z = float(z)
                    vec_x = float(vec_x)
                    vex_y = float(vec_y)
                    vec_z = float(vec_z)
                    # not using vectors for anything at the moment
                    # therefor we just drop them
                    atom = Atom(sym2no[symbol], [x, y, z])
                elif num_fields == 8:
                    symbol, x, y, z, q, vec_x, vec_y, vec_z = line.split()
                    x = float(x)
                    y = float(y)
                    z = float(z)
                    q = float(q)
                    vec_x = float(vec_x)
                    vex_y = float(vec_y)
                    vec_z = float(vec_z)
                    atom = Atom(sym2no[symbol], [x, y, z], charge=q)
                else:
                    print('Error reading frame.')
                atoms.append(atom)

            if num_frames == 1:
                return Molecule(atoms, name=name, filename=filename)

            mols.append(Molecule(atoms, name=name))
    return mols


def read_mol(infile):
    '''Read mol files

    Might work for sdf files. See: https://en.wikipedia.org/wiki/Chemical_table_file

    Args:
        infile (str): filename

    Returns:
        mol (molecule.Molecule): A molecule instance
    '''
    with open(infile, 'r') as inf:
        # Lines 1-3 are header.
        # Line 1 is formula. Line2 ???. Line 3 is empty.
        header1 = next(inf).strip()
        header2 = next(inf).strip()
        blank = next(inf)

        # Line 4 has number of atoms and bonds in first two fields
        num_atoms, num_bonds = next(inf).split()[:2]
        num_atoms = int(num_atoms)
        num_bonds = int(num_bonds)

        # The following lines have coordinates and element in first four fields
        # Coordinates are in angstrom.
        # Fields: x, y, z, element, ...other...
        atoms = []
        for i in range(num_atoms):
            x, y, z, symbol = next(inf).split()[:4]
            x = float(x)
            y = float(y)
            z = float(z)
            atom = Atom(sym2no[symbol], [x, y, z])
            atoms.append(atom)
        # We can return the molecule here or extract the remaining information
        return Molecule(atoms, name=header1)


def get_mol_by_name(molname, datadir=data_dir):
    '''Returns a molecule corresponding to the given identifier

    Searches in datadir for molecule and returns it if found.
    Otherwise downloads it, saves it to datadir and returns a molecule instance.

    Args:
        molname (str): The identifier string. Can be name, common name or CAS number
        datadir (str): The directory to search before downloading and to save
        downloaded data to.

    Returns:
        One of:
        mol (molecule.Molecule): A molecule instance
        or None (if molecule wasn't found)
    '''
    print('Searching in {:s}'.format(data_dir))
    hits = []
    for f in os.listdir(data_dir):
        name, ext = splitext(f)
        if name == molname and (ext in valid_extensions):
            hits.append(f)

    if len(hits) > 1:
        print('Found more than one candidate')
        mols = []
        for i in hits:
            mols.append(read_xyz(hits[i]))
    elif len(hits) == 1:
        return read_xyz(data_dir + '/' + hits[0])
    else:
        print('Downloading {:s} data'.format(molname))
        molstring = get_mol_from_url(molname)
        if molstring:
            filename = data_dir + '/' + molname + '.xyz'
            with open(filename, 'w') as outfile:
                outfile.write(molstring)
            print('Wrote file {:s}'.format(filename))
            return read_xyz(filename)
        return None


def get_mol_from_url(molname):
    '''Downloads a mol file in xyz format and returns it

    This function uses http://cactus.nci.nih.gov/ to search for the molecule
    using the given identifier.

    Args:
        molname (str): The identifier
    Returns:
        answer (str): The query result
        or
        None on failure.
    '''
    assert len(molname) >= 1

    url = 'http://cactus.nci.nih.gov/chemical/structure/' + molname + '/file?format=xyz&get3d=true'
    r = requests.get(url)
    if r.status_code == 200:
        return r.text
    elif r.status_code == 404:
        print('Unknow molecule: {:s}'.format(molname))
        return None
    elif r.status_code == 500:
        print('Invalid request string')
        return None
    else:
        print('Some error occcured while retrieving molecule data')
        return None


####################
# Write files
####################
def write_xyz(mol, file_object=None):
    '''Write molecule to file in xyz format

    Args:
        mol (molecule.Molecule): The molecule instance to write
        file_object (str): If given write data to this file.
        If None, print to stdout (default: None).

    Returns:
        None
    '''
    if file_object:
        with open(file_object, 'w') as outfile:
            outfile.write(mol.xyz)
        print('Wrote file {:s}'.format(file_object))
    else:
        print(record)


####################
# View molecules
####################
def view(mol, viewer='native'):
    '''Render the molecule

    The mayavi backend doesn't work under python 3. The native backend uses
    visvis to render the molecule. This is very slow.
    It's better to use the molecular viewer.

    Args:
        mol (molecule.Molecule): The molecule instance to render
        viewer: The backend to use. Valid choices are 'native', 'maya'
        and, 'avogadro' (default: native).

    Rturns:
        None
    '''
    # mayavi
    if viewer == 'maya':
        from mayavi import mlab
        for atom in mol.atoms:
            pts = mlab.points3d(atom.r[0], atom.r[1], atom.r[2],
                                scale_factor=0.75,
                                scale_mode='none',
                                resolution=20,
                                color=atom.color())
        for i, j in mol.bonds():
            mlab.plot3d([mol.atoms[i].r[0], mol.atoms[j].r[0]],
                        [mol.atoms[i].r[1], mol.atoms[j].r[1]],
                        [mol.atoms[i].r[2], mol.atoms[j].r[2]],
                        tube_radius=0.1,
                        tube_sides=20)
        mlab.show()

    # avogadro
    if viewer == 'avogadro':
        from subprocess import call
        write_xyz(mol, 'avogadro.xyz')
        call(['avogadro', 'avogadro.xyz'])
        call(['rm', 'avogadro.xyz'])

    # visvis
    if viewer == 'native':
        import visvis as vv

        for atom in mol.atoms:
            x, y, z = atom.r
            at = vv.solidSphere((x, y, z), atom.radius()*0.25)
            at.faceColor = atom.color()

        for bond in mol.bonds():
            pp = vv.Pointset(3)
            pp.append(bond.atoms[0].r)
            pp.append(bond.atoms[1].r)
            vv.solidLine(pp, radius=0.15, N=16)

        pp = vv.Pointset(3)
        pp.append([3, 3, 3])
        pp.append([4, 3, 3])
        x = vv.solidLine(pp, radius=0.05)
        x.faceColor = 'r'
        conex = vv.solidCone([4, 3, 3], scaling=[0.1, 0.1, 0.1], direction=[1, 0, 0])
        conex.faceColor = 'r'

        pp = vv.Pointset(3)
        pp.append([3, 3, 3])
        pp.append([3, 4, 3])
        y = vv.solidLine(pp, radius=0.05)
        y.faceColor = 'g'
        coney = vv.solidCone([3, 4, 3], scaling=[0.1, 0.1, 0.1], direction=[0, 1, 0])
        coney.faceColor = 'g'

        pp = vv.Pointset(3)
        pp.append([3, 3, 3])
        pp.append([3, 3, 4])
        z = vv.solidLine(pp, radius=0.05)
        z.faceColor = 'b'
        conez = vv.solidCone([3, 3, 4], scaling=[0.1, 0.1, 0.1], direction=[0, 0, 1])
        conez.faceColor = 'b'

        # Set axes settings
        axes = vv.gca()
        axes.SetLimits(rangeX=(-5, 5), rangeY=(-5, 5), rangeZ=(-5, 5))
        vv.axis('off')
        app = vv.use()
        app.Run()


if __name__ == '__main__':
    mol = get_mol_by_name('acetic acid')
    print(mol.mass())
    view(mol)
