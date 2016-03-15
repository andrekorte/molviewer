testmol = '''
14
C3H8O3
C    1.2492   -0.6653    0.2948
C    0.0000    0.0339   -0.2451
C   -1.2492   -0.6653    0.2948
O   -2.4131   -0.0873   -0.2998
O   -0.0000    1.4000    0.1742
O    2.4131   -0.0873   -0.2998
H    1.2952   -0.5428    1.3769
H    1.2056   -1.7267    0.0508
H    0.0000   -0.0119   -1.3341
H   -1.2056   -1.7267    0.0508
H   -1.2952   -0.5428    1.3769
H   -3.2441   -0.4828   -0.0029
H   -0.0000    1.5168    1.1342
H    3.2441   -0.4828   -0.0029
'''


class BaseParser(object):
    def __init__(self):
        pass

class XyzParser(object):
    def __init__(self, source):
        if isinstance(source, str):
            self.filename = source
        if hasattr(source, 'read'):
            print('Banana')
    
    def __str__(self):
        '''Return a string representation of the object.'''
        return 'XYZ format file {:s}'.format(self.filename)

    def __repr__(self):
        '''Return a representation of the object.'''
        return 'XYZ format("{:s}")'.format(self.filename)

if __name__ == '__main__':
    test2 = open('1,2,3-propantriol.xyz')
    parser = XyzParser(test2)
    test2.close()







#def count_frames(infile):
#    '''
#    Returns the number of frames in a xyz files
#    All frames must be the same molecule
#    '''
#    with open(infile, 'r') as inf:
#        num_atoms = int(inf.readline().strip())
#        inf.seek(0)
#        for i, l in enumerate(inf):
#            pass
#    num_frames = int((i + 1) / (num_atoms + 2))
#    return num_frames


#def read_xyz(infile):
#    '''
#    Read frames from a xyz format file
#    Returns Molecule if only one frame in file,
#    else retruns a list of Molecules.
#    '''
#    num_frames = count_frames(infile)
#    print('Found {:d} frame(s) in file {:s}'.format(num_frames, os.path.realpath(infile)))

#    mols = []
#    with open(infile, 'r') as inf:
#        filename = basename(infile)
#        print('Reading file: {:s}'.format(infile))
#        for i in range(num_frames):
#            print('Reading frame {:d}'.format(i+1))
#            num_atoms = int(next(inf).strip())
#            name = next(inf).strip()
#            atoms = []
#            for i in range(num_atoms):
#                line = next(inf)
#                # 4 fields --> symbol x y z
#                # 5 fields --> symbol x y z charge
#                # 7 fields --> symbol x y z vec_x vec_y vec_z
#                # 8 fields --> symbol x y z charge vec_x vec_y vec_z
#                num_fields = len(line.split())
#                if num_fields == 4:
#                    symbol, x, y, z = line.split()
#                    x = float(x)
#                    y = float(y)
#                    z = float(z)
#                    atom = Atom(sym2no[symbol], [x, y, z])
#                elif num_fields == 5:
#                    symbol, x, y, z, q = line.split()
#                    x = float(x)
#                    y = float(y)
#                    z = float(z)
#                    q = float(q)
#                    atom = Atom(sym2no[symbol], [x, y, z], charge=q)
#                elif num_fields == 7:
#                    symbol, x, y, z, vec_y, vec_y, vec_z = line.split()
#                    x = float(x)
#                    y = float(y)
#                    z = float(z)
#                    vec_x = float(vec_x)
#                    vex_y = float(vec_y)
#                    vec_z = float(vec_z)
#                    # not using vectors for anything at the moment
#                    # therefor we just drop them
#                    atom = Atom(sym2no[symbol], [x, y, z])
#                elif num_fields == 8:
#                    symbol, x, y, z, q, vec_x, vec_y, vec_z = line.split()
#                    x = float(x)
#                    y = float(y)
#                    z = float(z)
#                    q = float(q)
#                    vec_x = float(vec_x)
#                    vex_y = float(vec_y)
#                    vec_z = float(vec_z)
#                    atom = Atom(sym2no[symbol], [x, y, z], charge=q)
#                else:
#                    print('Error reading frame.')
#                atoms.append(atom)

#            if num_frames == 1:
#                return Molecule(atoms, name=name, filename=filename)

#            mols.append(Molecule(atoms, name=name))
#    return mols
