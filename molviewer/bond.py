# -*- coding: utf-8 -*-
'''
Bonds exist within a molecule.

'''


class Bond(object):
    '''The bond object.

    :param i: The index of the first atom in the bond.
    :type i: int
    :param j: The index of the second atom in the bond.
    :type j: int
    :param ati: Atom instance of the first atom.
    :type ati: atom.Atom
    :param atj: Atom instance of the second atom.
    :type atj: atom.Atom
    :param ishbond: Wether the bond is a hydrogen bond (optional).
    :type ishbond: bool
    :param bondtype: The type of bond (optional).
    :type bondtype: str
    :param p1: The position of the first atom.
    :type p1: numpy.ndarray
    :param p2: The position of the second atom.
    :type p2: numpy.ndarray

    '''

    def __init__(self, i, j, ati, atj, **kwargs):
        self.type = None  # One of single, double, triple, arom
        self.ishbond = kwargs.get('ishbond', False)
        self.index = [i, j]
        self.atoms = [ati, atj]
        self.p1 = ati.r
        self.p2 = atj.r

    def __repr__(self):
        s = '{0:s}{1:d}\t{2:s}{3:d}\t{4:12.6f}'.format(self.atoms[0].symbol, self.index[0],
                                                       self.atoms[1].symbol, self.index[1],
                                                       self.length()
                                                       )
        return s

    def asdict(self):
        '''Returns the bond information as a dict.

        :returns: dict -- The bond information dictionary.

        '''
        d = {
            'type': self.type,
            'ishbond': self.ishbond,
            'index': (self.index[0], self.index[1]),
            'atoms': (self.atoms[0], self.atoms[1]),
            'p1': self.p1,
            'p2': self.p2
        }
        return d

    def length(self):
        '''Return length of bond.

        :returns: float -- The euclidian distance.
        
        '''
        return self.atoms[0].distance(self.atoms[1])
