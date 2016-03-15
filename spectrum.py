# -*- coding: utf-8 -*-
'''

'''
from random import randint
from moldata import spectrumcolors
from jcamp import JCAMP_reader

class Spectrum(object):
    def __init__(self, filename):
        dictionary = JCAMP_reader(filename)
        for k, v in dictionary.items():
            setattr(self, k.replace(' ', ''), v)

        self.color = None

if __name__ == '__main__':
    s = Spectrum('chloroform.dx')
