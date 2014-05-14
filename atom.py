# atom.py
# Justin Clark - May 2014

import numpy as np

class Atom:
    def __init__(self, aid, name, residue, x, y, z, element):
        self.aid = aid
        self.name = name
        self.residue = residue
        self.x = x
        self.y = y
        self.z = z
        self.element = element
        self.coords = np.array([x, y, z])

    def get_coords(self):
        return self.coords

    def set_coords(self, coords):
        self.coords = np.array(coords)
        self.x, self.y, self.z = coords
