import re
from dataclasses import dataclass
import numpy as np
from ..molecule import *
from symel_generators import *

class PointGroup():
    def __init__(self, s, family, n, subfamily):
        self.str = s
        self.family = family
        self.n = n
        self.subfamily = subfamily
    
    @classmethod
    def from_string(cls, s):
        regex = r"([A-Z]+)(\d+)?([a-z]+)?"
        m = re.match(regex, s)
        family, n, subfamily = m.groups()
        if n is not None:
            n = int(n)
        if subfamily is not None:
            subfamily = str(subfamily)
        family = str(family)
        return cls(s, family, n, subfamily)

    def __str__(self):
        return "Family: "+self.family+"\nn: "+str(self.n)+"\nSubfamily: "+self.subfamily

@dataclass
class Symel():
    symbol:str
    rrep:np.array

class CharTable():
    def __init__(self) -> None:
        self.name = 0
        self.irreps = 0
        self.classes = 0
        self.class_orders = 0
        self.characters = 0
        self.irrep_dims = 0

class Symtext():
    def __init__(self) -> None:
        self.pg = 0
        self.symels = 0
        self.chartable = 0
        self.class_map = 0
        self.atom_map = 0
        self.mult_table = 0
        self.order = 0

def reduce(n, i):
    g = gcd(n, i)
    return n//g, i//g

def gcd(A, B):
    # A quick implementation of the Euclid algorithm for finding the greatest common divisor
    a = max(A,B)
    b = min(A,B)
    if a == 0:
        return b
    elif b == 0:
        return a
    else:
        r = a % b
        gcd(b, r)

if __name__ == "__main__":
    a = PointGroup.from_string("Ci")
    print(a)