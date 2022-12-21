import re
from dataclasses import dataclass
import numpy as np
from san_diego.molecule import *
#from symel_generators import *

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
    def __str__(self) -> str:
        return f"\nSymbol: {self.symbol}, Reducible Rep.:\n{self.rrep}"
    def __repr__(self) -> str:
        return self.__str__()

class CharTable():
    def __init__(self, pg, irreps, classes, class_orders, chars, irrep_dims) -> None:
        self.name = pg
        self.irreps = irreps
        self.classes = classes
        self.class_orders = class_orders
        self.characters = chars
        self.irrep_dims = irrep_dims
    def __repr__(self) -> str:
        return f"Character Table for {self.name}\nIrreps: {self.irreps}\nClasses: {self.classes}\nCharacters:\n{self.characters}\n"

class Symtext():
    def __init__(self, pg, symels, chartable, class_map, atom_map, mult_table) -> None:
        self.pg = pg
        self.symels = symels
        self.chartable = chartable
        self.class_map = class_map
        self.atom_map = atom_map
        self.mult_table = mult_table
        self.order = len(symels)

def reduce(n, i):
    g = gcd(n, i)
    return n//g, i//g # floor divide to get an int, there should never be a remainder since we are dividing by the gcd

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
        return gcd(b, r)

if __name__ == "__main__":
    a = PointGroup.from_string("Ci")
    print(a)