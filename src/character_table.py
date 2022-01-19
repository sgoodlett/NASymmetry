import numpy as np
from dataclasses import dataclass
import re

@dataclass
class SymmetryOperation:
    name: str
    order: int
    power: int
    axis: np.array
    cla: int

@dataclass
class Irrep:
    name: str
    characters: np.array
    idx: int

@dataclass
class DPResult:
    irreps: list
    mult: np.array

class CharacterTable():
    def __init__(self, name, irreps, symops, rep_symops):
        self.name = name
        self.irreps = irreps
        self.symops = symops
        self.rep_symops = rep_symops
        self.nirreps = len(irreps)
        self.table = np.zeros((self.nirreps, self.nirreps))
        self.h = 0
        for i, irrep in enumerate(self.irreps):
            self.h += irrep.characters[0]
            self.table[i,:] = irrep.characters
        self.class_orders = np.zeros((self.nirreps, ))
        for s, symop in enumerate(self.symops):
            self.class_orders[symop.cla] += 1
    
    @classmethod
    def from_libmsym(cls, ctx, pg):
        # Convert libmsym objects to our objects
        symops = []
        rgx = re.compile(r"^libmsym.SymmetryOperation\(\s*([^,\s]*)\s*")
        for s, symop in enumerate(ctx.symmetry_operations):
            m = re.match(rgx, str(symop))
            if m:
                name = m.group(1)
            else:
                raise RuntimeError("No name found for libmsym symmetry operation!")
            symops.append(SymmetryOperation(name, symop.order, symop.power, np.array(symop.vector), symop.conjugacy_class))
        irreps = []
        for i, irrep in enumerate(ctx.subrepresentation_spaces):
            irreps.append(Irrep(ctx.character_table.symmetry_species[i].name, ctx.character_table.table[i,:], i))
        rsymops = []
        for i in ctx.character_table.symmetry_operations:
            m = re.match(rgx, str(i))
            if m:
                rsymops.append(m.group(1))
            else:
                raise RuntimeError("No name found for libmsym symmetry operation!")
        
        return cls(pg, irreps, symops, rsymops)

    def __str__(self):
        strang = f"\nCharacter Table for {self.name} point group\n"
        for i, rsym in enumerate(self.rep_symops):
            if self.class_orders[i] == 1:
                n = ""
            else:
                n = int(self.class_orders[i])
            strang += f"\t{n}{rsym}"
        for i, irrep in enumerate(self.irreps):
            strang += f"\n{irrep.name}\t"
            for j in irrep.characters:
                strang += f"{int(j):2d}\t"
        return strang

    def __repr__(self):
        return self.__str__()

    def full_direct_product(self, a, b, *args):
        chars = a.characters * b.characters
        for arg in args:
            chars *= arg.characters
        m = chars * self.class_orders
        return self.find_irrep(m)

    def dp_contains(self, irrep, a, b, *args):
        chars = a.characters * b.characters
        for arg in args:
            chars *= arg.characters
        s = sum(chars * self.class_orders * irrep.characters)
        n = s // self.h
        if n > 0:
            return True
        return False

    def idx_from_irrep(self, irrep):
        for i, selfirreps in enumerate(self.irreps):
            if irrep == selfirreps:
                return i
        return RuntimeError("No match found in idx_from_irrep")

    def _direct_product(self, a, b):
        m = a.characters * b.characters * self.class_orders
        return self.find_irrep(m)

    def find_irrep(self, characters):
        for idx, irrep in enumerate(self.irreps):
            if np.array_equal(characters, irrep.characters):
                return DPResult([irrep], [1])
        return self.decompose_rep(characters)

    def decompose_rep(self, characters):
        constituent_irreps = []
        multiplicities = []
        for idx, irrep in enumerate(self.irreps):
            s = sum(characters * irrep.characters)
            n = s // self.order
            n = int(n)
            if n != 0:
                constituent_irreps.append(irrep)
                multiplicities.append(n)
        return DPResult(constituent_irreps, np.array(multiplicities))