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
    idx: int

@dataclass
class Irrep:
    name: str
    characters: np.array
    idx: int
    orbital_idxs: list

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
            self.h += irrep.characters[0]**2
            self.table[i,:] = irrep.characters
        self.h = int(self.h)
        self.class_orders = np.zeros((self.nirreps, ))
        for s, symop in enumerate(self.symops):
            self.class_orders[symop.cla] += 1
        self.build_mult_table()
    
    @classmethod
    def from_libmsym(cls, ctx, pg, irrep_idx_lengths):
        # Convert libmsym objects to our objects
        symops = []
        rgx = re.compile(r"^libmsym.SymmetryOperation\(\s*([^,\s]*)\s*")
        for s, symop in enumerate(ctx.symmetry_operations):
            m = re.match(rgx, str(symop))
            if m:
                name = m.group(1)
            else:
                raise RuntimeError("No name found for libmsym symmetry operation!")
            symops.append(SymmetryOperation(name, symop.order, symop.power, np.array(symop.vector), symop.conjugacy_class, s))
        irreps = []
        counter = 0
        for i, irrep in enumerate(ctx.subrepresentation_spaces):
            if irrep_idx_lengths == 0:
                idxs = []
            else:
                idxs = list(range(counter, counter+irrep_idx_lengths[i]))
                counter += irrep_idx_lengths[i]
            irreps.append(Irrep(ctx.character_table.symmetry_species[i].name, ctx.character_table.table[i,:], i, idxs))
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
        n = s / self.h
        if n > 0:
            return True
        return False

    def dp_symmetric(self, a, b, *args):
        pass

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

    def Cn(self, order, power, axis):
        ang = 2*np.pi*power / order
        sintheta = np.sin(ang)
        costheta = np.cos(ang)
        v = axis / np.linalg.norm(axis)
        a = [0,1,2]
        O = np.zeros((3,3))
        O += (1-costheta)
        for i in a:
            for j in a:
                if i == j:
                    O[i,i] *= v[i]**2
                    O[i,i] += costheta
                else:
                    O[i,j] *= v[i]*v[j]
                    b = [i,j]
                    C = [x for x in a if x not in b][0]
                    if i < j:
                        O[i,j] += (-1)**(i+j) * v[C]*sintheta
                    else:
                        O[i,j] += (-1)**(i+j-1) * v[C]*sintheta
        return O

    def sigma(self, axis):
        O = np.zeros((3,3))
        v = axis / np.linalg.norm(axis)
        for i in range(3):
            for j in range(i,3):
                if i == j:
                    O[i,i] += 1 - 2*v[i]**2
                else:
                    O[i,j] += -2 * v[i]*v[j]
                    O[j,i] += O[i,j]
        return O

    def Sn(self, order, power, axis):
        return np.dot(self.Cn(order, power, axis), self.sigma(axis))

    def get_matrix_rep(self, I):
        if I.name == "E":
            return np.identity(3)
        elif I.name == "i":
            return -1*np.identity(3)
        elif I.name[0] == "C":
            return self.Cn(I.order, I.power, I.axis)
        elif I.name[0] == "\u03C3":
            return self.sigma(I.axis)
        elif I.name[0] == "S":
            return self.Sn(I.order, I.power, I.axis)
        else:
            raise ValueError("Unfamiliar symmetry operation!")

    def abidx(self, A, B):
        C = np.dot(A,B)
        for i, r in enumerate(self.rreps):
            if np.allclose(C,r,rtol=1e-5,atol=1e-2):
                return i

    def build_mult_table(self):
        rreps = []
        for symop in self.symops:
            rreps.append(self.get_matrix_rep(symop))
        self.rreps = rreps
        mtab = np.zeros((self.h, self.h), dtype=int)
        for i in range(self.h):
            mtab[0,i] = i
            mtab[i,0] = i
        for i in range(1, self.h):
            for j in range(1, self.h):
                mtab[i,j] = self.abidx(self.rreps[i], self.rreps[j])
        self.mtab = mtab

    def multiply_symops(self, a, b):
        r = self.mtab[a.idx, b.idx]
        return self.symops[r]