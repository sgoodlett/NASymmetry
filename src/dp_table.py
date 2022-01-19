from typing import List
import numpy as np
from numpy.lib.financial import irr
from dataclasses import dataclass

@dataclass
class Irrep:
    name: str
    characters: np.array
    idx: int

@dataclass
class DPResult:
    irreps: List
    mult: np.array

class DPTable():
    def __init__(self, pg, irreps, class_orders) -> None:
        self.pg = pg
        self.irreps = irreps
        self.nirreps = len(irreps)
        self.irreps_byidx = range(self.nirreps)
        self.class_orders = class_orders
        self.order = sum(class_orders)
        self.table = np.empty((self.nirreps, self.nirreps), dtype=DPResult)
        self.build_table()

    def __repr__(self) -> str:
        strang = ""
        for i in range(self.nirreps):
            strang += f"\t{self.irreps[i].name}"
        for row in range(self.nirreps):
            strang += f"\n{self.irreps[row].name}"
            for column in range(self.nirreps):
                strang += "\t"
                result = self._direct_product(self.irreps[row], self.irreps[column])
                first_iter = True
                for r, irrep in enumerate(result.irreps):
                    if first_iter:
                        first_iter = False
                    else:
                        strang += "+"
                    strang += f"{result.mult[r]}{self.irreps[irrep.idx].name}"
        return strang

    def build_table(self):
        for i, row_irrep in enumerate(self.irreps):
            for j, column_irrep in enumerate(self.irreps):
                if i >= j:
                    self.table[i,j] = self._direct_product(row_irrep, column_irrep)
                    self.table[j,i] = self.table[i,j]

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
        n = s // self.order
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

if __name__ == "__main__":
    pg = "Ih"
    if pg == "C2v":
        A1 = Irrep("A1", np.asarray([1.0,  1.0,  1.0,  1.0]))
        A2 = Irrep("A2", np.asarray([1.0,  1.0, -1.0, -1.0]))
        B1 = Irrep("B1", np.asarray([1.0, -1.0,  1.0, -1.0]))
        B2 = Irrep("B2", np.asarray([1.0, -1.0, -1.0,  1.0]))
        class_orders = [1,1,1,1]
        irreps = [A1, A2, B1, B2]
    elif pg == "C3v":
        A1 = Irrep("A1", np.asarray([1.0,  1.0,  1.0]), 0)
        A2 = Irrep("A2", np.asarray([1.0,  1.0, -1.0]), 1)
        E  = Irrep("E" , np.asarray([2.0, -1.0,  0.0]), 2)
        irreps = [A1, A2, E]
        class_orders = [1,2,3]
    elif pg == "Ih":
        pr5 = 0.5*(1 + np.sqrt(5))
        mr5 = 0.5*(1 - np.sqrt(5))
        Ag  = Irrep("Ag",  np.asarray([1.0,  1.0,  1.0,  1.0,  1.0,  1.0,  1.0,  1.0,  1.0,  1.0]), 0)
        T1g = Irrep("T1g", np.asarray([3.0,  pr5,  mr5,  0.0, -1.0,  3.0,  mr5,  pr5,  0.0, -1.0]), 1)
        T2g = Irrep("T2g", np.asarray([3.0,  mr5,  pr5,  0.0, -1.0,  3.0,  pr5,  mr5,  0.0, -1.0]), 2)
        Gg  = Irrep("Gg",  np.asarray([4.0, -1.0, -1.0,  1.0,  0.0,  4.0, -1.0, -1.0,  1.0,  0.0]), 3)
        Hg  = Irrep("Hg",  np.asarray([5.0,  0.0,  0.0, -1.0,  1.0,  5.0,  0.0,  0.0, -1.0,  1.0]), 4)
        Au  = Irrep("Au",  np.asarray([1.0,  1.0,  1.0,  1.0,  1.0, -1.0, -1.0, -1.0, -1.0, -1.0]), 5)
        T1u = Irrep("T1u", np.asarray([3.0,  pr5,  mr5,  0.0, -1.0, -3.0, -mr5, -pr5,  0.0,  1.0]), 6)
        T2u = Irrep("T2u", np.asarray([3.0,  mr5,  pr5,  0.0, -1.0, -3.0, -pr5, -mr5,  0.0,  1.0]), 7)
        Gu  = Irrep("Gu",  np.asarray([4.0, -1.0, -1.0,  1.0,  0.0, -4.0,  1.0,  1.0, -1.0,  0.0]), 8)
        Hu  = Irrep("Hu",  np.asarray([5.0,  0.0,  0.0, -1.0,  1.0, -5.0,  0.0,  0.0,  1.0, -1.0]), 9)
        irreps = [Ag, T1g, T2g, Gg, Hg, Au, T1u, T2u, Gu, Hu]
        class_orders = [1, 12, 12, 20, 15, 1, 12, 12, 20, 15]
    dptab = DPTable(pg, irreps, class_orders)
    #print(dptab)
    for i in range(1000000):
        dptab.dp_contains(Ag, T1g, T2g, Gg, Hg, Au, T1u, T2u, Gu, Hu)
    #for i in range(1000000):
    #    dptab.full_direct_product(Ag, T1g, T2g, Gg, Hg, Au, T1u, T2u, Gu, Hu)
