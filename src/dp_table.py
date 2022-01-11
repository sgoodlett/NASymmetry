from typing import List
import numpy as np
from numpy.lib.financial import irr
from dataclasses import dataclass

@dataclass
class Irrep:
    name: str
    characters: np.array

@dataclass
class DPResult:
    irreps: List
    mult: List

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
                    strang += f"{result.mult[r]}{self.irreps[irrep].name}"
        return strang

    def build_table(self):
        for i, row_irrep in enumerate(self.irreps):
            for j, column_irrep in enumerate(self.irreps):
                if i >= j:
                    self.table[i,j] = self._direct_product(row_irrep, column_irrep)
                    self.table[j,i] = self.table[i,j]

    def direct_product(self, *args):
        print("AHHHHHHHHHHHHHHHHHH")

    def dp_from_table(self, *args):
        pass

    def _direct_product(self, a, b):
        m = a.characters * b.characters * self.class_orders
        return self.find_irrep(m)

    def find_irrep(self, characters):
        for idx, irrep in enumerate(self.irreps):
            if np.array_equal(characters, irrep.characters):
                return DPResult([idx], [1])
        return self.decompose_rep(characters)

    def decompose_rep(self, characters):
        constituent_irreps = []
        multiplicities = []
        for idx, irrep in enumerate(self.irreps):
            s = sum(characters * irrep.characters)
            n = s // self.order
            n = int(n)
            if n != 0:
                constituent_irreps.append(idx)
                multiplicities.append(n)
        return DPResult(constituent_irreps, multiplicities)

if __name__ == "__main__":
    #pg = "C2v"
    #A1 = Irrep("A1", np.asarray([1.0,  1.0,  1.0,  1.0]))
    #A2 = Irrep("A2", np.asarray([1.0,  1.0, -1.0, -1.0]))
    #B1 = Irrep("B1", np.asarray([1.0, -1.0,  1.0, -1.0]))
    #B2 = Irrep("B2", np.asarray([1.0, -1.0, -1.0,  1.0]))
    #class_orders = [1,1,1,1]
    #irreps = [A1, A2, B1, B2]
    pg = "C3v"
    A1 = Irrep("A1", np.asarray([1.0,  1.0,  1.0]))
    A2 = Irrep("A2", np.asarray([1.0,  1.0, -1.0]))
    E  = Irrep("E" , np.asarray([2.0, -1.0,  0.0]))
    irreps = [A1, A2, E]
    class_orders = [1,2,3]
    #pg = "Ih"
    #pr5 = 0.5*(1 + np.sqrt(5))
    #mr5 = 0.5*(1 - np.sqrt(5))
    #Ag  = Irrep("Ag",  np.asarray([1.0,  1.0,  1.0,  1.0,  1.0,  1.0,  1.0,  1.0,  1.0,  1.0]))
    #T1g = Irrep("T1g", np.asarray([3.0,  pr5,  mr5,  0.0, -1.0,  3.0,  mr5,  pr5,  0.0, -1.0]))
    #T2g = Irrep("T2g", np.asarray([3.0,  mr5,  pr5,  0.0, -1.0,  3.0,  pr5,  mr5,  0.0, -1.0]))
    #Gg  = Irrep("Gg",  np.asarray([4.0, -1.0, -1.0,  1.0,  0.0,  4.0, -1.0, -1.0,  1.0,  0.0]))
    #Hg  = Irrep("Hg",  np.asarray([5.0,  0.0,  0.0, -1.0,  1.0,  5.0,  0.0,  0.0, -1.0,  1.0]))
    #Au  = Irrep("Au",  np.asarray([1.0,  1.0,  1.0,  1.0,  1.0, -1.0, -1.0, -1.0, -1.0, -1.0]))
    #T1u = Irrep("T1u", np.asarray([3.0,  pr5,  mr5,  0.0, -1.0, -3.0, -mr5, -pr5,  0.0,  1.0]))
    #T2u = Irrep("T2u", np.asarray([3.0,  mr5,  pr5,  0.0, -1.0, -3.0, -pr5, -mr5,  0.0,  1.0]))
    #Gu  = Irrep("Gu",  np.asarray([4.0, -1.0, -1.0,  1.0,  0.0, -4.0,  1.0,  1.0, -1.0,  0.0]))
    #Hu  = Irrep("Hu",  np.asarray([5.0,  0.0,  0.0, -1.0,  1.0, -5.0,  0.0,  0.0,  1.0, -1.0]))
    #irreps = [Ag, T1g, T2g, Gg, Hg]#, Au, T1u, T2u, Gu, Hu]
    #class_orders = [1, 12, 12, 20, 15, 1, 12, 12, 20, 15]
    dptab = DPTable(pg, irreps, class_orders)
    print(dptab)
