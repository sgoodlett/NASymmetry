import numpy as np
from numpy.lib.financial import irr

class Irrep():
    def __init__(self, name, characters) -> None:
        self.name = name
        self.characters = characters
    def __repr__(self) -> str:
        return f"{self.name}:{*self.characters,}"
class DPTable():
    def __init__(self, pg, irreps, class_orders) -> None:
        self.pg = pg
        self.irreps = irreps
        self.nirreps = len(irreps)
        self.class_orders = class_orders
        self.order = sum(class_orders)
        self.table = np.zeros((self.nirreps, self.nirreps))
    def __repr__(self) -> str:
        strang = ""
        for i in range(self.nirreps):
            strang += f"\t{self.irreps[i].name}"
        for row in range(self.nirreps):
            strang += f"\n{self.irreps[row].name}"
            for column in range(self.nirreps):
                strang += "\t"
                result = self.direct_product(self.irreps[row], self.irreps[column])
                for r in result:
                    strang += f"{r.name} "
        return strang
    def direct_product(self, a, b):
        m = a.characters * b.characters * self.class_orders
        return self.find_irrep(m)
    def find_irrep(self, characters):
        for i in self.irreps:
            if np.array_equal(characters, i.characters):
                return [i]
        return self.decompose_rep(characters)
    def decompose_rep(self, characters):
        constituent_irreps = []
        for i in self.irreps:
            s = sum(characters * i.characters)
            n = s // self.order
            n = int(n)
            if n != 0:
                for mult in range(n):
                    constituent_irreps.append(i)
        return constituent_irreps

if __name__ == "__main__":
    pg = "C2v"
    #A1 = Irrep("A1", np.asarray([1.0,  1.0,  1.0,  1.0]))
    #A2 = Irrep("A2", np.asarray([1.0,  1.0, -1.0, -1.0]))
    #B1 = Irrep("B1", np.asarray([1.0, -1.0,  1.0, -1.0]))
    #B2 = Irrep("B2", np.asarray([1.0, -1.0, -1.0,  1.0]))
    #class_orders = [1,1,1,1]
    #irreps = [A1, A2, B1, B2]
    #pg = "C3v"
    #A1 = Irrep("A1", np.asarray([1.0,  1.0,  1.0]))
    #A2 = Irrep("A2", np.asarray([1.0,  1.0, -1.0]))
    #E  = Irrep("E" , np.asarray([2.0, -1.0,  0.0]))
    #irreps = [A1, A2, E]
    #class_orders = [1,2,3]
    pg = "Ih"
    pr5 = 0.5*(1 + np.sqrt(5))
    mr5 = 0.5*(1 - np.sqrt(5))
    Ag  = Irrep("Ag",  np.asarray([1.0,  1.0,  1.0,  1.0,  1.0,  1.0,  1.0,  1.0,  1.0,  1.0]))
    T1g = Irrep("T1g", np.asarray([3.0,  pr5,  mr5,  0.0, -1.0,  3.0,  mr5,  pr5,  0.0, -1.0]))
    T2g = Irrep("T2g", np.asarray([3.0,  mr5,  pr5,  0.0, -1.0,  3.0,  pr5,  mr5,  0.0, -1.0]))
    Gg  = Irrep("Gg",  np.asarray([4.0, -1.0, -1.0,  1.0,  0.0,  4.0, -1.0, -1.0,  1.0,  0.0]))
    Hg  = Irrep("Hg",  np.asarray([5.0,  0.0,  0.0, -1.0,  1.0,  5.0,  0.0,  0.0, -1.0,  1.0]))
    Au  = Irrep("Au",  np.asarray([1.0,  1.0,  1.0,  1.0,  1.0, -1.0, -1.0, -1.0, -1.0, -1.0]))
    T1u = Irrep("T1u", np.asarray([3.0,  pr5,  mr5,  0.0, -1.0, -3.0, -mr5, -pr5,  0.0,  1.0]))
    T2u = Irrep("T2u", np.asarray([3.0,  mr5,  pr5,  0.0, -1.0, -3.0, -pr5, -mr5,  0.0,  1.0]))
    Gu  = Irrep("Gu",  np.asarray([4.0, -1.0, -1.0,  1.0,  0.0, -4.0,  1.0,  1.0, -1.0,  0.0]))
    Hu  = Irrep("Hu",  np.asarray([5.0,  0.0,  0.0, -1.0,  1.0, -5.0,  0.0,  0.0,  1.0, -1.0]))
    irreps = [Ag, T1g, T2g, Gg, Hg]#, Au, T1u, T2u, Gu, Hu]
    class_orders = [1, 12, 12, 20, 15, 1, 12, 12, 20, 15]
    dptab = DPTable(pg, irreps, class_orders)
    print(dptab)
