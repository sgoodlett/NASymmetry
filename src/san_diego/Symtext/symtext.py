
class PointGroup():
    def __init__(self) -> None:
        self.str = 0
        self.family = 0
        self.n = 0
        self.subfamily = 0
    
    @classmethod
    def from_string(cls):
        pass

class Symel():
    def __init__(self) -> None:
        self.symbol = 0
        self.rrep = 0

    def __mul__(self):
        pass

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
