from libmsym_wrap import LibmsymWrapper
from character_table import CharacterTable, SymmetryOperation
import numpy as np
import psi4

from ammonia import Settings
molecule = psi4.geometry(Settings["molecule"])
basis = psi4.core.BasisSet.build(molecule, "BASIS", Settings["basis"])
def get_basis(molecule):
    num = molecule.natom()
    molecule_basis = []
    for x in range(0, num):
        atom_basis = []
        molecule_basis.append(str(molecule.symbol(x)))
        for y in range(0, basis.nshell_on_center(x)):
            atom_basis.append(basis.shell(y).am)
        molecule_basis.append(atom_basis)
    return molecule_basis

molecule_basis = get_basis(molecule)

exec_libmsym = LibmsymWrapper(molecule, molecule_basis)

exec_libmsym.run()
ctab = exec_libmsym.ctab
print(ctab.multiply_symops(ctab.symops[1], ctab.symops[2]))
#S = SymmetryOperation("Bob", 3, 2, np.array([1,2,3]), 0)
#exec_libmsym.ctab.get_matrix_rep(S)
#print(exec_libmsym.ctab)


