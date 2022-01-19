import libmsym as msym
import numpy as np
import psi4
import pprint
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

def find_bf_idx(target_bf, bf_map):
    for i,bf in bf_map:
        if target_bf == bf:
            return i

def set_basis(element, z):
    basis_functions = []
    check = []
    z = z*2 
    bf_on_atoms = []
    am_vals = molecule_basis[z + 1]
    iterator = 0
    for l in am_vals:
        n = l + 1 + iterator
        for m in range(0,l+1):
            if m == 0:
                bf = msym.RealSphericalHarmonic(element = element, n = n, l = l, m = m, name =  "whocares")
                basis_functions.append(bf)
                check.append([l, m])
            else:
                bf = msym.RealSphericalHarmonic(element = element, n = n, l = l, m = m, name =  "whocares")
                basis_functions.append(bf)
                check.append([l, m])
                bf = msym.RealSphericalHarmonic(element = element, n = n, l = l, m = -m, name =  "whocares")
                basis_functions.append(bf)
                check.append([l, -m])
        iterator +=1
    bf_on_atoms.append(check)
    element.basis_functions = basis_functions
    return basis_functions

def gen_salcs(mol):
    (coord, masses, atoms, *garbage) = mol.to_arrays()
    elements = []
    for i in range(mol.natom()):
        elements.append(msym.Element(name=atoms[i], mass=masses[i], coordinates=coord[i]))
    

    basis_functions = []
    for z, e in enumerate(elements):
        bf = set_basis(e, z)
        for bfi in bf:
            basis_functions.append(bfi)
    aos_map = []
    idx = 0
    for element in elements:
        for bf in element.basis_functions:
            aos_map.append((idx, bf))
            idx += 1
    with msym.Context(elements=elements, basis_functions=basis_functions) as ctx:
        pg = ctx.find_symmetry()
        selements = ctx.symmetrize_elements()
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
        
        ctab = CharacterTable(pg, irreps, symops, rsymops)
        
        # Grab SALCs from libmysm
        nbfxns = len(basis_functions)
        super_irrep_block = []
        for srsidx, srs in enumerate(ctx.subrepresentation_spaces):
            salcidx = 0
            nsalcs_per_irrep = sum([len(salc.partner_functions) for salc in srs.salcs])
            irrep_block = np.zeros((nbfxns,nsalcs_per_irrep))
            for salc in srs.salcs:
                bfidxs = []
                for bf in salc.basis_functions:
                    bfidxs.append(find_bf_idx(bf, aos_map))
                for pf in salc.partner_functions:
                    irrep_block[bfidxs,salcidx] = pf
                    salcidx += 1
            super_irrep_block.append(irrep_block)
        separate_pieces = []
        for srsidx in range(len(ctx.subrepresentation_spaces)):
            #print(ctx.character_table.symmetry_species[srsidx].name, super_irrep_block[srsidx]) # Will print aotoso matrix WITH irreps
            separate_pieces.append(super_irrep_block[srsidx])
        return super_irrep_block, separate_pieces, ctab

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


if __name__ == "__main__":
    pp = pprint.PrettyPrinter(indent=4)
    np.set_printoptions(precision=10, linewidth=140, suppress=True)
    
    # Load in the required options, print out
    from input import Settings
    pp.pprint(Settings)
    
    
    molecule = psi4.geometry(Settings['molecule'])
   
    #Puream MUST be set to True for the code to work, as only the transformation of
    #pure angular momentum (spherical harmonics) under the group operations is supported
    
    basis = psi4.core.BasisSet.build(molecule, 'BASIS', Settings['basis'], puream=True)
    mints = psi4.core.MintsHelper(basis)

    wfn = psi4.core.Wavefunction.build(molecule, basis)

    molecule_basis = get_basis(molecule)

    # print out psi4's AO -> SO transformation matrix, which is the current target of this script
    #print('psi4 ao->so')
    #print(wfn.aotoso().nph)
    
    salcs, blocks, ctab = gen_salcs(molecule)
 
    #print('libmsym ao->so')
    #print(salcs)
    

