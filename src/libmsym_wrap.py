import os
from character_table import CharacterTable

import libmsym as msym
import numpy as np
import pprint

"""
This class is used to generate Symmetry-Adapted Linear Combinations of Atomic Orbitals (SALCS)
for the use of symmetrizing wavefunctions in electronic structure methods.
"""

class LibmsymWrapper(object):
    def __init__(self, molecule, molecule_basis):
        self.molecule = molecule
        self.molecule_basis = molecule_basis
    
    def run(self):
        def find_bf_idx(target_bf, bf_map):
            for i,bf in bf_map:
                if target_bf == bf:
                    return i
        
        def set_basis(element, z, molecule_basis):
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
                bf = set_basis(e, z, self.molecule_basis)
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
                piece_sizes = []
                for srsidx, srs in enumerate(ctx.subrepresentation_spaces):
                    #print(ctx.character_table.symmetry_species[srsidx].name, super_irrep_block[srsidx]) # Will print aotoso matrix WITH irreps
                    separate_pieces.append(super_irrep_block[srsidx])
                    piece_sizes.append(np.shape(super_irrep_block[srsidx])[1])
                ctab = CharacterTable.from_libmsym(ctx, pg, piece_sizes)
                return super_irrep_block, ctab

        self.salcs, self.ctab = gen_salcs(self.molecule)

