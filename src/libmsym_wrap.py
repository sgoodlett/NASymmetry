import libmsym as msym
import numpy as np
import psi4

def find_bf_idx(target_bf, bf_map):
    for i,bf in bf_map:
        if target_bf == bf:
            return i

def set_basis(element):
    if element.name == "O":
        bf_1s = msym.RealSphericalHarmonic(element=element, n=1, l=0, m=0, name="1s")
        bf_2s = msym.RealSphericalHarmonic(element=element, n=2, l=0, m=0, name="2s")
        bf_2py = msym.RealSphericalHarmonic(element=element, n=2, l=1, m=0, name="2pzonder")
        bf_2pz = msym.RealSphericalHarmonic(element=element, n=2, l=1, m=-1, name="2pxonder")
        bf_2px = msym.RealSphericalHarmonic(element=element, n=2, l=1, m=1, name="2pyonder")
        basis_functions = [bf_1s, bf_2s, bf_2pz, bf_2px, bf_2py]
        element.basis_functions = basis_functions
        return basis_functions
    else:
        bf_1s = msym.RealSphericalHarmonic(element=element, n=1, l=0, m=0, name="1s")
        basis_functions = [bf_1s]
        element.basis_functions = basis_functions
        return basis_functions

def gen_salcs(mol):
    (coord, masses, atoms, *garbage) = mol.to_arrays()
    elements = []
    for i in range(mol.natom()):
        elements.append(msym.Element(name=atoms[i], mass=masses[i], coordinates=coord[i]))
    basis_functions = []
    for e in elements:
        bf = set_basis(e)
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
        for srs in ctx.subrepresentation_spaces:
            salcidx = 0
            nsalcs_per_irrep = len(srs.salcs)
            irrep_block = np.zeros((nbfxns,nsalcs_per_irrep))
            for salc in srs.salcs:
                bfidxs = []
                for bf in salc.basis_functions:
                    bfidxs.append(find_bf_idx(bf, aos_map))
                irrep_block[bfidxs,salcidx] = salc.partner_functions
                salcidx += 1
            super_irrep_block.append(irrep_block)
        for srsidx in range(len(ctx.subrepresentation_spaces)):
            continue
            print(ctx.character_table.symmetry_species[srsidx].name, super_irrep_block[srsidx])
        for e in elements:
            continue
            print(f"{e.name}:{e.coordinates}")
        return super_irrep_block
        #print(ctx.character_table.table)
        #print(ctx.character_table.symmetry_operations)
#        for i in range(len(super_irrep_block)):
#            if i == 0:
#                a = super_irrep_block[i]
#            else:
#                a = np.concatenate((a, super_irrep_block[i]), axis=1)
#        np.set_printoptions(formatter={"float": "{: 2.1f}".format})




if __name__ == "__main__":
    h2o = psi4.geometry("O \nH 1 0.96\nH 1 0.96 2 104.5")
    #psi4.set_options({"basis":"sto-3g"})
    #print(h2o.basis_on_atom(0))
    wfn = psi4.core.Wavefunction.build(h2o, "sto-3g")
    #print(h2o.to_arrays())
    print(wfn.aotoso().nph)
    salcs = gen_salcs(h2o)
    print(type(salcs[0][0][0]))

