from psi4.driver.molutil import geometry
import libmsym as msym
import numpy as np
from sys import argv
import psi4

def set_basis(element):
    bf_1s = msym.RealSphericalHarmonic(element=element, n=1, l=0, m=0, name="1s")
    bf_2s = msym.RealSphericalHarmonic(element=element, n=2, l=0, m=0, name="2s")
    bf_2pz = msym.RealSphericalHarmonic(element=element, n=2, l=1, m=0, name="2pz")
    bf_2px = msym.RealSphericalHarmonic(element=element, n=2, l=1, m=-1, name="2px")
    bf_2py = msym.RealSphericalHarmonic(element=element, n=2, l=1, m=1, name="2py")
    basis_functions = [bf_1s, bf_2s, bf_2pz, bf_2px, bf_2py]
    element.basis_functions = basis_functions
    return basis_functions

def gen_salcs(mol):
    (coord, masses, atoms, *garbage) = mol.to_arrays()
    elements = []
    for i in range(mol.natom()):
        elements.append(msym.Element(name = atoms[i], coordinates = coord[i]))
    basis_functions = list(np.asarray([set_basis(e) for e in elements]).flatten())
    #print(len(basis_functions))
    with msym.Context(elements=elements, basis_functions=basis_functions) as ctx:
        pg = ctx.find_symmetry()
        selements = ctx.symmetrize_elements()
        salcs = ctx.subrepresentation_spaces[0].salcs
        for srs in ctx.subrepresentation_spaces:
            salcidx = 0
            for salc in srs.salcs:
                for basis_function in salc.basis_functions:
                    print(f"{ctx.character_table.symmetry_species[srs.symmetry_species].name}:{salcidx}, {basis_function.element.name} {basis_function.n}{basis_function.l}{basis_function.m:2}\n")
                print(f"{salc.partner_functions}")
                salcidx += 1
        #print([salc.basis_functions[0].element.name for salc in salcs])

if __name__ == "__main__":
    h2o = psi4.geometry("O \nH 1 0.96\nH 1 0.96 2 104.5")
    salcs = gen_salcs(h2o)

