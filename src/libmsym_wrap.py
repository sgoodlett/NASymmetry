import libmsym as msym
import numpy as np
import psi4
import pprint
import scipy
from scipy.linalg import fractional_matrix_power
from scipy.linalg import block_diag
from scipy import linalg as LA

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
        nbfxns = len(basis_functions)
        super_irrep_block = []
        for srs in ctx.subrepresentation_spaces:
            salcidx = 0
            nsalcs_per_irrep = len(srs.salcs)
            irrep_block = np.zeros((nbfxns,nsalcs_per_irrep))
            for salc in srs.salcs:
                bfidxs = []
                for bf in salc.basis_functions:
                    #print(salc.basis_functions)
                    bfidxs.append(find_bf_idx(bf, aos_map))
                irrep_block[bfidxs,salcidx] = salc.partner_functions
                salcidx += 1
            super_irrep_block.append(irrep_block)
        separate_pieces = []
        for srsidx in range(len(ctx.subrepresentation_spaces)):
            #continue
            #print(ctx.character_table.symmetry_species[srsidx].name, super_irrep_block[srsidx])
            print(super_irrep_block[srsidx])
            print(type(super_irrep_block[srsidx]))
            separate_pieces.append(super_irrep_block[srsidx])
        for e in elements:
            continue
            print(f"{e.name}:{e.coordinates}")
        return super_irrep_block, separate_pieces
        #print(ctx.character_table.table)
        #print(ctx.character_table.symmetry_operations)
#        for i in range(len(super_irrep_block)):
#            if i == 0:
#                a = super_irrep_block[i]
#            else:
#                a = np.concatenate((a, super_irrep_block[i]), axis=1)
#        np.set_printoptions(formatter={"float": "{: 2.1f}".format})

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

def symm_oei(oei, salcs):
    oei_blocks = []
    for salc in salcs:
        oei_block = np.einsum('vj,uv,ui->ij', salc, oei, salc)
        #print('oei_block')
        #print(oei_block)
        #print(type(oei_block))
        oei_blocks.append(oei_block)
    blocked_oei = block_diag(oei_blocks)
    print('blocked_oei')
    print(blocked_oei)
    return blocked_oei, oei_blocks

if __name__ == "__main__":
    pp = pprint.PrettyPrinter(indent=4)
    np.set_printoptions(precision=14, linewidth=140, suppress=True)
    
    # Load in the required options, print out
    from input import Settings
    pp.pprint(Settings)
    
    
    molecule = psi4.geometry(Settings['molecule'])
    molecule.print_out() 
    #Puream MUST be set to True for the code to work, as only the transformation of
    #pure angular momentum (spherical harmonics) under the group operations is supported
    
    basis = psi4.core.BasisSet.build(molecule, 'BASIS', Settings['basis'], puream=True)

    #wfn = psi4.core.Wavefunction.build(molecule, basis)

    molecule_basis = get_basis(molecule)

    # print out psi4's AO -> SO transformation matrix, which is the current target of this script
    #print('psi4 ao->so')
    #print(wfn.aotoso().nph)
    
    salcs, blocks = gen_salcs(molecule)
    print('libmsym ao->so')
    print('salcs')
    print(salcs)
    #print('blocks')
    #print(blocks)



#    psi4.set_options({'basis': 'sto-3g', 
#                      'scf_type' : 'pk', 
#                      'e_convergence': 10e-10,
#                      'reference' : 'rhf'})
#    energy, wfn = psi4.energy('scf', return_wfn = True)
#    mints = psi4.core.MintsHelper(basis)
#    F = wfn.Fa_subset('SO')
#    F_sym = psi4.core.Matrix.to_array(F, copy=True, dense=False)
#    print(F_sym)
#    ndocc = Settings['nalpha'] #Settings['nbeta']
#    scf_max_iter = Settings['scf_max_iter']
#    symT = mints.so_kinetic()
#    symT.print_out()
#    
#    S = mints.ao_overlap()
#    T = mints.ao_kinetic()
#    V = mints.ao_potential()
#      
#    S_sym, S_sym_blocks = symm_oei(S, salcs) 
#    print('succ')
#    print(block_diag(*S_sym_blocks))
#    blocked_S = block_diag(*S_sym_blocks)
#
#    T_sym, T_sym_blocks = symm_oei(T, salcs) 
#    V_sym, V_sym_blocks = symm_oei(V, salcs) 
#    print('H_sym')
#    H_sym = T_sym + V_sym
#    print(H_sym)
#    I_sym = []
#    I = mints.ao_eri()
#    for salc in salcs:
#        I_block  = np.einsum('PQRS,Pp,Qq,Rr,Ss->pqrs',I, salc, salc, salc, salc, optimize = 'optimal')
#        I_sym.append(I_block)
#    print('I_sym') 
#    print(I_sym) 
#    print('now construct the fock matrix in the AO basis')
#    print('loop over the irreps')
#    for i, x in enumerate(S_sym_blocks):
#        T = T_sym[0][i]
#        V = V_sym[0][i]
#        S = S_sym[0][i]
#        I = I_sym[i]
#        #print('I')
#        #print(I)
#        H = T + V
#        #print('kinetic')
#        #print(T)
#        #print(i, x)
#        print('construct the orthonormalizer from the symmertrized overlap matrix blocks')
#        if S.size ==0:
#            continue
#        else:
#            A = fractional_matrix_power(S, -0.5)
#            print('number of doubly occupied orbitals')
#            print(ndocc)
#            #print('A.shape')
#            #print(A.shape)
#            #print(H.shape) 
#            #print(H_sym[i]) 
#            #construct the fock matrix
#            Ft = A.dot(H).dot(A)
#            #print('F')
#            #print(Ft)
#            _, C = np.linalg.eigh(Ft) #underscore b/c ingonring eigenvalues
#
#            C = A.dot(C)
#            #print(C)
#            
#            Cocc = C[:, :ndocc] #rows atomic, columns, molecular. only up to ndocc columns
#            print('printing c occupied')
#            print(Cocc)
#            D = np.einsum('pi,qi->pq', Cocc, Cocc)
#            convergence_criteria = 1e-10
#            En = 0
#            for i in range(1,51): #submitted assignment with range(1,50)
#                #print(I.shape)
#                #print(D.shape)
#                Ji = np.einsum('pqrs,rs->pq',I,D)
#            
#                Ki = np.einsum('prqs,rs->pq',I,D) #permute the indices here, b/c permuted in 2-e integral?
#            
#                Fn = H + 2*Ji - Ki #new fock matrix w/2-electron integrals (I)
#                #print(Fn)
#                E = np.einsum('rs,rs->',D,H) + np.einsum('rs,rs->',D,Fn) # arrow for scalar
#            
#                En_n = E + molecule.nuclear_repulsion_energy() #total energy
#                
#                Fn = A.dot(Fn).dot(A) #transform new fock matrix
#            
#                evals, C = np.linalg.eigh(Fn) #diagonalize fock matrix
#               
#                C = A.dot(C) 
#               
#                Cocc = C[:, :ndocc]
#                #Cocc = C
#                Dn = np.einsum('pi,qi->pq', Cocc, Cocc)
#                w = scipy.linalg.eigh(Dn, eigvals_only=True)
#                print('density')
#                print(w)
#                delta_E = np.format_float_scientific(np.absolute(En_n) - np.absolute(En), unique=False, precision=15) #Formats delta_E to scientific notation, takes difference of iteration energies
#                
#                delta_D = Dn - D
#                rms = np.format_float_scientific(np.linalg.norm(x=delta_D,ord='fro'))
#                
#                print('RHF iteration' + ' ' + str(i) + ':' + ' energy' + str(En_n) + ' dE ' + str(delta_E) + ' RMS ' + str(rms))#prints iteration status
#                
#                En = En_n #sets energy value for next iteration, equal to the previous iteration
#                
#                D = Dn #sets Density matrix for next iteration, equal to previous iteration
#                
#                if (float(delta_E) < convergence_criteria) and (float(rms) < convergence_criteria):
#                    print('fock')
#                    print(Fn)
#                    break             


    #for salc in salcs:
    #    print(salc)
    #    S_block = np.einsum('vj,uv,ui->ij', salc, S, salc)
    #    print(S_block)

    #symI = mints.so_eri()
    #symI.print_out()

