
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

