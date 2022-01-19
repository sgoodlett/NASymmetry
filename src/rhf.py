import psi4
import numpy as np
import pprint
import scipy
from scipy.linalg import fractional_matrix_power
from scipy.linalg import block_diag
from scipy import linalg as LA
from collections import Counter
from libmsym_wrap import LibmsymWrapper 
import collections


from input import Settings

print(Settings)

np.set_printoptions(suppress=True, linewidth=120, precision=14)


molecule = psi4.geometry(Settings['molecule'])
molecule.update_geometry()

ndocc = Settings['nalpha'] #Settings['nbeta']

scf_max_iter = Settings['scf_max_iter']

print(molecule.nuclear_repulsion_energy())

basis = psi4.core.BasisSet.build(molecule, 'BASIS', Settings['basis'])
mints = psi4.core.MintsHelper(basis)


psi4.set_options({'basis': 'sto-3g',
                  'scf_type': 'pk',
                  'e_convergence': 1e-10,
                 'reference': 'rhf',
                 'print' : 5,
                 'guess' : 'core'})


psi4.energy('scf')


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
salcs = exec_libmsym.salcs


def oei_transform(salcs, oei):
    S_oei = []
    for i, salc in enumerate(salcs):
        if len(salc) == 0:
            s_oei.appen(np.array([]))
            continue
        else:
            s_oei = np.einsum('vj,uv,ui->ij', salc, oei, salc)
            S_oei.append(s_oei)
    return S_oei

def TEI_transform(salcs, tei):
    Tei = []
    for i, salc in enumerate(salcs):
        if len(salc) == 0:
            Tei.appen(np.array([]))
            continue
        else:
            tei = np.einsum('PQRS,Pp,Qq,Rr,Ss->pqrs',I, salc, salc, salc, salc, optimize = 'optimal') 
            Tei.append(tei)
    return Tei

S = mints.ao_overlap().np
T = mints.ao_kinetic().np
V = mints.ao_potential().np
I = mints.ao_eri().np

S_sym = oei_transform(salcs, S)
T_sym = oei_transform(salcs, T)
V_sym = oei_transform(salcs, V)
I_sym = TEI_transform(salcs, I)

def build_H(T_sym, V_sym):
    H_sym = []
    for i, t in enumerate(T_sym):
        if len(t) == 0:
            H_sym.append(np.array([]))
            continue
        else:
            h = t + V_sym[i]
            H_sym.append(h)
    return H_sym


def get_norm(num):
    norm = (1/np.sqrt(num))
    return norm

def normalize_overlap(S_sym):
#    NS_sym = []
    A = []
    for i, s in enumerate(S_sym):
        if len(s) == 0:
#            NS_sym.append(np.array([]))
            A.append(np.array([]))
            continue
        else:
            normlist = []
            for i in range(len(s)):
                norm1 = get_norm(s[i,i])
                normlist.append(norm1)
                for j in range(len(s)):
                    norm2 = get_norm(s[j,j])
                    s[i,j] = s[i,j] * norm1 * norm2
            eigval, U = np.linalg.eigh(s)
            Us = U.copy()
            print('eigvecs')
            print(U)
            for i in range(len(eigval)):
                Us[:,i] = U[:,i] * 1.0/np.sqrt(eigval[i])
            print('half transformed?')
            print(Us)
            
            for i in range(len(eigval)):
                Us[i,:] = Us[i,:] * (normlist[i]) 
            anti = np.dot(Us, U.T)
            A.append(anti)
    return A

def init_fock(A, H_sym):
    Fn = []
    C = []
    Energy = []
    for i, a in enumerate(A):
        if len(a) == 0:
            Fn.append(np.array([]))
            C.append(np.array([]))
            Energy.append(np.array([]))
            continue
        else:
            fn = (a.T).dot(H_sym[i]).dot(a)
            energy, c = np.linalg.eigh(fn)
            Fn.append(fn)
            coeff = a.dot(c)
            C.append(coeff)
            Energy.append(energy)
    
    return Fn, C, Energy

H_sym = build_H(T_sym, V_sym)
A = normalize_overlap(S_sym)
print('antisymm with normalize')
print(A)

#A = antisymmetrize(NS_sym)



Fn, C, Energy = init_fock(A, H_sym)
print('init fock')
print(Fn)
print('coefficients')
print(C)



energies = []
aufbau = []
for x, e in enumerate(Energy):
    for y in e:
        aufbau.append(x)
        energies.append(y)
energies = np.array(energies)

aufbau = np.array(aufbau)
aufbau = aufbau[np.argsort(energies)]

energies = energies[np.argsort(energies)]
aufbau = aufbau[:ndocc]
collect = collections.Counter(aufbau)

def occupy(C, collect):
    Cocc = []
    D = []
    for i, c in enumerate(C):
        if len(c) == 0:
            Cocc.append(np.array([]))
            D.append(np.array([]))
            continue
        else:
            irr_occ = collect[i]
            cocc = c[:, :irr_occ]
            d = np.einsum('pi,qi->pq', cocc, cocc)
            Cocc.append(cocc)
            D.append(d)

    return Cocc, D

Cocc, D = occupy(C, collect)

print('init density')
print(D)




#def antisymmetrize(NS_sym):
#    A = []
#    for i, n in enumerate(NS_sym):
#        if len(n) == 0:
#            A.append(np.array([]))
#            continue
#        else:
#            #a = fractional_matrix_power(n, -0.5)
#            #A.append(a)
#            
#            tmp = n
#            u, v = np.linalg.eigh(n)
#            print('eigvecs')
#            print(v)
#            #print('anti?')
#            for i in range(len(v)):
#                #print('diag of eigvec?')
#                v[:,i] = v[:,i] * 1.0/np.sqrt(u[i])
#            print('half transformed?')
#            print(v)
#                #norm1 = get_norm(v[i,i])
#                #print('half tranformed?')
#                #print(v[:,i]* norm1)    
#            #norm3 = get_norm
#            #v[i,:] = v[i,:] * norm3
#            anti = np.dot(v, v.T)
#            print('the real anti?')
#            print(anti)
#            A.append(anti)
#    return A
