import psi4
import numpy as np
from scipy.linalg import fractional_matrix_power
from scipy.linalg import block_diag
from scipy import linalg as LA
from libmsym_wrap import LibmsymWrapper 
from input import Settings
from bdmats import BDMatrix
from rhf_helper import get_basis, build_fock_sparse, SparseERI, build_fock_sparse_sym

def aotoso(A, salcs):
    """
    AO->SO transformation for one electron integrals
    """
    B = []
    for i, irrep in enumerate(salcs):
        B.append(np.einsum("vj,uv,ui->ij", irrep, A, irrep))
    return BDMatrix(B)

def aotoso_2(ERI, salcs):
    """
    AO->SO transformation for two electron integrals
    """
    first = True
    for i, salc in enumerate(salcs):
        if first:
            s = salc
            first = False
        else:
            s = np.concatenate((s,salc), axis=1)
    E = np.einsum("PQRS,Pp,Qq,Rr,Ss->pqrs", ERI, s, s, s, s)
    return E

def get_norm(i):
    return 1/np.sqrt(i)

def normalize_S(S):
    """
    Nate's function for normalizing the antisymmetrizing matrix in the SO basis.
    Not sure why we call it normalize_S...
    Ripped straight from PSI4.
    """
    A = []
    for i, s in enumerate(S):
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
            for i in range(len(eigval)):
                Us[:,i] = U[:,i] * 1.0/np.sqrt(eigval[i])
            
            for i in range(len(eigval)):
                Us[i,:] = Us[i,:] * (normlist[i]) 
            anti = np.dot(Us, U.T)
            A.append(anti)
    return BDMatrix(A)

def order_eigval(eigval):
    """
    Sort the eigenvalues and eigenvectors from smallest to largest eigenvalue.
    Returns a list of lists containing eigenvalue and indices relating it's postion w.r.t. our irrep conventions.
    """
    B = []
    idx = 0
    for i, irrep in enumerate(eigval):
        for j, e in enumerate(irrep):
            B.append([e,i,j]) # i is the irrep index, j is the orbital index within irrep i, and idx is the index of the eigenvalue w.r.t. the full matrix
            idx += 1
    t = sorted(B, key=lambda x: x[0]) # sort based on eigenvalue
    return t

def build_D(C, eps, ndocc):
    """
    Builds the density matrix. There are two paths depending on whether D is a full matrix or a block diagonal matrix.
    """
    if isinstance(C, BDMatrix):
        T = order_eigval(eps)
        blocks = []
        for h, Cirrep in enumerate(C.blocks):
            Cidxs = []
            for t in T[:ndocc]:
                if t[1] == h:
                    Cidxs.append(t[2])
            irrep_size = np.shape(Cirrep)[0]
            Dblock = np.zeros(np.shape(Cirrep))
            for m in range(irrep_size):
                for n in range(irrep_size):
                    for cidx in Cidxs:
                        Dblock[m,n] += Cirrep[m,cidx]*Cirrep[n,cidx]
            blocks.append(Dblock)
        return BDMatrix(blocks)
        
    else:
        Cf = C.full_mat()
        Cocc = Cf[:,:ndocc]
        Df = np.einsum("pi,qi->pq", Cocc, Cocc)
        D = []
        sizes = [4,0,2,1]
        ctr = 0
        for i in sizes:
            if i == 0:
                D.append(np.empty((0,)))
            D.append(Df[ctr:ctr+i,ctr:ctr+i]) 
            ctr += i
        return BDMatrix(D)

def build_fock(H, D, ERI, nbfxns):
    """
    Build the Fock matrix given full matrices for H, D, and ERI. nbfxns is a useless parameter so I can change between this
    function and the build_fock_sparse quicker for testing
    """
    J = np.einsum("pqrs,rs->pq", ERI, D)
    K = np.einsum("prqs,rs->pq", ERI, D)
    if isinstance(H, BDMatrix):
        H = H.full_mat()
    F = H + 2*J - K
    return F

def rhf_energy(D, H, F):
    """
    Calculate HF energy
    """
    if isinstance(D, BDMatrix):
        #E = (D*(H+F)).sum()
        Df = D.full_mat()
        Hf = H.full_mat()
        Ff = F.full_mat()
        E = sum(sum(np.multiply(Df,(Hf+Ff))))
    else:
        E = sum(sum(np.multiply(D,(H+F))))
    return E

def myRHF(molecule, ints, enuc, basis):
    """
    RHF procedure.
    Starts by defining symmetry objects, then symmetrizes AO integrals.
    Then we build the initial matrices and loop over the SCF steps fifty times.
    """
    molecule_basis = get_basis(molecule, basis)
    
    # Being lazy
    nbfxns = 7
    ndocc = 5
    
    # AO integrals
    S = ints.ao_overlap().np
    T = ints.ao_kinetic().np
    V = ints.ao_potential().np
    bigERI = ints.ao_eri().np
    
    # Defining symmetry objects
    exec_libmsym = LibmsymWrapper(molecule, molecule_basis)
    exec_libmsym.run()
    salcs = exec_libmsym.salcs
    ctab = exec_libmsym.ctab
    
    # SO integrals
    S = aotoso(S, salcs)
    T = aotoso(T, salcs)
    V = aotoso(V, salcs)
    bigERI = aotoso_2(bigERI, salcs)
    
    # Building sparse ERI
    ERI = SparseERI(nbfxns)
    ERI.remove_asym(ctab)
    ERI.get_ints_dummy(bigERI)
    
    # Initializing RHF
    H = T + V
    A = normalize_S(S.blocks)
    Fs = A.transpose().dot(H.dot(A))
    eps, Cs = Fs.eigh()
    C = A.dot(Cs)
    D = build_D(C, eps, ndocc)
    for i in range(50):
        F = build_fock_sparse_sym(H, D, ERI, nbfxns)
        E = rhf_energy(D, H, F) + enuc
        Fs = A.transpose().dot(F.dot(A))
        eps, Cs = Fs.eigh()
        C = A.dot(Cs)
        D = build_D(C, eps, ndocc)
    print(E)
if __name__ == "__main__":
    np.set_printoptions(suppress=True, linewidth=120, precision=14)
    molecule = psi4.geometry(Settings['molecule'])
    molecule.update_geometry()
    ndocc = Settings['nalpha'] #Settings['nbeta']
    scf_max_iter = Settings['scf_max_iter']
    Enuc = molecule.nuclear_repulsion_energy()
    basis = psi4.core.BasisSet.build(molecule, 'BASIS', Settings['basis'])
    mints = psi4.core.MintsHelper(basis)
    psi4.set_options({'basis': 'sto-3g',
                      'scf_type': 'pk',
                      'e_convergence': 1e-10,
                     'reference': 'rhf',
                     'print' : 5,
                     'guess' : 'core'})
    #psi4.energy('scf')
    #-74.96466253910498
    E = myRHF(molecule, mints, Enuc, basis)
