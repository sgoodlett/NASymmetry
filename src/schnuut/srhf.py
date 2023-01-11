import psi4
import numpy as np
from scipy.linalg import fractional_matrix_power
from scipy.linalg import block_diag
from scipy import linalg as LA
from libmsym_wrap import LibmsymWrapper 
from bdmats import BDMatrix
from rhf_helper import get_basis, build_fock_sparse, SparseERI, build_fock_sparse_sym
import time
from copy import deepcopy

def aotoso(A, salcs):
    """
    AO->SO transformation for one electron integrals
    """
    B = []
    for i, irrep in enumerate(salcs):
        #B.append(np.einsum("vj,uv,ui->ij", irrep, A, irrep))
        #temp = np.einsum("vj,vu,ui->ji", irrep, A, irrep)
        temp1 = np.einsum('uv,ui->iv', A, irrep)
        temp  = np.einsum('iv,vj->ij', temp1, irrep)
        B.append(temp)
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
    normlists = []
    for i, s in enumerate(S):
        over = np.zeros((len(s),len(s)))
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
                    over[i,j] = s[i,j] * norm1 * norm2
            eigval, U = np.linalg.eigh(over)
            Us = deepcopy(U)
            for i in range(len(eigval)):
                Us[:,i] = U[:,i] * 1.0/np.sqrt(eigval[i])
            
            for i in range(len(eigval)):
                Us[i,:] = Us[i,:] * (normlist[i]) 
            anti = np.dot(Us, U.T)
            A.append(anti)
            normlists.append(normlist)
    print("Normlist")
    print(normlists)
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

def myRHF(molecule, ints, enuc, basis, ndocc, nbfxns, paotoso):
    """
    RHF procedure.
    Starts by defining symmetry objects, then symmetrizes AO integrals.
    Then we build the initial matrices and loop over the SCF steps fifty times.
    """
    print("Thank you for choosing our RHF code, now beginning. Please keep your hands and feet inside the ride at all times.")
    start = time.time()
    molecule_basis = get_basis(molecule, basis)
    # AO integrals
    S = ints.ao_overlap().np
    T = ints.ao_kinetic().np
    V = ints.ao_potential().np
    bigERI = ints.ao_eri().np
    now = time.time()
    print(f"AO integrals are calculated! {now-start:6.3f}")
    before = time.time()
    # Defining symmetry objects
    exec_libmsym = LibmsymWrapper(molecule, molecule_basis)
    exec_libmsym.run()
    salcs = exec_libmsym.salcs
    #salcs[0][:,[3,2]] = salcs[0][:,[2,3]]
    #salcs = paotoso
    print("Printing SALCs")
    print(salcs)
    ctab = exec_libmsym.ctab
    now = time.time()
    print(f"Point group is: {ctab.name}. {now-before:6.3f}")
    before = time.time()
    # SO integrals
    print("Printing AO Overlap")
    print(S)
    S = aotoso(S, salcs)
    print("Printing SO Overlap")
    print(S.full_mat())
    T = aotoso(T, salcs)
    V = aotoso(V, salcs)
    bigERI = aotoso_2(bigERI, salcs)
    now = time.time()
    print(f"SO integrals are formed! {now-before:6.3f}")
    before = time.time()
    # Building sparse ERI
    ERI = SparseERI(nbfxns)
    ERI.remove_asym(ctab)
    ERI.get_ints_dummy(bigERI)
    
    # Initializing RHF
    now = time.time()
    print(f"Building initial things! {now-before:6.3f}")
    before = time.time()
    H = T + V
    A = normalize_S(S.blocks)
    Fs = A.transpose().dot(H.dot(A))
    eps, Cs = Fs.eigh()
    C = A.dot(Cs)
    D = build_D(C, eps, ndocc)
    now = time.time()
    print(f"Starting iterations! {now-before:6.3f}")
    before = time.time()
    for i in range(50):
        F = build_fock_sparse_sym(H, D, ERI, nbfxns)
        E = rhf_energy(D, H, F) + enuc
        Fs = A.transpose().dot(F.dot(A))
        eps, Cs = Fs.eigh()
        C = A.dot(Cs)
        D = build_D(C, eps, ndocc)
        now = time.time()
        print(f"Step {i:2d}: RHF Energy = {E:14.10f}   {now-before:6.3f}")
        before = time.time()
    print(f"Energy is converged probably! E = {E:14.10f}. Total time = {now-start:6.3f}")
    return E

if __name__ == "__main__":
    #from input import Settings
    from ammonia import Settings
    np.set_printoptions(suppress=True, linewidth=12000, precision=7)
    molecule = psi4.geometry(Settings['molecule'])
    molecule.update_geometry()
    ndocc = Settings['nalpha'] #Settings['nbeta']
    scf_max_iter = Settings['scf_max_iter']
    Enuc = molecule.nuclear_repulsion_energy()
    basis = psi4.core.BasisSet.build(molecule, 'BASIS', Settings['basis'], puream=True)
    mints = psi4.core.MintsHelper(basis)
    psi4.set_options({'basis': Settings["basis"],
                      'scf_type': 'pk',
                      'e_convergence': 1e-10,
                     'reference': 'rhf',
                     'guess' : 'core',
                     "puream": True,
                     "print": 1},)
    ndocc = Settings["nalpha"]
    nbfxns = psi4.core.BasisSet.nbf(basis)
    pe, wfn = psi4.energy('scf', return_wfn=True)
    paotoso = wfn.aotoso().nph
    #-74.96466253910498
    #-55.44441714586791
    print("Here")
    E = myRHF(molecule, mints, Enuc, basis, ndocc, nbfxns, paotoso)
    print(f"Difference between us and PSI4: {abs(E-pe)}")
