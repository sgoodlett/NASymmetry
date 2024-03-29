# Because I'm really creative and can copy almost everything Gustavo did...
from copy import deepcopy
import numpy as np
from itertools import combinations_with_replacement
import psi4
from character_table import CharacterTable
from input import Settings
from libmsym_wrap import LibmsymWrapper

class SparseERI():
    """
    Object for creating and storing sparse ERI vectors.
    Generates unique indices, reduces indices by symmetry arguments, and grabs integral values from full matrix.
    """
    def __init__(self, nbfxns):
        self.nbfxns = nbfxns
        self.idxs = self.unique_idx()

    def unique_idx(self):
        d = []
        for i in range(self.nbfxns):
            for j in range(i,self.nbfxns):
                for k in range(self.nbfxns):
                    for l in range(k, self.nbfxns):
                        if idx2(i,j) < idx2(k,l):
                            continue
                        d.append((i,j,k,l))
        return d

    def remove_asym(self, ctab:CharacterTable):
        def dothing(i):
            for irrep in ctab.irreps:
                if i in irrep.orbital_idxs:
                    return irrep
        for j, idx in enumerate(self.idxs):
            irl = []
            for i in idx:
                irl.append(dothing(i))
            if ctab.dp_contains(ctab.irreps[0], irl[0], irl[1], irl[2], irl[3]):
                continue
            else:
                del self.idxs[j]

    def get_ints_dummy(self, bigERI):
        v = []
        for i,idx in enumerate(self.idxs):
            v.append(bigERI[idx])
        self.v = v

def get_basis(molecule, basis):
    num = molecule.natom()
    molecule_basis = []
    counter = 0
    for x in range(0, num):
        atom_basis = []
        molecule_basis.append(str(molecule.symbol(x)))
        for y in range(0, basis.nshell_on_center(x)):
            atom_basis.append(basis.shell(y+counter).am)
        counter += basis.nshell_on_center(x)
        molecule_basis.append(atom_basis)
    return molecule_basis
    
def idx2(i:int,j:int):
    if i < j:
        return ((j*(j+1)) >> 1) + i
    else:
        return ((i*(i+1)) >> 1) + j

def build_fock_sparse(H, D, sparse_eri, nbfxns):
    Fout = np.zeros(np.shape(H))
    Fout += H
    Ft = np.zeros(np.shape(Fout))
    for z, idx in enumerate(sparse_eri.idxs):
        i,j,k,l = idx
        v = sparse_eri.v[z]
        ij = idx2(i,j)
        kl = idx2(k,l)
        yij = i != j
        ykl = k != l
        yab = ij != kl

        Xik = 2.0 if i==k else 1.0
        Xjk = 2.0 if j==k else 1.0
        Xil = 2.0 if i==l else 1.0
        Xjl = 2.0 if j==l else 1.0

        if yij & ykl & yab:
            #J
            Ft[i,j] += 4.0*D[k,l]*v
            Ft[k,l] += 4.0*D[i,j]*v

            # K
            Ft[i,k] -= Xik*D[j,l]*v
            Ft[j,k] -= Xjk*D[i,l]*v
            Ft[i,l] -= Xil*D[j,k]*v
            Ft[j,l] -= Xjl*D[i,k]*v

        elif ykl & yab:
            # J
            Ft[i,j] += 4.0*D[k,l]*v
            Ft[k,l] += 2.0*D[i,j]*v

            # K
            Ft[i,k] -= Xik*D[j,l]*v
            Ft[i,l] -= Xil*D[j,k]*v
        elif yij & yab:
            # J
            Ft[i,j] += 2.0*D[k,l]*v
            Ft[k,l] += 4.0*D[i,j]*v

            # K
            Ft[i,k] -= Xik*D[j,l]*v
            Ft[j,k] -= Xjk*D[i,l]*v

        elif yij & ykl:
            # Only possible if i = k and j = l
            # and i < j ⇒ i < l

            # J
            Ft[i,j] += 4.0*D[k,l]*v

            # K
            Ft[i,k] -= D[j,l]*v
            Ft[i,l] -= D[j,k]*v
            Ft[j,l] -= D[i,k]*v
        elif yab:
            # J
            Ft[i,j] += 2.0*D[k,l]*v
            Ft[k,l] += 2.0*D[i,j]*v
            # K
            Ft[i,k] -= Xik*D[j,l]*v
        else:
            Ft[i,j] += 2.0*D[k,l]*v
            Ft[i,k] -= D[j,l]*v
    for i in range(nbfxns):
        Fout[i,i] += Ft[i,i]
        for j in range(i+1,nbfxns):
            Fout[i,j] += Ft[i,j] + Ft[j,i]
            Fout[j,i] = Fout[i,j]
    return Fout

def build_fock_sparse_sym(H, Dp, sparse_eri, nbfxns):
    Fout = np.zeros(np.shape(H))
    Fout = deepcopy(H)
    Ft = np.zeros((nbfxns, nbfxns))
    D = Dp.full_mat()
    for z, idx in enumerate(sparse_eri.idxs):
        i,j,k,l = idx
        v = sparse_eri.v[z]
        ij = idx2(i,j)
        kl = idx2(k,l)
        yij = i != j
        ykl = k != l
        yab = ij != kl

        Xik = 2.0 if i==k else 1.0
        Xjk = 2.0 if j==k else 1.0
        Xil = 2.0 if i==l else 1.0
        Xjl = 2.0 if j==l else 1.0

        if yij & ykl & yab:
            #J
            Ft[i,j] += 4.0*D[k,l]*v
            Ft[k,l] += 4.0*D[i,j]*v

            # K
            Ft[i,k] -= Xik*D[j,l]*v
            Ft[j,k] -= Xjk*D[i,l]*v
            Ft[i,l] -= Xil*D[j,k]*v
            Ft[j,l] -= Xjl*D[i,k]*v

        elif ykl & yab:
            # J
            Ft[i,j] += 4.0*D[k,l]*v
            Ft[k,l] += 2.0*D[i,j]*v

            # K
            Ft[i,k] -= Xik*D[j,l]*v
            Ft[i,l] -= Xil*D[j,k]*v
        elif yij & yab:
            # J
            Ft[i,j] += 2.0*D[k,l]*v
            Ft[k,l] += 4.0*D[i,j]*v

            # K
            Ft[i,k] -= Xik*D[j,l]*v
            Ft[j,k] -= Xjk*D[i,l]*v

        elif yij & ykl:
            # Only possible if i = k and j = l
            # and i < j ⇒ i < l

            # J
            Ft[i,j] += 4.0*D[k,l]*v

            # K
            Ft[i,k] -= D[j,l]*v
            Ft[i,l] -= D[j,k]*v
            Ft[j,l] -= D[i,k]*v
        elif yab:
            # J
            Ft[i,j] += 2.0*D[k,l]*v
            Ft[k,l] += 2.0*D[i,j]*v
            # K
            Ft[i,k] -= Xik*D[j,l]*v
        else:
            Ft[i,j] += 2.0*D[k,l]*v
            Ft[i,k] -= D[j,l]*v
    offset = 0
    for h, irrep in enumerate(Fout.blocks):
        hlen = np.shape(irrep)[0]
        for i in range(hlen):
            irrep[i,i] += Ft[i+offset,i+offset]
            for j in range(i+1,hlen):
                irrep[i,j] += Ft[i+offset,j+offset] + Ft[j+offset,i+offset]
                irrep[j,i] = irrep[i,j]
        offset += hlen
    return Fout

if __name__ == "__main__":
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
    molecule_basis = get_basis(molecule, basis)
    libmsym_o = LibmsymWrapper(molecule, molecule_basis)
    libmsym_o.run()
    salcs = libmsym_o.salcs
    ctab = libmsym_o.ctab
    eri = mints.ao_eri().np
    nbfxns = 7
    seri = SparseERI(nbfxns)
    #seri.remove_asym(ctab)
    seri.get_ints_dummy(eri)
    dif = [(0, 0, 0, 2), (0, 0, 0, 3), (0, 0, 0, 4), (0, 0, 1, 2), (0, 0, 1, 3), (0, 0, 1, 4), (0, 0, 2, 3), (0, 0, 2, 4), (0, 0, 3, 4), (0, 1, 0, 2), (0, 1, 0, 3), (0, 1, 0, 4), (0, 1, 1, 2), (0, 1, 1, 3), (0, 1, 1, 4), (0, 1, 2, 3), (0, 1, 2, 4), (0, 1, 3, 4), (0, 2, 0, 3), (0, 2, 0, 4), (0, 2, 1, 1), (0, 2, 1, 3), (0, 2, 1, 4), (0, 2, 2, 2), (0, 2, 2, 3), (0, 2, 2, 4), (0, 2, 3, 3), (0, 2, 3, 4), (0, 2, 4, 4), (0, 3, 0, 4), (0, 3, 1, 1), (0, 3, 1, 2), (0, 3, 1, 4), (0, 3, 2, 2), (0, 3, 2, 3), (0, 3, 2, 4), (0, 3, 3, 3), (0, 3, 3, 4), (0, 3, 4, 4), (0, 4, 1, 1), (0, 4, 1, 2), (0, 4, 1, 3), (0, 4, 2, 2), (0, 4, 2, 3), (0, 4, 2, 4), (0, 4, 3, 3), (0, 4, 3, 4), (0, 4, 4, 4), (0, 5, 1, 1), (0, 5, 1, 2), (0, 5, 1, 3), (0, 5, 1, 4), (0, 5, 2, 2), (0, 5, 2, 3), (0, 5, 2, 4), (0, 5, 3, 3), (0, 5, 3, 4), (0, 5, 4, 4), (0, 6, 1, 1), (0, 6, 1, 2), (0, 6, 1, 3), (0, 6, 1, 4), (0, 6, 1, 5), (0, 6, 2, 2), (0, 6, 2, 3), (0, 6, 2, 4), (0, 6, 2, 5), (0, 6, 3, 3), (0, 6, 3, 4), (0, 6, 3, 5), (0, 6, 4, 4), (0, 6, 4, 5), (0, 6, 5, 5), (1, 1, 1, 2), (1, 1, 1, 3), (1, 1, 1, 4), (1, 1, 2, 3), (1, 1, 2, 4), (1, 1, 3, 4), (1, 2, 1, 3), (1, 2, 1, 4), (1, 2, 2, 2), (1, 2, 2, 3), (1, 2, 2, 4), (1, 2, 3, 3), (1, 2, 3, 4), (1, 2, 4, 4), (1, 3, 1, 4), (1, 3, 2, 2), (1, 3, 2, 3), (1, 3, 2, 4), (1, 3, 3, 3), (1, 3, 3, 4), (1, 3, 4, 4), (1, 4, 2, 2), (1, 4, 2, 3), (1, 4, 2, 4), (1, 4, 3, 3), (1, 4, 3, 4), (1, 4, 4, 4), (1, 5, 2, 2), (1, 5, 2, 3), (1, 5, 2, 4), (1, 5, 3, 3), (1, 5, 3, 4), (1, 5, 4, 4), (1, 6, 2, 2), (1, 6, 2, 3), (1, 6, 2, 4), (1, 6, 2, 5), (1, 6, 3, 3), (1, 6, 3, 4), (1, 6, 3, 5), (1, 6, 4, 4), (1, 6, 4, 5), (1, 6, 5, 5), (2, 2, 2, 3), (2, 2, 2, 4), (2, 2, 3, 4), (2, 3, 2, 4), (2, 3, 3, 3), (2, 3, 3, 4), (2, 3, 4, 4), (2, 4, 3, 3), (2, 4, 3, 4), (2, 4, 4, 4), (2, 5, 3, 3), (2, 5, 3, 4), (2, 5, 4, 4), (2, 6, 3, 3), (2, 6, 3, 4), (2, 6, 3, 5), (2, 6, 4, 4), (2, 6, 4, 5), (2, 6, 5, 5), (3, 3, 3, 4), (3, 4, 4, 4), (3, 5, 4, 4), (3, 6, 4, 4), (3, 6, 4, 5), (3, 6, 5, 5), (4, 6, 5, 5)]
    mine = [(0, 0, 0, 0), (0, 0, 0, 1), (0, 0, 0, 2), (0, 0, 0, 3), (0, 0, 0, 4), (0, 0, 0, 5), (0, 0, 0, 6), (0, 0, 1, 1), (0, 0, 1, 2), (0, 0, 1, 3), (0, 0, 1, 4), (0, 0, 1, 5), (0, 0, 1, 6), (0, 0, 2, 2), (0, 0, 2, 3), (0, 0, 2, 4), (0, 0, 2, 5), (0, 0, 2, 6), (0, 0, 3, 3), (0, 0, 3, 4), (0, 0, 3, 5), (0, 0, 3, 6), (0, 0, 4, 4), (0, 0, 4, 5), (0, 0, 4, 6), (0, 0, 5, 5), (0, 0, 5, 6), (0, 0, 6, 6), (0, 1, 0, 1), (0, 1, 0, 2), (0, 1, 0, 3), (0, 1, 0, 4), (0, 1, 0, 5), (0, 1, 0, 6), (0, 1, 1, 1), (0, 1, 1, 2), (0, 1, 1, 3), (0, 1, 1, 4), (0, 1, 1, 5), (0, 1, 1, 6), (0, 1, 2, 2), (0, 1, 2, 3), (0, 1, 2, 4), (0, 1, 2, 5), (0, 1, 2, 6), (0, 1, 3, 3), (0, 1, 3, 4), (0, 1, 3, 5), (0, 1, 3, 6), (0, 1, 4, 4), (0, 1, 4, 5), (0, 1, 4, 6), (0, 1, 5, 5), (0, 1, 5, 6), (0, 1, 6, 6), (0, 2, 0, 2), (0, 2, 0, 3), (0, 2, 0, 4), (0, 2, 0, 5), (0, 2, 0, 6), (0, 2, 1, 1), (0, 2, 1, 2), (0, 2, 1, 3), (0, 2, 1, 4), (0, 2, 1, 5), (0, 2, 1, 6), (0, 2, 2, 2), (0, 2, 2, 3), (0, 2, 2, 4), (0, 2, 2, 5), (0, 2, 2, 6), (0, 2, 3, 3), (0, 2, 3, 4), (0, 2, 3, 5), (0, 2, 3, 6), (0, 2, 4, 4), (0, 2, 4, 5), (0, 2, 4, 6), (0, 2, 5, 5), (0, 2, 5, 6), (0, 2, 6, 6), (0, 3, 0, 3), (0, 3, 0, 4), (0, 3, 0, 5), (0, 3, 0, 6), (0, 3, 1, 1), (0, 3, 1, 2), (0, 3, 1, 3), (0, 3, 1, 4), (0, 3, 1, 5), (0, 3, 1, 6), (0, 3, 2, 2), (0, 3, 2, 3), (0, 3, 2, 4), (0, 3, 2, 5), (0, 3, 2, 6), (0, 3, 3, 3), (0, 3, 3, 4), (0, 3, 3, 5), (0, 3, 3, 6), (0, 3, 4, 4), (0, 3, 4, 5), (0, 3, 4, 6), (0, 3, 5, 5), (0, 3, 5, 6), (0, 3, 6, 6), (0, 4, 0, 4), (0, 4, 0, 5), (0, 4, 0, 6), (0, 4, 1, 1), (0, 4, 1, 2), (0, 4, 1, 3), (0, 4, 1, 4), (0, 4, 1, 5), (0, 4, 1, 6), (0, 4, 2, 2), (0, 4, 2, 3), (0, 4, 2, 4), (0, 4, 2, 5), (0, 4, 2, 6), (0, 4, 3, 3), (0, 4, 3, 4), (0, 4, 3, 5), (0, 4, 3, 6), (0, 4, 4, 4), (0, 4, 4, 5), (0, 4, 4, 6), (0, 4, 5, 5), (0, 4, 5, 6), (0, 4, 6, 6), (0, 5, 0, 5), (0, 5, 0, 6), (0, 5, 1, 1), (0, 5, 1, 2), (0, 5, 1, 3), (0, 5, 1, 4), (0, 5, 1, 5), (0, 5, 1, 6), (0, 5, 2, 2), (0, 5, 2, 3), (0, 5, 2, 4), (0, 5, 2, 5), (0, 5, 2, 6), (0, 5, 3, 3), (0, 5, 3, 4), (0, 5, 3, 5), (0, 5, 3, 6), (0, 5, 4, 4), (0, 5, 4, 5), (0, 5, 4, 6), (0, 5, 5, 5), (0, 5, 5, 6), (0, 5, 6, 6), (0, 6, 0, 6), (0, 6, 1, 1), (0, 6, 1, 2), (0, 6, 1, 3), (0, 6, 1, 4), (0, 6, 1, 5), (0, 6, 1, 6), (0, 6, 2, 2), (0, 6, 2, 3), (0, 6, 2, 4), (0, 6, 2, 5), (0, 6, 2, 6), (0, 6, 3, 3), (0, 6, 3, 4), (0, 6, 3, 5), (0, 6, 3, 6), (0, 6, 4, 4), (0, 6, 4, 5), (0, 6, 4, 6), (0, 6, 5, 5), (0, 6, 5, 6), (0, 6, 6, 6), (1, 1, 1, 1), (1, 1, 1, 2), (1, 1, 1, 3), (1, 1, 1, 4), (1, 1, 1, 5), (1, 1, 1, 6), (1, 1, 2, 2), (1, 1, 2, 3), (1, 1, 2, 4), (1, 1, 2, 5), (1, 1, 2, 6), (1, 1, 3, 3), (1, 1, 3, 4), (1, 1, 3, 5), (1, 1, 3, 6), (1, 1, 4, 4), (1, 1, 4, 5), (1, 1, 4, 6), (1, 1, 5, 5), (1, 1, 5, 6), (1, 1, 6, 6), (1, 2, 1, 2), (1, 2, 1, 3), (1, 2, 1, 4), (1, 2, 1, 5), (1, 2, 1, 6), (1, 2, 2, 2), (1, 2, 2, 3), (1, 2, 2, 4), (1, 2, 2, 5), (1, 2, 2, 6), (1, 2, 3, 3), (1, 2, 3, 4), (1, 2, 3, 5), (1, 2, 3, 6), (1, 2, 4, 4), (1, 2, 4, 5), (1, 2, 4, 6), (1, 2, 5, 5), (1, 2, 5, 6), (1, 2, 6, 6), (1, 3, 1, 3), (1, 3, 1, 4), (1, 3, 1, 5), (1, 3, 1, 6), (1, 3, 2, 2), (1, 3, 2, 3), (1, 3, 2, 4), (1, 3, 2, 5), (1, 3, 2, 6), (1, 3, 3, 3), (1, 3, 3, 4), (1, 3, 3, 5), (1, 3, 3, 6), (1, 3, 4, 4), (1, 3, 4, 5), (1, 3, 4, 6), (1, 3, 5, 5), (1, 3, 5, 6), (1, 3, 6, 6), (1, 4, 1, 4), (1, 4, 1, 5), (1, 4, 1, 6), (1, 4, 2, 2), (1, 4, 2, 3), (1, 4, 2, 4), (1, 4, 2, 5), (1, 4, 2, 6), (1, 4, 3, 3), (1, 4, 3, 4), (1, 4, 3, 5), (1, 4, 3, 6), (1, 4, 4, 4), (1, 4, 4, 5), (1, 4, 4, 6), (1, 4, 5, 5), (1, 4, 5, 6), (1, 4, 6, 6), (1, 5, 1, 5), (1, 5, 1, 6), (1, 5, 2, 2), (1, 5, 2, 3), (1, 5, 2, 4), (1, 5, 2, 5), (1, 5, 2, 6), (1, 5, 3, 3), (1, 5, 3, 4), (1, 5, 3, 5), (1, 5, 3, 6), (1, 5, 4, 4), (1, 5, 4, 5), (1, 5, 4, 6), (1, 5, 5, 5), (1, 5, 5, 6), (1, 5, 6, 6), (1, 6, 1, 6), (1, 6, 2, 2), (1, 6, 2, 3), (1, 6, 2, 4), (1, 6, 2, 5), (1, 6, 2, 6), (1, 6, 3, 3), (1, 6, 3, 4), (1, 6, 3, 5), (1, 6, 3, 6), (1, 6, 4, 4), (1, 6, 4, 5), (1, 6, 4, 6), (1, 6, 5, 5), (1, 6, 5, 6), (1, 6, 6, 6), (2, 2, 2, 2), (2, 2, 2, 3), (2, 2, 2, 4), (2, 2, 2, 5), (2, 2, 2, 6), (2, 2, 3, 3), (2, 2, 3, 4), (2, 2, 3, 5), (2, 2, 3, 6), (2, 2, 4, 4), (2, 2, 4, 5), (2, 2, 4, 6), (2, 2, 5, 5), (2, 2, 5, 6), (2, 2, 6, 6), (2, 3, 2, 3), (2, 3, 2, 4), (2, 3, 2, 5), (2, 3, 2, 6), (2, 3, 3, 3), (2, 3, 3, 4), (2, 3, 3, 5), (2, 3, 3, 6), (2, 3, 4, 4), (2, 3, 4, 5), (2, 3, 4, 6), (2, 3, 5, 5), (2, 3, 5, 6), (2, 3, 6, 6), (2, 4, 2, 4), (2, 4, 2, 5), (2, 4, 2, 6), (2, 4, 3, 3), (2, 4, 3, 4), (2, 4, 3, 5), (2, 4, 3, 6), (2, 4, 4, 4), (2, 4, 4, 5), (2, 4, 4, 6), (2, 4, 5, 5), (2, 4, 5, 6), (2, 4, 6, 6), (2, 5, 2, 5), (2, 5, 2, 6), (2, 5, 3, 3), (2, 5, 3, 4), (2, 5, 3, 5), (2, 5, 3, 6), (2, 5, 4, 4), (2, 5, 4, 5), (2, 5, 4, 6), (2, 5, 5, 5), (2, 5, 5, 6), (2, 5, 6, 6), (2, 6, 2, 6), (2, 6, 3, 3), (2, 6, 3, 4), (2, 6, 3, 5), (2, 6, 3, 6), (2, 6, 4, 4), (2, 6, 4, 5), (2, 6, 4, 6), (2, 6, 5, 5), (2, 6, 5, 6), (2, 6, 6, 6), (3, 3, 3, 3), (3, 3, 3, 4), (3, 3, 3, 5), (3, 3, 3, 6), (3, 3, 4, 4), (3, 3, 4, 5), (3, 3, 4, 6), (3, 3, 5, 5), (3, 3, 5, 6), (3, 3, 6, 6), (3, 4, 3, 4), (3, 4, 3, 5), (3, 4, 3, 6), (3, 4, 4, 4), (3, 4, 4, 5), (3, 4, 4, 6), (3, 4, 5, 5), (3, 4, 5, 6), (3, 4, 6, 6), (3, 5, 3, 5), (3, 5, 3, 6), (3, 5, 4, 4), (3, 5, 4, 5), (3, 5, 4, 6), (3, 5, 5, 5), (3, 5, 5, 6), (3, 5, 6, 6), (3, 6, 3, 6), (3, 6, 4, 4), (3, 6, 4, 5), (3, 6, 4, 6), (3, 6, 5, 5), (3, 6, 5, 6), (3, 6, 6, 6), (4, 4, 4, 4), (4, 4, 4, 5), (4, 4, 4, 6), (4, 4, 5, 5), (4, 4, 5, 6), (4, 4, 6, 6), (4, 5, 4, 5), (4, 5, 4, 6), (4, 5, 5, 5), (4, 5, 5, 6), (4, 5, 6, 6), (4, 6, 4, 6), (4, 6, 5, 5), (4, 6, 5, 6), (4, 6, 6, 6), (5, 5, 5, 5), (5, 5, 5, 6), (5, 5, 6, 6), (5, 6, 5, 6), (5, 6, 6, 6), (6, 6, 6, 6)]
    gus = [(0, 0, 0, 0), (0, 0, 0, 1), (0, 1, 0, 1), (0, 0, 1, 1), (0, 1, 1, 1), (1, 1, 1, 1), (0, 2, 0, 2), (0, 2, 1, 2), (1, 2, 1, 2), (0, 0, 2, 2), (0, 1, 2, 2), (1, 1, 2, 2), (2, 2, 2, 2), (0, 3, 0, 3), (0, 3, 1, 3), (1, 3, 1, 3), (2, 3, 2, 3), (0, 0, 3, 3), (0, 1, 3, 3), (1, 1, 3, 3), (2, 2, 3, 3), (3, 3, 3, 3), (0, 4, 0, 4), (0, 4, 1, 4), (1, 4, 1, 4), (2, 4, 2, 4), (3, 4, 3, 4), (0, 0, 4, 4), (0, 1, 4, 4), (1, 1, 4, 4), (2, 2, 4, 4), (3, 3, 4, 4), (4, 4, 4, 4), (0, 0, 0, 5), (0, 1, 0, 5), (1, 1, 0, 5), (0, 2, 0, 5), (1, 2, 0, 5), (2, 2, 0, 5), (0, 3, 0, 5), (1, 3, 0, 5), (2, 3, 0, 5), (3, 3, 0, 5), (0, 4, 0, 5), (1, 4, 0, 5), (2, 4, 0, 5), (3, 4, 0, 5), (4, 4, 0, 5), (0, 5, 0, 5), (0, 0, 1, 5), (0, 1, 1, 5), (1, 1, 1, 5), (0, 2, 1, 5), (1, 2, 1, 5), (2, 2, 1, 5), (0, 3, 1, 5), (1, 3, 1, 5), (2, 3, 1, 5), (3, 3, 1, 5), (0, 4, 1, 5), (1, 4, 1, 5), (2, 4, 1, 5), (3, 4, 1, 5), (4, 4, 1, 5), (0, 5, 1, 5), (1, 5, 1, 5), (0, 0, 2, 5), (0, 1, 2, 5), (1, 1, 2, 5), (0, 2, 2, 5), (1, 2, 2, 5), (2, 2, 2, 5), (0, 3, 2, 5), (1, 3, 2, 5), (2, 3, 2, 5), (3, 3, 2, 5), (0, 4, 2, 5), (1, 4, 2, 5), (2, 4, 2, 5), (3, 4, 2, 5), (4, 4, 2, 5), (0, 5, 2, 5), (1, 5, 2, 5), (2, 5, 2, 5), (0, 0, 3, 5), (0, 1, 3, 5), (1, 1, 3, 5), (0, 2, 3, 5), (1, 2, 3, 5), (2, 2, 3, 5), (0, 3, 3, 5), (1, 3, 3, 5), (2, 3, 3, 5), (3, 3, 3, 5), (0, 4, 3, 5), (1, 4, 3, 5), (2, 4, 3, 5), (3, 4, 3, 5), (4, 4, 3, 5), (0, 5, 3, 5), (1, 5, 3, 5), (2, 5, 3, 5), (3, 5, 3, 5), (0, 0, 4, 5), (0, 1, 4, 5), (1, 1, 4, 5), (0, 2, 4, 5), (1, 2, 4, 5), (2, 2, 4, 5), (0, 3, 4, 5), (1, 3, 4, 5), (2, 3, 4, 5), (3, 3, 4, 5), (0, 4, 4, 5), (1, 4, 4, 5), (2, 4, 4, 5), (3, 4, 4, 5), (4, 4, 4, 5), (0, 5, 4, 5), (1, 5, 4, 5), (2, 5, 4, 5), (3, 5, 4, 5), (4, 5, 4, 5), (0, 0, 5, 5), (0, 1, 5, 5), (1, 1, 5, 5), (0, 2, 5, 5), (1, 2, 5, 5), (2, 2, 5, 5), (0, 3, 5, 5), (1, 3, 5, 5), (2, 3, 5, 5), (3, 3, 5, 5), (0, 4, 5, 5), (1, 4, 5, 5), (2, 4, 5, 5), (3, 4, 5, 5), (4, 4, 5, 5), (0, 5, 5, 5), (1, 5, 5, 5), (2, 5, 5, 5), (3, 5, 5, 5), (4, 5, 5, 5), (5, 5, 5, 5), (0, 0, 0, 6), (0, 1, 0, 6), (1, 1, 0, 6), (0, 2, 0, 6), (1, 2, 0, 6), (2, 2, 0, 6), (0, 3, 0, 6), (1, 3, 0, 6), (2, 3, 0, 6), (3, 3, 0, 6), (0, 4, 0, 6), (1, 4, 0, 6), (2, 4, 0, 6), (3, 4, 0, 6), (4, 4, 0, 6), (0, 5, 0, 6), (1, 5, 0, 6), (2, 5, 0, 6), (3, 5, 0, 6), (4, 5, 0, 6), (5, 5, 0, 6), (0, 6, 0, 6), (0, 0, 1, 6), (0, 1, 1, 6), (1, 1, 1, 6), (0, 2, 1, 6), (1, 2, 1, 6), (2, 2, 1, 6), (0, 3, 1, 6), (1, 3, 1, 6), (2, 3, 1, 6), (3, 3, 1, 6), (0, 4, 1, 6), (1, 4, 1, 6), (2, 4, 1, 6), (3, 4, 1, 6), (4, 4, 1, 6), (0, 5, 1, 6), (1, 5, 1, 6), (2, 5, 1, 6), (3, 5, 1, 6), (4, 5, 1, 6), (5, 5, 1, 6), (0, 6, 1, 6), (1, 6, 1, 6), (0, 0, 2, 6), (0, 1, 2, 6), (1, 1, 2, 6), (0, 2, 2, 6), (1, 2, 2, 6), (2, 2, 2, 6), (0, 3, 2, 6), (1, 3, 2, 6), (2, 3, 2, 6), (3, 3, 2, 6), (0, 4, 2, 6), (1, 4, 2, 6), (2, 4, 2, 6), (3, 4, 2, 6), (4, 4, 2, 6), (0, 5, 2, 6), (1, 5, 2, 6), (2, 5, 2, 6), (3, 5, 2, 6), (4, 5, 2, 6), (5, 5, 2, 6), (0, 6, 2, 6), (1, 6, 2, 6), (2, 6, 2, 6), (0, 0, 3, 6), (0, 1, 3, 6), (1, 1, 3, 6), (0, 2, 3, 6), (1, 2, 3, 6), (2, 2, 3, 6), (0, 3, 3, 6), (1, 3, 3, 6), (2, 3, 3, 6), (3, 3, 3, 6), (0, 4, 3, 6), (1, 4, 3, 6), (2, 4, 3, 6), (3, 4, 3, 6), (4, 4, 3, 6), (0, 5, 3, 6), (1, 5, 3, 6), (2, 5, 3, 6), (3, 5, 3, 6), (4, 5, 3, 6), (5, 5, 3, 6), (0, 6, 3, 6), (1, 6, 3, 6), (2, 6, 3, 6), (3, 6, 3, 6), (0, 0, 4, 6), (0, 1, 4, 6), (1, 1, 4, 6), (0, 2, 4, 6), (1, 2, 4, 6), (2, 2, 4, 6), (0, 3, 4, 6), (1, 3, 4, 6), (2, 3, 4, 6), (3, 3, 4, 6), (0, 4, 4, 6), (1, 4, 4, 6), (2, 4, 4, 6), (3, 4, 4, 6), (4, 4, 4, 6), (0, 5, 4, 6), (1, 5, 4, 6), (2, 5, 4, 6), (3, 5, 4, 6), (4, 5, 4, 6), (5, 5, 4, 6), (0, 6, 4, 6), (1, 6, 4, 6), (2, 6, 4, 6), (3, 6, 4, 6), (4, 6, 4, 6), (0, 0, 5, 6), (0, 1, 5, 6), (1, 1, 5, 6), (0, 2, 5, 6), (1, 2, 5, 6), (2, 2, 5, 6), (0, 3, 5, 6), (1, 3, 5, 6), (2, 3, 5, 6), (3, 3, 5, 6), (0, 4, 5, 6), (1, 4, 5, 6), (2, 4, 5, 6), (3, 4, 5, 6), (4, 4, 5, 6), (0, 5, 5, 6), (1, 5, 5, 6), (2, 5, 5, 6), (3, 5, 5, 6), (4, 5, 5, 6), (5, 5, 5, 6), (0, 6, 5, 6), (1, 6, 5, 6), (2, 6, 5, 6), (3, 6, 5, 6), (4, 6, 5, 6), (5, 6, 5, 6), (0, 0, 6, 6), (0, 1, 6, 6), (1, 1, 6, 6), (0, 2, 6, 6), (1, 2, 6, 6), (2, 2, 6, 6), (0, 3, 6, 6), (1, 3, 6, 6), (2, 3, 6, 6), (3, 3, 6, 6), (0, 4, 6, 6), (1, 4, 6, 6), (2, 4, 6, 6), (3, 4, 6, 6), (4, 4, 6, 6), (0, 5, 6, 6), (1, 5, 6, 6), (2, 5, 6, 6), (3, 5, 6, 6), (4, 5, 6, 6), (5, 5, 6, 6), (0, 6, 6, 6), (1, 6, 6, 6), (2, 6, 6, 6), (3, 6, 6, 6), (4, 6, 6, 6), (5, 6, 6, 6), (6, 6, 6, 6)]
    summine = 0
    sumgus = 0
    print(seri.idxs)
    for i in seri.idxs:
        if eri[i] < 1E-14:
            summine += 1
    print(summine)
    for i in gus:
        if eri[i] < 1E-14:
            sumgus += 1
    print(sumgus)
    print(len(seri.idxs), len(gus))