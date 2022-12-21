import numpy as np
from math import isclose
from san_diego.molecule import *

class RotationElement():
    def __init__(self, axis, order):
        self.axis = axis
        self.order = order
    def __eq__(self, other):
        if isinstance(other, RotationElement):
            return issame_axis(self.axis, other.axis) and self.order == other.order

def normalize(a):
    n = np.linalg.norm(a)
    if n <= tol:
        return None
    return a / np.linalg.norm(a)

def issame_axis(a, b):
    A = normalize(a)
    B = normalize(b)
    d = abs(np.dot(A,B))
    return isclose(d, 1.0, abs_tol=tol)

def intersect(a, b):
    out = []
    for re_a in a:
        for re_b in b:
            if re_a == re_b:
                out.append(re_a)
                break
    return out

def rotation_set_intersection(rotation_set):
    out = rotation_set[0]
    if len(rotation_set) > 1:
        for i in range(len(rotation_set)):
            out = intersect(out, rotation_set[i])
    return out

def find_rotation_sets(mol, SEAs):
    out_all_SEAs = []
    for sea in SEAs:
        length = len(sea.subset)
        out_per_SEA = []
        if length < 2:
            sea.label = "Single Atom"
        elif length == 2:
            sea.label = "Linear"
            sea.axis = normalize(mol[sea.subset[0]].coords)
        else:
            sea_mol = mol[sea.subset]
            sea_mol.translate(sea_mol.find_com())
            evals, evecs = np.linalg.eigh(calcmoit(mol[sea.subset]))
            idx = evals.argsort()
            Ia, Ib, Ic = evals[idx]
            Iav, Ibv, Icv = [evecs[:,i] for i in idx]
            if np.isclose(Ia, Ib, atol=tol) and np.isclose(Ia, Ic, atol=tol):
                sea.label = "Spherical"
            elif np.isclose(Ia+Ib, Ic, atol=tol):
                axis = Icv
                sea.axis = axis
                if np.isclose(Ia, Ib, atol=tol):
                    sea.label = "Regular Polygon"
                    for i in range(2,length+1):
                        if isfactor(length,i):
                            re = RotationElement(axis,i)
                            out_per_SEA.append(re)
                else:
                    sea.label = "Irregular Polygon"
                    for i in range(2,length):
                        if isfactor(length, i):
                            re = RotationElement(axis, i)
                            out_per_SEA.append(re)
            else:
                if not (np.isclose(Ia, Ib, atol=tol) or np.isclose(Ib, Ic, atol=tol)):
                    sea.label = "Asymmetric Rotor"
                    for i in [Iav, Ibv, Icv]:
                        re = RotationElement(i, 2)
                        out_per_SEA.append(re)
                else:
                    if np.isclose(Ia, Ib, atol=tol):
                        sea.label = "Oblate Symmetric Top"
                        axis = Icv
                        sea.axis = Icv
                    else:
                        sea.label = "Prolate Symmetric Top"
                        axis = Iav
                        sea.axis = Iav
                    k = length//2
                    for i in range(2,k+1):
                        if isfactor(k,i):
                            re = RotationElement(axis, i)
                            out_per_SEA.append(re)
            if len(out_per_SEA) > 0:
                out_all_SEAs.append(out_per_SEA)
    return out_all_SEAs

def isfactor(n,a):
    if n % a == 0:
        return True
    else:
        return False

def find_rotations(mol, rotation_set):
    if len(rotation_set) < 1:
        return []
    molmoit = calcmoit(mol)
    evals = np.sort(np.linalg.eigh(molmoit)[0])
    if evals[0] == 0.0 and np.isclose(evals[1], evals[2], atol=tol):
        axis = normalize(mol.coords[0,:])
        re = rotation_matrix(axis, 0)
        return [re]
    rsi = rotation_set_intersection(rotation_set)
    out = []
    for i in rsi:
        rmat = Cn(i.axis, i.order)
        molB = mol.transform(rmat)
        if isequivalent(mol, molB):
            out.append(i)
    return out

def find_a_c2(mol, SEAs):
    for sea in SEAs:
        a = c2a(mol, sea)
        if a is not None:
            return a
        else:
            b = c2b(mol, sea)
            if b is not None:
                return b
            else:
                if sea.label == "Linear":
                    for sea2 in SEAs:
                        if sea == sea2:
                            continue
                        elif sea2.label == "Linear":
                            c = c2c(mol, sea, sea2)
                            if c is not None:
                                return c
    return None

def is_there_ortho_c2(mol, SEAs, paxis):
    for sea in SEAs:
        b = c2b(mol, sea, axis=paxis)
        if b is not None:
            return True, b
        else:
            a = c2a(mol, sea, axis=paxis)
            if a is not None:
                return True, a
            else:
                if sea.label == "Linear":
                    for sea2 in SEAs:
                        if sea == sea2:
                            continue
                        elif sea2.label == "Linear":
                            c = c2c(mol, sea, sea2, axis=paxis)
                            if c is not None:
                                return True, c
    return False, None

def num_C2(mol, SEAs):
    axes = []
    for sea in SEAs:
        a = c2a(mol, sea, all=True)
        if a is not None:
            for i in a:
                axes.append(i)
        b = c2b(mol, sea, all=True)
        if b is not None:
            for i in b:
                axes.append(i)
    if len(axes) < 1:
        return None
    unique_axes = [axes[0]]
    for i in axes:
        check = True
        for j in unique_axes:
            if issame_axis(i,j):
                check = False
                break
        if check:
            unique_axes.append(i)
    return len(unique_axes), unique_axes

def c2a(mol, sea, axis=None, all=False):
    length = len(sea.subset)
    out = []
    for i in range(length):
        for j in range(i+1,length):
            midpoint = mol.coords[sea.subset[i],:] + mol.coords[sea.subset[j],:]
            if np.isclose(midpoint, [0,0,0], atol=tol).all():
                continue
            else:
                midpoint = normalize(midpoint)
                if axis is not None and issame_axis(midpoint, axis):
                    continue
                c2 = Cn(midpoint, 2)
                molB = mol.transform(c2)
                if isequivalent(mol, molB):
                    if all:
                        out.append(midpoint)
                    else:
                        return midpoint
    if len(out) < 1:
        return None
    return out

def c2b(mol, sea, axis=None, all=False):
    length = len(sea.subset)
    out = []
    for i in range(length):
        c2_axis = normalize(mol.coords[sea.subset[i],:])
        if axis is not None and issame_axis(c2_axis, axis):
            continue
        c2 = Cn(c2_axis, 2)
        molB = mol.transform(c2)
        if isequivalent(mol, molB):
            if all:
                out.append(c2_axis)
            else:
                return c2_axis
    if len(out) < 1:
        return None
    return out

def c2c(mol, sea1, sea2, axis=None):
    rij = mol.coords[sea1.subset[0],:] - mol.coords[sea1.subset[1],:]
    rkl = mol.coords[sea2.subset[0],:] - mol.coords[sea2.subset[1],:]
    c2_axis = normalize(np.cross(rij, rkl))
    if axis is not None and issame_axis(c2_axis, axis):
        return None
    c2 = Cn(c2_axis,2)
    molB = mol.transform(c2)
    if isequivalent(mol, molB):
        return c2_axis
    return None

def highest_order_axis(rotations):
    ns = []
    for i in range(len(rotations)):
        ns.append(rotations[i].order)
    return np.sort(ns)[-1]

def is_there_sigmah(mol, paxis):
    sigmah = reflection_matrix(paxis)
    molB = mol.transform(sigmah)
    return isequivalent(mol, molB)

def is_there_sigmav(mol, SEAs, paxis):
    axes = []
    for sea in SEAs:
        length = len(sea.subset)
        if length < 2:
            continue
        A = sea.subset[0]
        for i in range(1,length):
            B = sea.subset[i]
            #n = normalize(mol[A].xyz - mol[B].xyz)
            n = normalize(mol.coords[A,:] - mol.coords[B,:])
            sigma = reflection_matrix(n)
            molB = mol.transform(sigma)
            if isequivalent(mol, molB):
                axes.append(n)
            else:
                continue
    if len(axes) < 1:
        if mol_is_planar(mol):
            return True, planar_mol_axis(mol)
        else:
            return False, None
    unique_axes = [axes[0]]
    for i in axes:
        check = True
        for j in unique_axes:
            if issame_axis(i,j):
                check = False
                break
        if check:
            unique_axes.append(i)
    for i in unique_axes:
        if issame_axis(i, paxis):
            continue
        else:
            return True, i
    return False, None

def mol_is_planar(mol):
    rank = np.linalg.matrix_rank(mol.coords)
    if rank < 3:
        return True
    return False

def planar_mol_axis(mol):
    for i in range(mol.natoms):
        for j in range(i,mol.natoms):
            a = normalize(mol.coords[i,:])
            b = normalize(mol.coords[j,:])
            chk = np.dot(a,b)
            if not np.isclose(chk, 1.0, atol=tol):
                return normalize(np.cross(a,b))
    return None