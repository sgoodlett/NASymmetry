import numpy as np
import psi4
import json

class Molecule():
    def __init__(self, squeema) -> None:
        self.atoms = squeema["elem"]
        self.natoms = len(self.atoms)
        self.coords = np.reshape(squeema["geom"], (self.natoms,3))
        self.masses = squeema["mass"]
        self.translate([10,10,10])
        print(self.coords)
        self.com = self.find_com()
        self.translate(self.com)
        self.is_at_com()

    def find_com(self):
        com = np.zeros(3)
        for i in range(self.natoms):
            com += self.masses[i]*self.coords[i]
        return com / sum(self.masses)

    def translate(self, r):
        for i in range(self.natoms):
            self.coords[i] -= r
        
    def is_at_com(self):
        if sum(abs(self.find_com())) < 1e-8:
            return True
        else:
            return False

    def transform(self, M):
        self.coords = np.dot(M,self.coords)

    def rotation_matrix(self, axis, theta):
        cos_t = np.cos(theta)
        sin_t = np.sin(theta)
        # NOT NORMALIZING AXIS!!!
        M = np.zeros((3,3))
        M += 1 - cos_t
        for i in range(3):
            for j in range(3):
                M[i,j] *= axis[i]*axis[j]
        M += np.asarray([[cos_t, -axis[2]*sin_t, axis[1]*sin_t],[axis[2]*sin_t, cos_t, -axis[0]*sin_t],[-axis[1]*sin_t, axis[0]*sin_t, cos_t]])
        return M
    
    def reflection_matrix(self, axis):
        M = np.zeros((0,0))
        for i in range(3):
            for j in range(i,3):
                if i == j:
                    M[i,i] = 1 - 2*(axis[i]**2)
                else:
                    M[i,j] = -2 * axis[i] * axis[j]
                    M[j,i] = M[i,j]
        return M

    def inversion_matrix(self):
        M = np.zeros((3,3))
        for i in range(3):
            M[i,i] = -1
        return M

if __name__ == "__main__":
    from input import Settings
    mol = psi4.geometry(Settings["molecule"])
    beebus = mol.to_schema("psi4")
    mul = Molecule(beebus)
    print(mul.is_at_com())
