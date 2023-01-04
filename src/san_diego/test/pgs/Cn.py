import numpy as np
from numpy.linalg import matrix_power
from san_diego.molecule import Cn, reflection_matrix, inversion_matrix, Sn
from san_diego.symtext.symtext import Symel

x = np.array([1,0,0])
y = np.array([0,1,0])
z = np.array([0,0,1])

C1s = [Symel("E",np.identity(3))]
Cis = [Symel("E",np.identity(3)), Symel("i",inversion_matrix())]
Css = [Symel("E",np.identity(3)), Symel("sigma_h",reflection_matrix(z))]

C2s = [Symel("E",np.identity(3)), Symel("C_2^1",Cn(z,2))]
C3s = [Symel("E",np.identity(3)), Symel("C_3^1",Cn(z,3)),Symel("C_3^2",matrix_power(Cn(z,3),2))]
C4s = [Symel("E",np.identity(3)), Symel("C_4^1",Cn(z,4)),Symel("C_2",Cn(z,2)),Symel("C_4^3",matrix_power(Cn(z,4),3))]
C5s = [Symel("E",np.identity(3)), Symel("C_5^1",Cn(z,5)),Symel("C_5^2",matrix_power(Cn(z,5),2)),Symel("C_5^3",matrix_power(Cn(z,5),3)),Symel("C_5^4",matrix_power(Cn(z,5),4))]
C6s = [Symel("E",np.identity(3)), Symel("C_6^1",Cn(z,6)),Symel("C_3",Cn(z,3)),Symel("C_2",Cn(z,2)),Symel("C_3^2",matrix_power(Cn(z,3),2)),Symel("C_6^5",matrix_power(Cn(z,6),5))]

