import numpy as np
from numpy.linalg import matrix_power
from san_diego.symtext.symtext import Symel
from san_diego.molecule import Cn, reflection_matrix, inversion_matrix, Sn

x = np.array([1,0,0])
y = np.array([0,1,0])
z = np.array([0,0,1])

C2hs = [Symel("E",np.identity(3)), Symel("sigma_h",reflection_matrix(z)), Symel("i",inversion_matrix()), Symel("C_2^1",Cn(z,2))]
C3hs = [Symel("E",np.identity(3)), Symel("sigma_h",reflection_matrix(z)),
    Symel("C_3^1",Cn(z,3)), Symel("C_3^2",matrix_power(Cn(z,3),2)), Symel("S_3^1",Cn(z,3)), Symel("S_3^5",matrix_power(Sn(z,3),5))]
C4hs = [Symel("E",np.identity(3)), Symel("sigma_h",reflection_matrix(z)), Symel("i",inversion_matrix()),
    Symel("C_4^1",Cn(z,4)), Symel("C_2^1",Cn(z,2)), Symel("C_4^3",matrix_power(Cn(z,4),3)), Symel("S_4^1",Sn(z,4)), Symel("S_4^3",matrix_power(Sn(z,4),3))]
C5hs = [Symel("E",np.identity(3)), Symel("sigma_h",reflection_matrix(z)),
    Symel("C_5^1",Cn(z,5)), Symel("C_5^2",matrix_power(Cn(z,5),2)), Symel("C_5^3",matrix_power(Cn(z,5),3)), Symel("C_5^4",matrix_power(Cn(z,5),4)),
    Symel("S_5^1",Sn(z,5)), Symel("S_5^7",matrix_power(Sn(z,5),7)), Symel("S_5^3",matrix_power(Sn(z,5),3)), Symel("S_5^9",matrix_power(Sn(z,5),9))]
C6hs = [Symel("E",np.identity(3)), Symel("sigma_h",reflection_matrix(z)), Symel("i",inversion_matrix()),
    Symel("C_6^1",Cn(z,6)), Symel("C_3^1",Cn(z,3)), Symel("C_2^1",Cn(z,2)), Symel("C_3^2",matrix_power(Cn(z,3),2)), Symel("C_6^5",matrix_power(Cn(z,6),5)),
    Symel("S_6^1",Sn(z,6)), Symel("S_3^1",Sn(z,3)), Symel("S_3^2",matrix_power(Sn(z,3),2)), Symel("S_6^5",matrix_power(Sn(z,6),5))]
