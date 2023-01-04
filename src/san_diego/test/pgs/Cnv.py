import numpy as np
from numpy.linalg import matrix_power
from san_diego.molecule import Cn, reflection_matrix, inversion_matrix, Sn
from san_diego.symtext.symtext import Symel

x = np.array([1,0,0])
y = np.array([0,1,0])
z = np.array([0,0,1])

C2s = [Symel("E",np.identity(3)), Symel("C_2^1",Cn(z,2)), Symel("sigma_v(1)",reflection_matrix(y)), Symel("sigma_d(1)",reflection_matrix(x))]
C3s = [Symel("E",np.identity(3)), Symel("C_3^1",Cn(z,3)),Symel("C_3^2",matrix_power(Cn(z,3),2)),
    Symel("sigma_v(1)",reflection_matrix(y)), Symel("sigma_v(2)",reflection_matrix(np.dot(y,Cn(z,3).transpose()))), Symel("sigma_v(3)",reflection_matrix(np.dot(y,matrix_power(Cn(z,3),2).transpose())))]
C4s = [Symel("E",np.identity(3)), Symel("C_4^1",Cn(z,4)),Symel("C_2",Cn(z,2)),Symel("C_4^3",matrix_power(Cn(z,4),3))]
C5s = [Symel("E",np.identity(3)), Symel("C_5^1",Cn(z,5)),Symel("C_5^2",matrix_power(Cn(z,5),2)),Symel("C_5^3",matrix_power(Cn(z,5),3)),Symel("C_5^4",matrix_power(Cn(z,5),4))]
C6s = [Symel("E",np.identity(3)), Symel("C_6^1",Cn(z,6)),Symel("C_3",Cn(z,3)),Symel("C_2",Cn(z,2)),Symel("C_3^2",matrix_power(Cn(z,3),2)),Symel("C_6^5",matrix_power(Cn(z,6),5))]

