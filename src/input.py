#Settings = dict()
#
#Settings["basis"] = "sto-3g"
#Settings["molecule"] = """
#  0 1
#  O
#  H 1 R
#  H 1 R 2 A
#  R = 1.0
#  A = 104.5
#  symmetry c2v
#"""
#Settings["nalpha"] = 5
#Settings["nbeta"] = 5
#Settings["scf_max_iter"] = 50

Settings = dict()

Settings["basis"] = "sto-3g"
Settings["molecule"] = """
  0 1
  N     0.0000000       0.0000000       0.1173470  
  H     0.0000000       0.9326490       -0.2738090 
  H     0.8076980       -0.4663250      -0.2738090 
  H     -0.8076980      -0.4663250      -0.2738090 
"""
Settings["nalpha"] = 5
Settings["nbeta"] = 5
Settings["scf_max_iter"] = 50

