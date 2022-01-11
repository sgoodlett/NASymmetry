Settings = dict()

Settings["basis"] = "sto-3g"
Settings["molecule"] = """
  0 1
  N     0.0000000       0.0000000       0.1173470  
  H     0.0000000       0.9326490       -0.2738090 
  H     0.8076980       -0.4663250      -0.2738090 
  H     -0.8076980      -0.4663250      -0.2738090 
  symmetry cs
"""
Settings["nalpha"] = 5
Settings["nbeta"] = 5
Settings["scf_max_iter"] = 50

