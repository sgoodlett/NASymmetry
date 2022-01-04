# NASymmetry
Non-Abelian Symmetry

A rough start to implementing non-abelian symmetry into electronic structure code.
Make sure to have Psi4 and Psi4numpy as importable modules?

# Setup
```
  git clone https://github.com/sgoodlett/NASymmetry.git
  cd NASymmetry
  cd libmsym
  mkdir build
  cd build
  cmake -DBUILD_SHARED_LIBS:BOOL=ON -DMSYM_BUILD_EXAMPLES:BOOL=ON ../.
  make
  sudo make install
  cd ../bindings/python
  python setup.py install --user
  cd ../../../
  ./exp_var.sh
```
