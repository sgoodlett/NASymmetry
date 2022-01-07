# NASymmetry
Non-Abelian Symmetry

A rough start to implementing non-abelian symmetry into electronic structure code.
Make sure to have Psi4!

# Setup
These are the steps to build libmsym and its Python bindings and also add the Python bindings to PYTHONPATH. Most of this is taken from the libmsym page. Hopefully I am not violating any licenses.
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
