# Tutorial codes for molecular dynamics.
[Aleks Donev](mailto:donev@courant.nyu.edu), Courant Institute, February 2023

These codes illustrate the (positional) Verlet algorithm for Molecular Dynamics and is meant for students, not for research. The parameters are set up to start in a gas state but then a small droplet should nucleate. Increase the number of atoms for larger droplets (but you may need to wait longer).

The main file is mddemo.py, and has some options at the top that you can change. To run this, execute (once):

_f2py3 -c -m ljlib ljlib.f90_ # If you have f2py and UseF2PY=True

_gfortran -fPIC --shared ljlib_new.f90 -o ljlib_new.so_ # If UseF2PY=False

and then run with _python3 mddemo.py_. In _ljlib_new.f90_ I modernized the Fortran code to use Fortran 2003 for Interoperability with C, together with [ctypes](https://docs.python.org/3/library/ctypes.html); the older code _ljlib.f90_ uses Fortran 90 and [_f2py_](https://numpy.org/doc/stable/f2py/). This makes the python caller code ugly. **If someone makes this work with [PyBind11](https://pybind11.readthedocs.io/en/stable/) please let me know**.

At present direct (runtime) visualization is broken. The original code used the visual package in python2, which has now become VPython in python3 (_pip3 install --user vpython_). I could not get it to work; **if someone gets visualization to work please let me know**. The code will write _animnew.pdb_ which contains samples from the trajectory and can be opened with [chimerax](https://www.cgl.ucsf.edu/chimerax/). There is a demo movie in _Movies/DemoMD.mp4_; here are instructions on how to make a movie using _chimerax_ commands:

Select the framerate both here and the last command:

open animnew.pdb coordset true

graphics rate maxFrameRate 2

movie record format png

To start playing and recording the movie:

coordset #1; movie record

To save the movie up to the current frame:

movie encode framerate 2 quality higher output Movies/DemoMD.mp4

