# Tutorial codes for molecular dynamics.
[Aleks Donev](mailto:donev@courant.nyu.edu), Courant Institute, February 2023

These codes illustrate the (positional) Verlet algorithm for Molecular Dynamics and is meant for students, not for research. The parameters are set up to start in a gas state but then a small droplet should nucleate. Increase the number of atoms for larger droplets (but you may need to wait longer).

At present direct (runtime) visualization is broken. The original code used the visual package in python2, which has now become VPython in python3 (pip3 install --user vpython). I could not get it to work; **if someone gets visualization to work please let me know**. The code will write _animnew.pdb_ which contains samples from the trajectory and can be opened with [chimera](https://www.cgl.ucsf.edu/chimera/) or [chimerax](https://www.cgl.ucsf.edu/chimerax/). I could not figure out how to make a movie though; **if someone makes a movie successfully please send it to me**.

The main file is mddemo.py, and has some options at the top that you can change. To run this, execute (once):

_f2py3 -c -m ljlib ljlib.f90_ # If you have f2py and UseF2PY=True
_gfortran -fPIC --shared ljlib_new.f90 -o ljlib_new.so_ # If UseF2PY=False

and then run with _python3 mddemo.py_. I modernized the Fortran code to use Fortran 2003 for Interoperability with C, together with [ctypes](https://docs.python.org/3/library/ctypes.html). This makes the python caller code ugly; **if someone makes this work with [PyBind11](https://pybind11.readthedocs.io/en/stable/) please let me know.** I assume the code will be easier to read but not really sure.




