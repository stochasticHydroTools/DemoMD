# A. Donev upgraded the code to python3
# Also made it possible to use ctypes+ISO_C_BINDING instead of f2py

#import python modules
import numpy as np
import time

#import pdb writing
import atomwrite

#=======================================    
# Some user options, see routine RunTest below for more

# whether or not to save the trajectory as a PdbFile to render in chimera(x)
SaveTrajectory = True

# distance cutoff for pairwise Lenard-Jones interactions
Cut = 3.0

# Simulate the liquid state (True) or 
# a phase transition to a drop of something glassy (False)
Liquid = False

#whether or not to use the 3d visualization
UseVisual = False # This requires vpython and I could not get it to work
if UseVisual:
    #import custom visualization library
    import atomvis

# How do we call Fortran from python?
UseF2PY = False # Use f2py/f2py3 or ISO_C_BINDING+ctypes
if UseF2PY:
    #import compiled Fortran library
    import ljlib
else:
    import ctypes as ct
    fortlib = ct.CDLL('./ljlib_new.so')

    EF = fortlib.EnergyForces # Lenard-Jones force field
    # void EnergyForces(double Pos[], double L, double rc, double *PEnergy, 
    #                   double Forces[], int Dim, int NAtom)
    EF.argtypes = [ct.POINTER(ct.c_double), ct.c_double, ct.c_double, ct.POINTER(ct.c_double), ct.POINTER(ct.c_double), ct.c_int, ct.c_int]
    
    VV = fortlib.VVIntegrate # Velocity Verlet
    # void VVIntegrate(double Pos[], double Vel[], double Accel[], 
    # double L, double rc, double dt, double *Kenergy, double *Penergy)
    VV.argtypes = [ct.POINTER(ct.c_double), ct.POINTER(ct.c_double), ct.POINTER(ct.c_double), ct.c_double, ct.c_double, ct.c_double, ct.POINTER(ct.c_double), ct.POINTER(ct.c_double)] 
   
#=======================================    

#NOTE:
#everything below assumes unit atomic masses,
#such that forces = accelerations.

def InitPositions(N, L):
    """Returns an array of initial positions of each atom,
placed on a cubic lattice for convenience.
Input:
    N: number of atoms
    L: box length
Output:
    Pos: (N,3) array of positions
"""
    #make the position array
    Pos = np.zeros((N,3), float)
    #compute integer grid # of locations for cubic lattice
    NLat = int(N**(1./3.) + 1.)
    #make an array of lattice sites
    r = L * (np.arange(NLat, dtype=float)/NLat - 0.5)
    
    #loop through x, y, z positions in lattice until done
    #for every atom in the system
    i = 0
    for x in r:
        for y in r:
            for z in r:
                Pos[i] = np.array([x,y,z], float)
                i += 1
                #if done placing atoms, return
                if i >= N:
                    return Pos
    return Pos


def RescaleVelocities(Vel, T):
    """Rescales velocities in the system to the target temperature.
Input:
    Vel: (N,3) array of atomic velocities
    T: target temperature
Output:
    Vel: same as above
"""
    #recenter to zero net momentum (assuming all masses same)
    Vel = Vel - Vel.mean(axis=0)
    #find the total kinetic energy
    KE = 0.5 * np.sum(Vel * Vel)
    #find velocity scale factor from ratios of kinetic energy
    VScale = np.sqrt(1.5 * len(Vel) * T / KE)
    Vel = Vel * VScale
    return Vel  


def InitVelocities(Pos, N, T):
    """Returns an initial random velocity set.
Input:
    Pos: (N,3) array of atomic positions
    N: number of atoms
    T: target temperature
Output:
    Vel: (N,3) array of atomic velocities
"""
    Vel = np.zeros_like(Pos) # A. Donev changed this
    # A. Donev changed this to randn instead of rand
    Vel[:,:] = np.random.randn(N, 3)
    Vel = RescaleVelocities(Vel, T)
    return Vel


def InitAccel(Pos,N,L):
    """Returns the initial acceleration array.
Input:
    Pos: (N,3) array of atomic positions
    N: number of atoms
    L: simulation box length
Output:
    Accel: (N,3) array of acceleration vectors
"""
    Accel = np.zeros_like(Pos)
    #get the acceleration from the forces
    if UseF2PY:
        PEnergy, Accel = ljlib.energyforces(Pos, L, Cut, Accel)
    else:
       # void EnergyForces(double Pos[], double L, double rc, double *PEnergy, 
       #                   double Forces[], int Dim, int NAtom)
        Pos_F   = np.zeros((N,3), float, order="F") # Fortran order
        Pos_F   = np.transpose(Pos)
        p_ptr = Pos_F.ctypes.data_as(ct.POINTER(ct.c_double))
        Accel_F = np.zeros_like(Pos_F)    
        Accel_F = np.transpose(Accel)
        a_ptr = Accel_F.ctypes.data_as(ct.POINTER(ct.c_double))
        
        PEnergy_F = ct.c_double(0.0)
        EF(p_ptr, ct.c_double(L), ct.c_double(Cut), ct.byref(PEnergy_F), a_ptr, ct.c_int(3), ct.c_int(N))
        Pos = np.transpose(Pos_F)
        Accel = np.transpose(Accel_F)
        PEnergy = np.double(PEnergy_F)
        #print("PEnergy=", PEnergy)
            
    return Accel  



def RunTest():

    if Liquid: # Liquid
       rho = 0.05 # density -- liquid
       InitTemp = 1.0 # Initial temperature
       Temp = 1.0 # Target temperature
       dt = 0.01 # timestep 0.001-0.01
    else: # Perhaps a solid, likely a glass
       rho = 0.15 # density -- solid
       InitTemp = 0.0 # Initial temperature
       Temp = 0.1 # Target temperature
       dt = 0.01 # timestep 0.001-0.01

    #set the init box width, number of particles, temperature, and timestep
    N = 125 # number of atoms = 5^3
    L = (N / rho)**(1./3.) # periodic box size
    
    print("packing density = %11.6g" % (N*(4/3*3.14159)/L**3))
 
    #set the max number of md steps; 0 for infinite loop
    MaxSteps = 10000

    # Set the frequency in md steps to write statistics to file
    WriteStatsSteps=10
    #set the frequency in md steps to rescale velocities
    #RescaleSteps = 1000 # Never rescale if <=0
    RescaleSteps = 100 # Never rescale if <=0
    #set the frequency in md steps to write coordinates to a file for visualization
    WriteConfSteps = 100 # Never if <=0

    #set the random number seed; useful for debugging
    np.random.seed = 342324

    #get the initial positions, velocities, and acceleration (forces)    
    Pos = InitPositions(N, L)
    Vel = InitVelocities(Pos, N, Temp)
    Accel = InitAccel(Pos, N, L)

    #set up the 3D display, if used    
    if UseVisual:
        atomvis.Init(Pos, L)
    
    # Initialize files
    if SaveTrajectory:
        Pdb = atomwrite.pdbfile("animnew.pdb", L)
    with open("log.dat",'w') as f: 
               f.write("# time energy potential kinetic\n")
    with open("position_first.dat",'w') as f:
       f.write("#atom 1: t x y z") 
    with open("position_last.dat",'w') as f:
       f.write("#atom N: t x y z") 
    with open("position_all.dat",'w') as f:
       # A. Donev added this to save initial and final configuration
       f.write("#atom N: x y z\n")
    
    # A. Donev: Write initial configuration
    with open("position_initial.dat",'a') as f:
        f.write("#t=%11.6g\n" % 0.0)
        for atom in range(0, N-1):                  
            f.write("%11.6g  %11.6g  %11.6g\n" % (Pos[atom,0], Pos[atom,1], Pos[atom,2]))

    #MD steps
    StartTime = time.time()
    i = 0
    # Begin time loop
    #atomwrite.write_csv_file(Pos, 1.0, "configuration", 0, L) # CSV file
    while i < MaxSteps or MaxSteps <= 0:
    
        #do one step of the integration by calling the Fortran libraries
        if UseF2PY:
            # Calling sequence: VVIntegrate(Pos, Vel, Accel, L, CutSq, dt, KEnergy, PEnergy, Dim, NAtom)
            Pos, Vel, Accel, KEnergy, PEnergy = ljlib.vvintegrate(Pos, Vel, Accel, L, Cut, dt)
        else:
        # void VVIntegrate(double Pos[], double Vel[], double Accel[], 
        # double L, double rc, double dt, double *Kenergy, double *Penergy)
            Pos_F   = np.zeros((N,3), float, order="F") # Fortran order
            Pos_F   = np.transpose(Pos)
            p_ptr = Pos_F.ctypes.data_as(ct.POINTER(ct.c_double))
            Vel_F = np.zeros_like(Pos_F)    
            Vel_F = np.transpose(Vel)
            v_ptr = Vel_F.ctypes.data_as(ct.POINTER(ct.c_double))
            Accel_F = np.zeros_like(Pos_F)    
            Accel_F = np.transpose(Accel)
            a_ptr = Accel_F.ctypes.data_as(ct.POINTER(ct.c_double))
        
            KEnergy_F = ct.c_double(0.0)
            PEnergy_F = ct.c_double(0.0)
            VV(p_ptr, v_ptr, a_ptr, ct.c_double(L), ct.c_double(Cut), ct.c_double(dt), ct.byref(KEnergy_F), ct.byref(PEnergy_F), ct.c_int(3), ct.c_int(N))
            Pos = np.transpose(Pos_F)
            Accel = np.transpose(Accel_F)  
            KEnergy = np.double(KEnergy_F)
            PEnergy = np.double(PEnergy_F)          
            #print("PEnergy=", PEnergy, " KEnergy=", KEnergy)
        
        
        #check if we need to output the positions 
        if SaveTrajectory and WriteConfSteps > 0 and i % WriteConfSteps == 0:
            Pdb.write(Pos) # PDB file
            #atomwrite.write_csv_file(Pos, 1.0, "trajectory", i, L) # CSV file
        
        # Write positions of just two atoms to visualize trajectory of an atom:
        if WriteConfSteps > 0 and i % WriteConfSteps == 0:    
            with open("position_first.dat",'a') as f:
               f.write("%11.6g %11.6g  %11.6g  %11.6g\n" % (i*dt, Pos[0,0], Pos[0,1], Pos[0,2])) 
            with open("position_last.dat",'a') as f:
               f.write("%11.6g %11.6g  %11.6g  %11.6g\n" % (i*dt, Pos[N-1,0], Pos[N-1,1], Pos[N-1,2]))
                                             
        #check if we need to update the display            
        if i % WriteStatsSteps == 0:   
            print("t,E,Ep,Ek = %11.6g  %11.6g  %11.6g  %11.6g" % (i*dt, PEnergy + KEnergy, PEnergy, KEnergy))
            with open("log.dat",'a') as f: 
               print("%11.6g  %11.6g  %11.6g  %11.6g\n" % (i*dt, PEnergy + KEnergy, PEnergy, KEnergy), file=f)
                        
        #update the 3D visualization
        if UseVisual:
            atomvis.Update(Pos, L)

        #check if we need to rescale the velocities 
        if RescaleSteps >0 and i % RescaleSteps == 0:
            Vel = RescaleVelocities(Vel, Temp)
            print("Thermostat: Rescaling velocities to T=%11.6g" % Temp)

        i += 1
    # end while loop
    #atomwrite.write_csv_file(Pos, 1.0, "configuration", MaxSteps, L) # CSV file
    
    # A. Donev: Write final configuration
    Pos = atomwrite.MinImage(Pos, L)    
    with open("position_final.dat",'a') as f:
        f.write("\n#t=%11.6g\n" % (MaxSteps*dt))
        for atom in range(0, N-1):                  
            f.write("%11.6g  %11.6g  %11.6g\n" % (Pos[atom,0], Pos[atom,1], Pos[atom,2]))
    
    if SaveTrajectory:
        Pdb.close()            

    #do one last update of the display         
    if UseVisual:
        atomvis.Update(Pos, L, Force = True)

    StopTime = time.time()
    print("Total time: %.1f s" % (StopTime - StartTime))


#check to see if we were run at the command line
if __name__ == '__main__':
    #run the test simulation
    RunTest()

    
