#import python modules
import numpy as np
import time

#import compiled Fortran library
import ljlib

#import pdb writing
import atomwrite


#whether or not to use the 3d visualization
UseVisual = False
if UseVisual:
    #import custom visualization library
    import atomvis
    
#whether or not to save the trajectory as a PdbFile
SaveTrajectory = True

#distance cutoff for pairwise interactions
Cut = 3.0


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


def InitVelocities(N, T):
    """Returns an initial random velocity set.
Input:
    N: number of atoms
    T: target temperature
Output:
    Vel: (N,3) array of atomic velocities
"""
    Vel = np.random.rand(N, 3)
    Vel = RescaleVelocities(Vel, T)
    return Vel


def InitAccel(Pos, L):
    """Returns the initial acceleration array.
Input:
    Pos: (N,3) array of atomic positions
    L: simulation box length
Output:
    Accel: (N,3) array of acceleration vectors
"""
    Accel = np.zeros_like(Pos)
    #get the acceleration from the forces
    PEnergy, Accel = ljlib.energyforces(Pos, L, Cut, Accel)
    return Accel  



def RunTest():

    if True: # Liquid
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
    
    print "packing density = %11.4f" % (N*(4/3*3.14159)/L**3)
 
    #set the max number of md steps; 0 for infinite loop
    MaxSteps = 10000

    # Set the frequency in md steps to write statistics to file
    WriteStatsSteps=100
    #set the frequency in md steps to rescale velocities
    #RescaleSteps = 1000 # Never rescale if <=0
    RescaleSteps = 0 # Never rescale if <=0
    #set the frequency in md steps to write coordinates to a file for visualization
    WriteConfSteps = 100 # Never if <=0

    #set the random number seed; useful for debugging
    np.random.seed = 342324

    #get the initial positions, velocities, and acceleration (forces)    
    Pos = InitPositions(N, L)
    Vel = InitVelocities(N, Temp)
    Accel = InitAccel(Pos, L)

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

    #MD steps
    StartTime = time.time()
    i = 0
    # Begin time loop
    #atomwrite.write_csv_file(Pos, 1.0, "configuration", 0, L) # CSV file
    while i < MaxSteps or MaxSteps <= 0:
        #do one step of the integration by calling the Fortran libraries
        # Calling sequence: VVIntegrate(Pos, Vel, Accel, L, CutSq, dt, KEnergy, PEnergy, Dim, NAtom)
        Pos, Vel, Accel, KEnergy, PEnergy = ljlib.vvintegrate(Pos, Vel, Accel, L, Cut, dt)
        
        #check if we need to rescale the velocities 
        if RescaleSteps >0 and i % RescaleSteps == 0:
            Vel = RescaleVelocities(Vel, Temp)

        #check if we need to output the positions 
        if SaveTrajectory and WriteConfSteps > 0 and i % WriteConfSteps == 0:
            Pdb.write(Pos) # PDB file
            #atomwrite.write_csv_file(Pos, 1.0, "trajectory", i, L) # CSV file
        
        # Write positions of just two atoms to visualize trajectory of an atom:
        if WriteConfSteps > 0 and i % WriteConfSteps == 0:    
            with open("position_first.dat",'a') as f:
               f.write("%11.4f %11.4f  %11.4f  %11.4f\n" % (i*dt, Pos[0,0], Pos[0,1], Pos[0,2])) 
            with open("position_last.dat",'a') as f:
               f.write("%11.4f %11.4f  %11.4f  %11.4f\n" % (i*dt, Pos[N-1,0], Pos[N-1,1], Pos[N-1,2]))
               
        #check if we need to update the display            
        if i % WriteStatsSteps == 0:   
            print "%11.4f  %11.4f  %11.4f  %11.4f" % (i*dt, PEnergy + KEnergy, PEnergy, KEnergy)
            with open("log.dat",'a') as f: 
               f.write("%11.4f  %11.4f  %11.4f  %11.4f\n" % (i*dt, PEnergy + KEnergy, PEnergy, KEnergy))
                        
        #update the 3D visualization
        if UseVisual:
            atomvis.Update(Pos, L)

        i += 1
    # end while loop
    #atomwrite.write_csv_file(Pos, 1.0, "configuration", MaxSteps, L) # CSV file
    
    if SaveTrajectory:
        Pdb.close()            

    #do one last update of the display         
    if UseVisual:
        atomvis.Update(Pos, L, Force = True)
    print "%d  %11.4f  %11.4f  %11.4f" % (i, PEnergy + KEnergy, PEnergy, KEnergy)

    StopTime = time.time()
    print "Total time: %.1f s" % (StopTime - StartTime)


#check to see if we were run at the command line
if __name__ == '__main__':
    #run the test simulation
    RunTest()

    
