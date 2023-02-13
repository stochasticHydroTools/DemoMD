! A. Donev improved/modernized this code (February 2023)
! It requires Fortran 2003 or at least Fortran 95 + ISO_C_BINDING
module LJlib
    use ISO_C_BINDING ! Make code callable from C
    implicit none
   
    integer, parameter :: wp = C_DOUBLE; ! Working precision
    ! single=C_FLOAT, double=C_DOUBLE, 
    ! quad=C_FLOAT128 (in gfortran, matches _Float128 in C)

contains

! Arrays here were declared to be of size (NAtom,Dim)
! which is in order to follow the python code
! However, the efficient thing is to make the size (Dim,NAtom)
! So I have changed the code
! I also made arrays start at 1 as is common in Fortran

! Compute forces and total energy for the Lenard-Jones potential
! C prototype:
! void EnergyForces(double Pos[], double L, double rc, double *PEnergy, 
!                   double Forces[], int Dim, int NAtom)
subroutine EnergyForces(Pos, L, rc, PEnergy, Forces, Dim, NAtom) bind(C, name="EnergyForces")
    implicit none
    integer(c_int), intent(in), value :: Dim, NAtom
    real(wp), intent(in), dimension(1:Dim, 1:NAtom) :: Pos
    real(wp), intent(in), value :: L, rc 
    real(wp), intent(out) :: PEnergy
    real(wp), intent(out), dimension(1:Dim, 1:NAtom) :: Forces

    ! Local arrays
    real(wp), dimension(Dim) :: rij, Fij, Posi
    real(wp) :: d2, id2, id6, id12
    real(wp) :: rc2, Shift
    integer :: i, j
    
    !write(*,*) "L, rc, Dim, NAtom =", L, rc, Dim, NAtom
    
    PEnergy = 0.0_wp
    Forces = 0.0_wp
    Shift = -4.0_wp * (rc**(-12) - rc**(-6))
    rc2 = rc * rc
    do i = 1, NAtom
        !store Pos(i,:) in a temporary array for faster access in j loop
        Posi = Pos(:,i)
        !write(*,*) "Fortran i,r_i=",i,Posi
        do j = i + 1, NAtom ! Loop over all other particles (upper triangle)
            rij = Pos(:,j) - Posi
            ! Implement periodic BCs                
            rij = rij - L * nint(rij / L) ! Minimum image convention
            
            !compute only the squared distance and compare to squared cut
            d2 = sum(rij * rij)
            if (d2 > rc2) then 
                ! Skip pairs further than cutoff (can be done faster)
                cycle
            end if
            !write(*,*) "Fortran: d=", sqrt(d2)
           
            id2 = 1.0_wp / d2        !inverse squared distance
            id6 = id2 * id2 * id2    !inverse sixth distance
            id12 = id6 * id6         !inverse twelvth distance
            PEnergy = PEnergy + 4.0_wp * (id12 - id6) + Shift
            Fij = rij * ((-48.0_wp * id12 + 24.0_wp * id6) * id2)
            
            Forces(:,i) = Forces(:,i) + Fij
            Forces(:,j) = Forces(:,j) - Fij
        enddo
    enddo
    !write(*,*) "Fortran: PEnergy=", PEnergy
end subroutine

! Verlet integrator
! C prototype
! void VVIntegrate(double Pos[], double Vel[], double Accel[], 
!   double L, double rc, double dt, double *Kenergy, double *Penergy)
subroutine VVIntegrate(Pos, Vel, Accel, L, CutSq, dt, KEnergy, PEnergy, Dim, NAtom) bind(C,name="VVIntegrate")
    implicit none
    integer, intent(in), value :: Dim, NAtom
    real(wp), intent(in), value :: L, CutSq, dt
    real(wp), intent(inout), dimension(1:Dim, 1:NAtom) :: Pos, Vel
    real(wp), intent(out), dimension(1:Dim, 1:NAtom) :: Accel
    real(wp), intent(out) :: KEnergy, PEnergy
    
    ! Positional Verlet integrator:
    Pos = Pos + dt * Vel + 0.5 * dt*dt * Accel ! 2nd order Taylor for positions
    ! Explicit trapezoidal method for velocity:
    Vel = Vel + 0.5 * dt * Accel ! First half velocity update
    ! mass=1, so acceleration=force
    call EnergyForces(Pos, L, CutSq, PEnergy, Accel, Dim, NAtom) ! Expensive call to force field
    Vel = Vel + 0.5 * dt * Accel ! Second half velocity update
    KEnergy = 0.5 * sum(Vel*Vel)
end subroutine

end module LJlib

