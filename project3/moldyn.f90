program mol_dyn
implicit none
integer :: Natoms,read_Natoms,i,j,nstep,l
integer, parameter :: n=3                         ! number of dimension (xyz)
integer, parameter :: unit_id=11                  ! id number for the unit input
integer, parameter :: out_id=12                   ! id number for the unit output
double precision, parameter :: dt=0.02                     ! time step (ps)
character(len=20) :: file_name,out_file
double precision :: e_tot,V,T,eps,sig,pot,kin
! declaring array 
double precision, allocatable :: xyz(:,:)
double precision, allocatable :: dist(:,:)
double precision, allocatable :: vel(:,:)
double precision, allocatable :: acc(:,:)
double precision, allocatable :: mass(:)

! check allocation status
integer :: i_stat

! preparing for reading the input file
write(*,*)"Write the name of the input file"
read(*,*)file_name
open(unit=unit_id, file=trim(file_name), status='old')

Natoms=read_Natoms(unit_id)                       ! reading atoms number

! allocation 2D array
allocate(xyz(Natoms,n) , stat=i_stat)             ! allocation coordinate array
if(i_stat /= 0) then
  print *, "MEMORY ALLOCATION FAILED"
  stop
end if
allocate(dist(Natoms,Natoms) , stat=i_stat)       ! allocation distance array
if(i_stat /= 0) then
  print *, "MEMORY ALLOCATION FAILED"
  stop
end if
allocate(vel(Natoms,n) , stat=i_stat)             ! allocation velocity array      
if(i_stat /= 0) then
  print *, "MEMORY ALLOCATION FAILED"
  stop
end if
allocate(acc(Natoms,n) , stat=i_stat)             ! allocation accelaration array      
if(i_stat /= 0) then
  print *, "MEMORY ALLOCATION FAILED"
  stop
end if
allocate(mass(Natoms) , stat=i_stat)             ! allocation mass vector
if(i_stat /= 0) then
  print *, "MEMORY ALLOCATION FAILED"
  stop
end if

read(unit_id,*)                                   ! skipping blank line
call read_mol(unit_id, Natoms, xyz, mass)

! inizializing the velocity array v0=0
do i=1,n
 do j=1,Natoms
  vel(j,i)=0.D0
 enddo
enddo

write(*,*)"Write the number of step for yout MD simulation"
read(*,*)nstep


write(*,*) "------------------------------------------------------------------------"
read(unit_id,*)eps
write(*,*)"Reading sigma for LJ potential:",eps,"kJ/mol" 
read(unit_id,*)sig
write(*,*)"Reading epsilon for LJ potential:",sig,"Angstrom"
sig=sig*10.D0  ! conversion to nm


write(*,*) "------------------------------------------------------------------------"
!write(*,*)"STARTING COORDINATES AND MASS OF EACH ATOM"
write(*, '(A6, 3A15, A18)') &
    adjustr("Atom"), &
    adjustr("x (Å)"), &
    adjustr("y (Å)"), &
    adjustr("z (Å)"), &
    adjustr("mass (g/mol)")

write(*,*) "------------------------------------------------------------------------"
do i=1,Natoms
 write(*,'(I6, 3F15.8, F14.8)')i,(xyz(i,j)/10.D0, j=1, 3),mass(i)
enddo

write(*,*) "------------------------------------------------------------------------"

call compute_distances(Natoms, xyz, dist)

kin=T(Natoms,vel,mass)
write(*,*)"INITIAL KINETIC ENERGY =",kin,"kJ/mol"
pot=V(eps,sig,Natoms,dist)
write(*,*)"INITIAL POTENTIAL ENERGY =",pot,"kJ/mol"
e_tot=kin+pot
write(*,*)"TOTAL ENERGY =",e_tot,"kJ/mol"

write(*,*) "------------------------------------------------------------------------"

out_file=trim(file_name)//'_traj.xyz'
open(unit=out_id, file=out_file) 

do j=1,nstep/10
 do i=1,10                                                                   ! progressing the system dynamic
  call verlet_update(Natoms, xyz, mass, dist, eps, sig, vel, dt)
 enddo
 write(out_id,*)Natoms
 write(out_id,*)"T=",T(Natoms,vel,mass),"kJ/mol","    V=",V(eps,sig,Natoms,dist),"kJ/mol","    Nstep=",j
 do l=1,Natoms 
  write(out_id,*)"Ar",xyz(l,1)/10.D0,xyz(l,2)/10.D0,xyz(l,3)/10.D0
 enddo 
enddo


kin=T(Natoms,vel,mass)
write(*,*)"FINAL KINETIC ENERGY =",kin,"kJ/mol"
pot=V(eps,sig,Natoms,dist)
write(*,*)"FINAL POTENTIAL ENERGY =",pot,"kJ/mol"
e_tot=kin+pot
write(*,*)"TOTAL ENERGY =",e_tot,"kJ/mol"

write(*,*) "------------------------------------------------------------------------"

! deallocate the array
deallocate(xyz)
deallocate(dist)
deallocate(vel)
deallocate(acc)
end program mol_dyn




! reading atom function
integer function read_Natoms(input_file) result(Natoms)
implicit none
integer, intent(in) :: input_file
read(input_file,*)Natoms
write(6,*)"ATOMS NUMBER =", Natoms
end function read_Natoms 


! distances calculations
subroutine compute_distances(Natoms, coord, dist)
implicit none
integer, intent(in) :: Natoms
double precision, intent(in) :: coord(Natoms,3)
double precision, intent(out) :: dist(Natoms,Natoms)
integer :: i,j
double precision :: xdist,ydist,zdist

do i=1, Natoms
 do j=i, Natoms
  xdist=coord(i,1)-coord(j,1)
  ydist=coord(i,2)-coord(j,2)
  zdist=coord(i,3)-coord(j,3)
  dist(i,j)=sqrt(xdist**2+ydist**2+zdist**2)          ! calculating the distance between two atom
  dist(j,i)=dist(i,j)                                 ! filling the symmetric elements
 enddo
enddo
end subroutine compute_distances



! reading atomic data
subroutine read_mol(input_file, n, c, m)
implicit none
integer, intent(in) :: input_file
integer :: i
integer, intent(in) :: n
doubleprecision, intent(out) :: c(n,3)
doubleprecision, intent(out) :: m(n)
do i=1,n
 read(input_file,*) c(i,1), c(i,2), c(i,3), m(i)
 c(i,1)=c(i,1)*10.D0                                    ! conversion A->nm
 c(i,2)=c(i,2)*10.D0
 c(i,3)=c(i,3)*10.D0
end do
end subroutine read_mol

double precision function V(eps, sig, n, dis)
        implicit none
        integer ::  i, t
        double precision, intent(in) :: eps, sig
        integer, intent(in) :: n
        double precision, intent(in) :: dis(n,n)
        double precision :: d,l,p
        
v=0.D0
do i=1,n
 do t=i+1,n
  d=4.D0*eps
  l=(sig/dis(i,t))**12
  p=(sig/dis(i,t))**6
  v=v+d*(l-p)
!  write(*,*)"l=",l,"p=",p,"d=",d,"v=",v,"           i=",i,"t=",t     debug line for checking value of the component of potential energy
 end do
end do
V=v
end function V


! calculating the kinetic energy
double precision function T(Natoms, vel, mass)
        implicit none
        integer, intent(in) :: Natoms
        double precision, intent(in) :: vel(Natoms,3)
        double precision, intent(in) :: mass(Natoms)
        double precision :: xvel, yvel, zvel
        integer :: i
t=0.D0
do i=1, Natoms
 xvel=vel(i,1)**2                       ! calculating the square of each component of the velocity
 yvel=vel(i,2)**2                       !        
 zvel=vel(i,3)**2                       !
 t=t+mass(i)*(xvel+yvel+zvel)           ! updating the total kinetic energy
enddo
T=0.5D0*t
end function T



! calculating the acceleration for each atom
subroutine compute_acc(Natoms, coord, mass, dist, acc, eps, sig)
        implicit none
        integer, intent(in) :: Natoms
        integer :: i,j
        double precision, intent(in) :: coord(Natoms,3)
        double precision, intent(in) :: dist(Natoms,Natoms)
        double precision, intent(in) :: mass(Natoms)
        double precision, intent(out) :: acc(Natoms,3)
        double precision :: eps,sig,ax,ay,az,rij,u,dpot

do i=1, Natoms
 ax=0.D0
 ay=0.D0
 az=0.D0
 do j=1, Natoms
  if (i==j) then
   ax=ax+0.D0
   ay=ay+0.D0
   az=az+0.D0
  else 
   u=dpot(eps, sig, dist(i,j))
   rij=dist(i,j)
   ax=ax+u*(coord(i,1)-coord(j,1))/rij
   ay=ay+u*(coord(i,2)-coord(j,2))/rij
   az=az+u*(coord(i,3)-coord(j,3))/rij
  endif
 enddo 
 acc(i,1)=-(1.D0/mass(i))*ax
 acc(i,2)=-(1.D0/mass(i))*ay
 acc(i,3)=-(1.D0/mass(i))*az
enddo




end subroutine compute_acc









! function calculating the energy potential derivate
double precision function dpot(eps, sig, r) result(u)
        implicit none
        double precision, intent(in) :: eps, sig
        double precision, intent(in) :: r
        double precision :: d,l,p

d=24.D0*eps/r
l=2.D0*(sig/r)**12
p=(sig/r)**6
u=d*(p-l)
end function dpot



subroutine verlet_update(Natoms, coord, mass, dist, eps, sig, vel, dt)
        implicit none
        integer, intent(in) :: Natoms
        integer :: i,j
        double precision :: dist(Natoms,Natoms)
        double precision, intent(in) :: mass(Natoms)
        double precision :: coord(Natoms,3)
        double precision :: acc(Natoms,3)
        double precision :: vel(Natoms,3)
        double precision :: dt, tempacc(3), eps, sig

        
do i=1,Natoms       
 call compute_acc(Natoms, coord, mass, dist, acc, eps, sig)
 do j=1,3
  tempacc(j)=acc(i,j)                                                          !save the acceleration for the step n
  coord(i,j)=coord(i,j)+vel(i,j)*dt+(tempacc(j)*dt**2) /2.D0
 enddo
 call compute_distances(Natoms, coord, dist)
 call compute_acc(Natoms, coord, mass, dist, acc, eps, sig)                  !calculating the new accelaration after the change of the posistions
 do j=1,3
  vel(i,j)=vel(i,j)+(acc(i,j)+tempacc(j))*dt/2.D0
 enddo
enddo

end subroutine verlet_update



















































































