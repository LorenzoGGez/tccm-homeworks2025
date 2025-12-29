program main
use sparse_mod
implicit none
integer :: ndim1, ndim2, i 
integer, allocatable :: R1(:), C1(:)
double precision, allocatable :: V1(:)
integer, allocatable :: R2(:), C2(:)
double precision, allocatable :: V2(:)
character(len=50) :: filename


! first matrix entering
write(*,*)"Insert the name of the matrix file"
read(*,*)filename
open(unit=10,file=trim(filename), status='old')
call read_matrix(10, ndim1, R1, C1, V1)
close(10)

! second matrix entering
write(*,*)"Insert the name of the matrix file"
read(*,*)filename
open(unit=11,file=trim(filename), status='old')
call read_matrix(11, ndim2, R2, C2, V2)
close(11)







deallocate(R1, C1, V1, R2, C2, V2)
end program main


