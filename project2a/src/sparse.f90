program main
use sparse_mod
implicit none
integer :: ndim1, ndim2, i, j, n_mul, n, i_stat 
double precision :: y
integer, allocatable :: R1(:), C1(:)
double precision, allocatable :: V1(:)
integer, allocatable :: R2(:), C2(:)
double precision, allocatable :: V2(:), C(:,:)
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



! check the dimension of the matrices for the multiplication
if(ndim1 /= ndim2) stop "Mismatch in matrices dimension"
n=ndim1

allocate(C(n,n),  stat= i_stat)
if(i_stat /= 0) then
  print *, "MEMORY ALLOCATION FAILED"
  stop
end if


call sparse_multiplication(R1, C1, V1, R2, C2, V2, C, n_mul, n)


write(*,*)"Matrix product"
do i=1, n
 write(*,'(*(F12.6,1X))')C(i,:)
enddo

write(*,*)"Number of multiplication=",n_mul




deallocate(R1, C1, V1, R2, C2, V2, C)
end program main


