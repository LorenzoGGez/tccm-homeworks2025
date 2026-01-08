program main
use sparse_mod
implicit none
integer :: ndim1, ndim2, i, j, n_mul, n, i_stat, maxiter, option
double precision :: y
integer, allocatable :: R1(:), C1(:)
double precision, allocatable :: V1(:)
integer, allocatable :: R2(:), C2(:)
double precision, allocatable :: V2(:), C(:,:), mat1(:,:), mat2(:,:)
double precision  :: t_start, t_end
character(len=50) :: filename


! first matrix entering
write(*,*)"Insert the name of the matrix file"
read(*,*)filename
open(unit=10,file=trim(filename), status='old')
call read_matrix(10, ndim1, R1, C1, V1)
close(10)
write(*,*)

! second matrix entering
write(*,*)"Insert the name of the matrix file"
read(*,*)filename
open(unit=11,file=trim(filename), status='old')
call read_matrix(11, ndim2, R2, C2, V2)
close(11)
write(*,*)

! check the dimension of the matrices for the multiplication
if(ndim1 /= ndim2) stop "Mismatch in matrices dimension"
n=ndim1

! chosing the maximum amout of iteratio of the matrix multiplication
write(*,*)"Insert the number of iteration for the matrix multiplication"
read(*,*)maxiter
write(*,*)
! chosing if print the matrix
write(*,*)"Do you want to print the matrices? (1 for yes, 0 for no)"
read(*,*)option
write(*,*)


! allocating memory to the new matrix product
allocate(C(n,n),  stat= i_stat)
if(i_stat /= 0) then
  print *, "MEMORY ALLOCATION FAILED"
  stop
end if

write(*,*)"---------------------------------------------------------------------------------"
write(*,*)"Matrix product from sparse matrix"
write(*,*)"---------------------------------------------------------------------------------"

call cpu_time(t_start)
do i=1, maxiter
 call sparse_multiplication(R1, C1, V1, R2, C2, V2, C, n_mul, n)
enddo
call cpu_time(t_end)

if (option == 1) then
 do i=1, n
  write(*,'(*(F12.6,1X))')C(i,:)
 enddo
 write(*,*)"---------------------------------------------------------------------------------"
endif

write(*,'(A, 1X, F14.8)')"Time for the operation:",t_end-t_start
write(*,*)"---------------------------------------------------------------------------------"
write(*,*)"Number of multiplication:",n_mul
write(*,*)"---------------------------------------------------------------------------------"
write(*,*)


! FIRST MATRIX
! allocating the non-sparse format matrix
allocate(mat1(ndim1,ndim1),  stat= i_stat)
if(i_stat /= 0) then
  print *, "MEMORY ALLOCATION FAILED"
  stop
end if

! turning spare matrix in non-sparse format
call desparsification(R1, C1, V1, mat1, ndim1)

write(*,*)"---------------------------------------------------------------------------------"
write(*,*)"Transforming the first matrix in the non-sparse format"
write(*,*)"---------------------------------------------------------------------------------"

if (option == 1) then
 do i=1, ndim1
  write(*,'(*(F12.6,1X))')mat1(i,:)                    
 enddo
 write(*,*)"---------------------------------------------------------------------------------"
endif

! SECOND MATRIX
! allocating the non-sparse format matrix
allocate(mat2(ndim2,ndim2),  stat= i_stat)
if(i_stat /= 0) then
  print *, "MEMORY ALLOCATION FAILED"
  stop
end if

! turning spare matrix in non-sparse format
call desparsification(R2, C2, V2, mat2, ndim2)

write(*,*)
write(*,*)"---------------------------------------------------------------------------------"
write(*,*)"Transforming the second matrix in the non-sparse format"
write(*,*)"---------------------------------------------------------------------------------"

if (option == 1) then
 do i=1, ndim2
  write(*,'(*(F12.6,1X))')mat2(i,:)
 enddo
 write(*,*)"---------------------------------------------------------------------------------"
endif

write(*,*)
write(*,*)"---------------------------------------------------------------------------------"
write(*,*)"Matrix product with Fortran subroutine"
write(*,*)"---------------------------------------------------------------------------------"

! matrix multiplication with non-sparse format
call cpu_time(t_start)
do i=1,maxiter
 call matrices_prod(mat1, mat2, C, n)
enddo
call cpu_time(t_end)

if (option == 1) then
 do i=1, n
  write(*,'(*(F12.6,1X))')C(i,:)
 enddo
 write(*,*)"---------------------------------------------------------------------------------"
endif

write(*,'(A, 1X, F14.8)')"Time for the operation:",t_end-t_start
write(*,*)"---------------------------------------------------------------------------------"
write(*,*)

write(*,*)"---------------------------------------------------------------------------------"
write(*,*)"Matrix product with DGEMM"
write(*,*)"---------------------------------------------------------------------------------"

! matrix multiplication with non-sparse format with DGEMM
call cpu_time(t_start)
do i=1,maxiter
 call dgemm('N','N',n,n,n,1.D0,mat1,n,mat2,n,0.D0,C,n)
enddo
call cpu_time(t_end)

if (option == 1) then
 do i=1, n
  write(*,'(*(F12.6,1X))')C(i,:)
 enddo
 write(*,*)"---------------------------------------------------------------------------------"
endif

write(*,'(A, 1X, F14.8)')"Time for the operation:",t_end-t_start
write(*,*)"---------------------------------------------------------------------------------"



deallocate(R1, C1, V1, R2, C2, V2, C, mat1, mat2)
end program main


