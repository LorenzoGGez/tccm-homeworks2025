program main
        use omp_lib

        implicit none
        integer :: n,i,j,nit,seed,k,l,nmax
        real*8,allocatable :: a(:,:),b(:),b0(:),c(:)
        real*8 :: lamb,tcpustart,tcpuend,tgpustart,tgpuend,dd,tmp,rootdd,eps,err

        !Request size, allocation and filler of the matrix
        
        write(*,*) 'matrix size'
        read(*,*) n
        write(*,*) n
        allocate(a(n:n))
        a=0.0d0
        do i=1,n-1
           j=i+1
           a(i,j)=1.0d0
           a(j,i)=1.0d0
        end do
        write(*,*)

        
