program main
        use omp_lib

        implicit none
        integer :: n,i,j,nit,seed,k,l,nmax
        real*8,allocatable :: a(:,:),b(:),b0(:),c(:)
        real*8 :: lamb,tcpustart,tcpuend,tgpustart,tgpuend,dd,tmp,rootdd,eps,error

        !Requesting size, allocation and filler of the matrix
        
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

       !Requesting treshold's convergence
       write(*,*) 'max iteration number'
       read(*,*) nmax
       write(*,*) nmax

       !Random vector
       allocate(b(n),b0(n))
       call random_number(b)

       b=b/dsqrt(dot_product(b,b))
       b0=b
       write(*,*)
       write(*,*) 'CPU FORTRAN FUNCTION'

       allocate(c(n))
       nit=0
       error=1
       tcpustart=omp_get_wtime()
       do while(error.gt.eps)
          c=matmul(a,b)
          c=c/dsqrt(dot_product(c,c))
          error=dot_product(c-b,c-b)
          b=c
          nit=nit+1
          if (nit.gt.nmax) then
                  write(*,*) 'No convergence'
                  exit
          end if
       end do
       tcpuend=omp_get_wtime()
       lamb=0.0d0
       lamb=dot_product(b,matmul(a,b))/dot_product(b,b)
       write(*,*) 'Tot. iterations', nit
       write(*,*) 'lambda', lamb
       write(*,*) 'CPU execution time', tcpuend-tcpustart

       write(*,*)
       write(*,*) 'GPU'

       b=b0
       nit=0
       error=1
       tgpustart=omp_get_wtime()
       do while(error.gt.eps)
          do k=1,n
             c(k)=0.d0
          end do
          dd=0.d0
          error=0.d0

          do k=1,n
           do l=1,n
              c(k)=c(k)+a(k,l)*b(l)
           end do
          end do

          do k=1,n
             dd=dd+c(k)*c(k)
          end do

          rootdd=dsqrt(dd)

          do k=1,n
             c(k)=c(k)/rootdd
          end do

          do k=1,n
             error=error+(c(k)-b(k))*(c(k)-b(k))
          end do

          do k=1,n
             b(k)=c(k)
          end do
          nit=nit+1
          if (nit.gt.nmax) then
                  write(*,*) 'No convergence'
                  exit
          end if
       end do

       tgpuend=omp_get_wtime()
       lamb=0.d0
       lamb=dot_product(b,matmul(a,b))/dot_product(b,b)
       write(*,*) 'Tot. iterations', nit
       write(*,*) 'lamb', lamb
       write(*,*) 'GPU execution time', tgpuend-tgpustart

       write(*,*)
