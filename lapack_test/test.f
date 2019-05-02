      program simple_test

      implicit none            
      
      integer n,m
      parameter (n=5)
      real matA(n,n),matV(n,n)
      real wr(n),wi(n)

      complex cA(n,n),cV(n,n)
      complex w(n)

      character outfmt*32
      character outfmt2*32

      integer seed

      integer i,j

!      write(outfmt,'(A32)') (' ', j=1,32)
!      write(outfmt2,'(A32)') (' ', j=1,32)

      write(outfmt,'(A1,I2.2,A13)') '(',n,'(E15.8E2,1x))'

      write(outfmt2,'(A1,I2.2,A27)') '(',n,
     $                '(E15.8E2,1x,E15.8E2,A1,1x))'

      write(6,*) outfmt
      write(6,*) outfmt2

      m = n*n            
      seed = 86456
      call srand(seed)
      do i=1,m
        matA(i,1) = rand()
        cA(i,1)   = complex(matA(i,1),0)
      enddo

      do i=1,n
        write(6,outfmt) (matA(i,j), j=1,n)
      enddo

      write(6,*) ' '
      do i=1,n
        write(6,outfmt2) (real(cA(i,j)),imag(cA(i,j)),'i',j=1,n)
      enddo
      
      call wrp_cschur(cA,n,n,cV,n,w) 

      do i=1,n
        write(6,outfmt2) (real(cA(i,j)),imag(cA(i,j)),'i',j=1,n)
      enddo

      end program simple_test
