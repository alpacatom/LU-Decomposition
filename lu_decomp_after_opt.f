      program main 
       include 'mpif.h'
c       include 'lu.inc'
       common /mpienv/myid,numprocs
     
       integer MATRIX, NN, count
       parameter (MATRIX = 1,NN = 7680, count = 100)
       double precision EPS
       parameter(EPS = 2.220446e-16)

       double precision A(NN, NN)
       double precision b(NN)
       double precision x(NN)

       double precision t0, t1, t2, t_w, sum_t_w, sum_mflops
       double precision d_mflops, dtemp, dtemp2, dtemp_t
       integer ierr
       integer i, j
       integer ii      
       integer ib


       call MPI_INIT(ierr)
       call MPI_COMM_RANK(MPI_COMM_WORLD, myid, ierr)
       call MPI_COMM_SIZE(MPI_COMM_WORLD, numprocs, ierr)

!     cache
c      --- matrix generation --------------------------
       if (MATRIX .eq. 1) then
         do j=1, NN
           ii = 0
           do i=j, NN
             A(j, i) = dble(NN-j-ii+1)
             A(i, j) = A(j, i)
             ii = ii + 1
           enddo
         enddo

       else 
         call RANDOM_SEED
         do j=1, NN
           do i=1, NN
             call RANDOM_NUMBER(dtemp)
             A(i, j) = dtemp
           enddo
         enddo  
      
       endif
c      --- end of matrix generation --------------------------

c      --- set vector b  -----------------------------------
       do i=1, NN
         b(i) = 0.0d0
         do j=1, NN
           b(i) = b(i) + A(j,i)
         enddo
       enddo
c      -----------------------------------------------------


c      --- Start of LU routine ----------------------------
       call MPI_BARRIER(MPI_COMM_WORLD, ierr)
       t1 = MPI_WTIME(ierr)

       call MyLUsolve(A, b, x, NN)

       call MPI_BARRIER(MPI_COMM_WORLD, ierr)
       t2 = MPI_WTIME(ierr)

       t0 =  t2 - t1
       call MPI_REDUCE(t0, t_w, 1, MPI_DOUBLE_PRECISION,
     &            MPI_MAX, 0, MPI_COMM_WORLD, ierr)
c      --- End of LU routine --------------------------- */

       if (myid .eq. 0) then
c         print *, " NN  = ", NN
c         print *, " LU solve time[sec.] = ",t_w

         d_mflops = 2.0d0/3.0d0*dble(NN)*dble(NN)*dble(NN)
         d_mflops = d_mflops + 7.0d0/2.0d0*dble(NN)*dble(NN)
         d_mflops = d_mflops + 4.0d0/3.0d0*dble(NN)
         d_mflops = d_mflops/t_w
         d_mflops = d_mflops * 1.0d-6
c         print *, " MFLOPS = ", d_mflops
       endif

c      --- Verification routine ----------------- 
       ib = NN / numprocs
       dtemp_t = 0.0d0
       do j=myid*ib+1, (myid+1)*ib
         dtemp2 = x(j) - 1.0d0
         dtemp_t = dtemp_t + dtemp2 * dtemp2
       enddo
       dtemp_t = dsqrt(dtemp_t)
c      ------------------------------------------

       call MPI_REDUCE(dtemp_t, dtemp, 1, MPI_DOUBLE_PRECISION, 
     &        MPI_SUM, 0, MPI_COMM_WORLD, ierr)

c      --- Do not modify follows. -------------------- 
       if (myid .eq. 0) then

         if (MATRIX .eq. 1) then
           dtemp2 = dble(NN)*dble(NN)*dble(NN)
         else 
           dtemp2 = dble(NN)*dble(NN)
         endif
         dtemp_t = EPS * dble(NN) * dtemp2
c         print *, "Pass value: ", dtemp_t
c         print *, "Calculated value: ", dtemp
         if (dtemp .gt. dtemp_t) then
c           print *, "Error! Test is falled."
           stop
         endif
c         print *, " OK! Test is passed."

       endif
c      ----------------------------------------- 

!     start
       sum_t_w = 0.0d0
       sum_mflops = 0.0d0
       do k = 1, count, 1
c      --- matrix generation --------------------------
       if (MATRIX .eq. 1) then
         do j=1, NN
           ii = 0
           do i=j, NN
             A(j, i) = dble(NN-j-ii+1)
             A(i, j) = A(j, i)
             ii = ii + 1
           enddo
         enddo

       else 

         call RANDOM_SEED
         do j=1, NN
           do i=1, NN
             call RANDOM_NUMBER(dtemp)
             A(i, j) = dtemp
           enddo
         enddo  
      
       endif
c      --- end of matrix generation --------------------------

c      --- set vector b  -----------------------------------
       do i=1, NN
         b(i) = 0.0d0
         do j=1, NN
           b(i) = b(i) + A(j,i)
         enddo
       enddo
c      -----------------------------------------------------


c      --- Start of LU routine ----------------------------
       call MPI_BARRIER(MPI_COMM_WORLD, ierr)
       t1 = MPI_WTIME(ierr)

       call MyLUsolve(A, b, x, NN)

       call MPI_BARRIER(MPI_COMM_WORLD, ierr)
       t2 = MPI_WTIME(ierr)

       t0 =  t2 - t1
       call MPI_REDUCE(t0, t_w, 1, MPI_DOUBLE_PRECISION,
     &            MPI_MAX, 0, MPI_COMM_WORLD, ierr)
c      --- End of LU routine --------------------------- */

       if (myid .eq. 0) then
c         print *, " NN  = ", NN
c         print *, " LU solve time[sec.] = ",t_w
         d_mflops = 2.0d0/3.0d0*dble(NN)*dble(NN)*dble(NN)
         d_mflops = d_mflops + 7.0d0/2.0d0*dble(NN)*dble(NN)
         d_mflops = d_mflops + 4.0d0/3.0d0*dble(NN)
         d_mflops = d_mflops/t_w
         d_mflops = d_mflops * 1.0d-6
c         print *, " MFLOPS = ", d_mflops
         sum_t_w = sum_t_w + t_w
         sum_mflops = sum_mflops + d_mflops
       endif

c      --- Verification routine ----------------- 
       ib = NN / numprocs
       dtemp_t = 0.0d0
       do j=myid*ib+1, (myid+1)*ib
         dtemp2 = x(j) - 1.0d0
         dtemp_t = dtemp_t + dtemp2 * dtemp2
       enddo
       dtemp_t = dsqrt(dtemp_t)
c      ------------------------------------------

       call MPI_REDUCE(dtemp_t, dtemp, 1, MPI_DOUBLE_PRECISION, 
     &        MPI_SUM, 0, MPI_COMM_WORLD, ierr)

c      --- Do not modify follows. -------------------- 
       if (myid .eq. 0) then

         if (MATRIX .eq. 1) then
           dtemp2 = dble(NN)*dble(NN)*dble(NN)
         else 
           dtemp2 = dble(NN)*dble(NN)
         endif
         dtemp_t = EPS * dble(NN) * dtemp2
c         print *, "Pass value: ", dtemp_t
c         print *, "Calculated value: ", dtemp
         if (dtemp .gt. dtemp_t) then
c           print *, "Error! Test is falled."
           stop
         endif
c         print *, " OK! Test is passed."
       endif
c      ----------------------------------------- 
       enddo

       if(myid .eq. 0) then

          print *, "procs: ", numprocs, " time-mean = ", 
&                  sum_t_w / count, "MFLOPS-mean = ", 
&                  sum_mflops / count 

       endif

       call MPI_FINALIZE(ierr)

       stop
       end


       subroutine MyLUsolve(A, b, x, n)
       include 'mpif.h'
       common /mpienv/myid,numprocs
       double precision A(n, n)
       double precision b(n)
       double precision x(n)
       double precision c(n) 

       double precision ctmp
       integer n
       integer size, s, e, w, judge, colnum, dif
 
       integer i, j, k, tags, h 
       double precision dtemp

       double precision SBUF, RBUF

       size = n / numprocs
       s =  1 + myid * size
       e =  (myid+1) * size
       colnum = n - 1
       tags = 0
       dif = e - s + 1

       if(myid .eq. numprocs-1 .and. e .ne. n) then
          e = n
       endif 

!     1  
       do k = 0 , n - 1
          w = k / size 
 
          if (w .eq. myid) then
             dtemp = 1.0d0 / A(k+1, k+1)
             do i = k + 2, n
                A(i, k+1) = A(i, k+1) * dtemp
                c(i) = A(i, k+1) 
             enddo
          endif

          call MPI_BCAST( c(k + 2), n - k - 1 , MPI_DOUBLE_PRECISION, w,
&                     MPI_COMM_WORLD, ierr )
          if(myid .eq. 0) then
             judge = s + 1 + k
          else
             if(k+2 < s) then
                judge = s
             else
                judge = k + 2
             endif
          endif
 
          do j = judge, e, 1
             dtemp = A(k+1,j)
            do i = k + 2, n
                A(i,j) = A(i,j) - dtemp * c(i)
             enddo
          enddo
      enddo

      dif = e - s + 1
      do k = 0, numprocs - 1, 1                              
         call MPI_BCAST(A(2:n, 1+k*dif:(k+1)*dif), (n-1)*dif,
&                       MPI_DOUBLE_PRECISION,k,MPI_COMM_WORLD,ierr)
      enddo 

!     2     
       do k = 1, n, 1
          c(k) = 0.0d0
          x(k) = 0.0d0
       enddo

       do k = 1, n, size
          if (k >= s) then             
             w = (k-1) / size
c             print *, "Worker is ", w 
             if(myid .ne. 0) then
                call MPI_RECV(c(k), size, MPI_DOUBLE_PRECISION,
&                             myid-1, k, MPI_COMM_WORLD, istatus,
&                             ierr)
c                print *, "MyID:", myid, " Recv C:",c(k)
             endif
!     worker
             if(myid .eq. w) then
                do j = 1, size, 1
                   c(k + j - 1) = b(k + j - 1) + c(k + j - 1)
                   do i = s, s + j - 2
                      c(k + j - 1) = c(k + j - 1) - A(k + j - 1, i) 
&                        * c(i)
                   enddo
                enddo
!     sub worker
             else
                do j = 1, size
                   do i = s, e
                      c(k + j - 1) = c (k + j - 1) - A(k + j - 1, i)
&                        * c(i)
                   enddo
                enddo
                if(myid .ne. numprocs-1) then
                      call MPI_SEND(c(k), size, MPI_DOUBLE_PRECISION,
&                             myid+1, k, MPI_COMM_WORLD, ierr)
                endif
             endif
          endif
       enddo
            
!     3      
      h = size
      if(myid .eq. numprocs-1) then         
         x(n) = c(n) / A(n,n) 
      endif
      k = n-1
      do
         if(k .eq. n-1) then
            h = size-1
         else
            h = size
         endif

         if (k <= e) then
            w = (k-1) / size
             if(myid .ne. numprocs-1) then
                call MPI_RECV(x(k-size+1), size, MPI_DOUBLE_PRECISION,
&                             myid+1, k, MPI_COMM_WORLD, istatus,
&                             ierr)
             endif

             if(myid .eq. w) then 
                do j =  1, h, 1
                   x(k - j + 1) = c(k - j + 1) + x(k - j + 1)
                   if(h .eq. size-1) then
                   do i = e, e-j+1, -1
                      x(k - j + 1) = x(k - j + 1) - A(k - j + 1, i) 
&                        * x(i)
                   enddo
                   else 
                      do i = e, e-j+2, -1
                         x(k - j + 1) = x(k - j + 1) - A(k - j + 1, i) 
&                        * x(i)
                   enddo     
                   endif
                 x(k-j+1) = x(k-j+1) / A(k-j+1,k-j+1)
                enddo
             else
                do j = 1, h, 1
                   do i = s, e, 1
                      x(k - j + 1) = x (k - j + 1) - A(k - j + 1, i)
&                        * x(i)
                   enddo
                enddo
                if(myid .ne. 0) then
                      call MPI_SEND(x(k-size+1), size,
&                                  MPI_DOUBLE_PRECISION, myid-1,
&                                  k, MPI_COMM_WORLD, ierr)
                endif           ! myid .ne. 0
             endif              ! myid . eq. w 
          endif 
          k = k - h
          if(k < 1) then
             exit
          endif
       enddo

      return
      end

