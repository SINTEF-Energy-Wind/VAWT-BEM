!-------------------------------------------------------------------------------



module wind
! Contains subroutines to calculate turbulent wind velocity at
! specified points in space.
!
!
! Revisions:
!     Date          Programmer          Notes
!  ----------      ------------         ----------------------------------------
!  20.10.2011        K. Merz            Original code.
!  21.10.2013        K. Merz            Updated the expression for alpha*epsilon^(2/3)
!                                       to better match isotropic von Karman.
!
!

   implicit none

   contains




   subroutine windField (seed,Np,Nbox,dL,Vinf,ti,Lu,fid,ierr)
   !  This subroutine calculates velocities at points over a box
   !  of space, using the method of Mann, "Wind Field Simulation",
   !  Probabilistic Engineering Mechanics 13 (1998) 269-282.
   !
   !  Currently an isotropic turbulence spectrum is used, not 
   !  Mann's extension to sheared boundary layer flow.
   !
   !  Inputs:
   !  -------
   !  seed            : Seed for random number generator.
   !  Np              : = Nbox(1)*Nbox(2)*Nbox(3)
   !  Nbox            : Number of points along x,y,z dimensions.
   !  dL              : Spatial separation of x,y,z points.
   !  Vinf            : Remote average windspeed.
   !  ti              : Turbulence intensity.
   !  Lu              : Length scale of turbulence.
   !  fname,fid       : Name and ID of output file.
   !
   !  Outputs:
   !  --------
   !  [To file:]
   !  u               : A vector of turbulent velocity components
   !                    in space.  Used sequentially for ux, uy, uz.
   !
   !  Local variables:
   !  ----------------
   !  L               : Length of sides of spatial box.
   !  dk              : Wave number increment for each dimension.
   !  k               : Wave number.
   !  kk              : Sum of squares of wave number components.
   !  theta           : Phase angle.
   !  KE              : Turbulent kinetic energy.
   !  C               : Mann's C matrix.
   !  Cnx,y,z         : Mann's coefficients Cij times random phase nj.

      use fourier

      implicit none

      double precision, parameter :: PI = 3.141592653589793d0

      integer, parameter :: ikind = SELECTED_REAL_KIND(p=4,r=2)

      integer,                            intent(in) :: seed,Np,fid
      integer,          dimension(3),     intent(in) :: Nbox
      double precision, dimension(3),     intent(in) :: dL
      double precision,                   intent(in) :: Vinf,ti,Lu

      integer,                         intent(inout) :: ierr

      character(LEN=51)                              :: fname
      integer                                        :: i,j,div,ik,ik1,ik2,ik3, &
                                                        ikre,ikim,irec,irecl,icomp
      double precision                               :: twoPI,sqdk3,            &
                                                        Luk,aeps23,KE,          &
                                                        val,sigu,kk,root2,temp
      double precision, dimension(3)                 :: L,dk,k
      double precision, dimension(3,3)               :: C

      !  Reduce memory usage by switching to lower precision.
      real(kind=ikind)                               :: temp1,temp2
      real(kind=ikind), dimension(6)                 :: gaus
      real(kind=ikind), dimension(2*Np)              :: Cnxyz,u

      !  Verify Nbox dimensions are an integer power of 2.
      do i = 1,3
         div = 1
         do
            if (div .eq. Nbox(i)) exit
            if (div .gt. Nbox(i)) then
print *,'Nbox(',i,') not a power of 2.'
               ierr = -1
               return
            end if
            div = div*2
         end do
      end do

      twoPI = 2.d0*PI
      sigu = Vinf*ti
      root2 = sqrt(2.d0)

      !  Wave number increment is based on size of the spatial
      !  box.
      do i = 1,3
         L(i) = Nbox(i)*dL(i)
         dk(i) = twoPI/L(i)
      end do

      sqdk3 = sqrt(dk(1)*dk(2)*dk(3))

      !  Generate all the random numbers I will need for a
      !  run, compute the inverse Gaussian distribution, and
      !  write to file.  In this way, I can avoid storing all
      !  velocity components in memory at once.  Rather, I can
      !  work with one velocity component at a time, and link
      !  the velocity components together by selecting the
      !  appropriate Gaussian numbers from the file.
print *,'Computing Gaussian variables and writing to file.'
print *,'Progress:'
      val = random (-seed)  !  Seed the random number generator.
      inquire (IOLENGTH=irecl) gaus  !  6 reals.
      fname = 'gaus.bin'
      open (UNIT=fid,FILE=fname,ACCESS='DIRECT', &
            FORM='UNFORMATTED',STATUS='REPLACE',RECL=irecl)
      irec = 0
      do ik3 = 1,Nbox(3)
print *,'   ',ik3,' of ',Nbox(3)
         do ik2 = 1,Nbox(2)
            do ik1 = 1,Nbox(1)
               do j = 1,6
                  !  Simulate an independent random variable with 
                  !  unit variance.  (Avoid exactly 1 or 0.)
                  val = 0.9999998d0*random(1) + 0.0000001d0
                  call invCDF (val,temp,ierr)
                  if (ierr .ne. 0) then
print *,'Error simulating Gaussian variables.'
                     return
                  end if
                  gaus(j) = real(temp,ikind)
               end do
               irec = irec + 1
               write (fid,REC=irec) gaus
            end do
         end do
      end do
      close(fid)

      do icomp = 1,3  !  X, Y, and Z components.
print *,'Velocity component ',icomp
         !  Open the Gauss file for reading.
         inquire (IOLENGTH=irecl) gaus
         fname = 'gaus.bin'
         open (UNIT=fid,FILE=fname,ACCESS='DIRECT',STATUS='OLD', &
               ACTION='READ',RECL=irecl)
         irec = 0

         !  Calculate components of (Cij nj) for all k and x,y,z
         !  velocity components.
print *,'   Generating Fourier coefficients with random phase.'
print *,'   Progress:'
         ik = 0
         do ik3 = 1,Nbox(3)
print *,'      ',ik3,' of ',Nbox(3)
            if (ik3 .le. Nbox(3)/2 + 1) then            
               !  Positive frequencies from 0 up to Nbox/2.
               k(3) = dble(ik3-1)*dk(3)
            else
               !  Negative frequencies from -(Nbox/2)+1 to -1.
               k(3) = dble(-Nbox(3)+ik3-1)*dk(3)
            end if
            do ik2 = 1,Nbox(2)
               if (ik2 .le. Nbox(2)/2 + 1) then            
                  k(2) = dble(ik2-1)*dk(2)
               else
                  k(2) = dble(-Nbox(2)+ik2-1)*dk(2)
               end if
               do ik1 = 1,Nbox(1)
                  if (ik1 .le. Nbox(1)/2 + 1) then            
                     k(1) = dble(ik1-1)*dk(1)
                  else
                     k(1) = dble(-Nbox(1)+ik1-1)*dk(1)
                  end if

                  ik = ik + 1  !  One complex value of (Cij nj)
                               !  for each k.
                  ikre = 2*ik-1
                  ikim = 2*ik

                  !  k(scalar) = |k(vector)|, when used in the energy expression.
                  kk = sqrt(k(1)**2 + k(2)**2 + k(3)**2)
                  Luk = Lu*kk

                  !  Turbulent kinetic energy.  Select (alpha*
                  !  epsilon)^(2/3) such that the nominal von
                  !  Karman spectrum from Burton is recovered.
                  if ((ik1 .eq. 1) .and. (ik2 .eq. 1) .and. (ik3 .eq. 1)) then
                     !  kk is zero, and C will come out as zero too.
                     aeps23 = 0.d0
                     KE = 0.d0
                     val = 0.d0
                  else
                     !  These give rough approximations to the von Karman spectrum in Burton.
!                     aeps23 = 1.7d0*(sigu**2/Lu**(2.d0/3.d0))
!                     aeps23 = 1.8d0*(sigu**2/Lu**(2.d0/3.d0))

                     !  Here is a much better calibration.  The input length scale 
                     !  must be adjusted also.  For Burton Lu = 160 m, use here Lu = 210 m.
                     !  This is valid over the entire operating windspeed range.
                     aeps23 = 1.5d0*(sigu**2/Lu**(2.d0/3.d0))

                     KE = aeps23*(Lu**(5.d0/3.d0))*(Luk**4) &
                        / ((1.d0 + Luk**2)**(17.d0/6.d0))

                     val = sqdk3*sqrt(KE/(4.d0*PI))/(kk**2)
                  end if

                  !  Mann Equations 50 and 46 give Cij.
                  C(1,1) = 0.d0
                  C(1,2) = val*k(3)
                  C(1,3) = -val*k(2)
                  C(2,1) = -val*k(3)
                  C(2,2) = 0.d0
                  C(2,3) = val*k(1)
                  C(3,1) = val*k(2)
                  C(3,2) = -val*k(1)
                  C(3,3) = 0.d0

                  !  Read the Gaussian variables.
                  irec = irec + 1
                  read(fid,REC=irec) gaus

                  !  Contribution to real and imaginary parts
                  !  of Cij nj.
                  if (icomp .eq. 1) then  !  X
                     Cnxyz(ikre) = real(C(1,2),ikind)*gaus(2) &
                                 + real(C(1,3),ikind)*gaus(3)
                     Cnxyz(ikim) = real(C(1,2),ikind)*gaus(5) &
                                 + real(C(1,3),ikind)*gaus(6)
                  else if (icomp .eq. 2) then  !  Y
                     Cnxyz(ikre) = real(C(2,1),ikind)*gaus(1) &
                                 + real(C(2,3),ikind)*gaus(3)
                     Cnxyz(ikim) = real(C(2,1),ikind)*gaus(4) &
                                 + real(C(2,3),ikind)*gaus(6)
                  else if (icomp .eq. 3) then  !  Z
                     Cnxyz(ikre) = real(C(3,1),ikind)*gaus(1) &
                                 + real(C(3,2),ikind)*gaus(2)
                     Cnxyz(ikim) = real(C(3,1),ikind)*gaus(4) &
                                 + real(C(3,2),ikind)*gaus(5)
                  end if
               end do  !  k1
            end do  !  k2
         end do  !  k3   

         !  Close the Gaussian file.
         close (fid)   

         !  FFTs to get velocity components, write the real part of each
         !  to file.
print *,'   FFT to get u.'
         call rFFTn (Cnxyz,Nbox,3,-1,u,ierr)

         !  Write to a binary file consisting of a list of double
         !  precision numbers.  Get the appropriate record length.
         !  See Chapman page 676.
print *,'   Writing u to velocity file.'
print *,'   Progress:'
         inquire (IOLENGTH=irecl) u(1),u(1),u(1)
         fname = 'u.bin'
         if (icomp .eq. 1) then  !  X
            open (UNIT=fid,FILE=fname,ACCESS='DIRECT', &
                  FORM='UNFORMATTED',STATUS='REPLACE',RECL=irecl)
            irec = 1
            write (fid,REC=irec) Nbox(1)
            irec = 2
            write (fid,REC=irec) Nbox(2)
            irec = 3
            write (fid,REC=irec) Nbox(3)
            irec = 4
            write (fid,REC=irec) dL(1)
            irec = 5
            write (fid,REC=irec) dL(2)
            irec = 6
            write (fid,REC=irec) dL(3)
            do i = 1,Np
if (mod(i,1000000) .eq. 1) then
print *,'      ',i,' of ',Np
end if
               irec = irec + 1
               write (fid,REC=irec) u(2*i-1)
            end do
         else
            open (UNIT=fid,FILE=fname,ACCESS='DIRECT',STATUS='OLD', &
                  ACTION='READWRITE',RECL=irecl)
            irec = 6
            if (icomp .eq. 2) then
               do i = 1,Np
if (mod(i,1000000) .eq. 1) then
print *,'      ',i,' of ',Np
end if
                  irec = irec + 1
                  read (fid,REC=irec) temp1
                  write (fid,REC=irec) temp1,u(2*i-1)
               end do
            else if (icomp .eq. 3) then
               do i = 1,Np
if (mod(i,1000000) .eq. 1) then
print *,'      ',i,' of ',Np
end if
                  irec = irec + 1
                  read (fid,REC=irec) temp1,temp2
                  write (fid,REC=irec) temp1,temp2,u(2*i-1)
               end do
            end if
         end if
         close(fid)

!if (icomp .eq. 1) then
!  Redefine k to mean x for debug output.
i = 0
do ik3 = 1,Nbox(3)
if (ik3 .le. Nbox(3)/2 + 1) then
k(3) = dble(ik3-1)*dL(3)
else
k(3) = dble(-Nbox(3)+ik3-1)*dL(3)
end if
do ik2 = 1,Nbox(2)
if (ik2 .le. Nbox(2)/2 + 1) then            
k(2) = dble(ik2-1)*dL(2)
else
k(2) = dble(-Nbox(2)+ik2-1)*dL(2)
end if
val = 0.d0
do ik1 = 1,Nbox(1)
if (ik1 .le. Nbox(1)/2 + 1) then            
k(1) = dble(ik1-1)*dL(1)
else
k(1) = dble(-Nbox(1)+ik1-1)*dL(1)
end if
i = i + 1
!write (41,'(1X,5ES13.4)') k(1),k(2),k(3),u(2*i-1),u(2*i)
if ((ik3 .eq. 1) .and. (ik2 .eq. 1)) then
if (ik1 .eq. 1) then
sigu = 0.d0
end if
write (41,'(1X,2ES13.4)') val,u(2*i-1)
val = val + dL(1)
sigu = sigu + u(2*i-1)**2
end if
end do
end do
end do
sigu = sqrt(sigu/Nbox(1))
print *,'stdev = ',sigu
!end if

      end do  !  icomp: X,Y,Z component

   end subroutine windField










   !  For g95 compiler.
   function random (k)
   !  Press et al. p. 195.
   !  This subroutine must be called once with a negative integer
   !  before beginning to use it.

      implicit none

      integer :: j,k,junk1,junk2
      integer :: seed
      real, dimension(97), save :: v
      real, save :: y
      real :: temp
      double precision :: random

      if (k .lt. 0) then
         seed = abs(k)
         temp = rand(seed)
         do j = 1,97
            temp = rand(0)
         end do
         do j = 1,97
            v(j) = rand(0)
         end do
         y = rand(0)
      end if

      j = 1 + int(97.0*y)
      if (j .gt. 97) j = 97
      if (j .lt. 1) j = 1
      y = v(j)
      random = dble(y)
      v(j) = rand(0)

   end function random


   !  For Silverfrost Plato compiler.
!   function random (k)
   !  Press et al. p. 195.
   !  This subroutine must be called once with a negative integer
   !  before beginning to use it.

!      implicit none

!      integer :: j,k,seed
!      double precision, dimension(97), save :: v
!      double precision, save :: y
!      double precision :: random,temp

!      if (k .lt. 0) then
!         seed = abs(k)
!         call random_seed(put=seed)
!         do j = 1,97
!            call random_number(temp)
!         end do
!         do j = 1,97
!            call random_number(v(j))
!         end do
!         call random_number(y)
!      end if

!      j = 1 + int(97.d0*y)
!      if (j .gt. 97) j = 97
!      if (j .lt. 1) j = 1
!      y = v(j)
!      random = y
!      call random_number(v(j))

!   end function random





   subroutine invCDF (y,x,ierr)
   !  Returns the inverse of the cumulative density
   !  function of a normal distribution.

      implicit none

      double precision, intent(in) :: y
      double precision, intent(out) :: x
      integer, intent(inout) :: ierr

      double precision, parameter :: PI = 3.14159265358979d0

      double precision :: ycalc,logy,A,B,C,D

      if ((y .le. 0.d0) .or. (y .ge. 1.d0)) then
         ierr = -1
         return
      end if

      if (y .le. 0.5d0) then
         ycalc = y
      else
         ycalc = 1.d0 - y
      end if

      if (ycalc .lt. 1.d-9) then
         ycalc = 1.d-9  !  Truncate to valid range of 
                        !  equations below.
      end if

      logy = log10(ycalc)

      if (logy .le. -0.521790760d0) then
         A = -16.47336308d0
         B = -16.11480089d0
         C = 0.169768608d0
         D = -1.628363803d0
         x = (-PI - asin((logy - A)/B) - D)/C
!write (41,'(1X,6ES18.9)') y,ycalc,logy,(logy - A)/B,asin((logy - A)/B),x
      else
         A = 2.815859327d0
         B = 4.58184272d0
         C = 1.12410157d0
         x = C + logy*(B + A*logy)
      end if

      if (y .gt. 0.5d0) then
         x = -x
      end if

   end subroutine invCDF





end module wind 