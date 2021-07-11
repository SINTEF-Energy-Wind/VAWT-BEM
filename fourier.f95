!-------------------------------------------------------------------------------



module fourier

contains


   subroutine FFT (data,nn,isign,dataout,ierr)
   !  Press et al. "Numerical Recipes" Cambridge University Press, 1986
   !  pages 394-395
   !  "[Computes the] discrete Fourier transform [of data], if isign is input as 1;
   !  or [computes] ... its inverse discrete Fourier transform, if isign is
   !  input as -1.  data is a complex array of length nn or, equivalently, a real
   !  array of length 2nn.  nn must be an integer power of 2 ..." (p 394)

      implicit none

      integer, intent(in) :: nn, isign
      double precision, dimension(2*nn), intent(in) :: data
      double precision, dimension(2*nn), intent(out) :: dataout

      integer :: ierr

      integer :: i,j,m,n,mmax,istep,div
      double precision :: wr,wi,wpr,wpi,wtemp,theta,tempi,tempr

      dataout = data

      !  Verify nn is an integer power of 2.
      div = 1

      do
         if (div .eq. nn) exit
         if (div .gt. nn) then
print *,'nn not a power of 2.'
            ierr = -1
            return
         end if
         div = div*2
      end do

      n = 2*nn
      j = 1

      !  This is the bit-reversal section of the routine.
      do i = 1,n,2
         if (j .gt. i) then
            !  Exchange the two complex numbers.
            tempr = dataout(j)
            tempi = dataout(j+1)
            dataout(j) = dataout(i)
            dataout(j+1) = dataout(i+1)
            dataout(i) = tempr
            dataout(i+1) = tempi
         end if
         m=n/2
         do
            if ((m .lt. 2) .or. (j .le. m)) exit
            j = j - m
            m = m/2
         end do
         j = j + m
      end do

      !  Here begins the Danielson-Lanczos section of the routine.
      mmax = 2

      !  Outer loop executed log_2 nn times.
      do

         if (n .le. mmax) exit

         istep = 2*mmax
         !  Initialize for the trigonometric recurrence.
         theta = 6.28318530717959d0/dble(isign*mmax)
         wpr = -2.d0*dsin(0.5d0*theta)**2
         wpi = dsin(theta)
         wr = 1.d0
         wi = 0.d0

         do m = 1,mmax,2
            do i = m,n,istep
               !  This is the Danielson-Lanczos formula:
               j = i + mmax
               tempr = wr*dataout(j) - wi*dataout(j+1)
               tempi = wr*dataout(j+1) + wi*dataout(j)
               dataout(j) = dataout(i) - tempr
               dataout(j+1) = dataout(i+1) - tempi
               dataout(i) = dataout(i) + tempr
               dataout(i+1) = dataout(i+1) + tempi
            end do
            !  Trigonometric recurrence.
            wtemp = wr
            wr = wr*wpr - wi*wpi + wr
            wi = wi*wpr + wtemp*wpi + wi
         end do

         mmax = istep

      end do

      !  If this is an inverse transform, multiply by the normalizing factor 1/nn.
!  CANCEL THAT:  SEE AUG 22, 2009 NOTES.
!      if (isign .eq. -1) then
!         dataout = dataout/dble(nn)
!      end if

   end subroutine FFT







   subroutine FFTn (datain,NN,Ndim,isign,data,ierr)
   !  Press et al. "Numerical Recipes" Cambridge University Press, 1986
   !  pages 451-453

      implicit none

      integer,                           intent(in)  :: Ndim,isign
      integer,          dimension(Ndim), intent(in)  :: NN
      double precision, dimension(:),    intent(in)  :: datain
      double precision, dimension(:),    intent(out) :: data

      integer, intent(inout) :: ierr

      integer :: Ntot,idim,Nprev,N,Nrem,ip1,ip2,ip3,i2rev,i1,i2,i3, &
                 i3rev,ibit,ifp1,ifp2,k1,k2,div
      double precision :: tempr,tempi,theta,wpr,wpi,wr,wi,wtemp

      !  Verify NN's are an integer power of 2.
      do idim = 1,Ndim
         div = 1

         do
            if (div .eq. NN(idim)) exit
            if (div .gt. NN(idim)) then
print *,'NN(',idim,') not a power of 2.'
               ierr = -1
               return
            end if
            div = div*2
         end do
      end do

      data = datain

      Ntot = 1
      do idim = 1,Ndim
         Ntot = Ntot*NN(idim)
      end do

      Nprev = 1
      do idim = 1,Ndim
         N = NN(idim)
         Nrem = Ntot/(N*Nprev)
         ip1 = 2*Nprev
         ip2 = ip1*N
         ip3 = ip2*Nrem
         i2rev = 1
         do i2 = 1,ip2,ip1
            if (i2 .lt. i2rev) then
               do i1 = i2,i2+ip1-2,2
                  do i3 = i1,ip3,ip2
                     i3rev = i2rev + i3 - i2
                     tempr = data(i3)
                     tempi = data(i3+1)
                     data(i3) = data(i3rev)
                     data(i3+1) = data(i3rev+1)
                     data(i3rev) = tempr
                     data(i3rev+1) = tempi
                  end do
               end do
            end if
            ibit = ip2/2
            do
               if (.not. ((ibit .ge. ip1) .and. (i2rev .gt. ibit))) then
                  exit
               else
                  i2rev = i2rev - ibit
                  ibit = ibit/2
               end if
            end do
            i2rev = i2rev + ibit
         end do
         ifp1 = ip1
         do
            if (.not. (ifp1 .lt. ip2)) then
               exit
            else
               ifp2 = 2*ifp1
               theta = isign*6.28318530717959d0/(ifp2/ip1)
               wpr = -2.d0*sin(0.5d0*theta)**2
               wpi = sin(theta)
               wr = 1.d0
               wi = 0.d0
               do i3 = 1,ifp1,ip1
                  do i1 = i3,i3+ip1-2,2
                     do i2 = i1,ip3,ifp2
                        k1 = i2
                        k2 = k1 + ifp1
                        tempr = wr*data(k2) - wi*data(k2+1)
                        tempi = wr*data(k2+1) + wi*data(k2)
                        data(k2) = data(k1) - tempr
                        data(k2+1) = data(k1+1) - tempi
                        data(k1) = data(k1) + tempr
                        data(k1+1) = data(k1+1) + tempi
                     end do
                  end do
                  wtemp = wr
                  wr = wr*wpr - wi*wpi + wr
                  wi = wi*wpr + wtemp*wpi + wi
               end do
               ifp1 = ifp2
            end if
         end do
         Nprev = N*Nprev
      end do



   end subroutine FFTn








   subroutine rFFTn (datain,NN,Ndim,isign,data,ierr)
   !  Press et al. "Numerical Recipes" Cambridge University Press, 1986
   !  pages 451-453

      implicit none

      integer, parameter :: ikind = SELECTED_REAL_KIND(p=4,r=2)
      integer,                           intent(in)  :: Ndim,isign
      integer,          dimension(Ndim), intent(in)  :: NN
      real(kind=ikind), dimension(:),    intent(in)  :: datain
      real(kind=ikind), dimension(:),    intent(out) :: data

      integer, intent(inout) :: ierr

      integer :: Ntot,idim,Nprev,N,Nrem,ip1,ip2,ip3,i2rev,i1,i2,i3, &
                 i3rev,ibit,ifp1,ifp2,k1,k2,div
      real(kind=ikind) :: tempr,tempi,theta,wpr,wpi,wr,wi,wtemp

      !  Verify NN's are an integer power of 2.
      do idim = 1,Ndim
         div = 1

         do
            if (div .eq. NN(idim)) exit
            if (div .gt. NN(idim)) then
print *,'NN(',idim,') not a power of 2.'
               ierr = -1
               return
            end if
            div = div*2
         end do
      end do

      data = datain

      Ntot = 1
      do idim = 1,Ndim
         Ntot = Ntot*NN(idim)
      end do

      Nprev = 1
      do idim = 1,Ndim
         N = NN(idim)
         Nrem = Ntot/(N*Nprev)
         ip1 = 2*Nprev
         ip2 = ip1*N
         ip3 = ip2*Nrem
         i2rev = 1
         do i2 = 1,ip2,ip1
            if (i2 .lt. i2rev) then
               do i1 = i2,i2+ip1-2,2
                  do i3 = i1,ip3,ip2
                     i3rev = i2rev + i3 - i2
                     tempr = data(i3)
                     tempi = data(i3+1)
                     data(i3) = data(i3rev)
                     data(i3+1) = data(i3rev+1)
                     data(i3rev) = tempr
                     data(i3rev+1) = tempi
                  end do
               end do
            end if
            ibit = ip2/2
            do
               if (.not. ((ibit .ge. ip1) .and. (i2rev .gt. ibit))) then
                  exit
               else
                  i2rev = i2rev - ibit
                  ibit = ibit/2
               end if
            end do
            i2rev = i2rev + ibit
         end do
         ifp1 = ip1
         do
            if (.not. (ifp1 .lt. ip2)) then
               exit
            else
               ifp2 = 2*ifp1
               theta = isign*6.28318530717959/(ifp2/ip1)
               wpr = -2.0*sin(0.5*theta)**2
               wpi = sin(theta)
               wr = 1.
               wi = 0.
               do i3 = 1,ifp1,ip1
                  do i1 = i3,i3+ip1-2,2
                     do i2 = i1,ip3,ifp2
                        k1 = i2
                        k2 = k1 + ifp1
                        tempr = wr*data(k2) - wi*data(k2+1)
                        tempi = wr*data(k2+1) + wi*data(k2)
                        data(k2) = data(k1) - tempr
                        data(k2+1) = data(k1+1) - tempi
                        data(k1) = data(k1) + tempr
                        data(k1+1) = data(k1+1) + tempi
                     end do
                  end do
                  wtemp = wr
                  wr = wr*wpr - wi*wpi + wr
                  wi = wi*wpr + wtemp*wpi + wi
               end do
               ifp1 = ifp2
            end if
         end do
         Nprev = N*Nprev
      end do



   end subroutine rFFTn











   subroutine twoFFT (nn,data1,data2,fft1,fft2,ierr)

      implicit none

      integer :: j,n2
      integer, intent(in) :: nn
      integer, intent(inout) :: ierr
      double precision, dimension(nn), intent(in) :: data1,data2
      double precision, dimension(2*nn), intent(out) :: fft1,fft2
      double precision, dimension(2*nn) :: qr,wr
      complex(KIND=KIND(1.d0)) :: h1,h2,c1,c2
      complex(KIND=KIND(1.d0)), dimension(nn) :: w1,w2

      c1 = cmplx(0.5d0,0.d0,KIND(1.d0))
      c2 = cmplx(0.d0,-0.5d0,KIND(1.d0))

      do j = 1,nn
         qr(2*j-1) = data1(j)
         qr(2*j) = data2(j)
      end do

      call FFT (qr,nn,1,wr,ierr)
      if (ierr .ne. 0) then
         return
      end if

      do j = 1,nn
         w1(j) = cmplx(wr(2*j-1),wr(2*j),KIND(1.d0))
      end do

      w2(1) = cmplx(aimag(w1(1)),0.d0,KIND(1.d0))
      w1(1) = cmplx(real(w1(1)),0.d0,KIND(1.d0))
      n2 = nn + 2
      do j = 2,nn/2+1
         h1 = c1*(w1(j) + conjg(w1(n2-j)))
         h2 = c2*(w1(j) - conjg(w1(n2-j)))
         w1(j) = h1
         w1(n2-j) = conjg(h1)
         w2(j) = h2
         w2(n2-j) = conjg(h2)
      end do

      do j = 1,nn
         fft1(2*j-1) = real(w1(j))
         fft1(2*j) = aimag(w1(j))
         fft2(2*j-1) = real(w2(j))
         fft2(2*j) = aimag(w2(j))
      end do

   end subroutine twoFFT





end module fourier