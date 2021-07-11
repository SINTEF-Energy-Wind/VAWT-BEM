include "fourier.f95"
include "wind.f95"

program testMann

   use wind

   implicit none

   integer                        :: fid,seed,Np
   integer,          dimension(3) :: Nbox
   double precision, dimension(3) :: dL
   double precision               :: Vinf,ti,Lu

   integer :: ierr

   open (UNIT=41, FILE='wind.txt', STATUS='REPLACE')

   ierr = 0

   fid = 49
   Nbox(1) = 4096
   Nbox(2) = 64
   Nbox(3) = 64
   Np = Nbox(1)*Nbox(2)*Nbox(3)

   dL(1) = 2.d0
   dL(2) = 1.d0
   dL(3) = 1.d0

   Vinf = 8.d0
   ti = 0.230d0
   Lu = 36.5d0
!   Lu = 210.d0  !  Calibrated, corresponds to Burton Lu = 160.d0.
                !  Valid at least for 7 to 25 m/s.

   seed = 356

   call windField (seed,Np,Nbox,dL,Vinf,ti,Lu,fid,ierr)
print *,'ierr = ',ierr

   close (41)

end program testMann