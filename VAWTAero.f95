!-------------------------------------------------------------------------------



module VAWTBEM
! Contains subroutines to perform a steady blade element-momentum calculation.
!
!
! Revisions:
!     Date          Programmer          Notes
!  ----------      ------------         ----------------------------------------
!  21.03.2009        K. Merz            Original code, quasi-steady.
!  15.06.2011        K. Merz            Rewrote code for time-domain analysis,
!                                       eliminating iteration on induced velocity.
!                                       Condensed code into a single VAWTBEM module,
!                                       with all stored variables passed through
!                                       subroutine arguments.
!  25.09.2011        K. Merz            Restructured code to be suited for numerical
!                                       integration.
!  08.11.2011        K. Merz            Corrected errors in Re definition for dynamic
!                                       stall initialization, matching of upwind and
!                                       downwind induced velocities, and velocity used
!                                       for force-to-coefficient-to-force transformation
!                                       during dynamic stall calculation (had to omit 
!                                       the spanwise component to be consistent).
!                                       
!
! 

   implicit none

   integer, parameter :: MAXSTRLEN = 80
   integer, parameter :: Maoa = 361 ! 180
   integer, parameter :: MRe = 7
   integer, parameter :: nErrorTypes = 53
!   integer, parameter :: coeffFlag = 2  !  1: normal.  2: fixed -180 to 180, 1 deg



   integer, parameter :: coeffFlag = 1




   double precision, parameter :: SMALLV = 0.1d0
   double precision, parameter :: PI = 3.1415926535897932384626433832795d0

   integer, dimension(nErrorTypes), save :: errorLog
   
   contains













   subroutine aeroTimestep (debug,QSflag0,x,dViqsdt,psi0h_g,omega, &
                            Nb,Nrows,Ncols,Ndof,t,                 &
                            Vinf_g,Vh_g,turbFid,Nbox,dL,Tg_w,xw0,  &
                            twist0,node,profile,                   &
                            Lblade,Ro,xs_r,n_r,density,viscosity,  &
                            chord,Lelem,Zelem,                     &
                            NRe,Naoa,Rea,aoaa,Cla,Cda,Cma,Tsr_s,   &
                            aoaz,aoafs,dClda_a,pitchAng,           &
                            QSflag,dxdt,aoa,fqs,Viqs,              &
                            V0_r,V0avg_g,Fb_r,F_g,ierr)
   !  This subroutine prepares the input for the VAWTBEMTimestep
   !  subroutine, calls it, and processes the output to get F_g.
   !
   !  Inputs:
   !  -------
   !  QSflag0          : The input dynamic stall flag.
   !  x                : The aerodynamic state vector: f, Vi, Viint.
   !  dViqsdt          : Rate of change of quasi-steady induced
   !                     velocity, used in the dynamic inflow 
   !                     calculation.
   !  psi0h_g          : Angle of the rotor relative to the global
   !                     coordinate system.
   !  omega            : rotational speed of the rotor.
   !  Nb,Nrows,Ncols   : Number blades, rows and cols of surface 
   !                     elements.
   !  Ndof             : Length of aerodynamic state variable vector.
   !  t                : The present time (for turbulence look-up).
   !  Vinf_g           : remote velocity, including wind shear,
   !                     in the global coordinate system.  Defined
   !                     at each surface element centroid, assuming
   !                     that the rotor is oriented in the global
   !                     coordinate system.
   !  Vh_g             : rotor translational velocity, global CS.
   !                     Defined at each surface element centroid.
   !  turbFid          : File ID of turbulence file.
   !  Nbox             : Number of turbulence grid points in the Xw,
   !                     Yw, and Zw directions (wind coordinates).
   !  dL               : Separation between grid points (Xw,Yw,Zw).
   !  Tg_w             : Transform matrix from global to wind coords.
   !  xw0              : Position of the rotor center in wind coords.
   !  twist0           : Blade twist (rad).
   !  node             : Nodal r-z-offset.
   !  profile          : Element center r-z-offset.
   !  Lblade           : Total blade length.
   !  Ro               : Outer radius.
   !  xs_r             : Surface element centroid coordinates.
   !  n_r              : Surface element normal vectors.
   !  density,viscosity
   !  chord            : Blade element chord.
   !  Lelem            : Element length in span direction.
   !  Zelem            : Along-span distance from equator.
   !  NRe ... Cma      : Airfoil coefficient data.
   !  Tsr_s            : Transform from surface to rotor coordinates.
   !  aoaz,aoafs,      : Zero-lift AOA, AOA at full stall (f = 1), and
   !  dClda_a            attached-flow lift coefficient slope.  Used in
   !                     the dynamic stall calculation.
   !  pitchAng         : Blade pitch angle.
   !  
   !  Outputs:
   !  --------
   !  QSflag           : Updated QSflag.
   !  dxdt             : Rate of change of f, Vi, and Viint.
   !  aoa,fqs          : Variables for dynamic stall logic.
   !  Viqs             : Quasi-steady induced velocity, for dynamic inflow.
   !  V0_r             : Incoming velocity, output in order to truncate Vi.
   !  Fb_r             : Local aerodynamic forces on each blade element.
   !  F_g              : Aerodynamic forces in global coordinates.

      implicit none

      integer,                                      intent(in) :: Nb,Nrows,Ncols,Ndof, &
                                                                  debug

      integer,          dimension(Nrows,Ncols),     intent(in) :: QSflag0
      double precision, dimension(Ndof),            intent(in) :: x
      double precision, dimension(3,Nrows,Ncols),   intent(in) :: dViqsdt

      double precision,                             intent(in) :: Lblade,Ro,t,psi0h_g, &
                                                                  density,viscosity,   &
                                                                  omega
      double precision, dimension(3),               intent(in) :: Vh_g
      double precision, dimension(Nrows),           intent(in) :: twist0,chord,Lelem,  &
                                                                  Zelem
      double precision, dimension(3,Nrows),         intent(in) :: profile
      double precision, dimension(3,Nrows+1),       intent(in) :: node
      double precision, dimension(Nrows,Ncols),     intent(in) :: aoaz,pitchAng
      double precision, dimension(2,Nrows,Ncols),   intent(in) :: aoafs,dClda_a
      double precision, dimension(3,Nrows,Ncols/2), intent(in) :: Vinf_g
      double precision, dimension(3,Nrows,Ncols),   intent(in) :: xs_r,n_r
      double precision, dimension(3,3,Nrows,Ncols), intent(in) :: Tsr_s

      !                             Airfoil properties.
      integer,          dimension(Nrows),           intent(in) :: NRe
      integer,          dimension(MRe,Nrows),       intent(in) :: Naoa
      double precision, dimension(MRe,Nrows),       intent(in) :: Rea
      double precision, dimension(Maoa,MRe,Nrows),  intent(in) :: aoaa,Cla,Cda,Cma

      !                             Turbulence.
      integer,                                      intent(in) :: turbFid
      integer,          dimension(3),               intent(in) :: Nbox
      double precision, dimension(3),               intent(in) :: dL,xw0
      double precision, dimension(3,3),             intent(in) :: Tg_w

      integer,          dimension(Nrows,Ncols),    intent(out) :: QSflag
      double precision, dimension(Ndof),           intent(out) :: dxdt
      double precision, dimension(Nrows,Ncols),    intent(out) :: aoa,fqs
      double precision, dimension(3),              intent(out) :: V0avg_g
      double precision, dimension(3,Nrows,Ncols),  intent(out) :: Viqs,V0_r
      double precision, dimension(6,Nrows,Ncols),  intent(out) :: Fb_r
      double precision, dimension(6),              intent(out) :: F_g

      integer,                                   intent(inout) :: ierr

      !  Local variables:
      integer                                       :: k,icol,irow,ib,icolup

      double precision                              :: Asum,r,area,dpsi,psi,cpsi,spsi,  &
                                                       anorm,theta,ct,st,psi0_g,psi0_r
      double precision, dimension(3)                :: vec2,VAsum,vel,Vh_r,Vi_r
      double precision, dimension(6)                :: F_r
      double precision, dimension(3,3)              :: Tg_r,Tr_g,Tss_r
      double precision, dimension(3,Nrows,Ncols/2)  :: Vinf_r
      double precision, dimension(3,Nrows,Ncols)    :: u_g,u_r,Vb_r
      double precision, dimension(3,3,Nrows,Ncols)  :: Tbr_a

!write (41,'(1X,2ES13.4)') density,viscosity
!write (41,'(1X,A)') 'node1-3 profile1-3 Lelem Zelem'
!do irow = 1,Nrows+1
!if (irow .eq. 1) then
!write (41,'(I4,3F10.3)') irow-1,node(1,irow),node(2,irow),node(3,irow)
!else 
!write (41,'(I4,8F10.3)')                               &
!irow-1,node(1,irow),node(2,irow),node(3,irow),         &
!profile(1,irow-1),profile(2,irow-1),profile(3,irow-1), &
!Lelem(irow-1),Zelem(irow-1)
!end if
!end do
!write (41,'(1X,A)') 'chord twist0 xs_r1-3 n_r1-3'
!do irow = 1,Nrows
!write (41,'(I4,8F10.3)')                      &
!irow,chord(irow),twist0(irow),                &
!xs_r(1,irow,9),xs_r(2,irow,9),xs_r(3,irow,9), &
!n_r(1,irow,9),n_r(2,irow,9),n_r(3,irow,9)
!end do
!do icol = 1,Ncols
!do irow = 1,Nrows
!write (41,'(1X,2I4,5ES13.4,A)')                             &
!irow,icol,aoaz(irow,icol),dClda_a(1,irow,icol),             &
!dClda_a(2,irow,icol),aoafs(1,irow,icol),aoafs(2,irow,icol), &
!'   irow icol aoaz dClda_a+- aoafs+-'
!end do
!end do
!do icol = 1,Ncols/2
!do irow = 1,Nrows
!write (41,'(1X,2I4,3ES13.4,A)') irow,icol,Vinf_g(:,irow,icol), &
!'   irow icol Vinf_g'
!end do
!end do
! Verified.

      dpsi = 2.d0*PI/dble(Ncols)
      
      !  Get turbulence at time t corresponding to the current
      !  timestep.  Assume the rotor is aligned with global
      !  coordinates.  (I shouldn't input the turbulence 
      !  velocities because I need to recompute turbulence
      !  anyways once I find the instantaneous "upwind" and
      !  "downwind" orientations.)
      do icol = 1,Ncols
         do irow = 1,Nrows
            call readTurbulence (turbFid,dL,Nbox,xw0,xs_r(:,irow,icol), &
                                 Tg_w,u_g(:,irow,icol))
         end do
      end do





!u_g = 0.d0






!irow = 10
!icol = 9
!write (41,'(1X,2I4,9F10.4,A)')                    &
!irow,icol,xw0,xs_r(:,irow,icol),u_g(:,irow,icol), &
!'   irow icol xw0_x,y,z xs_r_x,y,z u_g_x,y,z'
! Verified for the case where Tg_w is the identity matrix.

!if (t .gt. 30.d0) then
!u_g(1,:,:) = 2.d0*sin(2.d0*PI*0.05*(t-30.d0))
!else
!u_g(1,:,:) = 0.d0
!end if
!u_g(2,:,:) = 0.d0
!u_g(3,:,:) = 0.d0

      !  Calculate the area-weighted average velocity vector in
      !  global coordinates.  Take advantage of symmetry: assume 
      !  during this calculation that the rotor is aligned with
      !  the global coordinate system.  This means that n_g can
      !  be taken as n_r, and same with xs_g => xs_r.
      VAsum = 0.d0  !  3-vector.
      Asum = 0.d0   !  Scalar.
      do icol = 1,Ncols/2
         do irow = 1,Nrows
            r = sqrt(xs_r(1,irow,icol)**2 + xs_r(2,irow,icol)**2)  !  Note xs_r<--!
            area = Lelem(irow)*r*dpsi
            vel = Vinf_g(:,irow,icol) + u_g(:,irow,icol)
            anorm = area*abs(dot(vel,n_r(:,irow,icol)))/magnitude(vel)  !  Note n_r<--!
            VAsum = VAsum + vel*anorm
            Asum = Asum + anorm
!write (41,'(1X,2I4,10ES13.4,A)') irow,icol,r,area,vel,anorm,VAsum,Asum, &
!'   irow icol r area (Vinf+u)1-3 anorm VAsum Asum'
         end do
      end do
      V0avg_g = -Vh_g + VAsum/Asum
!write (41,'(1X,12ES13.4,A)') u_g(:,20,9),Vh_g,VAsum/Asum,V0avg_g, &
!'   u_g(20,9) Vh_g (Vinf+u)avg V0avg_g'
! Verified.

      !  Calculate the X^r direction and the Tg_r and Tr_g 
      !  transform pair.  Theta is the angle from the X^g
      !  to the X^r axis assuming no tilt.
      !
      !  (NOTE: IN THE EVENT OF LARGE TILT THIS WILL NEED
      !  TO BE UPDATED WITH A MORE ROBUST TRANSFORM.)
      if ((V0avg_g(1) .ge. SMALLV) .or. (V0avg_g(2) .ge. SMALLV)) then
         theta = atan2(V0avg_g(2),V0avg_g(1))
      else
         theta = 0.d0
      end if
      ct = cos(theta)
      st = sin(theta)
      Tg_r(1,1) = ct
      Tg_r(1,2) = st
      Tg_r(1,3) = 0.d0
      Tg_r(2,1) = -st
      Tg_r(2,2) = ct
      Tg_r(2,3) = 0.d0
      Tg_r(3,1) = 0.d0
      Tg_r(3,2) = 0.d0
      Tg_r(3,3) = 1.d0
      Tr_g(1,1) = ct
      Tr_g(1,2) = -st
      Tr_g(1,3) = 0.d0
      Tr_g(2,1) = st
      Tr_g(2,2) = ct
      Tr_g(2,3) = 0.d0
      Tr_g(3,1) = 0.d0
      Tr_g(3,2) = 0.d0
      Tr_g(3,3) = 1.d0

      !  Compute a normalized rotor head angle.
      psi0_g = mod(psi0h_g,2.d0*PI)
      psi0_r = psi0_g - theta

      !  Re-compute turbulence at the element centroids considering
      !  the actual rotor orientation.
      !  [SKIP FOR NOW, ONLY RELEVANT WHEN FLOW IS TURBULENT AND
      !  WHEN THE SUM OF TURBULENCE AND HEAD MOTION LEADS TO LARGE
      !  CHANGES IN THE INCOMING VELOCITY VECTORS.]
      u_r = u_g  !  Use the turbulence in the default rotor position.
!if (t .gt. 30.d0) then
!u_r(1,:,:) = 2.d0*sin(2.d0*PI*0.05*(t-30.d0))
!else
!u_r(1,:,:) = 0.d0
!end if
!u_r(2,:,:) = 0.d0
!u_r(3,:,:) = 0.d0

      !  Transform velocities into rotor coordinates.
      Vh_r = matvecmult(Tg_r,Vh_g,3)
      do icol = 1,Ncols/2
         do irow = 1,Nrows
            Vinf_r(:,irow,icol) = matvecmult(Tg_r,Vinf_g(:,irow,icol),3)
            V0_r(:,irow,icol) = Vinf_r(:,irow,icol) + u_r(:,irow,icol) - Vh_r
!write (41,'(1X,2I4,12F9.3,A)') &
!irow,icol,Vh_r,Vinf_r(:,irow,icol),u_r(:,irow,icol),V0_r(:,irow,icol), &
!'   irow icol Vh_r Vinf_r u_r V0_r'
         end do
      end do

      icolup = Ncols/2
      do icol = (Ncols/2)+1,Ncols
         do irow = 1,Nrows
            k = Nrows*Ncols + 3*(icolup-1)*Nrows + 3*(irow-1)
            Tss_r = transpose(Tsr_s(:,:,irow,icolup))
            !  With this k, x(k) is the array entry before the
            !  X^s component of upwind induced velocity.
            Vi_r = matvecmult(Tss_r,x(k+1:k+3),3)
            V0_r(1,irow,icol) = V0_r(1,irow,icolup) + 2.d0*Vi_r(1)
!  Wrong:
!            V0_r(2,irow,icol) = V0_r(2,irow,icolup)
!            V0_r(3,irow,icol) = V0_r(3,irow,icolup)
            !  Correct results are obtained with:
            V0_r(2,irow,icol) = V0_r(2,irow,icolup) + 2.d0*Vi_r(2)
            V0_r(3,irow,icol) = V0_r(3,irow,icolup) + 2.d0*Vi_r(3)


!if ((irow .eq. 20) .and. (icol .eq. 36)) then
!print '(1X,I4,3ES12.3)',icol,V0_r(1,irow,icolup),Vi_r(1),V0_r(1,irow,icol)
!end if

!write (41,'(1X,3I4,I7,9F9.3,A)') &
!irow,icol,icolup,k,Vi_r,V0_r(:,irow,icolup),V0_r(:,irow,icol), &
!'   irow icol icolup k Vi_r_up V0_r_up V0_r'
! Verified.

         end do
         icolup = icolup - 1
      end do

      call prepVAWTBEMTimestep (Nrows,Ncols,psi0_r,pitchAng,twist0,  &
                                node,profile,                        &
                                Tbr_a,ierr)

      !  Compute velocity of each blade element in rotor coordinates.
      do icol = 1,Ncols
         psi = psi0_r + dble(icol-1)*2.d0*PI/dble(Ncols)
         cpsi = cos(psi)
         spsi = sin(psi)
         do irow = 1,Nrows
            r = sqrt(xs_r(1,irow,icol)**2 + xs_r(2,irow,icol)**2)
            Vb_r(1,irow,icol) = -r*omega*spsi
            Vb_r(2,irow,icol) = r*omega*cpsi
            Vb_r(3,irow,icol) = 0.d0
         end do
      end do

      call VAWTBEMTimestepWithGhostBlades (                               &
                               debug,QSflag0,x,dViqsdt,                   &
                               Nb,Nrows,Ncols,Ndof,                       &
                               Lblade,Ro,xs_r,n_r,psi0_r,omega,           &
                               density,viscosity,chord,Lelem,Zelem,       &
                               NRe,Naoa,Rea,aoaa,Cla,Cda,Cma,             &
                               Tbr_a,Tsr_s,                               &
                               Vb_r,V0_r,                                 &
                               aoafs,aoaz,dClda_a,                        &
                               QSflag,dxdt,aoa,fqs,Viqs,Fb_r,             &
                               ierr)
      if (ierr .ne. 0) then
         return
      end if

      !  Calculate instantaneous rotor loads.
      F_r = 0.d0
      do ib = 1,Nb
         icol = 1 + (ib-1)*Ncols/Nb
         psi = psi0_r + dble(ib-1)*2.d0*PI/dble(Nb)
         cpsi = cos(psi)
         spsi = sin(psi)
         do irow = 1,Nrows
            vec2(1) = profile(1,irow)*cpsi &
                    - profile(3,irow)*spsi
            vec2(2) = profile(1,irow)*spsi &
                    + profile(3,irow)*cpsi
            vec2(3) = profile(2,irow)
            F_r(1:3) = F_r(1:3) + Fb_r(1:3,irow,icol)
            F_r(4:6) = F_r(4:6) + Fb_r(4:6,irow,icol) &
                     + cross(vec2,Fb_r(1:3,irow,icol))
         end do
      end do

      !  Transform back to global coordinates for output.
      F_g(1:3) = matvecmult(Tr_g,F_r(1:3),3)
      F_g(4:6) = matvecmult(Tr_g,F_r(4:6),3)

!write (41,'(1X,F10.4,12F12.1)') t,F_r,F_g
! Check torque.

   end subroutine aeroTimestep













   subroutine VAWTBEMTimestepWithGhostBlades                              &
                              (debug,QSflag0,x,dViqsdt,                   &
                               Nb,Nrows,Ncols,Ndof,                       &
                               Lblade,Ro,xs_r,n_r,psi0,omega,             &
                               density,viscosity,chord,Lelem,Zelem,       &
                               NRe,Naoa,Rea,aoaa,Cla,Cda,Cma,             &
                               Tbr_a,Tsr_s,                               &
                               Vb_r,V0_r,                                 &
                               aoafs,aoaz,dClda_a,                        &
                               QSflag,dxdt,aoa,fqs,Viqs_s,Fb_r,           &
                               ierr)
   !
   !  This subroutine executes a dynamic blade element momentum timestep
   !  for a vertical-axis wind turbine.  A double-multiple streamtube
   !  approach is used.
   !
   !  The difference between this and VAWTBEMTimestep is that here ghost
   !  blades are used to assign aerodynamic loads to each surface element.
   !  These surface element loads are used in computing the induced 
   !  velocity.  The ghost blades are exactly like real blades, except 
   !  there is one associated with each surface element about the 
   !  circumference of the VAWT.  This means that the induced velocity
   !  at each surface element begins to adapt immediately in response to
   !  changes in the incoming flow.  By contrast, using only real blades,
   !  updating the surface forces and induced velocity must wait until
   !  the next time one of the blades passes the element.  This delay can
   !  result in a numerical instability, which is prevented by the ghost
   !  blade approach.  The downside to ghost blades is that the speed
   !  of the calculation is slowed in proportion to the number of blades.
   !
   !  Inputs are fixed variables, and state variables for the current
   !  timestep.  Outputs are loads and such, and the rates of change of
   !  the aerodynamic state variables for the next timestep.  The aero-
   !  dynamic state variables are the separation-point position f and
   !  the induced velocity Vi (as well as the intermediate value Viint).
   !
   !  There are four coordinate systems related to the rotor which are 
   !  important.  One is a rectangular coordinate system in which the X_r
   !  axis is aligned with the rotor-average direction of the instantaneous
   !  remote wind Vinf, and the Z_r axis is the rotor's instantaneous axis
   !  of revolution.  This definition was determined to be the best 
   !  compromise between complexity of programming and accuracy of the
   !  methods.  It is valid when the mean direction of the wind changes on
   !  a timescale that is large in comparison with the period of rotation
   !  of the rotor.  
   !
   !  The second coordinate system is the blade coordinate system.  This
   !  is associated with the rotor coordinate system, but rotated by the
   !  azimuth angle about the Z_r axis.  That is, Z_b and Z_r are the same.
   !
   !  The third coordinate system is associated with the swept surface,
   !  divided into a number of discrete elements.  For each element, the
   !  Z_s axis is normal to the swept surface, and the X_s axis is
   !  orthogonal to the Z_r axis, pointing in the direction of rotation.
   !
   !  There is also an airfoil coordinate system.  The X_a axis points
   !  along the chordline, towards the tail.  The Y_a axis is orthogonal
   !  to the X_a axis, and is pointed towards the low-pressure surface.
   !  The Z_a axis is then spanwise, for a straight blade.
   !
   !  The rotor coordinate system is important in the calculation of total
   !  rotor loads, for application to the top of the support structure.
   !  The surface coordinate system is important in the calculation of
   !  induced velocity and momentum balance.
   !
   !  Inputs:
   !  -------
   !  ----- States ------
   !  QSflag0         : Flag for special dynamic stall behavior.
   !  x               : Aerodynamic state vector, see dxdt below.
   !  -- Other inputs ---
   !  dViqsdt         : Rate of change of quasi-steady induced velocity,
   !                    in surface coordinates.
   !  Nb              : Number of blades.
   !  Nrows           : Number of latitudinal divisions.
   !  Ncols           : Number of longitudinal divisions.
   !  Ndof            : Number of aerodynamic states.
   !  Lblade          : Length of blade.  Used in Prandtl factor.
   !  Ro              : Outer radius of rotor.  Used in calculating dynamic
   !                    inflow time constants.
   !  xs_r            : Surface element centroid coordinates, in the rotor 
   !                    coordinate system.
   !  n_r             : Surface element unit normal vectors, in the rotor 
   !                    coordinate system.
   !  psi0            : Azimuth angle in rotor coordinates of blade 1.  Zero
   !                    azimuth is directly downwind in the mean wind 
   !                    direction.
   !  omega           : Instantaneous rotational speed of the rotor.
   !  density,viscosity
   !  chord           : Chord length of each blade element.
   !  Lelem           : Spanwise length of each blade element.
   !  Zelem           : Spanwise coordinate of each element, measured from the
   !                    center span.
   !  NRe, ..., Cma   : Airfoil data arrays, coefficients as functions of 
   !                    aoa, Re.
   !  Tbr_a           : Transformation matrix which takes an XYZ vector from
   !                    rotor to airfoil (at the blade elements) coordinates.
   !  Tsr_s           : Transformation matrix which takes an XYZ vector from
   !                    rotor to surface element coordinates.
   !  Vb_r            : Velocity of the motion of each blade element, in
   !                    rotor coordinates.
   !  V0_r            : Remote velocity incoming to each surface element, in
   !                    rotor coordinates.  This differs from the reference
   !                    Vinf, because it includes turbulence.
   !  aoafs           : The angle-of-attack at full flow separation, for
   !                    positive (1) and negative (2) angles-of-attack.
   !  aoaz            : The angle-of-attack at zero lift.
   !  dClda_a         : The dCl/da slope when flow is attached, for positive (1)
   !                    and negative (2) angles-of-attack.
   !
   !  Outputs:
   !  --------
   !  ----- States ------
   !  QSflag          : Designates whether the AOA has recently reversed on the
   !                    airfoil, requiring special dynamic stall behavior.
   !  dxdt            : Rate of change of states.  Stored in a 7-by-Nrows-by-Ncols
   !                    vector.  Elements 1 through Nrows*Ncols are dfdt, and
   !                    elements Nrows*Ncols + 1 through 7*Nrows*Ncols are Vi_s.
   !                    dfdt: Rate of change of separation-point position.
   !                    dVidt: Rate of change of induced velocity.  A 6-vector.
   !                    The first 3 are dVi/dt and the next 3 are dViint/dt.  This 
   !                    is kept in surface coordinates because it is most natural
   !                    for the time-delay to be applied in surface coordinates.
   !  -- Other outputs --
   !  aoa, fqs        : Variables for dynamic stall logic.
   !  Viqs_s          : Quasi-steady induced velocity for dynamic inflow.
   !  Fb_r            : Forces and moments on each blade element, in rotor
   !                    coordinates.  These include ghost blades; simply extract
   !                    the array entries associated with real blades.
   !  ierr            : Set to a negative number if an error occurs.
   !  
   !  Local variables:
   !  ----------------
   !  ibrow,ibcol     : Row (latitude) and column (longitude) indices associated
   !                    with the blade coordinate system (azimuth rotated with
   !                    respect to the rotor coordinate system).
   !  isrow,iscol     : Row and column indices associated with surface elements.
   !  Tba_r,Tss_r     : Inverse tranforms from airfoil and surface to rotor 
   !                    coordinates.
   !  psi             : Azimuth angle of a blade element.
   !  Vb_a            : Velocity of the motion of a blade element, in airfoil
   !                    coordinates.
   !  V_r,V_a         : Total flow velocity incoming to a blade element, in
   !                    rotor and airfoil coordinates.
   !  V0_loc_r,       : Remote velocity incoming to a blade element, in rotor
   !  V0_loc_a          and airfoil coordinates.
   !  Vi_loc_a        : Induced velocity at a blade element, in airfoil
   !                    coordinates.
   !  V0_s            : Remote velocity incoming to a surface element, in surface
   !                    coordinates.
   !  Vi_s            : Induced velocity at a surface element, in surface 
   !                    coordinates.
   !  Ro              : External radius of the rotor.
   !  area            : Swept area of a given surface element.
   !  f               : Prandtl factor.
   !  aoarad          : Angle-of-attack in radians.
   !  Vmag            : Velocity magnitude.
   !  Re              : Reynolds number.
   !  taus            : Dynamic stall time constant.
   !  Cl,Clqs         : Dynamic and quasi-steady lift coefficients.
   !  Fl,Fd,M         : Local airfoil forces.
   !  Fs_r            : These are the aerodynamic forces at the surface
   !                    element centroids that are obtained by interpolating
   !                    from adjacent blade elements.  These surface forces
   !                    are knocked down from the raw blade element forces
   !                    by the ratio of real to total (real plus ghost)
   !                    blades, in this case Nb/Ncols.
   !  Vs_r            : Same as Vb_r, but for each surface element centroid,
   !                    rather than each blade, and including only rotational
   !                    velocity, not blade vibration.  Used in calculating
   !                    the Prandtl factor and momentum balance.

      implicit none

      !  Inputs:
      integer,                                      intent(in)    :: Nb,Nrows,Ncols,Ndof, &
                                                                     debug
      integer,          dimension(Nrows,Ncols),     intent(in)    :: QSflag0
      double precision, dimension(Ndof),            intent(in)    :: x
      double precision, dimension(3,Nrows,Ncols),   intent(in)    :: dViqsdt
      
      double precision,                             intent(in)    :: Lblade,Ro,psi0,omega, &
                                                                     density,viscosity
      double precision, dimension(Nrows),           intent(in)    :: chord,Lelem,Zelem
      double precision, dimension(3,Nrows,Ncols),   intent(in)    :: xs_r,n_r,Vb_r,V0_r
      double precision, dimension(3,3,Nrows,Ncols), intent(in)    :: Tbr_a,Tsr_s

      !                             Airfoil properties.
      integer,          dimension(Nrows),           intent(in)    :: NRe
      integer,          dimension(MRe,Nrows),       intent(in)    :: Naoa
      double precision, dimension(MRe,Nrows),       intent(in)    :: Rea
      double precision, dimension(Maoa,MRe,Nrows),  intent(in)    :: aoaa,Cla,Cda,Cma

      double precision, dimension(Nrows,Ncols),     intent(in)    :: aoaz
      double precision, dimension(2,Nrows,Ncols),   intent(in)    :: aoafs,dClda_a

      !  Outputs:
      integer,                                      intent(inout) :: ierr
      integer,          dimension(Nrows,Ncols),     intent(out)   :: QSflag
      double precision, dimension(Ndof),            intent(out)   :: dxdt
      double precision, dimension(Nrows,Ncols),     intent(out)   :: aoa,fqs
      double precision, dimension(3,Nrows,Ncols),   intent(out)   :: Viqs_s
      double precision, dimension(6,Nrows,Ncols),   intent(out)   :: Fb_r

      !  Local variables:
      integer                                                     :: ibrow,ibcol,           &
                                                                     isrow,iscol,           &
                                                                     icol1,icol2,i,j,k

      double precision                                            :: psi,dpsi,psiup,psilow, &
                                                                     area,f,aoarad,         &
                                                                     Vmag,val,taus,Cl,      &
                                                                     Fl,Fd,M,Clqs,R,        &
                                                                     term,sa,ca,comp
      double precision, dimension(3)                              :: Vb_a,V_a,              &
                                                                     V0_loc_r,V0_loc_a,     &
                                                                     Vi_loc_r,Vi_loc_a,     &
                                                                     V0_s,Vi_s,             &
                                                                     Fs_s
      double precision, dimension(6)                              :: Fb_a
      double precision, dimension(3,3)                            :: Tba_r,Tss_r
      double precision, dimension(3,Nrows,Ncols)                  :: Fs_r,Vi_r


      !  Zero forces associated with the swept surface.  These will be 
      !  overwritten based upon updated blade forces.
      Fs_r = 0.d0

      dpsi = 2.d0*PI/dble(Ncols)

      !  Load current induced velocity and separation point position
      !  into a more usable format.
      !  (Vi_r needed especially for interpolateLongitudinal function.)
      k = Nrows*Ncols
      do iscol = 1,Ncols
         do isrow = 1,Nrows
            Tss_r = transpose(Tsr_s(:,:,isrow,iscol))
            Vi_r(:,isrow,iscol) = matvecmult(Tss_r,x(k+1:k+3),3)
            k = k + 3
         end do
      end do

      k = 0
      do ibcol = 1,Ncols

         !  Normalize the azimuth angle of the present blade such
         !  that 0 corresponds to the angle from the Y_r axis
         !  instead of the X_r axis.  This makes the 0 reference
         !  equal to when the blade passes from the downwind to
         !  upwind halves of the swept surface.
         !
         !  The azimuth angle psi is thus associated with the
         !  current blade, out of the Ncols total (real plus
         !  ghost) blades.
         psi = psi0 - 0.5d0*PI + dble(ibcol-1)*dpsi
         if (psi .lt. 0.d0) then
            psi = psi + 2.d0*PI
         else if (psi .ge. 2.d0*PI) then
            psi = psi - 2.d0*PI
         end if

         !  Calculate the surface element column associated with
         !  the blade element column.  That is, the surface element
         !  centerline before which the blade X_b axis immediately
         !  lies.
         iscol = int((psi + 0.5d0*dpsi)/dpsi) + 1
         if (iscol .gt. Ncols) then
            !  Account for the case where 2PI-0.5dpsi < psi < 2PI.
            iscol = iscol - Ncols
         end if
         psiup = (0.5d0 + dble(iscol - 1))*dpsi
         psilow = psiup - dpsi
         if (psilow .lt. 0.d0) then
            psilow = psilow + 2.d0*PI
         end if
!if (ibcol .eq. 1) then
!write (41,'(1X,I4,2ES13.4,I4,2ES13.4,A)') &
!ibcol,psi0,psi,iscol,psiup,psilow,'   ibcol psi0 psi iscol psiup psilow'
!end if
!Verified.
         do ibrow = 1,Nrows

            k = k + 1  !  Used for x(k) = f below, in Oye.

            isrow = ibrow

            !  Calculate the local velocity vector due to blade motion in
            !  airfoil coordinates.
            Vb_a = matvecmult (Tbr_a(:,:,ibrow,ibcol),Vb_r(:,ibrow,ibcol),3)
!write (41,'(1X,2I4,6ES13.4,A)') &
!ibrow,ibcol,Vb_r(:,ibrow,ibcol),Vb_a,'   ibrow ibcol Vb_r Vb_a'
!Verified.
            !  Interpolate the local incoming windspeed at the blade element
            !  location.
            call interpolateLongitudinal (psi,psilow,dpsi,isrow,iscol,Ncols, &
                                          V0_r,V0_loc_r)

            !  Same with induced velocity.  First transform from surface
            !  to blade coordinates.
            call interpolateLongitudinal (psi,psilow,dpsi,isrow,iscol,Ncols, &
                                          Vi_r,Vi_loc_r)
!if ((ibrow .eq. 10) .and. (ibcol .eq. 1)) then
!write (41,'(1X,2I4,7ES13.4,A)') &
!ibrow,ibcol,psi*180.d0/PI,V0_loc_r,Vi_loc_r,'   ibrow ibcol psi V0_loc_r Vi_loc_r'
!end if
!Check under nonuniform flow.
            !  Transform to airfoil coordinates.
            V0_loc_a = matvecmult (Tbr_a(:,:,ibrow,ibcol),V0_loc_r,3)
            Vi_loc_a = matvecmult (Tbr_a(:,:,ibrow,ibcol),Vi_loc_r,3)

            !  Calculate the total velocity vector.  This is the vector sum
            !  of the incoming, induced, and blade motion velocities.
            V_a = V0_loc_a + Vi_loc_a - Vb_a
            Vmag = sqrt(V_a(1)**2 + V_a(2)**2)  !  Disregard local spanwise
                                                !  velocity for Re, Cl.
!if ((ibrow .eq. 10) .and. (ibcol .eq. 1)) then
!write (41,'(1X,2I4,13ES13.4,A)') &
!ibrow,ibcol,psi*180.d0/PI,V0_loc_a,Vi_loc_a,Vb_a,V_a,'   ibrow ibcol psi V0_loc_a Vi_loc_a Vb_a V_a'
!end if
!Verified.

            !  Calculate airfoil forces using the local velocity vector.
            call steadyAeroSectionLoads (NRe(ibrow),Naoa(:,ibrow),             &
                                         density,viscosity,chord(ibrow),       &
                                         V_a,Rea(:,ibrow),aoaa(:,:,ibrow),     &
                                         Cla(:,:,ibrow),Cda(:,:,ibrow),        &
                                         Cma(:,:,ibrow),aoafs(1,ibrow,ibcol),  &
                                         Fl,Fd,M,aoa(ibrow,ibcol),ierr)
            if (ierr .ne. 0) then
               print *,'Error calculating aero section loads.'
               return
            end if
!if ((ibrow .eq. 20) .and. (ibcol .eq. 19)) then
!print '(1X,7ES12.3)', &
!aoa(ibrow,ibcol),V0_loc_a(1),V0_loc_a(2),Vi_loc_a(1),Vi_loc_a(2),Vb_a(1),Vb_a(2)
!end if
            !  Calculate quasi-steady lift coefficient.
            val = 0.5d0*density*chord(ibrow)
            Clqs = Fl/(val*max(Vmag,SMALLV)**2)
!write (41,'(1X,2I4,5ES13.4,A)')                            &
!ibrow,ibcol,psi,aoa(ibrow,ibcol),                          &
!Fl/(val*max(Vmag,SMALLV)**2),Fd/(val*max(Vmag,SMALLV)**2), &
!M/(val*chord(ibrow)*max(Vmag,SMALLV)**2),'   ibrow ibcol psi aoa Clqs Cd Cm'
!Verified.
            !  Update the time constant (and prevent divide-by-zero).
            taus = 3.7d0*chord(ibrow)/max(Vmag,SMALLV)

            !  Øye dynamic stall model.  Note: x(k) contains f for the current
            !  timestep.  Importantly, f may need to be overwritten.  This is
            !  done according to QSflag0, in the main time integration loop.
            !  But here in oye is where QSflag is updated for the next 
            !  timestep.
            call oye (QSflag0(ibrow,ibcol),x(k),aoa(ibrow,ibcol),       &
                      Clqs,aoaz(ibrow,ibcol),                           &
                      dClda_a(:,ibrow,ibcol),aoafs(:,ibrow,ibcol),taus, &
                      Cl,fqs(ibrow,ibcol),QSflag(ibrow,ibcol),dxdt(k))
!if ((ibrow .eq. 10) .and. (ibcol .eq. 18)) then
!print '(1X,3ES12.3,2I4)', &
!x(k),fqs(ibrow,ibcol),aoa(ibrow,ibcol),QSflag0(ibrow,ibcol),QSflag(ibrow,ibcol)
!end if

!if ((ibrow .eq. 16) .and. (ibcol .eq. 1)) then
!write (41,                                                                &
!'(1X,2I4,I6,F10.3,2I3,2F8.4,F9.2,2F10.3,2F8.4,2F10.3,3F9.4,A)')           &
!ibrow,ibcol,k,psi*180.d0/PI,QSflag0(ibrow,ibcol),QSflag(ibrow,ibcol),     &
!x(k),fqs(ibrow,ibcol),dxdt(k),aoa(ibrow,ibcol),aoaz(ibrow,ibcol),         &
!dClda_a(1,ibrow,ibcol),dClda_a(2,ibrow,ibcol),                            &
!aoafs(1,ibrow,ibcol),aoafs(2,ibrow,ibcol),                                &
!taus,Clqs,Cl,                                                             &
!'   ibrow ibcol k psi QSflag0 QSflag f fqs dfdt aoa aoaz dClda_a+,- aoafs+,- taus Clqs Cl'
!end if
! Verified.

!if ((debug .eq. 1) .and. (ibrow .eq. 10) .and. (ibcol .eq. 1)) then
!write (42,                                                                &
!'(1X,2I4,I6,F10.3,2I3,2F10.4,F9.2,F10.3,3F10.4,A)')                       &
!ibrow,ibcol,k,psi*180.d0/PI,QSflag0(ibrow,ibcol),QSflag(ibrow,ibcol),     &
!x(k),fqs(ibrow,ibcol),dxdt(k),aoa(ibrow,ibcol),                           &
!Clqs,Cl,Fd/(val*max(Vmag,SMALLV)**2),                                     &
!'   ibrow ibcol k psi QSflag0 QSflag f fqs dfdt aoa Clqs Cl Cd'
!end if

            !  Overwrite the lift force with the dynamic value.
            Fl = Cl*val*max(Vmag,SMALLV)**2

            !  Transform to rotor coordinates.
            Tba_r = transpose(Tbr_a(:,:,ibrow,ibcol))
            aoarad = aoa(ibrow,ibcol)*PI/180.d0
            sa = sin(aoarad)
            ca = cos(aoarad)
            Fb_a(1) = (-Fl*sa + Fd*ca)*Lelem(ibrow)
            Fb_a(2) = (Fl*ca + Fd*sa)*Lelem(ibrow)
            Fb_a(3) = 0.d0
            Fb_a(4) = 0.d0
            Fb_a(5) = 0.d0
            Fb_a(6) = -M*Lelem(ibrow)
            Fb_r(1:3,ibrow,ibcol) = matvecmult(Tba_r,Fb_a(1:3),3)
            Fb_r(4:6,ibrow,ibcol) = matvecmult(Tba_r,Fb_a(4:6),3)

!write (41,'(1X,2I4,12ES13.4,A)') &
!ibrow,ibcol,Fb_a,Fb_r(:,ibrow,ibcol),'   ibrow ibcol Fb_a Fb_r'
!Verified.

            !  Add contribution of airfoil force to adjacent surface element
            !  centroids.  (Note that Fs_r was zeroed at the start of this
            !  subroutine.)  Scale the loads by the ratio of the number of
            !  real blades to the number of real plus ghost blades.  This
            !  is necessary because the blade forces were calculated as if
            !  there were only Nb blades, yet they were calculated at Ncols
            !  locations about the azimuth.
            if (psi .gt. 2.d0*PI - 0.5d0*dpsi) then
               icol1 = Ncols
               icol2 = 1
               val = (psi - 2.d0*PI + 0.5d0*dpsi)/dpsi
            else if (psi .lt. 0.5d0*dpsi) then
               icol1 = Ncols
               icol2 = 1
               val = (psi + 0.5d0*dpsi)/dpsi
            else
               icol1 = iscol - 1
               icol2 = iscol
               val = (psi - psilow)/dpsi
            end if
            !  Note that Fs_r has no moment associated with it.  Only the
            !  forces, not the moments, will be relevant for calculating
            !  the induction.
            Fs_r(:,isrow,icol2) = Fs_r(:,isrow,icol2) &
                                + val*Fb_r(1:3,ibrow,ibcol)*dble(Nb)/dble(Ncols)
            Fs_r(:,isrow,icol1) = Fs_r(:,isrow,icol1) &
                                + (1.d0 - val)*Fb_r(1:3,ibrow,ibcol)*dble(Nb)/dble(Ncols)

!write (41,'(1X,5I4,7ES13.4,A)')                                              &
!ibrow,ibcol,isrow,icol1,icol2,val,Fb_r(1:3,ibrow,ibcol),Fs_r(:,isrow,icol2), &
!' ibrow ibcol isrow icol1 icol2 val Fb_r Fs_r'
!Verified.

         end do  !  Rows.

      end do  !  Columns (blades).

      !  Apply momentum balance to update Vi.  Work in surface coordinates,
      !  which are most convenient for this portion of the calculation.
      j = Nrows*Ncols
      k = 4*Nrows*Ncols
      do iscol = 1,Ncols
         do isrow = 1,Nrows
            Fs_s = matvecmult (Tsr_s(:,:,isrow,iscol),Fs_r(:,isrow,iscol),3)
            V0_s = matvecmult (Tsr_s(:,:,isrow,iscol),V0_r(:,isrow,iscol),3)
            Vi_s = matvecmult (Tsr_s(:,:,isrow,iscol),Vi_r(:,isrow,iscol),3)
!write (41,'(1X,2I4,9ES13.4,A)') &
!isrow,iscol,V0_s,Vi_s,Fs_s,'   isrow iscol V0_s Vi_s Fs_s'
!Verified (except for Vi_s).

            !  Calculate the blade velocity at each surface element, neglecting
            !  local blade motion (vibration).
            r = sqrt(xs_r(1,isrow,iscol)**2 + xs_r(2,isrow,iscol)**2)
            !  ang differs from psi in that it is referred to rotor coordinates,
            !  same as psi0.
!            ang = atan2(xs_r(2,isrow,iscol),xs_r(1,isrow,iscol))
!            Vs_r(1) = -r*omega*sin(ang)
!            Vs_r(2) = r*omega*cos(ang)
!            Vs_r(3) = 0.d0

            !  Calculate the Prandtl factor.
!            V_r = V0_r(:,isrow,iscol) + Vi_r(:,isrow,iscol) - Vs_r
!            call prandtl (Nb,0.5d0*Lblade,Zelem(isrow), &
!                          V_r,n_r(:,isrow,iscol),f)
            f = 1.0d0  !  Deactivate Prandtl for numerical stability at the
                       !  outer rows, near the attachment point.
!write (41,'(1X,2I4,9ES13.4,A)')                                 &
!isrow,iscol,V_r,n_r(:,isrow,iscol),0.5d0*Lblade,Zelem(isrow),f, &
!'   isrow iscol V_r n_r Lblade/2 Zelem f'
!Verified.
            !  Update the quasi-steady induced velocity with momentum-balance.
            !  Work in surface coordinates.  The reason is that the surface-normal
            !  component of the induced velocity should dominate, because this
            !  is the direction in which the majority of lift acts.  Thus we can
            !  apply momentum-balance and dynamic inflow to each component
            !  independently, since the tangential components are of secondary
            !  importance.
            area = Lelem(isrow)*r*dpsi
            term = 2.d0*density*area*f*max(abs(V0_s(3) + f*Vi_s(3)),SMALLV)

            !  Momentum theory becomes invalid when the induction factor a
            !  exceeds 0.5, in the case of flow that is aligned with the
            !  surface element.  Under yawed conditions, such as elements
            !  towards the outside of a VAWT, the theory is not entirely
            !  clear on what is physically correct.  (In reality, the
            !  results from the momentum equation do not obey global
            !  conservation of mass, because the downwind streamtubes cross
            !  into adjacent streamtubes.)  It was observed, though, that
            !  implementing the raw momentum equation can result in 
            !  large induced velocities and physically unrealistic loads
            !  in some elements, especially those near the blade attachment
            !  points.
            !
            !  It is proposed to limit (Vi_s)Z <= 0.5 (V0_s)Z and
            !  (Vi_s)Z >= 0, for upwind elements.  (Signs are reversed
            !  for downwind elements because of the definition of the
            !  surface coordinate system.)
            do i = 1,3
               Viqs_s(i,isrow,iscol) = -Fs_s(i)/term
               comp = sign(max(abs(V0_s(i)),SMALLV),V0_s(i))
               val = Viqs_s(i,isrow,iscol)/comp
               if (val .lt. -0.5d0) then
                  val = -0.5d0  !  Do a is to i.e.  The Nodal the the the the the are not allow flow to "backwind" downstream.
               else if (val .gt. 0.d0) then
                  val = 0.d0    !  Do not allow blades to accelerate flow.
               end if
!if ((isrow .eq. 10) .and. (iscol .eq. 12)) then
!print '(1X,I4,5ES12.3)',i,Viqs_s(i,isrow,iscol),V0_s(i),comp,val,val*comp
!end if
               Viqs_s(i,isrow,iscol) = val*comp
            end do

!  Original logic applied truncation only to Z component.  Led
!  to numerical instabilty when integrating X and Y components.
!  The new version above gives about the same global load 
!  results and is more stable.
!            val = 0.5d0*abs(V0_s(3))
!            if (iscol .le. Ncols/2) then  !  Normal points upwind.
!               if (Viqs_s(3,isrow,iscol) .ge. val) then
!                  Viqs_s(3,isrow,iscol) = val
!               else if (Viqs_s(3,isrow,iscol) .lt. 0.d0) then
!                  Viqs_s(3,isrow,iscol) = 0.d0
!               end if
!            else  !  Normal points downwind.
!               if (-Viqs_s(3,isrow,iscol) .ge. val) then
!                  Viqs_s(3,isrow,iscol) = -val
!               else if (-Viqs_s(3,isrow,iscol) .lt. 0.d0) then
!                  Viqs_s(3,isrow,iscol) = 0.d0
!               end if
!            end if
!            Viqs_s(1,isrow,iscol) = -Fs_s(1)/term
!            Viqs_s(2,isrow,iscol) = -Fs_s(2)/term

            call inflowDelay (Viqs_s(:,isrow,iscol),dViqsdt(:,isrow,iscol), &
                              x(j+1:j+3),x(k+1:k+3),Ro,V0_s,                &
                              dxdt(j+1:j+3),dxdt(k+1:k+3))
!if ((isrow .eq. 10) .and. (iscol .eq. 9)) then
!write (41,'(1X,2I4,2I7,21F10.4,A)')                           &
!isrow,iscol,j,k,Viqs_s(:,isrow,iscol),dViqsdt(:,isrow,iscol), &
!x(j+1:j+3),x(k+1:k+3),V0_s,                                   &
!dxdt(j+1:j+3),dxdt(k+1:k+3),                                  &
!'   isrow iscol j k Viqs_s dViqs/dt Vi Viint V0_s dVi/dt dViint/dt'
!end if
! Check under timestepping.

            j = j + 3
            k = k + 3

         end do
      end do

   end subroutine VAWTBEMTimestepWithGhostBlades













   subroutine prepVAWTBEMTimestep (Nrows,Nb,psi0,pitchAng,twist0,  &
                                   node,profile,                   &
                                   Tbr_a,ierr)
   !  This subroutine prepares BEM analysis for a timestep.  The
   !  transform matrices between blade airfoil and rotor coordinates
   !  are updated based upon the azimuth angle and blade pitch.
   !  
   !  Inputs:
   !  -------
   !  Nrows,Nb         : The number of latitudinal divisions in the swept
   !                     surface, and number of blades.
   !  psi0             : The azimuth angle of the first blade, relative to 
   !                     the rotor coordinate system (0 degrees directly 
   !                     downwind).
   !  pitchAng         : The pitch angle of the blade at each blade element
   !                     (real and ghost blades).
   !  twist0           : Baseline (zero-pitch) blade twist.
   !  node             : Nodal coordinates of rotated profile, (r,z,offset).
   !  profile          : Blade element coords of rotated profile, (r,z,offset).
   !
   !  Outputs:
   !  --------
   !  Tbr_a            : The transform from rotor to airfoil coordinates,
   !                     for each element of the real and ghost blades.
   !
   !  Local variables:
   !  ----------------
   !

      implicit none

      !  Inputs:
      integer,                                      intent(in)  :: Nrows,Nb
      double precision,                             intent(in)  :: psi0
      double precision, dimension(Nrows),           intent(in)  :: twist0
      double precision, dimension(3,Nrows),         intent(in)  :: profile
      double precision, dimension(3,Nrows+1),       intent(in)  :: node
      double precision, dimension(Nrows,Nb),        intent(in)  :: pitchAng

      !  Outputs:
      double precision, dimension(3,3,Nrows,Nb),    intent(out) :: Tbr_a

      integer, intent(inout) :: ierr

      !  Local variables:
      integer :: irow,ib
      double precision :: cpsi,spsi,psi,dp,magns,ca,sa,ct,st,cp,sp, &
                          azi,tilt,pitch
      double precision, dimension(3) :: v1,v2,profrot,xb_r,nbv
      double precision, dimension(3,2) :: nodrot

      double precision, dimension(3,3,4) :: T

      !  To get from the rotor coordinate system to the airfoil 
      !  coordinate system associated with a blade element:
      !  (1)  Define Za parallel to Zr and Ya parallel to Xr and
      !       Xa opposite Yr.
      !  (2)  Rotate about Zr by the pitch angle.
      !  (3)  Rotate about Yr by the tilt angle.
      !  (4)  Rotate about Zr by the azimuth angle.
      !
      !  NOTE THAT THIS LOGIC WILL NOT WORK FOR PITCH ANGLE
      !  CHANGES IF THE BLADE IS TILTED OR TAPERED.
      do ib = 1,Nb
         psi = psi0 + dble(ib-1)*2.d0*PI/dble(Nb)
         cpsi = cos(psi)
         spsi = sin(psi)
         do irow = 1,Nrows
            nodrot(1,1) = node(1,irow)*cpsi &
                        - node(3,irow)*spsi
            nodrot(2,1) = node(1,irow)*spsi &
                        + node(3,irow)*cpsi
            nodrot(3,1) = node(2,irow)
            nodrot(1,2) = node(1,(irow+1))*cpsi &
                        - node(3,(irow+1))*spsi
            nodrot(2,2) = node(1,(irow+1))*spsi &
                        + node(3,(irow+1))*cpsi
            nodrot(3,2) = node(2,(irow+1))
            profrot(1) = profile(1,irow)*cpsi &
                       - profile(3,irow)*spsi
            profrot(2) = profile(1,irow)*spsi &
                       + profile(3,irow)*cpsi
            profrot(3) = profile(2,irow)
            v1 = nodrot(:,2) - nodrot(:,1)
            v2(1) = profrot(2)
            v2(2) = -profrot(1)
            v2(3) = 0.d0
 
            nbv = cross(v1,v2)

            !  Verify the sign.
            dp = dot(nbv,profrot)
            if (dp .lt. 0.d0) then
               nbv = -nbv
            end if
            magns = magnitude(nbv)
            if (magns .lt. 1.d-6) then
               ierr = -1
               call logError(53)
               return
            else
               nbv = nbv/magns
            end if

            xb_r = profrot

            azi = atan2(xb_r(2),xb_r(1))
            tilt = atan2(-nbv(3),sqrt(nbv(1)**2 + nbv(2)**2))
            pitch = pitchAng(irow,ib) + twist0(irow)

            ca = cos(azi)
            sa = sin(azi)
            ct = cos(tilt)
            st = sin(tilt)
            cp = cos(pitch)
            sp = sin(pitch)

            T(1,1,1) = 0.d0
            T(1,2,1) = -1.d0
            T(1,3,1) = 0.d0
            T(2,1,1) = 1.d0
            T(2,2,1) = 0.d0
            T(2,3,1) = 0.d0
            T(3,1,1) = 0.d0
            T(3,2,1) = 0.d0
            T(3,3,1) = 1.d0

            T(1,1,2) = cp
            T(1,2,2) = sp
            T(1,3,2) = 0.d0
            T(2,1,2) = -sp
            T(2,2,2) = cp
            T(2,3,2) = 0.d0
            T(3,1,2) = 0.d0
            T(3,2,2) = 0.d0
            T(3,3,2) = 1.d0

            T(1,1,3) = ct
            T(1,2,3) = 0.d0
            T(1,3,3) = -st
            T(2,1,3) = 0.d0
            T(2,2,3) = 1.d0
            T(2,3,3) = 0.d0
            T(3,1,3) = st
            T(3,2,3) = 0.d0
            T(3,3,3) = ct

            T(1,1,4) = ca
            T(1,2,4) = sa
            T(1,3,4) = 0.d0
            T(2,1,4) = -sa
            T(2,2,4) = ca
            T(2,3,4) = 0.d0
            T(3,1,4) = 0.d0
            T(3,2,4) = 0.d0
            T(3,3,4) = 1.d0

            Tbr_a(:,:,irow,ib) = matmatmult(T(:,:,1),T(:,:,2),3,3)
            Tbr_a(:,:,irow,ib) = matmatmult(Tbr_a(:,:,irow,ib),T(:,:,3),3,3)
            Tbr_a(:,:,irow,ib) = matmatmult(Tbr_a(:,:,irow,ib),T(:,:,4),3,3)

         end do

      end do


   end subroutine prepVAWTBEMTimestep













   subroutine interpolateLongitudinal (psi,psilow,dpsi,isrow,iscol,Ncols,Vin,Vout)
   !  This subroutine interpolates a vector between values given at surface
   !  element centroids.  Some logic is required to account for the jump in
   !  the index in passing from the downwind to upwind side of the rotor.
   !
   !  Note that magnitude and direction are interpolated separately.  This
   !  provides more accurate results for vectors that are misaligned.

      integer,                            intent(in)   :: isrow,iscol,Ncols
      double precision,                   intent(in)   :: psi,psilow,dpsi
      double precision, dimension(:,:,:), intent(in)   :: Vin

      double precision, dimension(:),     intent(out)  :: Vout

      integer                                          :: icol1,icol2
      double precision                                 :: mag1,mag2,mag,magdir,val
      double precision, dimension(3)                   :: dir

      if (psi .le. 0.5d0*dpsi) then
         icol1 = Ncols
         icol2 = 1
         val = psi + 0.5d0*dpsi
      else if (psi .gt. 2.d0*PI - 0.5d0*dpsi) then
         icol1 = Ncols
         icol2 = 1
         val = psi - 2.d0*PI + 0.5d0*dpsi
      else
         icol1 = iscol-1
         icol2 = iscol
         val = psi - psilow
      end if

      mag1 = sqrt(Vin(1,isrow,icol1)**2 &
           +      Vin(2,isrow,icol1)**2 &
           +      Vin(3,isrow,icol1)**2)
      mag2 = sqrt(Vin(1,isrow,icol2)**2 &
           +      Vin(2,isrow,icol2)**2 &
           +      Vin(3,isrow,icol2)**2)
      mag = (val/dpsi)*(mag2 - mag1) + mag1

      dir = (val/dpsi)                                &
          * (Vin(:,isrow,icol2) - Vin(:,isrow,icol1)) &
          + Vin(:,isrow,icol1)

      magdir = sqrt(dir(1)**2 + dir(2)**2 + dir(3)**2)

      if (magdir .gt. SMALLV) then
         Vout = dir*mag/magdir
      else
         !  Protect against divide-by-zero by using linear interpolation.
         Vout = dir
      end if

   end subroutine interpolateLongitudinal











   subroutine oye (QSflag0,f,aoa,Clqs,aoazi,dClda_ai,aoafsi,taus, &
                   Cl,fqs,QSflag,dfdt)
   !  This subroutine implements one timestep of the Øye 
   !  dynamic stall model.  Ref: Hansen 2008 p 96.
   !
   !  Inputs:
   !  -------
   !  QSflag0         : if AOA has switched signs, special behavior
   !                    may be required.  Must be initialized to
   !                    0 on the first timestep.
   !  f               : separation point position.
   !  aoa             : the instantaneous angle-of-attack.
   !  Clqs            : the quasi-steady lift coefficient.
   !  aoazi           : the zero-lift AOA.
   !  dClda_ai        : the slope of the linear portion of the curve.
   !  aoafsi          : AOA at full flow separation.
   !  taus            : time constant.
   !
   !
   !  Outputs:
   !  --------
   !  Cl              : the lift coefficient.
   !  fqs             : quasi-steady separation point position.
   !  QSflag          : if AOA has switched signs, special behavior
   !                    may be required.  Must be initialized to
   !                    0 on the first timestep.
   !  dfdt            : rate of change of state variable f.
   !
   !  Local variables:
   !  ----------------
   !
      
      implicit none

      !  Inputs:
      integer, intent(in) :: QSflag0
      double precision, intent(in) :: aoa,Clqs,aoazi,taus,f
      double precision, dimension(:), intent(in) :: dClda_ai,aoafsi

      !  Outputs:
      double precision, intent(out) :: Cl,fqs

      integer, intent(out) :: QSflag
      double precision, intent(out) :: dfdt

      !  Local variables:
      double precision :: slope,fcalc

      if (aoa .ge. aoazi) then
         slope = dClda_ai(1)
      else
         slope = dClda_ai(2)
      end if

      !  Convert the equivalent AOA into a physical
      !  separation-point location.
      if (aoa .gt. aoafsi(1)) then
         fqs = 1.d0
      else if (aoa .lt. aoafsi(2)) then
         fqs = -1.d0
      else if (abs(aoa - aoazi) .lt. 2.d0) then
         fqs = 0.d0
      else
         fqs = sign(1.d0,Clqs) &
             * (1.d0 - (2.d0*sqrt(Clqs/(slope*(aoa - aoazi))) - 1.d0)**2)
      end if

      dfdt = (fqs - f)/taus

      !  Let the sign on Cl follow the instantaneous AOA.  But special
      !  steps must be taken in the event that aoa and f have opposite
      !  signs.  The reason is that a very fast sign reversal of AOA 
      !  (like the case of a VAWT in high winds) should result in
      !  flow separation on the new low-pressure surface as vorticity
      !  starts up.  Without this special consideration, the value of
      !  f has to traverse near zero, which can result in anomalously
      !  high lift coefficients.
      ! 
      !  The way it works is this.  If it is detected that the
      !  instantaneous AOA is of opposite sign to f, then the 
      !  lift coefficient is made to follow its quasi-steady value
      !  until f "rounds the corner" at the trailing-edge.  (This
      !  time delay hopefully allows the new AOA to max out.)  Then,
      !  f jumps to fqs associated with the new AOA.  In the meantime,
      !  the quasi-steady coefficients are used.      
      !
      !  Note that under normal conditions, when f is oscillating
      !  about 0 with small changes in angle-of-attack, the fact
      !  that it jumps to fqs is insignificant, because fqs will be
      !  near zero too.
      fcalc = f
      if ((aoa .gt. aoazi) .and. (f .lt. 0.d0)) then
         QSflag = -1
      else if ((aoa .lt. aoazi) .and. (f .gt. 0.d0)) then
         QSflag = 1
      else if ((QSflag0 .eq. -1) .and. (aoa .lt. aoazi) .and. (f .lt. 0.d0)) then
         !  False alarm.
         QSflag = 0
      else if ((QSflag0 .eq. 1) .and. (aoa .gt. aoazi) .and. (f .gt. 0.d0)) then
         !  False alarm.
         QSflag = 0
      else if ((QSflag0 .eq. -1) .and. (aoa .gt. aoazi) .and. (f .gt. 0.d0)) then         
         !  This means that f has "rounded the corner" at the trailing
         !  edge.  Return to normal behavior.
         QSflag = 0
         fcalc = fqs
         !  f = fqs assignment done in the main integration routine.
      else if ((QSflag0 .eq. 1) .and. (aoa .lt. aoazi) .and. (f .lt. 0.d0)) then         
         !  This means that f has "rounded the corner" at the trailing
         !  edge.  Return to normal behavior.
         QSflag = 0
         fcalc = fqs
         !  f = fqs assignment done in the main integration routine.
      else
         QSflag = QSflag0
      end if

      !  Truncate fcalc to prevent attempting to take the square
      !  root of a negative number.
      fcalc = max(min(f,1.d0),-1.d0)

      if (QSflag0 .ne. 0) then
         Cl = Clqs
      else if ((aoa .le. aoafsi(1)) .and. (aoa .ge. aoafsi(2))) then
         !  Approach the Kirchoff-based curve.
         Cl = 0.25d0*slope*(aoa - aoazi)*(1.d0 + sqrt(1.d0 - abs(fcalc)))**2
      else
         !  In extreme cases (typically a retreating blade in high
         !  winds), |f| can be less than 1 while aoa is above 90 (or
         !  below -90).  In such a case, just set Cl = Clqs.  It will
         !  be continuous regardless, because as aoa crosses (about) 
         !  90 Clqs and Cl will both cross zero at the same time.
         if (Clqs*aoa .lt. 0.d0) then
            !  The signs on Clqs and aoa are opposite: this occurs only
            !  in the very deep-stall region above 90 or below -90 
            !  (assuming that aoa is normalized -180 to +180).
            Cl = Clqs
         else
            !  Approach the quasi-steady curve.
            Cl = Clqs*(1.d0 + sqrt(1.d0 - abs(fcalc)))**2
         end if
      end if

   end subroutine oye











   subroutine inflowDelay (Viqs,dViqsdt,Vi,Viint,Ro,Vinf,dVidt,dViintdt)
   !
   !  Inputs:
   !  -------
   !  Viqs           : Quasi-steady induced velocity for the current timestep.
   !  dViqsdt        : Rate of change of quasi-steady induced velocity.
   !  Vi             : Induced velocity from the previous timestep.
   !  Viint          : Intermediate value of induced velocity from the 
   !                   previous timestep.
   !  Ro             : outer radius of the rotor at the mid-height.
   !  Vinf           : remote incoming velocity.
   !  
   !  Outputs:
   !  --------
   !  dVidt          : time derivative of induced velocity.
   !  dViintdt       :  "      "          intermediate induced velocity.
   !  
   !  Local variables:
   !  ----------------

      implicit none

      !  Inputs:
      double precision,               intent(in)    :: Ro
      double precision, dimension(3), intent(in)    :: Vinf,Viqs,dViqsdt,Vi,Viint

      double precision, dimension(3), intent(out)   :: dVidt,dViintdt

      !  Local variables:
      integer                                       :: j
      double precision                              :: T1,T2,a,val

      !  Calculate time constants.
      val = max(sqrt(Vinf(1)**2 + Vinf(2)**2 + Vinf(3)**2),SMALLV)
      a = min(sqrt(Viqs(1)**2 + Viqs(2)**2 + Viqs(3)**2)/val, 0.5d0)
      T1 = (1.1d0/(1.d0 - 1.3d0*a))*(Ro/val)
      !  Original: T2 = (0.39d0 - 0.26d0*(r(iazi)/Ro)**2)*T1,
      !  but let r/Ro = 0.7 for a typical value that might apply
      !  for a VAWT.
      T2 = 0.263d0*T1

      do j = 1,3  !  X,Y,Z components
         dViintdt(j) = (Viqs(j) + 0.6*T1*dViqsdt(j) - Viint(j))/T1
         dVidt(j) = (Viint(j) - Vi(j))/T2
      end do

   end subroutine inflowDelay





   
            




   subroutine prandtl (Nb,Ro,r,V_r,n,f)

      implicit none

      integer, intent(in) :: Nb
      double precision, intent(in) :: Ro,r
      double precision, dimension(:), intent(in) :: V_r,n
      double precision, intent(out) :: f

      double precision :: sinphi,Vmag

      Vmag = magnitude(V_r)

      if (Vmag .gt. SMALLV) then
         !  Compute the angle that the flow relative to the blade 
         !  (including rotational velocity) makes with the swept 
         !  surface.
         sinphi = abs(dot(V_r,n)/Vmag)
         if (sinphi .gt. 1.d-6) then
            f = (2.d0/pi)*acos(exp(-dble(Nb)*(Ro - r)/(2.d0*r*sinphi)))
         else
            f = 1.0d0
         end if
      else
         f = 1.d0
      end if

      if (f .lt. 0.1d0) then
         !  Truncate for numerical stability.
         f = 0.1d0
      end if

   end subroutine prandtl









   subroutine readCoefficientLibrary (fid,fname,airfoilName,  &
                                      NElems,ar,              &
                                      NRe,Naoa,               &
                                      Re,aoa,Cl,Cd,Cm,        &
                                      ierr)
   ! Scans a library file of a given filename for a set of data associated 
   ! with a given airfoil ID.  The format of the library file must be:
   ! 
   ! Airfoil_ID , N_Re
   ! Re_1 , N_aoa
   ! aoa_1 , Cl_1 , Cd_1 , Cm_1
   ! ...
   ! aoa_N_aoa , Cl_N_aoa , Cd_N_aoa , Cm_N_aoa
   ! Re_2 , N_aoa_2
   ! ...
   ! Airfoil_ID_2 , N_Re_2
   ! ...
   !
   ! Inputs:
   ! -------
   ! fid           : Airfoil coefficient file ID number.
   ! fname         : Airfoil coefficient path and filename.
   ! airfoilName   : ID of the airfoil.
   ! NElems        : Number of airfoil sections that need look-up.
   ! ar            : aspect ratio of the blade.
   !
   ! Outputs:
   ! --------
   ! Reynolds number Re, angle-of-attack aoa, and lift, drag, and
   ! moment coefficients Cl, Cd, and Cm.
   ! Number of entries NRe, and Naoa for each Re.
   ! ierr          : -1 if terminal error.
   !
   ! Local Variables:
   ! ----------------
   !

      implicit none
         
      ! Inputs:
      integer,                                intent(in)  :: fid
      integer,                                intent(in)  :: Nelems
      double precision,                       intent(in)  :: ar
      character(len=MAXSTRLEN), dimension(:), intent(in)  :: airfoilName
      character(len=*),                       intent(in)  :: fname

      ! Outputs:
      integer,          dimension(:),         intent(out) :: NRe
      integer,          dimension(:,:),       intent(out) :: Naoa
      double precision, dimension(:,:),       intent(out) :: Re
      double precision, dimension(:,:,:),     intent(out) :: aoa, Cl, Cd, Cm
      integer,                              intent(inout) :: ierr

      ! Local variables: 
      character(len=MAXSTRLEN)                            :: str
      integer                                             :: iElem,iRe,iaoa, &
                                                             fstat,Nlow,Nhigh
      double precision, dimension(Maoa)                   :: tempaoa,tempCl,tempCd,tempCm
      logical                                             :: found

      integer,          dimension(Nelems)                 :: extendTable
      double precision, dimension(Nelems)                 :: rnose_c
      double precision, dimension(2,Nelems)               :: atail,anose
       
      !  Initialize output arrays.  This is necessary to avoid the compiler
      !  Bombing when undefined elements are passed around.
      NRe = 0
      Naoa = 0
      Re = 0.d0
      aoa = 0.d0
      Cl = 0.d0
      Cd = 0.d0
      Cm = 0.d0

      ! Open the requested library file.
      open (UNIT=fid, FILE=fname, STATUS='OLD', IOSTAT=fstat, &
            ACCESS='SEQUENTIAL', FORM='FORMATTED',            &
            ACTION='READ', POSITION='REWIND', DELIM='NONE')
      if (fstat .gt. 0) then
         call logError(10)
         ierr = -1
         return
      end if

      do iElem = 1,Nelems

         ! Scan the file until the appropriate airfoil ID is located.
         rewind(fid)
         found = .false.

         do while ((fstat .eq. 0) .and. (.not. found))

            read (UNIT=fid, FMT=*, IOSTAT=fstat) str

            if (trim(str) .eq. trim(airfoilName(iElem))) then
               found = .true.
            end if

         end do

         if (fstat .gt. 0) then
            call logError(11)
            ierr = -1
            return
         else if (.not. found) then
            call logError(12)
            ierr = -1
            return
         end if 

         ! The next line indicates whether the table should be extended
         ! to the deep stall regime with approximate equations,
         ! and, if so, contains necessary parameters.
         read (UNIT=fid, FMT=*, IOSTAT=fstat)                   &
              extendTable(iElem),atail(1,iElem),atail(2,iElem), &
              anose(1,iElem),anose(2,iElem),rnose_c(iElem)
         if (fstat .ne. 0) then
            call logError(11)
            ierr = -1
            return
         end if

         ! The next line should be the number of Reynolds numbers.
         read (UNIT=fid, FMT=*, IOSTAT=fstat) NRe(iElem)
         if (fstat .ne. 0) then
            call logError(11)
            ierr = -1
            return
         else if (NRe(iElem) .gt. MRe) then
            call logError(13)
            ierr = -1
            return
         end if

         ! Read data.  For each Reynolds number:
         do iRe = 1,NRe(iElem)
            ! Read the value of Re and the number of aoa entries.
            read (UNIT=fid, FMT=*, IOSTAT=fstat) Re(iRe,iElem), Naoa(iRe,iElem)
            if (fstat .ne. 0) then
               call logError(11)
               ierr = -1
               return
            else if (Re(iRe,iElem) .lt. 0.d0) then
               call logError(14)
               ierr = -1
               return
            else if (Naoa(iRe,iElem) .gt. Maoa) then
               call logError(15)
               ierr = -1
               return
            else if ((iRe .gt. 1) .and. (Re((iRe-1),iElem) .ge. Re(iRe,iElem))) then
               ! Re must be in ascending order.
               call logError(16)
               ierr = -1
               return
            end if
          
            do iaoa = 1,Naoa(iRe,iElem)
               ! Read data.
               read (UNIT=fid, FMT=*, IOSTAT=fstat) aoa(iaoa,iRe,iElem), &
                                                     Cl(iaoa,iRe,iElem), &
                                                     Cd(iaoa,iRe,iElem), &
                                                     Cm(iaoa,iRe,iElem)

               if (fstat .ne. 0) then
                  call logError(11)
                  ierr = -1
                  return
               else if ((iaoa .gt. 1) .and. &
                        (aoa((iaoa-1),iRe,iElem) .ge. aoa(iaoa,iRe,iElem))) then
                  ! Aoa must be in ascending order.
                  call logError(17)
                  ierr = -1
                  return
               end if
            end do

            !  If requested, fill out the airfoil coefficient table through high angles of attack.
            if (extendTable(iElem) .eq. 1) then
                  
               !  Store the old coefficients in temporary arrays.
               do iaoa = 1,Naoa(iRe,iElem)

                  tempaoa(iaoa) = aoa(iaoa,iRe,iElem)
                  tempCl(iaoa) = Cl(iaoa,iRe,iElem)
                  tempCd(iaoa) = Cd(iaoa,iRe,iElem)
                  tempCm(iaoa) = Cm(iaoa,iRe,iElem)

               end do

               !  Compute the number of data points between -90 degrees and the first
               !  existing value.  Create new values every five degrees.
               Nlow = int((180.d0 + aoa(1,iRe,iElem))/5.d0) - 1

               !  Compute the number of data points between the last existing value 
               !  and 90 degrees.
               Nhigh = int((180.d0 - aoa(Naoa(iRe,iElem),iRe,iElem))/5.d0) 

               if (Nlow+Naoa(iRe,iElem)+Nhigh .gt. Maoa) then
                  call logError(15)
                  ierr = -1
                  return
               end if

               !  Fill the lower part of the table.
               do iaoa = 1,Nlow
                  aoa(iaoa,iRe,iElem) = -180.d0 + 5.d0*dble(iaoa - 1)

                  call deepStall(aoa(iaoa,iRe,iElem),tempaoa(1),tempCl(1),tempCd(1),     & 
                                 ar,atail(2,iElem),anose(2,iElem),rnose_c(iElem),        &
                                 Cl(iaoa,iRe,iElem),                                     &
                                 Cd(iaoa,iRe,iElem),                                     &
                                 Cm(iaoa,iRe,iElem) )
               end do

               !  Insert the existing data.
               do iaoa = (Nlow+1),(Nlow+Naoa(iRe,iElem))
                  aoa(iaoa,iRe,iElem) = tempaoa(iaoa-Nlow)
                  Cl(iaoa,iRe,iElem) = tempCl(iaoa-Nlow)
                  Cd(iaoa,iRe,iElem) = tempCd(iaoa-Nlow)
                  Cm(iaoa,iRe,iElem) = tempCm(iaoa-Nlow)
               end do

               !  Fill the upper part of the table.
               do iaoa = (Nlow+Naoa(iRe,iElem)+1),(Nlow+Naoa(iRe,iElem)+Nhigh)
                  aoa(iaoa,iRe,iElem) = dble(180 + (iaoa - (Nlow+Naoa(iRe,iElem)+Nhigh))*5)
                  call deepStall(aoa(iaoa,iRe,iElem),tempaoa(Naoa(iRe,iElem)),                         &
                                 tempCl(Naoa(iRe,iElem)),tempCd(Naoa(iRe,iElem)),                      & 
                                 ar,atail(1,iElem),anose(1,iElem),rnose_c(iElem),                      &
                                 Cl(iaoa,iRe,iElem),                                                   &
                                 Cd(iaoa,iRe,iElem),                                                   &
                                 Cm(iaoa,iRe,iElem))
               end do

               !  Set new values for the number of entries in the table 
               !  (that is, the active dimensions of the arrays).
               Naoa(iRe,iElem) = (Nlow+Naoa(iRe,iElem)+Nhigh)

            end if
         end do
      end do

      close(fid)

   end subroutine readCoefficientLibrary

   subroutine deepStall (alphadeg,alpha_sdeg,Cl_s,Cd_s,ar,ataildeg,anosedeg,rnose_c,Cl,Cd,Cm)
   !
   !  This subroutine calculates airfoil coefficients during a condition
   !  of deep stall.  Of course, these can only be considered approximate.
   !
   !  Coefficient calculation comes from Lindenburg.
   !
   !  Inputs:
   !  -------
   !  alphadeg       : current angle-of-attack (degrees)
   !  alpha_sdeg,    : angle-of-attack (degrees) and lift, drag coefficients at which
   !  Cl_s,Cd_s        the data ends.
   !  ar             : aspect ratio of the airfoil section.
   !  ataildeg,      : if the side of the airfoil facing the wind were
   !  anosedeg         represented by a wedge, this is the angle between
   !                   a line perpendicular to the flow and the line
   !                   from the tip of the wedge. Units of degrees. 
   !  rnose_c        : the ratio of the nose radius to the chord of the airfoil.
   !  
   !
   !  Outputs:
   !  --------
   !  Cl,Cd,Cm       : lift, drag, and moment coefficients.
   !

      implicit none

      !  Inputs:
      double precision, intent(in) :: alphadeg,alpha_sdeg,Cl_s,Cd_s,ar,ataildeg,anosedeg,rnose_c

      !  Outputs:
      double precision, intent(out) :: Cl,Cd,Cm

      !  Local variables:
      double precision :: alpha,alpha_s,ca,sa,cas,sas,xcp,Cn,Cdmax,atail,anose

      alpha = alphadeg*PI/180.d0
      alpha_s = alpha_sdeg*PI/180.d0
      atail = ataildeg*PI/180.d0
      anose = anosedeg*PI/180.d0

      !  Compute cosines and sines for brevity.
      ca = cos(alpha)
      sa = sin(alpha)
      cas = cos(alpha_s)
      sas = sin(alpha_s)

      !  Compute maximum drag coefficient.
      Cdmax = 2.d0 - 0.82d0*(1.d0 - exp(-17.d0/ar))

      !  Lift coefficient:
      if (abs(alpha) .le. PI/2.d0) then
         Cl = (Cdmax/2.d0)*sin(2.d0*alpha) + ((Cl_s - Cdmax*sas*cas)*sas/(cas*cas))*ca*ca/sa
      else
         Cl = Cdmax*sa*ca
      end if

      !  Drag coefficient:
      if (abs(alpha) .le. PI/2.d0) then
         Cd = Cdmax*sa*sa + ((Cd_s - Cdmax*sas*sas)/cas)*ca
      else
         Cd = Cdmax*sa*sa
      end if

      !  Transform to find normal coefficient.
      Cn = Cl*ca + Cd*sa

      !  Compute the aerodynamic center at 90 degrees flow.
      xcp = 0.5d0 - 0.35d0*(atail*(0.2d0+0.08d0*atail)   &
          + (0.3d0 - anose*(0.2d0+0.08d0*anose))         &
          * (1.d0 - 1.8d0*sqrt(rnose_c))                 &
          - 0.3d0)

      !  Moment coefficient:
      Cm = -Cn*(xcp - 0.16d0*(1.d0 - 2.d0*abs(alpha)/PI) - 0.25d0)

   end subroutine deepStall










   subroutine steadyAeroSectionLoads(NRe,Naoa,density,viscosity,chord,    &
                                     Vrel,Rea,aoaa,Cla,Cda,Cma,aoafs,     &
                                     Fl,Fd,M,aoa,ierr)
   ! Computes aerodynamic loads on an airfoil section in flow that has been
   ! constant for a long time relative to c/V.
   !
   ! Inputs:
   ! -------
   ! MRe,Maoa       : max number for array dimensioning.
   ! NRe,Naoa       : number of entries in coefficient table arrays.
   ! Density, viscosity, chord
   ! Vrel           : incoming velocity in airfoil coordinates.
   ! Rea,aoaa,Cla,Cda,Cma: airfoil coefficient table.
   ! aoafs          : angle-of-attack at full-stall; used only to set
   !                  aoa in the case of near-zero flow velocity.
   !
   ! Outputs:
   ! --------
   ! Fl,Fd,M        : lift, drag, and moment forces per unit spanwise length.
   ! aoa            : angle of attack.
   ! ierr           : negative if an error has occurred.
   !
   ! Local variables:
   ! ----------------
   !  

      implicit none

      ! Inputs:
      integer,                             intent(in)  :: NRe
      integer,          dimension(:),      intent(in)  :: Naoa
      double precision,                    intent(in)  :: density,viscosity
      double precision,                    intent(in)  :: chord,aoafs
      double precision, dimension(:),      intent(in)  :: Vrel
      double precision, dimension(:),      intent(in)  :: Rea      
      double precision, dimension(:,:),    intent(in)  :: aoaa,Cla,Cda,Cma 

      ! Outputs:
      double precision,                    intent(out) :: Fl,Fd,M,aoa

      integer,                           intent(inout) :: ierr

      ! Local Variables:
      double precision                                 :: Re,Cl,Cd,Cm,magVrel

      !  Compute the magnitude of the incoming velocity, ignoring the
      !  spanwise component.
      magVrel = sqrt(Vrel(1)**2 + Vrel(2)**2)

      ! Compute angle-of-attack based upon input velocities.
      if (magVrel .lt. SMALLV) then
         ! Prevent divide-by-zero in the trivial case of no velocity.
         aoa = aoafs + 1.d0
      else
         ! The chord direction, from which aoa is measured, is (1,0,0)
         ! in airfoil coordinates.
         aoa = (180.d0/PI)*atan2(Vrel(2),Vrel(1))
      end if

      !  Compute the current Reynolds number.
      Re = density*magVrel*chord/viscosity

      call interpolateCoefficients (coeffFlag,Re,aoa,NRe,Naoa,Rea, &
                                    aoaa,Cla,Cda,Cma,              &
                                    Cl,Cd,Cm,ierr)
      if (ierr .lt. 0) then
         call logError(20)
         return
      end if

      ! Compute forces.
      Fl = 0.5d0*Cl*density*(magVrel**2)*chord
      Fd = 0.5d0*Cd*density*(magVrel**2)*chord
      M = 0.5d0*Cm*density*(magVrel**2)*(chord**2)

   end subroutine steadyAeroSectionLoads









   subroutine interpolateCoefficients (flag,Re,aoa,NRe,Naoa,     &
                                       Rea,aoaa,Cla,Cda,Cma,     &
                                       Cl,Cd,Cm,ierr)
   ! Performs a two-dimensonal linear interpolation in order to obtain
   ! coefficients at desired Re, aoa.
   ! 
   ! flag: Set to 1 for normal interpolation.
   !       Set to 2 if coefficient file is fixed for the entire
   !       analysis, and is defined each degree -180 < alpha < 180.
   !       

      implicit none

      ! Inputs:
      integer,                          intent(in)  :: NRe,flag
      integer,          dimension(:),   intent(in)  :: Naoa
      double precision,                 intent(in)  :: Re,aoa
      double precision, dimension(:),   intent(in)  :: Rea
      double precision, dimension(:,:), intent(in)  :: aoaa,Cla,Cda,Cma

      ! Outputs:
      double precision,                 intent(out) :: Cl,Cd,Cm

      integer,                        intent(inout) :: ierr

      ! Local Variables:
      integer :: irel,ireu,ial_rel,iau_rel,ial_reu,iau_reu
      double precision :: aoal,aoau,Cll,Clu,Cdl,Cdu,Cml,Cmu

      if (flag .eq. 2) then

         !  (Use variable names already declared for procedure below.)
         aoal = aoa + 181.d0
         irel = floor(aoal)
         ireu = ceiling(aoal)

         aoau = aoa - floor(aoa)
         Cl = aoau*(Cla(ireu,1) - Cla(irel,1)) + Cla(irel,1)
         Cd = aoau*(Cda(ireu,1) - Cda(irel,1)) + Cda(irel,1)
         Cm = aoau*(Cma(ireu,1) - Cma(irel,1)) + Cma(irel,1)

      else

         ! Scan the Reynolds array to bracket the desired Re value.
         call bisectSearch(Rea(1:NRe),NRe,Re,irel,ireu,ierr)
         if (ierr .eq. -1) then
            call logError(18)
            return
         else if (ierr .eq. -2) then
            call logError(19)
            return            
         end if

         ! Scan the aoa arrays corresponding to each of the bracketing
         ! Re numbers to obtain bracketing aoa values.
         call bisectSearch(aoaa(1:Naoa(irel),irel),Naoa(irel),aoa, &
                           ial_rel,iau_rel,ierr) 
         if (ierr .eq. -1) then
            call logError(18)
            return
         else if (ierr .eq. -2) then
            call logError(19)
            return            
         end if
         call bisectSearch(aoaa(1:Naoa(ireu),ireu),Naoa(ireu),aoa, &
                           ial_reu,iau_reu,ierr) 
         if (ierr .eq. -1) then
            call logError(18)
            return
         else if (ierr .eq. -2) then
            call logError(19)
            return            
         end if

         ! Interpolate aoa and coefficient values to obtain values at the
         ! desired Re.
         call interpolate (Rea(irel),aoaa(ial_rel,irel), &
                           Rea(ireu),aoaa(ial_reu,ireu), &
                           Re,aoal)
         call interpolate (Rea(irel),aoaa(iau_rel,irel), &
                           Rea(ireu),aoaa(iau_reu,ireu), &
                           Re,aoau)
         call interpolate (Rea(irel),Cla(ial_rel,irel), &
                           Rea(ireu),Cla(ial_reu,ireu), &
                           Re,Cll)
         call interpolate (Rea(irel),Cla(iau_rel,irel), &
                           Rea(ireu),Cla(iau_reu,ireu), &
                           Re,Clu)
         call interpolate (Rea(irel),Cda(ial_rel,irel), &
                           Rea(ireu),Cda(ial_reu,ireu), &
                           Re,Cdl)
         call interpolate (Rea(irel),Cda(iau_rel,irel), &
                           Rea(ireu),Cda(iau_reu,ireu), &
                           Re,Cdu)
         call interpolate (Rea(irel),Cma(ial_rel,irel), &
                           Rea(ireu),Cma(ial_reu,ireu), &
                           Re,Cml)
         call interpolate (Rea(irel),Cma(iau_rel,irel), &
                           Rea(ireu),Cma(iau_reu,ireu), &
                           Re,Cmu)

         ! Interpolate again, on the basis of aoa, to obtain coefficients
         ! at the desired Re and aoa.
         call interpolate (aoal,Cll,aoau,Clu,aoa,Cl)
         call interpolate (aoal,Cdl,aoau,Cdu,aoa,Cd)
         call interpolate (aoal,Cml,aoau,Cmu,aoa,Cm)

      end if  !  Short or long version of interpolation.

   end subroutine interpolateCoefficients

   subroutine bisectSearch (x,n,xreq, &
                            ilb,iub,ierr)
   ! Search an array for bracketing indices.
         
      implicit none
         
      ! Inputs:
      integer, intent(in) :: n
      double precision, dimension(n), intent(in) :: x
      double precision, intent(in) :: xreq
      ! Outputs:
      integer, intent(out) :: ilb, iub
      integer :: ierr

      ! Local Variables:
      integer :: i,is
      logical :: flag
      double precision :: s

      if (n .eq. 1) then
         ilb = 1
         iub = 1
         ierr = -2
         return 
      else if ((xreq .lt. x(1)) .or. (xreq .gt. x(n))) then
         ierr = -2
         return
      end if

      !  Initialize bounds and search point.
      ilb = 1
      iub = n
      is = n/2
      s = x(is)

      !  Perform a bisecting search.
      flag = .false.
      i = 0
      do
         !  Check convergence.
         if ((iub-ilb) .eq. 1) then
            flag = .true.
         else if (i .ge. 100) then
            flag = .true.
            ierr = -1
         end if
         if (flag) exit

         !  Move bounds based upon search point.
         if (s .gt. xreq) then
            !  Narrow search to lower half.
            iub = is
         else
            !  Narrow search to upper half.
            ilb = is
         end if

         is = (iub+ilb)/2
         s = x(is) 
         i = i + 1
      end do

   end subroutine bisectSearch

   subroutine interpolate (x1,y1,x2,y2,x,y)
      implicit none
      double precision :: x1,y1,x2,y2,x,y
      y = ((x-x1)/(x2-x1))*(y2-y1)+y1
   end subroutine interpolate












   subroutine initializeDynamicStall (Re,NRe,Naoa,           &
                                      Rea,aoaa,Cla,Cda,Cma,  &
                                      aoazi,dClda_ai,aoafsi, &
                                      ierr)
   !  This subroutine sets parameters that are relatively
   !  constant during the dynamic stall calculation, and
   !  thus do not need to be updated every timestep.
   !
   !  Inputs: 
   !  -------
   !  Re              : Reynolds number
   !  NRea ... Cma    : airfoil data array information.
   !
   !  Outputs:
   !  --------
   !  aoazi           : angle-of-attack at zero lift.
   !  dClda_ai        : maximum slope of the Cl-a curve in
   !                    the linear region.
   !  aoafsi          : angle-of-attack at full separation.
   !
   !  Local variables:
   !  ----------------
   !

      !  Inputs:
      integer, intent(in)                               :: NRe
      integer, dimension(:), intent(in)                 :: Naoa 
      double precision, intent(in)                      :: Re
      double precision, dimension(:), intent(in)        :: Rea      
      double precision, dimension(:,:), intent(in)      :: aoaa,Cla,Cda,Cma

      !  Outputs:
      double precision                                  :: aoazi
      double precision, dimension(2)                    :: dClda_ai,aoafsi

      integer, intent(inout) :: ierr

      !  Local variables:
      double precision :: aoaL,aoaU,temp1,temp2, &
                          ClL,ClU,Cl,tempCd,tempCm

      !  Locate the zero-lift angle-of-attack.
      aoaL = -5.d0
      aoaU = 5.d0
      call interpolateCoefficients (coeffFlag,Re,aoaL,NRe,Naoa,    &
                                    Rea,aoaa,Cla,Cda,Cma,          &
                                    ClL,temp1,temp2,ierr)
      if (ierr .ne. 0) then
         call logError(49)
         return
      end if
      call interpolateCoefficients (coeffFlag,Re,aoaU,NRe,Naoa,    &
                                    Rea,aoaa,Cla,Cda,Cma,          &
                                    ClU,temp1,temp2,ierr)
      if (ierr .ne. 0) then
         call logError(49)
         return
      end if
      do
         aoazi = 0.5d0*(aoaU + aoaL)
         if (aoaU - aoaL .lt. 0.001d0) then
            exit
         else if ((abs(aoaU) .gt. 20.d0) .or. (abs(aoaL) .gt. 20.d0)) then
            call logError(50)
            ierr = -1
            return
         end if

         if (ClU .lt. 0.d0) then
            aoaU = aoaU + 4.d0
         else if (ClL .gt. 0.d0) then
            aoaL = aoaL - 4.d0
         else
            call interpolateCoefficients (coeffFlag,Re,aoazi,NRe,Naoa,   &
                                          Rea,aoaa,Cla,Cda,Cma,          &
                                          Cl,temp1,temp2,ierr)   
            if (ierr .ne. 0) then
               call logError(49)
               return
            end if      

            if (Cl .gt. 0.d0) then
               !  Zero is to the left.
               aoaU = aoazi
               ClU = Cl
            else
               aoaL = aoazi
               ClL = Cl
            end if
         end if
      end do

      !  Locate the maximum attached-flow slopes.
      aoaU = aoazi + 2.d0
      dClda_ai(1) = 1.d-6
      do
         call interpolateCoefficients (coeffFlag,Re,aoaU,NRe,Naoa,    &
                                       Rea,aoaa,Cla,Cda,Cma,          &
                                       ClU,temp1,temp2,ierr)     
         if (ierr .ne. 0) then
            call logError(49)
            return
         end if    
         temp1 = ClU/(aoaU - aoazi)
         if (temp1 .gt. dClda_ai(1)) then
            dClda_ai(1) = temp1
         end if
         if (aoaU .gt. 10.d0) then
            exit
         else
            aoaU = aoaU + 0.5d0
         end if
      end do
      aoaU = aoazi - 2.d0
      dClda_ai(2) = 1.d-6
      do
         call interpolateCoefficients (coeffFlag,Re,aoaU,NRe,Naoa,    &
                                       Rea,aoaa,Cla,Cda,Cma,          &
                                       ClU,temp1,temp2,ierr)     
         if (ierr .ne. 0) then
            call logError(49)
            return
         end if    
         temp1 = ClU/(aoaU - aoazi)
         if (temp1 .gt. dClda_ai(2)) then
            dClda_ai(2) = temp1
         end if
         if (aoaU .lt. -10.d0) then
            exit
         else
            aoaU = aoaU - 0.5d0
         end if
      end do

      !  Locate the maximum attached-flow slopes and upper
      !  bounds, for positive and negative AOA.
      aoaL = 10.d0
      aoaU = 40.d0
      call interpolateCoefficients (coeffFlag,Re,aoaL,NRe,Naoa,    &
                                    Rea,aoaa,Cla,Cda,Cma,          &
                                    ClL,tempCd,tempCm,ierr)
      if (ierr .ne. 0) then
         call logError(49)
         return
      end if
      call interpolateCoefficients (coeffFlag,Re,aoaU,NRe,Naoa,    &
                                    Rea,aoaa,Cla,Cda,Cma,          &
                                    ClU,tempCd,tempCm,ierr)
      if (ierr .ne. 0) then
         call logError(49)
         return
      end if

      do
         aoafsi(1) = 0.5d0*(aoaU + aoaL)

         if (aoaU - aoaL .lt. 0.001d0) then
            exit
         else if ((aoaU .gt. 180.d0) .or. (aoaL .lt. -180.d0)) then
            call logError(53)
            ierr = -1
            return
         end if

         if (dClda_ai(1) .gt. 1.d-3) then
            
            if (aoaU .lt. 4.d0*ClU/dClda_ai(1) + aoazi) then
               aoaU = aoaU + 4.d0
            else if (aoaL .gt. 4.d0*ClL/dClda_ai(1) + aoazi) then
               aoaL = aoaL - 4.d0
            else
               call interpolateCoefficients (coeffFlag,Re,aoafsi(1),NRe,Naoa, &
                                             Rea,aoaa,Cla,Cda,Cma,            &
                                             Cl,tempCd,tempCm,ierr)   
               if (ierr .ne. 0) then
                  call logError(49)
                  return
               end if   

               if (aoafsi(1) .gt. 4.d0*Cl/dClda_ai(1) + aoazi) then
                  !  Zero is to the left.
                  aoaU = aoafsi(1)
                  ClU = Cl
               else
                  aoaL = aoafsi(1)
                  ClL = Cl
               end if
            end if
         else
            !  Probably the root cylinder.
            aoafsi(1) = 0.d0
            exit
         end if
      end do

      aoaL = -10.d0
      aoaU = -40.d0
      call interpolateCoefficients (coeffFlag,Re,aoaL,NRe,Naoa,   &
                                    Rea,aoaa,Cla,Cda,Cma,         &
                                    ClL,tempCd,tempCm,ierr)
      if (ierr .ne. 0) then
         call logError(49)
         return
      end if

      call interpolateCoefficients (coeffFlag,Re,aoaU,NRe,Naoa,  &
                                    Rea,aoaa,Cla,Cda,Cma,        &
                                    Clu,tempCd,tempCm,ierr)
      if (ierr .ne. 0) then
         call logError(49)
         return
      end if

      do
         aoafsi(2) = 0.5d0*(aoaU + aoaL)

         if (aoaU - aoaL .gt. -0.001d0) then
            exit
         else if ((aoaU .lt. -180.d0) .or. (aoaL .gt. 180.d0)) then
            call logError(53)
            ierr = -1
            return
         end if

         if (dClda_ai(2) .gt. 1.d-3) then

            if (aoaU .gt. 4.d0*ClU/dClda_ai(2) + aoazi) then
               aoaU = aoaU - 4.d0
            else if (aoaL .lt. 4.d0*ClL/dClda_ai(2) + aoazi) then
               aoaL = aoaL + 4.d0
            else
               call interpolateCoefficients (coeffFlag,Re,aoafsi(2),NRe,Naoa, &
                                             Rea,aoaa,Cla,Cda,Cma,            &
                                             Cl,tempCd,tempCm,ierr)   
               if (ierr .ne. 0) then
                  call logError(49)
                  return
               end if      

               if (aoafsi(2) .lt. 4.d0*Cl/dClda_ai(2) + aoazi) then
                  !  Zero is to the right.
                  aoaU = aoafsi(2)
                  ClU = Cl
               else
                  aoaL = aoafsi(2)
                  ClL = Cl
               end if
            end if
         else
            aoafsi(2) = 0.d0
            exit
         end if
      end do
!print '(1X,5ES13.4)',aoazi,aoafsi(1),aoafsi(2),dClda_ai(1),dClda_ai(2)
   end subroutine initializeDynamicStall








   subroutine readTurbulence (fid,dL,Nbox,x0_w,x_r,Tr_w,u_r)
   !  This subroutine finds the indices in a turbulence grid
   !  surrounding a given point in space, and interpolates
   !  the turbulence velocity components at this point.
   !
   !  Inputs: 
   !  -------
   !  fid             : ID of the turbulence direct-access file.
   !  dL              : Spacing of turbulence grid in Xw,Yw,Zw.
   !  Nbox            : Size of turbulence grid (number of points in
   !                    Xw, Yw, and Zw directions).
   !  x0_w            : Position of the rotor coordinate origin
   !                    in wind coordinates.
   !  x_r             : Coordinates of the point at which turbulence
   !                    is desired, in rotor coordinates, measured
   !                    from the rotor coordinate origin.
   !  Tr_w            : Transform from rotor to wind coordinates.
   !
   !  Outputs:
   !  --------
   !  u_r             : Turbulent velocity vector at the given point
   !                    in space, in rotor coordinates.
   !  
   !  Local variables:
   !  ----------------
   !

      implicit none

      integer, parameter :: ikind = SELECTED_REAL_KIND(p=4,r=2)
      integer, parameter :: noff = 6

      integer,                          intent(in)  :: fid
      integer,          dimension(3),   intent(in)  :: Nbox
      double precision, dimension(3),   intent(in)  :: dL,x0_w,x_r
      double precision, dimension(3,3), intent(in)  :: Tr_w

      double precision, dimension(3),   intent(out) :: u_r

      integer                                       :: j,k,irec
      integer,          dimension(3)                :: ig
      integer,          dimension(2,3)              :: i
      integer,          dimension(3,8)              :: p
      double precision                              :: val,A,B,C,D, &
                                                       xlow,xhigh,  &
                                                       ylow,yhigh,  &
                                                       zlow ! ,zhigh
      double precision, dimension(3)                :: dx_w,x_w,xp_int, &
                                                       uzlow,uzhigh,    &
                                                       u_w
      double precision, dimension(3,8)              :: ug
      double precision, dimension(3,3)              :: Tw_r

      real(kind=ikind), dimension(3)                :: utemp

      u_r = 0.d0

      !  Find x in wind coordinates.
      dx_w = matvecmult(Tr_w,x_r,3)
      x_w = x0_w + dx_w

      !  In each coordinate direction, find the grid indices 
      !  bracketing x_w.
      do j = 1,3
         val = x_w(j)/dL(j)
         if ((val .gt. -1.d-6) .and. (val .lt. 1.d-6)) then
            !  Prevent exactly zero, which messes up some of
            !  the logic below.
            val = 1.d-6
         end if
         i(1,j) = floor(val) + 1
         i(2,j) = ceiling(val) + 1
         if (i(1,j) .eq. i(2,j)) then
            !  Correct for deficiencies in the floor and 
            !  ceiling functions when val is exactly an
            !  integer.
            i(2,j) = i(1,j) + 1
         end if
         !  Take advantage of (assumed) periodicity in the
         !  turbulence grid and normalize the indices to the
         !  range 1 <= i <= Nbox.
         do k = 1,2
            do
               if (i(k,j) .lt. 1) then
                  i(k,j) = i(k,j) + Nbox(j)
               else if (i(k,j) .gt. Nbox(j)) then
                  i(k,j) = i(k,j) - Nbox(j)
               else
                  exit
               end if
            end do
         end do
      end do

!write (41,'(1X,6F10.2,6I6,A)')                    &
!x_w,dL,i(1,1),i(2,1),i(1,2),i(2,2),i(1,3),i(2,3), &
!'   xw_x,y,z dL_x,y,z ixlow ixhigh iylow iyhigh izlow izhigh'

      !  We now have eight points in the turbulence file
      !  to find:
      !       ix  iy  iz
      !  ----------------
      !   1    L   L   L
      !   2    H   L   L
      !   3    L   H   L
      !   4    H   H   L
      !   5    L   L   H
      !   6    H   L   H
      !   7    L   H   H
      !   8    H   H   H
      p(1,1) = 1
      p(1,2) = 2
      p(1,3) = 1
      p(1,4) = 2
      p(1,5) = 1
      p(1,6) = 2
      p(1,7) = 1
      p(1,8) = 2
      p(2,1) = 1
      p(2,2) = 1
      p(2,3) = 2
      p(2,4) = 2
      p(2,5) = 1
      p(2,6) = 1
      p(2,7) = 2
      p(2,8) = 2
      p(3,1) = 1
      p(3,2) = 1
      p(3,3) = 1
      p(3,4) = 1
      p(3,5) = 2
      p(3,6) = 2
      p(3,7) = 2
      p(3,8) = 2

      do k = 1,8

         do j = 1,3
            ig(j) = i(p(j,k),j)
         end do

         irec = noff + ig(1) + Nbox(1)*(ig(2)-1) + Nbox(1)*Nbox(2)*(ig(3)-1)
         read(fid,REC=irec) utemp
         do j = 1,3
            ug(j,k) = dble(utemp(j))
         end do   

!write (41,'(1X,4I3,4I6,4I10,3F9.4,A)')                                           &
!k,p(:,k),noff,ig,ig(1),Nbox(1)*(ig(2)-1),Nbox(1)*Nbox(2)*(ig(3)-1),irec,ug(:,k), &
!'   k px,y,z Noffset igx,y,z irec_x,y,z,tot ugx,y,z'

         !  Coordinates (for interpolation only, take account of
         !  periodicity) of the present grid point.
         if (k .eq. 1) then
            xlow = dble(ig(1)-1)*dL(1)
            ylow = dble(ig(2)-1)*dL(2)
            zlow = dble(ig(3)-1)*dL(3)
         else if (k .eq. 8) then
            !  If ig = 1 for xhigh, then it will be compared against
            !  xlow who has a coordinate of (Nbox-1)*dL.  In this case,
            !  it is desireable to let xhigh = Nbox*dL, rather than 0.
            if (ig(1) .eq. 1) then
               xhigh = Nbox(1)*dL(1)
            else  !  Valid in all other cases.
               xhigh = dble(ig(1)-1)*dL(1)
            end if

            !  Same goes for y.
            if (ig(2) .eq. 1) then
               yhigh = Nbox(2)*dL(2)
            else 
               yhigh = dble(ig(2)-1)*dL(2)
            end if

            ! zhigh is not needed for the formulas as implemented below.
         end if

      end do

!write (41,'(1X,5F10.2,A)')   &
!xlow,ylow,zlow,xhigh,yhigh,  &
!'   xlow ylow zlow xhigh yhigh'

      do j = 1,3
         !  Normalize x_w to positive to compare with grid point
         !  coordinates.
         xp_int(j) = x_w(j)
         do                                 
            val = dble(Nbox(j))*dL(j)
            if (xp_int(j) .lt. 0.d0) then
               xp_int(j) = xp_int(j) + val
            else if (xp_int(j) .gt. val) then
               xp_int(j) = xp_int(j) - val
            else
               exit
            end if
         end do
      end do

!write (41,'(1X,6F9.2,A)')                           &
!x_w(1),x_w(2),x_w(3),xp_int(1),xp_int(2),xp_int(3), &
!'   x_w xp_int'

      !  Spatial interpolation.
      val = dL(1)*dL(2)
      A = (xhigh - xp_int(1))*(yhigh - xp_int(2))/val
      B = (xp_int(1) - xlow)*(yhigh - xp_int(2))/val
      C = (xhigh - xp_int(1))*(xp_int(2) - ylow)/val
      D = (xp_int(1) - xlow)*(xp_int(2) - ylow)/val
      do j = 1,3
         uzlow(j) = A*ug(j,1) + B*ug(j,2) + C*ug(j,3) + D*ug(j,4)
         !  Take advantage of a regular grid, don't need to recompute
         !  A, B, C, D for high-Z case.
         uzhigh(j) = A*ug(j,5) + B*ug(j,6) + C*ug(j,7) + D*ug(j,8)

         u_w(j) = ((xp_int(3) - zlow)/dL(3))*(uzhigh(j) - uzlow(j)) + uzlow(j)
      end do

!write (41,'(1X,13F9.4,A)')  &
!A,B,C,D,uzlow,uzhigh,u_w,   &
!'   A B C D uzlow uzhigh uw'

      Tw_r = transpose(Tr_w)
      u_r = matvecmult(Tw_r,u_w,3)

   end subroutine readTurbulence











   function matmatmult(a, b, rowsOut, colsOut)
      implicit none
      integer rowsOut, colsOut  ! Necessary to dimension the function result.
      double precision, dimension(rowsOut,colsOut) :: matmatmult
      double precision, dimension(:,:), intent(in) :: a
      double precision, dimension(:,:), intent(in) :: b
      integer row, col, i

      if ((size(a,2) .eq. size(b,1)) .and. (size(a,1) .eq. rowsOut) .and. (size(b,2) .eq. colsOut)) then
         do row = 1, rowsOut
            do col = 1, colsOut
               matmatmult(row,col) = 0.d0
               do i = 1, size(a,2)
                  matmatmult(row,col) = matmatmult(row,col) + a(row,i)*b(i,col)
               end do
            end do
         end do
      else  ! Invalid matrix dimensions
         matmatmult = 0.d0
      end if
      
   end function matmatmult

   function matvecmult(a, b, rowsOut)
      implicit none
      integer rowsOut  ! Necessary to dimension the function result.
      double precision, dimension(rowsOut) :: matvecmult
      double precision, dimension(:,:), intent(in) :: a
      double precision, dimension(:), intent(in) :: b
      integer row, i

      if ((size(a,2) .eq. size(b)) .and. (size(a,1) .eq. rowsOut)) then
         do row = 1, rowsOut
            matvecmult(row) = 0.d0
            do i = 1, size(b)
               matvecmult(row) = matvecmult(row) + a(row,i)*b(i)
            end do
         end do
      else  ! Invalid matrix dimensions
         matvecmult = 0.d0
      end if
      
   end function matvecmult

   function dot(vec1, vec2)
      implicit none
      double precision :: dot
      double precision, dimension(:), intent(in) :: vec1
      double precision, dimension(:), intent(in) :: vec2
      integer i

      dot = 0.d0
      
      ! Verify that the result is defined.  For simplicity's sake, if the product
      ! is undefined, zero is returned.
      if ((size(vec1) .eq. size(vec2)) .and. (size(vec1) .gt. 0)) then
         do i = 1, size(vec1)
            dot = dot + vec1(i)*vec2(i)
         end do
      else
         call logError(2)
      end if
         
   end function dot

   function cross(vec1, vec2)
      implicit none
      ! Note that this is defined only for three dimensions.
      double precision, dimension(3) :: cross
      double precision, dimension(:), intent(in) :: vec1, vec2
      if ((size(vec1) .eq. 3) .and. (size(vec2) .eq. 3)) THEN
         cross(1) = vec1(2)*vec2(3) - vec1(3)*vec2(2)
         cross(2) = vec1(3)*vec2(1) - vec1(1)*vec2(3)
         cross(3) = vec1(1)*vec2(2) - vec1(2)*vec2(1)
      else
         cross = 0.d0
         call logError(2)
      end if
   end function cross

   function magnitude(vec)
      implicit none
      double precision :: magnitude
      double precision, dimension(:), intent(in) :: vec
      double precision :: sum
      integer i
      sum = 0.d0
      do i = 1, size(vec)
         sum = sum + vec(i)**2
      end do
      magnitude = sqrt(sum)
   end function magnitude











   subroutine initErrorLog()
      implicit none
      integer i
      do i = 1,nErrorTypes
         errorLog(i) = 0
      end do
   end subroutine

   subroutine logError(errorNum)
      implicit none
      integer errorNum
      errorlog(errorNum) = errorLog(errorNum) + 1
   end subroutine

   subroutine writeErrors(errorFile)
      implicit none
      integer errorFile
         
      character(len=70), dimension(nErrorTypes) :: errorText
      integer i
         
      errorText(1) = &
         'A matrix operation was performed with mismatched matrix sizes.'
      errorText(2) = &
         'A vector operation was performed with mismatched vector lengths.'
      errorText(3) = &
         'A numerical problem was encountered during matrix inversion.'
      errorText(4) = &
         'The magnitude of a coordinate system axis was zero.'
      errorText(5) = &
         'Bad input during blade element definition.'
      errorText(6) = &
         'File error during blade element definition.'
      errorText(7) = &
         'A negative or out-of-order Re number was in the aero coefficient file.'
      errorText(8) = &
         'An out-of-range or out-of-order AOA was in the aero coefficient file.'
      errorText(9) = &
         'The static BEM analysis did not converge.'
      errorText(10) = &
         'Error opening airfoil coefficient file.'
      errorText(11) = &
         'Error reading airfoil coefficient file.'
      errorText(12) = &
         'The airfoil name was not found.'
      errorText(13) = &
         'Too many Re numbers in airfoil coefficient library.'
      errorText(14) = &
         'Re number cannot be negative.'
      errorText(15) = &
         'Too many AOA points in airfoil coefficient library.'
      errorText(16) = &
         'Reynolds numbers must be listed in increasing order.'
      errorText(17) = &
         'Angles-of-Attack must be listed in increasing order.'
      errorText(18) = &
         'Library search did not converge.'
      errorText(19) = &
         'Requested value outside of library bounds.'
      errorText(20) = &
         'Error during coefficient interpolation.'
      errorText(21) = &
         'Warning, spanwise velocity is not negligible.'
      errorText(22) = &
         'Prandtl factor is unreasonably small.'
      errorText(23) = &
         'Thrust coefficient is negative.'
      errorText(24) = &
         'Axial induced velocity calculation did not converge.'
      errorText(25) = &
         'The yaw angle must be less than 90 degrees.'
      errorText(26) = &
         'Incorrect tangential induction option.'
      errorText(27) = &
         'Error encountered during axial induction calculation.'
      errorText(28) = &
         'Input error during induced velocity calculation.'
      errorText(29) = &
         'Error encountered during momentum balance calculation.'
      errorText(30) = &
         'Induced velocity cannot be calculated when incoming velocity is zero.'
      errorText(31) = &
         'Invalid option selected for induced velocity calculation.'
      errorText(32) = &
         'Unrealistic wake angle calculated.'
      errorText(33) = &
         'Invalid wake angle calculation option.'
      errorText(34) = &
         'Invalid input to steady BEM calculation.'
      errorText(35) = & 
         'Error encountered while computing steady aerodynamic loads.'
      errorText(36) = & 
         'Error encountered during induced velocity calculation.'
      errorText(37) = & 
         'Error encountered while calculating wake direction vector.'
      errorText(38) = &
         'Error in the input to the element length calculation.'
      errorText(39) = &
         'Too many data points: coefficient table cannot be extended.'
      errorText(40) = &
         'Error encountered during quasi-static BEM calculation.'
      errorText(41) = &
         'Invalid element length input.'
      errorText(42) = &
         'Too many airfoil elements.'
      errorText(43) = &
         'Error reading option flags.'
      errorText(44) = &
         'Error reading density and viscosity.'
      errorText(45) = &
         'Error reading radius and cone angle.'
      errorText(46) = &
         'Error reading number of blades and elements.'
      errorText(47) = &
         'Error reading airfoil geometry.'
      errorText(48) = &
         'Error opening input file.'
      errorText(49) = & 
         'Error encountered during dynamic stall analysis.'
      errorText(50) = &
         'Zero-lift AOA is outside reasonable bounds.'
      errorText(51) = &
         'Zero relative velocity at an element.'
      errorText(52) = &
         'Divided-by-zero when calculating thrust coefficient.'
      errorText(53) = &
         'Deep-stall AOA is outside reasonable bounds.'

      write (errorfile,'(/1x,a,a,/,/)')                     &
            'The aerodynamics module has encountered ',   &
            'the following errors:'
      do i = 1,nerrortypes
         if (errorlog(i) .ne. 0) then
            write (errorfile,'(/1x,a70,a3,i3,a12/)') &
                  errortext(i),' : ',errorlog(i),' occurrences'
         end if
      end do
   end subroutine writeErrors



end module VAWTBEM


























!-------------------------------------------------------------------------------



module setupVAWTBEM
! Contains subroutines to perform a steady blade element-momentum calculation.
!
! Note: in a real generalized dynamic analysis, surface elements must be
! redefined after each timestep, since a change in wind direction will
! change which elements are on the upwind side and which are on the 
! downwind side.
!
! Revisions:
!     Date          Programmer          Notes
!  ----------      ------------         ----------------------------------------
!  20.05.2009        K. Merz            Original code.
!  30.06.2011        K. Merz            Rewritten for the new dynamic VAWTBEM code.
!                                       Eliminated stored variables.
!
!

   implicit none

   contains











   subroutine readSetupParameters (setupFid,setupFile,Nrows,Ncols,    &
                                   density,viscosity,                 &
                                   node,profile,Lelem,Lblade,Zelem,   &
                                   chord,twist0,xs_r,n_r,Tsr_s,       &
                                   airfoilName,ierr)
   !  Reads an input file containing setup parameters.
   !
   !  Inputs:
   !  -------
   !  setupFile,     : Name and file ID of geometry input file.
   !  setupFid
   !  Nrows, Ncols   : Number of elements per blade and number of
   !                   azimuth elements.
   !
   !  Outputs:
   !  --------
   !  density, viscosity
   !  node           : Nodes of the profile that is rotated to form the
   !                   swept surface.  Components are (r,z,offset).
   !  profile        : Element centroid locations along the profile.
   !  Lelem          : Spanwise length of each element.
   !  Lblade         : Total spanwise length along the profile.
   !  Zelem          : Spanwise coordinate along the length, measured
   !                   from the center span.  Used for Prandtl factor
   !                   calculation.
   !  chord
   !  twist0         : Baseline blade twist (at zero pitch angle).
   !  xs_r           : Coordinates of each surface element in the rotor
   !                   coordinate system.
   !  n_r            : Normal to each surface element in the rotor
   !                   coordinate system.
   !  Tsr_s          : Transform from rotor to surface normal coordinates.
   !  airfoilName    : Airfoil ID for coefficient table look-up.
   !  
   !  Local variables:
   !  ----------------
   !
   !  

      use VAWTBEM

      implicit none

      !  Inputs:
      integer,                                    intent(in)   :: Nrows,Ncols
      character(len=*),                           intent(in)   :: setupFile
      integer,                                    intent(in)   :: setupFid

      !  Outputs:
      character(len=MAXSTRLEN),dimension(Nrows),  intent(out)  :: airfoilName
      double precision,                           intent(out)  :: density,viscosity,Lblade
      double precision, dimension(Nrows),         intent(out)  :: Lelem,Zelem,chord,twist0
      double precision, dimension(3,Nrows),       intent(out)  :: profile
      double precision, dimension(3,Nrows+1),     intent(out)  :: node
      double precision, dimension(3,Nrows,Ncols), intent(out)  :: xs_r,n_r
      double precision, dimension(3,3,Nrows,Ncols),intent(out) :: Tsr_s

      integer,                                   intent(inout) :: ierr

      !  Local variables:
      character(len=1) :: junk
      integer :: irow,icol,fstat
      double precision :: r1,z1,offset1,twistdeg,dp,magnr,dpsi, &
                          mag,magdir,z,cpsi,spsi,psi,tilt,      &
                          offsetang,ct,st,co,so
      double precision, dimension(3) :: dir,v1,v2,profrot
      double precision, dimension(3,2) :: nodrot
      double precision, dimension(3,3,4) :: T

      fstat = 0

      !  Open the file.
      open (UNIT=setupFid, FILE=setupFile, STATUS='OLD', IOSTAT=fstat, &
            ACCESS='SEQUENTIAL', FORM='FORMATTED',                     &
            ACTION='READ', POSITION='REWIND', DELIM='NONE')
      if (fstat .gt. 0) then
         ierr = -1
         call logError(48)
         return
      end if

      read (UNIT=setupFid, FMT=*, IOSTAT=fstat) &
      density,viscosity
      if (fstat .gt. 0) then
         ierr = -1
         call logError(44)
         return
      end if

      read (UNIT=setupFid, FMT=*, IOSTAT=fstat) junk  !  Already read array size.
      if (fstat .gt. 0) then
         ierr = -1
         call logError(46)
         return
      end if
      read (UNIT=setupFid, FMT=*, IOSTAT=fstat) &
      r1,z1,offset1
      if (fstat .gt. 0) then
         ierr = -1
         call logError(47)
         return
      end if

      node(1,1) = r1
      node(2,1) = z1
      node(3,1) = offset1

      do irow = 1,Nrows

         read (UNIT=setupFid, FMT=*, IOSTAT=fstat) &
         airfoilName(irow),r1,z1,offset1,chord(irow),twistdeg
         if (fstat .gt. 0) then
            ierr = -1
            call logError(47)
            return
         end if

         node(1,(irow+1)) = r1
         node(2,(irow+1)) = z1
         node(3,(irow+1)) = offset1

         !  Interpolate magnitude and direction to find
         !  the midpoint defining the blade profile.
         mag = 0.5d0*(magnitude(node(:,(irow+1))) &
                      + magnitude(node(:,irow)))
         dir = 0.5d0*(node(:,(irow+1)) + node(:,irow))
         magdir = magnitude(dir)
         Lelem(irow) = magnitude(node(:,(irow+1)) - node(:,irow))

         if (magdir .gt. 1.d-6) then
            profile(:,irow) = mag*(dir/magdir)     
         else
            !  Use straight linear interpolation.
            profile(:,irow) = dir
         end if

         twist0(irow) = twistdeg*PI/180.d0

      end do

      !  Compute element along-blade length position,
      !  and blade length.
      z = 0.d0
      do irow = 1,Nrows
         z = z + Lelem(irow)
      end do
      Lblade = z
      z = 0.d0
      do irow = 1,Nrows
         z = z + 0.5d0*Lelem(irow)
         if (z .le. Lblade/2.d0) then
            Zelem(irow) = Lblade/2.d0 - z
         else
            Zelem(irow) = z - Lblade/2.d0
         end if
         z = z + 0.5d0*Lelem(irow)
      end do

      !  Define the surface elements based upon the blade
      !  profile.  This involves taking the cross product
      !  between the line connecting adjacent nodes and
      !  the normal to the plane containing Z and one of
      !  the points in profile.
      !  First, though, rotate the profile curve about the
      !  Z axis such that it is aligned with the longitude
      !  of a column of surface elements.
      dpsi = 2.d0*PI/dble(Ncols)
      do icol = 1,Ncols
         psi = dpsi/2.d0 + dble(icol-1)*dpsi + PI/2.d0
         cpsi = cos(psi)
         spsi = sin(psi)
         do irow = 1,Nrows
            nodrot(1,1) = node(1,irow)*cpsi &
                        - node(3,irow)*spsi
            nodrot(2,1) = node(1,irow)*spsi &
                        + node(3,irow)*cpsi
            nodrot(3,1) = node(2,irow)
            nodrot(1,2) = node(1,(irow+1))*cpsi &
                        - node(3,(irow+1))*spsi
            nodrot(2,2) = node(1,(irow+1))*spsi &
                        + node(3,(irow+1))*cpsi
            nodrot(3,2) = node(2,(irow+1))
            profrot(1) = profile(1,irow)*cpsi &
                       - profile(3,irow)*spsi
            profrot(2) = profile(1,irow)*spsi &
                       + profile(3,irow)*cpsi
            profrot(3) = profile(2,irow)
            v1 = nodrot(:,2) - nodrot(:,1)
            v2(1) = profrot(2)
            v2(2) = -profrot(1)
            v2(3) = 0.d0

            n_r(:,irow,icol) = cross(v1,v2)

            !  Verify the sign.
            dp = dot(n_r(:,irow,icol),profrot)
            if (dp .lt. 0.d0) then
               n_r(:,irow,icol) = -n_r(:,irow,icol)
            end if
            magnr = magnitude(n_r(:,irow,icol))
            if (magnr .lt. 1.d-6) then
               call logError(53)
               ierr = -1
               return
            else
               n_r(:,irow,icol) = n_r(:,irow,icol)/magnr
            end if

            xs_r(:,irow,icol) = profrot

            !  Compute transformation matrices. 
            !
            !  The surface element coordinate system will be defined
            !  with the Zs axis normal to the surface, and the XsYs
            !  plane aligned with the Zr axis.
            !
            !  To get from the surface to the rotor coordinate system:
            !  (1)  Define Ys parallel to Zr and Zs parallel to Xr.
            !  (2)  Rotate about the Zr axis by the offset angle.
            !  (3)  Rotate about the Yr axis by the tilt of the blade
            !       element.
            !  (4)  Rotate about the Zr axis by the azimuth angle.
            tilt = atan2(-n_r(3,irow,icol), &
                         sqrt(n_r(1,irow,icol)**2 + n_r(2,irow,icol)**2))
            offsetang = atan2(profile(3,irow),profile(1,irow))

            ct = cos(tilt)
            st = sin(tilt)
            co = cos(offsetang)
            so = sin(offsetang)

            T(1,1,1) = 0.d0
            T(1,2,1) = 1.d0
            T(1,3,1) = 0.d0
            T(2,1,1) = 0.d0
            T(2,2,1) = 0.d0
            T(2,3,1) = 1.d0
            T(3,1,1) = 1.d0
            T(3,2,1) = 0.d0
            T(3,3,1) = 0.d0

            T(1,1,2) = co
            T(1,2,2) = so
            T(1,3,2) = 0.d0
            T(2,1,2) = -so
            T(2,2,2) = co
            T(2,3,2) = 0.d0
            T(3,1,2) = 0.d0
            T(3,2,2) = 0.d0
            T(3,3,2) = 1.d0

            T(1,1,3) = ct
            T(1,2,3) = 0.d0
            T(1,3,3) = -st
            T(2,1,3) = 0.d0
            T(2,2,3) = 1.d0
            T(2,3,3) = 0.d0
            T(3,1,3) = st
            T(3,2,3) = 0.d0
            T(3,3,3) = ct

            T(1,1,4) = cpsi
            T(1,2,4) = spsi
            T(1,3,4) = 0.d0
            T(2,1,4) = -spsi
            T(2,2,4) = cpsi
            T(2,3,4) = 0.d0
            T(3,1,4) = 0.d0
            T(3,2,4) = 0.d0
            T(3,3,4) = 1.d0

            Tsr_s(:,:,irow,icol) = matmatmult(T(:,:,1),T(:,:,2),3,3)
            Tsr_s(:,:,irow,icol) = matmatmult(Tsr_s(:,:,irow,icol),T(:,:,3),3,3)
            Tsr_s(:,:,irow,icol) = matmatmult(Tsr_s(:,:,irow,icol),T(:,:,4),3,3)
 
         end do
      end do


   end subroutine readSetupParameters


end module setupVAWTBEM






















!-------------------------------------------------------------------------------



module control
!  Contains subroutines related to the control systems.
!
! Revisions:
!     Date          Programmer          Notes
!  ----------      ------------         ----------------------------------------
!  07.10.2011        K. Merz            Original code.
!  02.11.2011        K. Merz            Implemented genT and estimateAvgSpeed
!                                       subroutines in order to attempt to 
!                                       reduce the fluctuating torque.
!

   implicit none

   contains





      subroutine genT (psi,w,T,ierr)

         implicit none

         double precision, intent(in) :: psi,w
         double precision, intent(out) :: T

         integer, intent(inout) :: ierr

         double precision :: wavg

         !  Accounting for the rotor azimuth relative to
         !  the wind, and the instantaneous rotational speed,
         !  estimate the revolution-average rotational speed.
!         call estimateAvgSpeed (psi,w,wavg,ierr)
!         if (ierr .ne. 0) then
!            return
!         end if

         !  Use the average speed to look up the (nominally
         !  uniform) torque.
         call Tlookup (w,T)

      end subroutine genT




      subroutine estimateAvgSpeed (psi,w,wavg,ierr)
      !  Uses false position to solve for average speed, given
      !  instantaneous speed.  (Press et al. p 249)
      !
      !  Values are particular to the turbine, in this case the
      !  Deepwind 5 MW 1st design.

         implicit none

         double precision, intent(in) :: psi,w
         double precision, intent(out) :: wavg

         integer, intent(inout) :: ierr

         integer :: iter
         double precision :: wL,wH,fL,fH,dw,del,f
       
         !  Bracket the root.  Amplitude larger than max of sine
         !  part of the function.
         dw = 0.002d0
         iter = 0
         do
            iter = iter + 1
            wL = w - dw
            call fspeed (wL,w,psi,fL)
!write (41,'(1X,I4,5ES12.3,A)')   &
!iter,psi,w,dw,wL,fL,             &
!'   iter psi w dw wL fL'
            if (fL .le. 0.d0) then
               exit
            else if (iter .ge. 50) then
print *,'Failed to bracket root.'
               ierr = -1
               return
            else
               dw = dw + 0.002d0
            end if
         end do
         dw = 0.002d0
         iter = 0
         do
            iter = iter + 1
            wH = w + dw
            call fspeed (wH,w,psi,fH)
!write (41,'(1X,I4,5ES12.3,A)')   &
!iter,psi,w,dw,wH,fH,             &
!'   iter psi w dw wH fH'
            if (fH .ge. 0.d0) then
               exit
            else if (iter .ge. 50) then
print *,'Failed to bracket root.'
               ierr = -1
               return
            else
               dw = dw + 0.002d0
            end if
         end do
        
         if (fL*fh .gt. 0.d0) then
print *,'Error, speed not bracketed.'
            ierr = -1
            return
         end if

         dw = wH - wL
         iter = 0
         do
            iter = iter + 1
            wavg = wL + dw*fL/(fL - fH)
            call fspeed (wavg,w,psi,f)
            if (f .lt. 0.d0) then
               del = wL - wavg
               wL = wavg
               fL = f
            else
               del = wH - wavg
               wH = wavg
               fH = f
            end if
            dw = wH - wL
!write (41,'(1X,I4,7ES12.3,A)') &
!iter,psi,w,wL,wH,fL,fH,wavg,   &
!'   iter psi w wL wH fL fH wavg'
            if (abs(del) .lt. 0.0001d0) then
               exit
            else if (iter .ge. 50) then
print *,'Calculation of average omega did not converge.'
               ierr = -1
               return
            end if
         end do

      end subroutine estimateAvgSpeed



      subroutine fspeed (wavg,w,psi,f)
         implicit none
         double precision, intent(in) :: psi,w,wavg
         double precision, intent(out) :: f
         double precision :: A,eps,ang
         if (wavg .le. 0.497d0) then
            A = ((wavg - 0.2145d0)/0.2825d0)*0.0119d0 + 0.006d0
            eps = ((wavg - 0.214d0)/0.283d0)*0.0124d0 - 0.0610d0
         else
            A = ((wavg - 0.497d0)/0.021d0)*0.0218d0 + 0.0179d0
            eps = ((wavg - 0.497d0)/0.027d0)*0.2017d0 - 0.0486d0
         end if
         ang = 2.d0*psi + eps
         f = wavg + A*sin(ang) - w
!write (41,'(1X,6ES12.3,A)')  &
!wavg,w,psi,A,eps,f,          &
!'   wavg w psi A eps f'
      end subroutine fspeed




      subroutine TLookup (wg,Tcontrol)
         implicit none
         double precision, intent(in) :: wg
         double precision, intent(out) :: Tcontrol
         Tcontrol = max((-0.01312d0 - 1.386d0*wg + 19.55d0*wg*wg), &
                        (-232.744d0 + 473.851d0*wg))*1.d6
!write (41,'(1X,2ES12.3,A)') wg,Tcontrol,'   wg Tcontrol'
      end subroutine TLookup


end module control




























!-------------------------------------------------------------------------------



module platformDynamics
!  Contains subroutines to analyze the Deepwind turbine dynamics.
!  Adaptive step-size Runge-Kutta comes from Press et al.
!
! Revisions:
!     Date          Programmer          Notes
!  ----------      ------------         ----------------------------------------
!  24.09.2011        K. Merz            Original code.
!
!

   use setupVAWTBEM
   use VAWTBEM
   use control

   implicit none

   contains







   subroutine simulateDeepwind (fid,setupFile,setupFid,          &
                                airfoilFile,airfoilFid,          &
                                turbFile,turbFid,                &
                                Nb,Nrows,Ncols,Ndof,Nadof,Nt,dt, &
                                wgset0,psisyn0,Tgset0,y0,ydot0,  &
                                psi0hin,omegah0,psi0gin,omegag0, &
                                xw0,Vref,z0,h0,xvane_g,          &
                                ierr)
   !  This is the top-level subroutine containing the time loop.
   !
   !  Inputs:
   !  -------
   !  fid             : Output file ID.
   !  setupFile,      : File name and ID of file describing rotor
   !  setupFid          geometry
   !  aeroFile,Fid    : Airfoil coefficient file.
   !  turbFile,Fid    : Turbulent windfield file.
   !  Nb,Nrows,Ncols  : Number blades, rows and columns of blade (or
   !                    surface) elements.
   !  Ndof,Nadof,     : Number of state variables, and aero state vars.
   !  Nt              : Number timesteps.
   !  dt              : timestep size.
   !  wgset0,psisyn0  : Initial generator synchronous speed command and
   !                    initial synchronous angle (global coordinates).
   !  Tg0             : Initial generator torque command.
   !  y0,ydot0        : Initial position and linear speed of the rotor
   !                    head.
   !  psi0hin,omegah0,: Initial azimuth and rotational speed of head
   !  psi0gin,omegag0   and generator.
   !  xw0             : Initial position of the rotor center in the 
   !                    turbulent windfield.
   !  Vref            : The reference mean windspeed for the analysis.
   !  z0,h0           : Height of rotor center and surface roughness.
   !  xvane_g         : Position of the wind vane measuring wind 
   !                    direction, in global coordinates.  (Though it
   !                    follows motion of the platform.)
   !
   !  Outputs:
   !  --------
   !  (Writes to file.)
   !  ierr            : set to a negative number if an error occurs.
   !
   !  Other variables:
   !  ----------------
   !  wystore         : An array of (omega,psi) points encountered over
   !                    the previous revolution.

!      use VAWTBEM

      implicit none

      character(len=*),                  intent(in)    :: setupFile,airfoilFile, &
                                                          turbFile
      integer,                           intent(in)    :: setupFid,airfoilFid, &
                                                          turbFid,fid
      integer,                           intent(in)    :: Nb,Nrows,Ncols,Ndof,Nadof,Nt
      double precision,                  intent(in)    :: dt,wgset0,psisyn0,     &
                                                          Tgset0,Vref,           &
                                                          psi0hin,omegah0,       &
                                                          psi0gin,omegag0,z0,h0
      double precision, dimension(3),    intent(in)    :: y0,ydot0,xw0,xvane_g

      integer,                           intent(inout) :: ierr

      integer                                          :: istep
      double precision                                 :: t
      double precision, dimension(Ndof)                :: x0,x,dxdt      
      integer,          dimension(Nrows,Ncols)         :: QSflag0,QSflag

      !  Geometry and airfoil parameters.
      character(len=MAXSTRLEN),dimension(Nrows)        :: airfoilName
      double precision                                 :: density,viscosity,Lblade,ab, &
                                                          ar,Ro
      double precision, dimension(Nrows)               :: Lelem,Zelem,chord,twist0
      double precision, dimension(3,Nrows)             :: profile
      double precision, dimension(3,Nrows+1)           :: node
      double precision, dimension(3,Nrows,Ncols)       :: xs_r,n_r
      double precision, dimension(3,3,Nrows,Ncols)     :: Tsr_s
      integer,          dimension(Nrows)               :: NRe
      integer,          dimension(MRe,Nrows)           :: Naoa
      double precision, dimension(MRe,Nrows)           :: Rea
      double precision, dimension(Maoa,MRe,Nrows)      :: aoaa,Cla,Cda,Cma
      double precision, dimension(Nrows,Ncols)         :: aoaz,pitchAng
      double precision, dimension(2,Nrows,Ncols)       :: aoafs,dClda_a

      !  Other DeepwindTimestep variables.
      double precision, dimension(3,Nrows,Ncols/2)     :: Vinf_g
      double precision, dimension(Nrows,Ncols)         :: aoa,fqs
      double precision, dimension(3,Nrows,Ncols)       :: Viqs,Viqsprev,dViqsdt,V0_r
      double precision, dimension(6,Nrows,Ncols)       :: Fb_r
      double precision, dimension(6)                   :: F_g
      double precision                                 :: Tg
      integer, parameter                               :: Nwy = 360
      double precision, dimension(Nwy)                 :: wystore
      integer,          dimension(3)                   :: Nbox
      double precision, dimension(3)                   :: dL,V0avg_g
      double precision, dimension(3,3)                 :: Tg_w

      !  Miscellaneous.
      integer                                          :: irow,icol,N,NrNc,j,k,irecl
      double precision                                 :: maxval,Re,val,lnz0h0,R,Vmag
      integer, parameter :: ikind = SELECTED_REAL_KIND(p=4,r=2)
      real(kind=ikind)                                 :: junk
double precision, dimension(6) :: Fb_a

      !  Initialize the analysis...

      !  At the moment there is no active aerodynamic control.  Set
      !  pitch angle to a constant zero.
      pitchAng = 0.d0  !  Radians.

      !  Read geometry from input file.
      call readSetupParameters (setupFid,setupFile,Nrows,Ncols,    &
                                density,viscosity,                 &
                                node,profile,Lelem,Lblade,Zelem,   &
                                chord,twist0,xs_r,n_r,Tsr_s,       &
                                airfoilName,ierr)

!write (fid,'(1X,2ES13.4)') density,viscosity
!write (fid,'(1X,A)') 'node1-3 profile1-3 Lelem Zelem'
!do irow = 1,Nrows+1
!   if (irow .eq. 1) then
!      write (fid,'(I4,3F10.3)') irow-1,node(1,irow),node(2,irow),node(3,irow)
!   else 
!      write (fid,'(I4,8F10.3)') irow-1,node(1,irow),node(2,irow),node(3,irow),         &
!                                profile(1,irow-1),profile(2,irow-1),profile(3,irow-1), &
!                                Lelem(irow-1),Zelem(irow-1)
!   end if
!end do
!write (fid,'(1X,A)') 'chord twist0 xs_r1-3 n_r1-3'
!do irow = 1,Nrows
!      write (fid,'(I4,8F10.3)') irow,chord(irow),twist0(irow),                &
!                                xs_r(1,irow,9),xs_r(2,irow,9),xs_r(3,irow,9), &
!                                n_r(1,irow,9),n_r(2,irow,9),n_r(3,irow,9)
!end do
!write (fid,'(1X,A)') 'airfoilName'
!do irow = 1,Nrows
!      write (fid,'(I4,A,A)') irow,' ',airfoilName(irow)
!end do
! Verified.

      !  Calculate aspect ratio and outer radius.  Aspect ratio is
      !  used in coefficient calculation, while outer radius is used
      !  in dynamic inflow, to set the time constant.
      ab = 0.d0
      maxval = 0.d0
      do irow = 1,Nrows
         ab = ab + chord(irow)*Lelem(irow)
         if (profile(1,irow) .gt. maxval) then
            maxval = profile(1,irow)
         end if
      end do
      ar = (Lblade**2)/ab
      Ro = maxval

      !  Read airfoil coefficients from file.
      call readCoefficientLibrary (airfoilFid,airfoilFile,airfoilName,  &
                                   Nrows,ar,                            &
                                   NRe,Naoa,                            &
                                   Rea,aoaa,Cla,Cda,Cma,                &
                                   ierr)

      !  Incoming mean windspeed is based upon the reference windspeed,
      !  adjusted to include wind shear.
      !
      !  By definition, the mean incoming wind will start pointing in
      !  the X^g direction.
      lnz0h0 = log(z0/h0)
      do irow = 1,Nrows
         val = Vref*log((z0 + xs_r(3,irow,1))/h0)/lnz0h0
         do icol = 1,Ncols/2
            Vinf_g(1,irow,icol) = val
            Vinf_g(2,irow,icol) = 0.d0
            Vinf_g(3,irow,icol) = 0.d0
         end do
      end do

      !  Initialize the state variables...
     
      !  First, QSflag begins at 0 (no special considerations for dynamic
      !  stall) for all elements.
      QSflag0 = 0

      !  The initial generator torque and speed commands are input.
      x0(1) = Tgset0
      x0(2) = wgset0
      x0(3) = psisyn0  !  Initial synchronous angle for synchronous generator.
                       !  (Meaningless for induction generator.)

      !  Initial position, velocity, and rotation speed are input.
      x0(4:6) = y0
      x0(7) = psi0hin
      x0(8) = psi0gin
      x0(9:11) = ydot0
      x0(12) = omegah0
      x0(13) = omegag0
      
      !  Initial position in the turbulent windfield.
      x0(14:16) = xw0
      !  Initial wind direction measured.
      x0(17) = 0.d0

      N = Ndof - Nadof
      NrNc = Nrows*Ncols

      !  Separation point position is initialized to zero
      !  (attached flow).
      x0(N+1:NrNc+N) = 0.d0

      !  Induced velocity can be initialized to zero.  (Might
      !  eventually be better to load this from file, in order
      !  to reduce transient time?)
      x0(NrNc+N+1:7*NrNc+N) = 0.d0
      Viqsprev = 0.d0
      dViqsdt = 0.d0

      !  Initialize airfoil parameters needed in the dynamic
      !  stall calculation.  Assign a characteristic Reynolds
      !  number and do not further update the coefficients and
      !  parameters as a function of Re.  [This can be changed
      !  later if needed.]
      do icol = 1,Ncols
         do irow = 1,Nrows
            R = sqrt(xs_r(1,irow,icol)**2 + xs_r(2,irow,icol)**2)
            Vmag = sqrt((R*omegah0)**2 + Vref**2)
            Re = density*chord(irow)*Vmag/viscosity
            call initializeDynamicStall (Re,NRe(irow),Naoa(:,irow),                  &
                                         Rea(:,irow),aoaa(:,:,irow),                 &
                                         Cla(:,:,irow),Cda(:,:,irow),Cma(:,:,irow),  &
                                         aoaz(irow,icol),dClda_a(:,irow,icol),       &
                                         aoafs(:,irow,icol),ierr)
            if (ierr .ne. 0) then
               return
            end if
!write (41,'(1X,2I4,5ES13.4,A)')                             &
!irow,icol,aoaz(irow,icol),dClda_a(1,irow,icol),             &
!dClda_a(2,irow,icol),aoafs(1,irow,icol),aoafs(2,irow,icol), &
!'   irow icol aoaz dClda_a+- aoafs+-'
! Verified.
         end do
      end do

      !  Initialize variables related to the turbulent wind field.
      !  These are:
      !  Nbox         : The number of grid points in the Xw, Yw,
      !                 and Zw directions (wind coordinate system).
      !  dL           : Separation between grid points (Xw,Yw,Zw).
      !  Tg_w         : Rotational transform matrix taking a vector
      !                 from global to wind coordinates.  Currently
      !                 hard-coded as the identity matrix.
      inquire (IOLENGTH=irecl) junk,junk,junk
      open (UNIT=turbfid,FILE=turbFile,ACCESS='DIRECT',STATUS='OLD', &
            ACTION='READ',RECL=irecl)
      read(turbfid,REC=1) Nbox(1)
      read(turbfid,REC=2) Nbox(2)
      read(turbfid,REC=3) Nbox(3)
      read(turbfid,REC=4) dL(1)
      read(turbfid,REC=5) dL(2)
      read(turbfid,REC=6) dL(3)
      Tg_w(1,1) = 1.d0
      Tg_w(1,2) = 0.d0
      Tg_w(1,3) = 0.d0
      Tg_w(2,1) = 0.d0
      Tg_w(2,2) = 1.d0
      Tg_w(2,3) = 0.d0
      Tg_w(3,1) = 0.d0
      Tg_w(3,2) = 0.d0
      Tg_w(3,3) = 1.d0

      t = 0.d0
      istep = 0
      do istep = 1,Nt

         if (mod(istep-1,100) .eq. 0) then
            print '(F12.4)',t
         end if

         !  Evaluate the state variable derivatives at the current timestep.
         !  The outputs from this call should be considered the "authoritative"
         !  outputs for the beginning of the timestep, at time t... that is,
         !  the loads that should be output to file.
         call DeepwindTimestep (1,QSflag0,x0,dViqsdt,                   &
                                Nb,Nrows,Ncols,Ndof,Nadof,              &
                                t,Vref,Vinf_g,                          &
                                turbFid,Nbox,dL,Tg_w,xvane_g,           &
                                twist0,node,profile,Lblade,Ro,xs_r,n_r, &
                                density,viscosity,                      &
                                chord,Zelem,Lelem,                      &
                                NRe,Naoa,Rea,aoaa,Cla,Cda,Cma,          &
                                Tsr_s,aoaz,aoafs,dClda_a,pitchAng,      &
                                QSflag,dxdt,aoa,fqs,Viqs,V0_r,V0avg_g,  &
                                Fb_r,F_g,Tg,                            &
                                ierr)
         if (ierr .ne. 0) then
            return
         end if

         !  Initialize induced velocities.
         if (istep .eq. 1) then
            Viqsprev = Viqs
            k = N + Nrows*Ncols
            j = N + 4*Nrows*Ncols
            do icol = 1,Ncols
               do irow = 1,Nrows
                  x0(k+1) = Viqs(1,irow,icol)
                  x0(k+2) = Viqs(2,irow,icol)
                  x0(k+3) = Viqs(3,irow,icol)
                  x0(j+1) = Viqs(1,irow,icol)
                  x0(j+2) = Viqs(2,irow,icol)
                  x0(j+3) = Viqs(3,irow,icol)
                  k = k + 3
                  j = j + 3
               end do
            end do
         end if

!print '(1X,I3,13ES12.3)',1,x0(1:13)

irow = 10
icol = 9
k = N + (icol+9)*Nrows + irow
j = N + Nrows*Ncols + 3*((icol-1)*Nrows + irow - 1)
write (41,'(1X,10F10.4,F10.3,F10.4,F10.2,F9.4,8ES12.3,A)')       &
t,V0avg_g,V0_r(:,irow,icol),x0(j+1),x0(j+2),x0(j+3),             &
x0(4),x0(9),mod(x0(8)*180.d0/PI,360.d0),x0(13),F_g,Tg,Tg*x0(13), &
'   t V0avg_g V0_r Vi X Xdot psi0g wg F_g(1:6) Tg Pg'

!write (44,'(1X,F9.3,3F10.4,F9.3,8ES12.3,A)')            &
!t,V0avg_g,mod(x0(8)*180.d0/PI,360.d0),F_g,Tg,Tg*x0(13), &
!'   t V0avg_g azi F_g(1:6) Tg Pg'

!irow = 20
!icol = 1
!k = N + Nrows*Ncols + 3*((icol-1)*Nrows + irow - 1)
!print '(1X,2I4,6ES12.3)',irow,icol,x0(k+1),x0(k+2),x0(k+3),dxdt(k+1),dxdt(k+2),dxdt(k+3)

!write (41,'(1X,8ES13.4)') &
!t,F_g,F_g(6)*x0(13)

write (42,'(1X,2ES13.4)') t,F_g(1)
write (43,'(1X,2ES13.4)') t,F_g(6)

k = N + irow
write (44,'(1X,F8.3,F8.2,F9.3,I4,2F8.4,2ES13.4,A)')                         &
 t,mod((x0(7) - PI/2.d0)*180.d0/PI,360.d0),aoa(irow,1),QSflag0(irow,1),     &
fqs(irow,1),x0(k),Fb_r(1,irow,1)*cos(x0(7)) + Fb_r(2,irow,1)*sin(x0(7)),    &
profile(1,irow)*(-Fb_r(1,irow,1)*sin(x0(7)) + Fb_r(2,irow,1)*cos(x0(7))),   &
'   t psi0 aoa QSflag0 fqs f Fn Torque'

         !  The following block of code is initimately connected with the
         !  first call to DeepwindTimestep, at the initial x for the 
         !  timestep.  This fixes the starting value of the separation
         !  point state variables if certain logic is met.
         k = N
         do icol = 1,Ncols
            do irow = 1,Nrows
               k = k + 1

               !  Special dynamic stall logic which requires artificially 
               !  setting f (state variable) outside the integration scheme.
               !
               !  See oye subroutine for a detailed description.  Note that
               !  the old (not updated) QSflag must be used in the logic here.
               if (((QSflag0(irow,icol) .eq. -1) .and.           &
                    (aoa(irow,icol) .gt. aoaz(irow,icol)) .and.  &
                    (x0(k) .gt. 0.d0)) .or.                      &
                   ((QSflag0(irow,icol) .eq. 1) .and.            &
                    (aoa(irow,icol) .lt. aoaz(irow,icol)) .and.  &
                    (x0(k) .lt. 0.d0))) then         
                  !  This means that f has "rounded the corner" at the trailing
                  !  edge.
                  x0(k) = fqs(irow,icol)
                  dxdt(k) = 0.d0
               end if
            end do
         end do

         !  Let the result of the initial evaluation of DeepwindTimestep
         !  set QSflag0 for the next evaluations in rk4 or modEuler.
         QSflag0 = QSflag

         !  Estimate the rate of change of quasi-static induced velocity
         !  based on a simple finite-difference approximation, for use in
         !  the next timestep.  This is used in the dynamic inflow
         !  calculation.  The timescale of dynamic inflow is slow, so
         !  it should not make a large difference that a coarse finite
         !  difference formula is used to estimate the rate of change.
         if (istep .gt. 1) then
            dViqsdt = (Viqs - Viqsprev)/dt
         else
            dViqsdt = 0.d0
         end if
         Viqsprev = Viqs

         !  Integrate to update state variables.
         call rk4 (x0,dxdt,t,dt,x,                         &
                   QSflag0,dViqsdt,                        &
                   Nb,Nrows,Ncols,Ndof,Nadof,              &
                   Vref,Vinf_g,                            &
                   turbFid,Nbox,dL,Tg_w,xvane_g,           &
                   twist0,node,profile,Lblade,Ro,xs_r,n_r, &
                   density,viscosity,                      &
                   chord,Zelem,Lelem,                      &
                   NRe,Naoa,Rea,aoaa,Cla,Cda,Cma,          &
                   Tsr_s,aoaz,aoafs,dClda_a,pitchAng,      &
                   ierr)
         if (ierr .ne. 0) then
            return
         end if

         !  Update value of x0 for the beginning of the next timestep.
         x0 = x
         t = t + dt

      end do

      close (turbFid)

   end subroutine simulateDeepwind









   subroutine rk4 (x0,dxdt,t,dt,x,                         &
                   QSflag0,dViqsdt,                        &
                   Nb,Nrows,Ncols,Ndof,Nadof,              &
                   Vref,Vinf_g,                            &
                   turbFid,Nbox,dL,Tg_w,xvane_g,           &
                   twist0,node,profile,Lblade,Ro,xs_r,n_r, &
                   density,viscosity,                      &
                   chord,Zelem,Lelem,                      &
                   NRe,Naoa,Rea,aoaa,Cla,Cda,Cma,          &
                   Tsr_s,aoaz,aoafs,dClda_a,pitchAng,      &
                   ierr)

      implicit none

      double precision,                             intent(in) :: t,dt

      integer,                                      intent(in) :: Nb,Nrows,Ncols,Ndof, &
                                                                  Nadof

      integer,          dimension(Nrows,Ncols),     intent(in) :: QSflag0
      double precision, dimension(Ndof),            intent(in) :: x0
      double precision, dimension(Ndof),            intent(in) :: dxdt
      double precision, dimension(3,Nrows,Ncols),   intent(in) :: dViqsdt

      double precision,                             intent(in) :: Lblade,Ro,          &
                                                                  density, viscosity
      double precision, dimension(Nrows),           intent(in) :: twist0,chord,Lelem,Zelem
      double precision, dimension(3,Nrows),         intent(in) :: profile
      double precision, dimension(3,Nrows+1),       intent(in) :: node
      double precision, dimension(Nrows,Ncols),     intent(in) :: aoaz,pitchAng
      double precision, dimension(2,Nrows,Ncols),   intent(in) :: aoafs,dClda_a
      double precision, dimension(3,Nrows,Ncols/2), intent(in) :: Vinf_g
      double precision, dimension(3,Nrows,Ncols),   intent(in) :: xs_r,n_r
      double precision, dimension(3,3,Nrows,Ncols), intent(in) :: Tsr_s

      integer,          dimension(Nrows),           intent(in) :: NRe
      integer,          dimension(MRe,Nrows),       intent(in) :: Naoa
      double precision, dimension(MRe,Nrows),       intent(in) :: Rea
      double precision, dimension(Maoa,MRe,Nrows),  intent(in) :: aoaa,Cla,Cda,Cma

      integer,                                      intent(in) :: turbFid
      integer,          dimension(3),               intent(in) :: Nbox
      double precision,                             intent(in) :: Vref
      double precision, dimension(3),               intent(in) :: dL,xvane_g
      double precision, dimension(3,3),             intent(in) :: Tg_w

      double precision, dimension(Ndof),           intent(out) :: x

      integer,                                   intent(inout) :: ierr

      !  (No need to save intermediate outputs from DeepwindTimestep.)
!integer :: i,j,k,irow,icol
      integer,          dimension(Nrows,Ncols)                 :: QSflag
      double precision, dimension(Nrows,Ncols)                 :: aoa,fqs
      double precision, dimension(3)                           :: V0avg_g
      double precision, dimension(3,Nrows,Ncols)               :: Viqs,V0_r
      double precision, dimension(6,Nrows,Ncols)               :: Fb_r
      double precision, dimension(6)                           :: F_g
      double precision                                         :: Tg

      double precision                                         :: dt2,dt6,t2
      double precision, dimension(Ndof)                        :: xt2,dxdt2,dxdtm

!  Need to verify against Press et al.; variable names changed.
      dt2 = 0.5d0*dt
      dt6 = dt/6.d0
      t2 = t + dt2
      xt2 = x0 + dt2*dxdt
      call DeepwindTimestep (0,QSflag0,xt2,dViqsdt,                  &
                             Nb,Nrows,Ncols,Ndof,Nadof,              &
                             t2,Vref,Vinf_g,                         &
                             turbFid,Nbox,dL,Tg_w,xvane_g,           &
                             twist0,node,profile,Lblade,Ro,xs_r,n_r, &
                             density,viscosity,                      &
                             chord,Zelem,Lelem,                      &
                             NRe,Naoa,Rea,aoaa,Cla,Cda,Cma,          &
                             Tsr_s,aoaz,aoafs,dClda_a,pitchAng,      &
                             QSflag,dxdt2,aoa,fqs,Viqs,V0_r,V0avg_g, &
                             Fb_r,F_g,Tg,                            &
                             ierr)
      if (ierr .ne. 0) then
         return
      end if
      xt2 = x0 + dt2*dxdt2
      call DeepwindTimestep (0,QSflag0,xt2,dViqsdt,                  &
                             Nb,Nrows,Ncols,Ndof,Nadof,              &
                             t2,Vref,Vinf_g,                         &
                             turbFid,Nbox,dL,Tg_w,xvane_g,           &
                             twist0,node,profile,Lblade,Ro,xs_r,n_r, &
                             density,viscosity,                      &
                             chord,Zelem,Lelem,                      &
                             NRe,Naoa,Rea,aoaa,Cla,Cda,Cma,          &
                             Tsr_s,aoaz,aoafs,dClda_a,pitchAng,      &
                             QSflag,dxdtm,aoa,fqs,Viqs,V0_r,V0avg_g, &
                             Fb_r,F_g,Tg,                            &
                             ierr)
      if (ierr .ne. 0) then
         return
      end if
      xt2 = x0 + dt*dxdtm
      dxdtm = dxdt2 + dxdtm
      call DeepwindTimestep (0,QSflag0,xt2,dViqsdt,                  &
                             Nb,Nrows,Ncols,Ndof,Nadof,              &
                             t+dt,Vref,Vinf_g,                       &
                             turbFid,Nbox,dL,Tg_w,xvane_g,           &
                             twist0,node,profile,Lblade,Ro,xs_r,n_r, &
                             density,viscosity,                      &
                             chord,Zelem,Lelem,                      &
                             NRe,Naoa,Rea,aoaa,Cla,Cda,Cma,          &
                             Tsr_s,aoaz,aoafs,dClda_a,pitchAng,      &
                             QSflag,dxdt2,aoa,fqs,Viqs,V0_r,V0avg_g, &
                             Fb_r,F_g,Tg,                            &
                             ierr)
      if (ierr .ne. 0) then
         return
      end if
      x = x0 + dt6*(dxdt + dxdt2 + 2.d0*dxdtm)

   end subroutine rk4










   subroutine DeepwindTimestep (debug,                                     &
                                QSflag0,x,dViqsdt,                         &
                                Nb,Nrows,Ncols,Ndof,Nadof,                 &
                                t,Vref,Vinf_g,                             &
                                turbFid,Nbox,dL,Tg_w,xvane_g,              &
                                twist0,node,profile,Lblade,Ro,xs_r,n_r,    &
                                density,viscosity,                         &
                                chord,Zelem,Lelem,                         &
                                NRe,Naoa,Rea,aoaa,Cla,Cda,Cma,             &
                                Tsr_s,aoaz,aoafs,dClda_a,pitchAng,         &
                                QSflag,dxdt,aoa,fqs,Viqs,V0_r,V0avg_g,     &
                                Fb_r,F_g,Tg,ierr)
   !  This subroutine calculates the state of the degrees-of-freedom
   !  describing the Deepwind aerodynamics, platform, and controller.
   !  The states of the platform are described by the usual mass-
   !  spring-damper equation.  However, the applied forces are nonlinear
   !  functions of the aerodynamics and the control system inputs.
   ! 
   !  States are contained in a vector x.  The states are:
   !  indices                   state
   !  -------                   -----
   !  1                         Generator torque command Tgset.
   !  2                         Generator (synchronous) speed set-point wgset.
   !  3                         Generator synchronous angle psisyn.
   !  [-]                       Pitch angle command [or other aero control]
   !  4 to 13                   y(3),psi0h,psi0g,ydot(3),omegah,omegag where:
   !                               y is the linear position of the center of the
   !                               rotor,
   !                               ydot is the linear velocity of the rotor,
   !                               psi0h is the azimuth of the rotor,
   !                               omegah is the actual rotor (head) speed,
   !                               psi0g is the generator rotor azimuth, and
   !                               omegag is the actual generator speed.
   !  14 to 16                  X,Y,Z coordinates in the wind coordinate system
   !                            (typically in the same orientation as the global
   !                            coordinate system, but with an X offset) of the
   !                            rotor center: the origin of the rotor coordinate
   !                            system.  In other words, where does the rotor 
   !                            center lie within the turbulent windfield?
   !  17                        Measured wind direction, in rotor coordinates.
   !  18 to Nr*Nc+17            f, the separation point position at each blade
   !                            element.
   !  Nr*Nc+18 to 4*Nr*Nc+17    Vi_s for each surface element.
   !  4*Nr*Nc+18 to 7*Nr*Nc+17  Viint_s  "    "    "    "    .
   !
   !  There is also an integer array QSflag that does not function as
   !  a continuous state, but is a discrete state which is used in the
   !  dynamic stall analysis to give a smooth output.  This includes
   !  a jump in f to a new value, although outputs remain smooth.  Thus
   !  f is to be overwritten (assigned a value fqs) before numerical
   !  integration, on timesteps when it jumps.  df/dt will naturally
   !  evaluate to zero on these timesteps, when f is set to fqs.
   !  
   !
   !  Inputs:
   !  -------
   !  QSflag0          : Flag for special dynamic stall behavior.
   !  x                : State vector, as described above.
   !  dViqsdt          : Rate of change of quasi-steady induced
   !                     velocity, used in the dynamic inflow 
   !                     calculation.
   !  Nb,Nrows,Ncols   : Number blades, rows and cols of surface 
   !                     elements.
   !  Ndof,Nadof       : Length of state variable vector.
   !  t                : The present time (for turbulence look-up).
   !  Vref             : Nominal incoming windspeed.
   !  Vinf_g           : remote velocity, including wind shear,
   !                     in the global coordinate system.  Defined
   !                     at each surface element centroid, assuming
   !                     that the rotor is oriented in the global
   !                     coordinate system.
   !  turbFid          : File ID of turbulence file.
   !  Nbox             : Number of turbulence grid points in the Xw,
   !                     Yw, and Zw directions (wind coordinates).
   !  dL               : Separation between grid points (Xw,Yw,Zw).
   !  Tg_w             : Transform matrix from global to wind coords.
   !  xvane_g          : Position of the wind vane measuring wind 
   !                     direction, in global coordinates. (Though the
   !                     vane location follows motion of the platform, 
   !                     it does not follow the mean wind direction, 
   !                     so it cannot be given in rotor coordinates.)
   !  twist0           : Blade twist (rad).
   !  node             : Nodal r-z-offset.
   !  profile          : Element center r-z-offset.
   !  Lblade           : Total blade length.
   !  Ro               : Outer radius.
   !  xs_r             : Surface element centroid coordinates.
   !  n_r              : Surface element normal vectors.
   !  density,viscosity
   !  chord            : Blade element chord.
   !  Lelem            : Element length in span direction.
   !  Zelem            : Along-span distance from equator.
   !  NRe ... Cma      : Airfoil coefficient data.
   !  Tsr_s            : Transform from surface to rotor coordinates.
   !  aoaz,aoafs,      : Zero-lift AOA, AOA at full stall (f = 1), and
   !  dClda_a            attached-flow lift coefficient slope.  Used in
   !                     the dynamic stall calculation.
   !  pitchAng         : Blade pitch angle.
   !  
   !
   !  Outputs:
   !  --------
   !  QSflag           : Updated QSflag.
   !  dxdt             : Rate of change of f, Vi, and Viint.
   !  aoa,fqs          : Variables for dynamic stall logic.
   !  Viqs             : Quasi-steady induced velocity, for dynamic inflow.
   !  V0_r             : Incoming velocity, for Vi truncation.
   !  V0avg_g          : Rotor-average velocity in global coordinates.
   !                     Includes Vinf, rotor movement, and turbulence.
   !  Fb_r             : Local aero forces on each blade element.
   !  F_g              : Aerodynamic forces in global coordinates.
   !  Tg               : Generator torque.
   !

      implicit none

      integer,                                      intent(in) :: Nb,Nrows,Ncols,Ndof, &
                                                                  Nadof,debug

      integer,          dimension(Nrows,Ncols),     intent(in) :: QSflag0
      double precision, dimension(Ndof),            intent(in) :: x
      double precision, dimension(3,Nrows,Ncols),   intent(in) :: dViqsdt

      double precision,                             intent(in) :: Lblade,Ro,t,   &
                                                                  density,viscosity
      double precision, dimension(Nrows),           intent(in) :: twist0,chord,Lelem,Zelem
      double precision, dimension(3,Nrows),         intent(in) :: profile
      double precision, dimension(3,Nrows+1),       intent(in) :: node
      double precision, dimension(Nrows,Ncols),     intent(in) :: aoaz,pitchAng
      double precision, dimension(2,Nrows,Ncols),   intent(in) :: aoafs,dClda_a
      double precision, dimension(3,Nrows,Ncols/2), intent(in) :: Vinf_g
      double precision, dimension(3,Nrows,Ncols),   intent(in) :: xs_r,n_r
      double precision, dimension(3,3,Nrows,Ncols), intent(in) :: Tsr_s

      !                             Airfoil properties.
      integer,          dimension(Nrows),           intent(in) :: NRe
      integer,          dimension(MRe,Nrows),       intent(in) :: Naoa
      double precision, dimension(MRe,Nrows),       intent(in) :: Rea
      double precision, dimension(Maoa,MRe,Nrows),  intent(in) :: aoaa,Cla,Cda,Cma

      !                             Turbulence.
      integer,                                      intent(in) :: turbFid
      integer,          dimension(3),               intent(in) :: Nbox
      double precision,                             intent(in) :: Vref
      double precision, dimension(3),               intent(in) :: dL,xvane_g
      double precision, dimension(3,3),             intent(in) :: Tg_w

      integer,          dimension(Nrows,Ncols),    intent(out) :: QSflag
      double precision, dimension(Ndof),           intent(out) :: dxdt
      double precision, dimension(Nrows,Ncols),    intent(out) :: aoa,fqs
      double precision, dimension(3),              intent(out) :: V0avg_g
      double precision, dimension(3,Nrows,Ncols),  intent(out) :: Viqs,V0_r
      double precision, dimension(6,Nrows,Ncols),  intent(out) :: Fb_r
      double precision, dimension(6),              intent(out) :: F_g
      double precision,                            intent(out) :: Tg

      integer,                                   intent(inout) :: ierr

integer :: irow,icol,i,j,k
      integer :: N
      double precision :: kgen,Jrotor,Jgen,ktower,dtower,Tcontrol, &
                          dirvane,psivr,psirnorm,mcart,kcart,ccart,wcart
      double precision, dimension(3) :: Vh_g,xw0,uvane
      double precision, dimension(Nadof) :: xaero,dxdtaero

      N = Ndof - Nadof

      Vh_g = x(9:11)
 
      !  Call the aerodynamic module to update aero loads and get
      !  aero state derivatives.
      xaero = x(N+1:Nadof+N)
      xw0 = x(14:16)
      call aeroTimestep (debug,QSflag0,xaero,dViqsdt,x(7),x(12),    &
                         Nb,Nrows,Ncols,Nadof,t,                    &
                         Vinf_g,Vh_g,turbFid,Nbox,dL,Tg_w,xw0,      &
                         twist0,node,profile,                       &
                         Lblade,Ro,xs_r,n_r,density,viscosity,      &
                         chord,Lelem,Zelem,                         &
                         NRe,Naoa,Rea,aoaa,Cla,Cda,Cma,Tsr_s,       &
                         aoaz,aoafs,dClda_a,pitchAng,               &
                         QSflag,dxdtaero,aoa,fqs,Viqs,              &
                         V0_r,V0avg_g,Fb_r,F_g,ierr)
      dxdt(N+1:Nadof+N) = dxdtaero

!irow = 10
!icol = 9
!k = (icol+9)*Nrows + irow
!j = Nrows*Ncols + 3*(icol*Nrows + irow - 1)
!i = 4*Nrows*Ncols + 3*(icol*Nrows + irow - 1)
!write (41,'(1X,A,F8.2,2I3,2F7.3,ES12.3,5F8.2,F8.3,6ES12.3,A)')  &
!'aero:',t,QSflag0(irow,icol),QSflag(irow,icol),                 &
!fqs(irow,icol+9),xaero(k),dxdtaero(k),                          &
!Viqs(3,irow,icol),dViqsdt(3,irow,icol),                         &
!xaero(j+3),dxdtaero(j+3),                                       &
!x(8)*180.d0/PI,x(13),                                           &
!F_g,                                                            &
!'   t QS0 QS fqs f dfdt Viqs dViqs/dt Vi dVi/dt psi0g wg F_g(6) Tg Pg'

      !  Simulate a wind vane on top of the turbine tower.
      !
      !  [At the moment, the vane location follows the translational motion
      !  of the platform, but does not follow the tilt.] 
      !
      !  uvane is the instantaneous measurement of wind by the vane.  This
      !  cannot be fed directly into the state variable, because the
      !  numerical integration scheme might not allow for instantaneous
      !  updating of state variables.  Therefore, implement a time lag
      !  of a short duration.
      call readTurbulence (turbFid,dL,Nbox,x(14:16),xvane_g,Tg_w,uvane)
      dirvane = atan2(uvane(2),Vref + uvane(1))
      dxdt(17) = (dirvane - x(17))/0.5d0
      psivr = x(7) - x(17)  !  Estimate of rotor azimuth relative to wind
                            !  (from vane measured direction to rotor
                            !  reference). 
      psivr = mod(psivr,2.d0*PI)

      !  Call the generator module to update generator torque based
      !  on current generator speed and commanded torque.
      !  
      !  [For now, let the generator behave as an induction generator
      !  with respect to a constant set-point.]
      !  Induction generator (Hau p 326)
!      kgen = (5.d6/x(2))/(0.01d0*x(2))  !  dT/dw = d(P/w)/dw, 1% slip.
!      Tg = kgen*(x(13) - x(2))
      !  Synchronous generator (Hau p 323)
!      kgen = (5.d6/x(2))/0.5236d0  !  dT/da = d(P/w)/da, 30-deg load angle.
!      Tg = kgen*(x(8) - x(3))
      !  Simple controlled generator: let Tg = Tcontrol with no delay.
      !  (Tg = Tcontrol is implemented below after lookup.)

      !  Call the control module to get control state derivatives.
      call genT (psivr,x(13),Tcontrol,ierr)
      if (ierr .ne. 0) then
print *,'Error during control torque calculation.'
         return
      end if
      Tg = Tcontrol

      !  Const speed, no dynamics, generator absorbs aero torque.
!      Tcontrol = F_g(6) 
!      Tg = Tcontrol

!if (debug .eq. 1) then
!write (41,'(1X,3F9.4,2F9.3,F10.3,F9.3,F9.4,ES12.3,A)')             &
!Vref+uvane(1),uvane(2),uvane(3),dirvane*180.d0/PI,x(17)*180.d0/PI, &
!x(7)*180.d0/PI,psivr*180.d0/PI,x(13),Tg,                           &
!'   uvane_x,y,z dirvane_qs dirvane azi psivr omega Tg'
!end if

      !  It is necessary to update the control variables using a
      !  proportional control scheme, because the numerical integration
      !  methods being used may not allow for discrete updating of
      !  state variables.  The idea here is to set the gain large enough
      !  (the timescale small enough) that the delay in the control signal
      !  does not significantly influence the dynamics.  As long as
      !  the speed control signal changes much more quickly than the speed,
      !  it should be OK.
      !
      !  [Note: I have bypassed this by setting Tg directly to Tcontrol.]
      dxdt(1) = (Tcontrol - x(1))/0.5d0
      !  No control, simple induction generator.
!      dxdt(1) = 0.d0
      !  wgset contains the estimate of revolution-average speed.
      !  Update with a short time constant.
      !  [Currently bypassed since input x(13) is used directly in the
      !  lookup table.]
      dxdt(2) = (x(13) - x(2))/0.5d0
      !  Constant wgset
!      dxdt(2) = 0.d0

      !  Call the platform module to get platform state derivatives.
!      call platform ()
      dxdt(3) = x(2)                !  d/dt(syn. angle) = wgset (assume instant response).
      dxdt(4) = x(9)                !  dy/dt = ydot
      dxdt(5) = x(10)               !   
      dxdt(6) = x(11)               !
      dxdt(7) = x(12)               !  dpsi/dt = omega
      dxdt(8) = x(13)               !

      !  Put the platform on a "cart with wheels".
      mcart = 5.d6
      kcart = 7.9d6
      ccart = 0.05d0*2.d0*sqrt(kcart*mcart)
!      dxdt(9) = (F_g(1) - ccart*x(9) - kcart*x(4))/mcart  !  Linear acceleration of platform.
!      wcart = 2.d0*PI*0.02d0
!      dxdt(9) = 1.d0*wcart*cos(wcart*t)
dxdt(9) = 0.d0
      dxdt(10) = 0.d0               !
      dxdt(11) = 0.d0               !

      Jrotor = 2.48d8               !  Svendsen, simplified model report
      Jgen = 2.17d7                 !  Jgen includes tower.
!  This code can be used for a flexible shaft.
!      ktower = 8.07d9
!      dtower = 2.d0*0.01d0*sqrt(ktower*Jgen)  !  d = 2 zeta sqrt(km).
!                                    !  Rotational acceleration of head.
!      dxdt(12) = (1.d0/Jrotor)                                               &
!               * ((-ktower*x(7) - dtower*x(12) + ktower*x(8) + dtower*x(13)) &
!               + F_g(6))
!                                    !  Rotational acceleration of generator.
!      dxdt(13) = (1.d0/Jgen)                                                &
!               * ((ktower*x(7) + dtower*x(12) - ktower*x(8) - dtower*x(13)) &
!               - Tg)


!  This code can be used for a rigid shaft.
      dxdt(13) = (F_g(6) - Tg)/(Jrotor + Jgen)
      dxdt(12) = dxdt(13)           !  Rotational acceleration of head.  Deactivate shaft.


!  This code can be used for a constant speed.
!      dxdt(13) = 0.d0
!      dxdt(12) = dxdt(13)           !  Rotational acceleration of head.  Deactivate shaft.


!if (debug .eq. 1) then
!write (41,'(1X,3ES12.3,A)') &
!F_g(6),Tg,dxdt(13),         &
!'   Taero Tg omegadot'
!end if

      !  Motion of the rotor center relative to the wind coordinate
      !  system is the sum of head motion and nominal incoming 
      !  velocity.
      dxdt(14:16) = matvecmult(Tg_w,x(9:11),3)  !  Contribution of rotor motion.
      dxdt(14) = dxdt(14) - Vref                !  Contribution of remote wind.

   end subroutine DeepwindTimestep


end module platformDynamics































program test

   use VAWTBEM
   use setupVAWTBEM
   use platformDynamics

   implicit none

   character(len=50)                             :: setupFile,airfoilFile,initFile, &
                                                    turbFile
   integer                                       :: setupFid,airfoilFid,initFid, &
                                                    turbFid

   character(len=1)                              :: junk
   integer                                       :: Nb,Nrows,Ncols,Ndof,Nadof,Nt, &
                                                    fstat,ierr
   double precision                              :: dt,wgset0,psisyn0,Tgset0, &
                                                    psi0hin,omegah0,          &
                                                    psi0gin,omegag0,Vref,z0,h0
   double precision, dimension(3)                :: y0,ydot0,xw0,xvane_g

   ierr = 0
   call initErrorLog()

   setupFile = 'in_Sandia17m.txt'
   setupFid = 49
   airfoilFile = 'AirfoilLibrary.dat'
   airfoilFid = 48
   initFile = 'init_Sandia17m.txt'
   initFid = 47
   turbFile = 'WindFiles\u08_I20.bin'  !  u15_I15.bin'
   turbFid = 46

   open (UNIT=41, FILE='out.txt', STATUS='REPLACE')
   open (UNIT=42, FILE='out2.txt', STATUS='REPLACE')
   open (UNIT=43, FILE='out3.txt', STATUS='REPLACE')
   open (UNIT=44, FILE='out4.txt', STATUS='REPLACE')

   open (UNIT=setupFid, FILE=setupFile, STATUS='OLD', IOSTAT=fstat, &
         ACCESS='SEQUENTIAL', FORM='FORMATTED',                     &
         ACTION='READ', POSITION='REWIND', DELIM='NONE')
   if (fstat .gt. 0) then
      ierr = -1
      call logError(48)
   end if
   read (UNIT=setupFid,FMT=*,IOSTAT=fstat) junk
   read (UNIT=setupFid,FMT=*,IOSTAT=fstat) Nb,Nrows,Ncols
   close(setupFID)

   Nadof = 7*Nrows*Ncols
   Ndof = 17 + Nadof

   open (UNIT=initFid, FILE=initFile, STATUS='OLD', IOSTAT=fstat, &
         ACCESS='SEQUENTIAL', FORM='FORMATTED',                   &
         ACTION='READ', POSITION='REWIND', DELIM='NONE')
   if (fstat .gt. 0) then
      ierr = -1
      call logError(48)
   end if
   read (UNIT=initFid,FMT=*,IOSTAT=fstat) Nt,dt
   read (UNIT=initFid,FMT=*,IOSTAT=fstat) Vref,z0,h0
   read (UNIT=initFid,FMT=*,IOSTAT=fstat) wgset0,psisyn0,Tgset0
   read (UNIT=initFid,FMT=*,IOSTAT=fstat) psi0hin,omegah0,psi0gin,omegag0
   read (UNIT=initFid,FMT=*,IOSTAT=fstat) y0(1),y0(2),y0(3)
   read (UNIT=initFid,FMT=*,IOSTAT=fstat) ydot0(1),ydot0(2),ydot0(3)
   read (UNIT=initFid,FMT=*,IOSTAT=fstat) xw0(1),xw0(2),xw0(3)
   read (UNIT=initFid,FMT=*,IOSTAT=fstat) xvane_g(1),xvane_g(2),xvane_g(3)
   close(initFid)

   call simulateDeepwind (41,setupFile,setupFid,           &
                          airfoilFile,airfoilFid,          &
                          turbFile,turbFid,                &
                          Nb,Nrows,Ncols,Ndof,Nadof,Nt,dt, &
                          wgset0,psisyn0,Tgset0,y0,ydot0,  &
                          psi0hin,omegah0,psi0gin,omegag0, &
                          xw0,Vref,z0,h0,xvane_g,ierr)

   close (41)
   close (42)
   close (43)
   close (44)

end program test


