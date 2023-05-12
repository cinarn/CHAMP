      subroutine sites(x,nelec,nsite)
! Written by Cyrus Umrigar
      use constants_mod
      use atom_mod
      use dim_mod
      use pseudo_mod
      use jel_sph2_mod
      use contrl_per_mod
      use orbpar_mod
      use wfsec_mod
      use dorb_mod
      use periodic_mod, only: rlatt_sim
      use gs_mod

      implicit real*8(a-h,o-z)

! Routine to put electrons down around centers for a VERY crude initial
! configuration if nothing else is available.  It is better to put them
! too close than to put them too far away because they equilibrate faster
! when they are too close.

      common /dot/ w0,we,bext,emag,emaglz,emagsz,glande,p1,p2,p3,p4,rring
      common /cyldot/ cyldot_v, cyldot_s, cyldot_rho !GO
      common /gndot/ gndot_v0, gndot_rho, gndot_s, gndot_k !GO
      common /wire/ wire_w,wire_length,wire_length2,wire_radius2, wire_potential_cutoff,wire_prefactor,wire_root1

      dimension x(3,*),nsite(*)

! Loop over spins and centers. If odd number of electrons on all
! atoms then the up-spins have an additional electron.
! So assumption is that system is not strongly polarized.

!     gauss()=dcos(two*pi*rannyu(0))*dsqrt(-two*dlog(rannyu(0)))
      pi=4*datan(1.d0)

      if(iperiodic.eq.3) then ! periodic solids
! Generate points uniformly in lattice coordinates and convert back to cartesian coordinates
        do ielec=1,nelec
          do k=1,ndim
            x(k,ielec)=0
            do i=1,ndim
              r_basis=rannyu(0)-0.5
              x(k,ielec)=x(k,ielec)+rlatt_sim(k,i)*r_basis
            enddo
          enddo
        enddo
        goto 20
      endif

      if((nloc.eq.-1).or.(nloc.eq.-5)) then ! parabolic quantum dot
        if(we.eq.0.d0) stop 'we should not be 0 in sites for quantum dots (nloc=-1)'
       elseif(nloc.eq.-3) then ! jellium RM
        if(zconst.eq.0.d0) stop 'zconst should not be 0 in sites for atoms in jellium (nloc=-3)'
      endif

      write(6, '(''nsite = '',x,*(i6))') (nsite(i),i=1,ncent)

      ielec=0
      do 10 ispin=1,2
        do 10 i=1,ncent
          if((nloc.eq.-1).or.(nloc.eq.-5)) then ! parabolic quantum dot
            znucc=dsqrt(we)
           elseif(nloc.eq.-3) then ! jellium RM
            znucc=zconst
           elseif(nloc.eq.-4) then ! quantum wire
            znucc=dsqrt(wire_w)
           elseif(nloc.eq.-6) then ! cylindrical quantum dot !GO
            znucc = cyldot_rho
           elseif(nloc.eq.-7) then ! gaussian quantum dot !GO
            znucc = gndot_rho 
           elseif(nloc.eq.-8) then !NC
            znucc = dexp(sum(gs_ncent*dlog(gs_rho))/gs_npot/sum(gs_ncent)) ! geometric mean of all site potentials
           else ! atoms and molecules
            if(znuc(iwctype(i)).eq.0.d0) stop 'znuc should not be 0 in sites for atoms and molecules'
            znucc=znuc(iwctype(i))
          endif
          ju=(nsite(i)+2-ispin)/2
          do 10 j=1,ju
            ielec=ielec+1
            if(ielec.gt.nelec) return
            if(nloc.eq.-1 .or. nloc.eq.-5 .or. nloc.eq.-4) then
              sitsca=1/znucc
             elseif(nloc.eq.-6 .or. nloc.eq.-7) then
              sitsca=znucc/2 !GO
             elseif(nloc.eq.-8) then !NC
              sitsca=znucc
             elseif(j.eq.1) then
              sitsca=1/max(znucc,1.d0)
             elseif(j.le.5) then
              sitsca=2/max(znucc-2,1.d0)
             elseif(j.le.9) then
              sitsca=3/max(znucc-10,1.d0)
             elseif(j.le.18) then
              sitsca=4/max(znucc-18,1.d0)
             else
              sitsca=5/max(znucc-36,1.d0)
            endif


! sample position from exponentials or gaussian around center
! A.D.Guclu 5/2008: need circular coo. for ring shaped quantum dots
            if((nloc.eq.-1 .or. nloc.eq.-5) .and. rring.gt.0.d0) then
              if(ibasis.eq.5) then
                site = (0.5d0 - rannyu(0))/dsqrt(we*oparm(3, iworbd(ielec,1), iwf))
                angle = (0.5d0 - rannyu(0))/dsqrt(oparm(4, iworbd(ielec,1), iwf))
                site = site + oparm(1, iworbd(ielec,1), iwf)
                angle = angle + oparm(2, iworbd(ielec,1), iwf)
!  Make sure electron is near the center of some gaussian - might not work
!     if there's more than 1 slater determinant
                x(1,ielec)=site*dcos(angle)
                x(2,ielec)=site*dsin(angle)
              else
!               This code sampled from a gaussian:
!                site=-dlog(rannyu(0))
!                site=dsqrt(site)
!                site=sign(site,(rannyu(0)-half))
!               This code samples from a smaller, uniform region:
!                site = 2.0d0*(0.5d0 - rannyu(0))
                site = (0.5d0 - rannyu(0))/dsqrt(we)
                angle=2.0d0*pi*(dble(ielec) - rannyu(0))/dble(nelec)
                x(1,ielec)=(site+rring)*dcos(angle)
                x(2,ielec)=(site+rring)*dsin(angle)
!               x(1,ielec)=(sitsca*site+rring)*dcos(angle)
!               x(2,ielec)=(sitsca*site+rring)*dsin(angle)
              endif
! sample position gaussian around center   
! G.Oztarhan 06/2021: for cylindrical and gaussian quantum dots, make sure that electrons are
!                     located around centers of dots within an effective radius, 
!                     for gaussian basis set ensure that electrons are within the width of the floating gaussians,
!                     might not work if there is more than 1 slater determinant
            elseif(nloc.eq.-6 .or. nloc.eq.-7) then
              site = dsqrt(-dlog(rannyu(0)))
              if(ibasis.eq.4) then
                 site = site*min(sitsca,1.d0/dsqrt(oparm(3, iworbd(ielec,1), iwf)))
              else
                 site = site*sitsca
              endif
              angle = 2.0d0*pi*rannyu(0)
              x(1,ielec) = site*dcos(angle) + cent(1,iworbd(i,1))
              x(2,ielec) = site*dsin(angle) + cent(2,iworbd(i,1))
              
            elseif (nloc .eq. -8) then !NC
              do k=1, ndim
                x(k,ielec) = (0.5d0-rannyu(0))*sitsca
              end do
              
            else
               do 5 k=1,ndim
! sample position from exponentials or gaussian around center
! a.d.guclu: for wires distribute electrons linearly in y direction
! a.c.mehta: unless floating gaussians, then make sure electrons
!             are close to centers of gaussians
!  Warning:  this might not work if we have multiple slater determinants
                 site=-dlog(rannyu(0))
                 if(nloc.eq.-1 .or. nloc.eq.-4 .or. nloc.eq.-5) site=dsqrt(site)
                 site=sign(site,(rannyu(0)-half))

                 if(nloc.eq.-4) then
                   if (ibasis.eq.6 .or. ibasis.eq.7) then
                     site = (0.5d0 - rannyu(0))/dsqrt(we*oparm(k+2, iworbd(ielec,1), iwf))
!  Make sure electron is near the center of some gaussian - might not work
!     if there's more than 1 slater determinant
                     x(k,ielec) = site + oparm(k, iworbd(ielec,1), iwf)
                   else
                     if(k.eq.2) then
                       x(k,ielec)=sitsca*site+cent(k,i)
                     elseif(iperiodic.eq.0) then
                       x(k,ielec)=wire_length*(0.5d0-rannyu(0))
                     else
                       x(k,ielec)=wire_length*rannyu(0)
                     endif
                   endif
                 else                      ! molecules (mostly at least)
                   x(k,ielec)=sitsca*site+cent(k,i)
                 endif ! nloc.eq.-4
   5           enddo
            endif

   10     continue
      if(ielec.lt.nelec) stop 'bad input to sites'

!     write(6,*)
   20 write(6,'(a,i3,a)') '1 configuration for',nelec,' electrons has been generated by routine sites.'
      write(6,'(''sites:'',1000f12.6)') ((x(k,i),k=1,ndim),i=1,nelec)

      return
      end
