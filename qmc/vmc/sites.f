      subroutine sites(x,nelec,nsite)
c Written by Cyrus Umrigar
      use atom_mod
      implicit real*8(a-h,o-z)
!JT      include 'vmc.h'
!JT      include 'force.h'
!JT      include 'pseudo.h'

!JT      parameter(half=0.5d0)

c Routine to put electrons down around centers for a VERY crude initial
c configuration if nothing else is available.  It is better to put them
c too close than to put them too far away because they equilibrate faster
c when they are too close.

      common /dim/ ndim
      common /pseudo/ vps(MELEC,MCENT,MPS_L),vpso(MELEC,MCENT,MPS_L,MFORCE)
     &,npotd(MCTYPE),lpotp1(MCTYPE),nloc
!JT      common /atom/ znuc(MCTYPE),cent(3,MCENT),pecent
!JT     &,iwctype(MCENT),nctype,ncent
      common /dot/ w0,we,bext,emag,emaglz,emagsz,glande,p1,p2,p3,p4,rring
       common /wire/ wire_w,wire_length,wire_length2,wire_radius2, wire_potential_cutoff,wire_prefactor,wire_root1
      common /jel_sph2/ zconst  ! RM

      dimension x(3,*),nsite(*)

c Loop over spins and centers. If odd number of electrons on all
c atoms then the up-spins have an additional electron.
c So assumption is that system is not strongly polarized.

c     gauss()=dcos(two*pi*rannyu(0))*dsqrt(-two*dlog(rannyu(0)))
      pi=4*datan(1.d0)

      if((nloc.eq.-1).or.(nloc.eq.-5)) then ! parabolic quantum dot
        if(we.eq.0.d0) stop 'we should not be 0 in sites for quantum dots (nloc=-1)'
       elseif(nloc.eq.-3) then ! jellium RM
        if(zconst.eq.0.d0) stop 'zconst should not be 0 in sites for atoms in jellium (nloc=-3)'
      endif

      ielec=0
      do 10 ispin=1,2
        do 10 i=1,ncent
          if((nloc.eq.-1).or.(nloc.eq.-5)) then ! parabolic quantum dot
            znucc=dsqrt(we)
           elseif(nloc.eq.-3) then ! jellium RM
            znucc=zconst
           elseif(nloc.eq.-4) then ! quantum wire
            znucc=dsqrt(wire_w)
           else ! atoms and molecules
            if(znuc(iwctype(i)).eq.0.d0) 
     &stop 'znuc should not be 0 in sites for atoms and molecules'
            znucc=znuc(iwctype(i))
          endif
          ju=(nsite(i)+2-ispin)/2
          do 10 j=1,ju
            ielec=ielec+1
            if(ielec.gt.nelec) return
            if(nloc.eq.-1 .or. nloc.eq.-5) then
              sitsca=1/znucc
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


c sample position from exponentials or gaussian around center
c A.D.Guclu 5/2008: need circular coo. for ring shaped quantum dots            
            if((nloc.eq.-1 .or. nloc.eq.-5) .and. rring.gt.0.d0) then
              site=dsqrt(site)
              site=sign(site,(rannyu(0)-half))
              angle=2*pi*rannyu(0)
              x(1,ielec)=(sitsca*site+rring)*dcos(angle)
              x(2,ielec)=(sitsca*site+rring)*dsin(angle)
             else
               do 5 k=1,ndim
c sample position from exponentials or gaussian around center
c a.d.guclu: for wires distribute electrons linearly in y direction  
                 site=-dlog(rannyu(0))
                 if(nloc.eq.-1 .or. nloc.eq.-4 .or. nloc.eq.-5) site=dsqrt(site)
                 site=sign(site,(rannyu(0)-half))

                 if(nloc.eq.-4) then 
                   if(k.eq.2) then
                     x(k,ielec)=sitsca*site+cent(k,i)
                   else
                     x(k,ielec)=wire_length*(0.5d0-rannyu(0))
                   endif
                 else
                   x(k,ielec)=sitsca*site+cent(k,i)
                 endif
   5           enddo
            endif

   10     continue


!      write(6,*)
      write(6,'(a,i3,a)') '1 configuration for',ielec,' electrons has been generated by routine sites.'
      write(6,'(''sites:'',100d12.4)') ((x(k,i),k=1,ndim),i=1,nelec)

      if(ielec.lt.nelec) stop 'bad input to sites'
      return
      end

