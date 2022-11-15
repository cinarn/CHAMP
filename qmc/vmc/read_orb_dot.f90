      subroutine read_orb_dot
! Written by A.D.Guclu, Feb 2004.
! Reads in 2-dimensional basis fns info for circular quantum dots.

      use all_tools_mod
      use control_mod
      use contr2_mod
      use contrl_per_mod
      use dim_mod
      use numbas_mod
      use basis2_mod
      use coefs_mod
      use forcepar_mod
      use pseudo_mod
      implicit real*8(a-h,o-z)

! do some debugging. not sure if all these are necessary to check
! but take no chance for now
      if(ndim.ne.2) stop 'ndim must be 2 for quantum dots'
      if(nforce.ne.1) stop 'nforce must be 1 for quantum dots'
!     if(nloc.ne.-1) stop 'nloc must be -1 for quantum dots'
      if(numr.ne.0) stop 'numr must be 0 in read_orb_dot'
      if(inum_orb.ne.0) stop 'inum_orb must be 0 for quantum dots'

! For a dot or a ring there is just one center, so set for consistency with molecules the foll:
      call alloc ('nrbas_analytical', nrbas_analytical, 1)
      call alloc ('nrbas_numerical', nrbas_numerical, 1)
      call alloc ('nrbas', nrbas, 1)
      call alloc ('', ictype_basis, nbasis)
!     call alloc ('iwrwf2', iwrwf2, nbasis)
!     allocate iwrwf here in the case when it is not allocated and read in before (new style format input)
!     call alloc ('iwrwf', iwrwf, mbasis_ctype ,nctype)

      nrbas_analytical(1)=nbasis
      nrbas_numerical(1)=0
      nrbas(1)=nbasis
      do 10 ib=1,nbasis
   10   ictype_basis(ib)=1

      if(ibasis.eq.3) then
        call read_orb_dot_fd
      elseif(ibasis.ge.4 .or. ibasis.le.7) then
        call read_orb_dot_gauss
      else
        stop 'In read_orb_dot: only ibasis=3,4,5,6,7 allowed'
      endif
      return
      end
!-----------------------------------------------------------------------

      subroutine read_orb_dot_fd
! Written by A.D.Guclu, Feb 2004.
! Reads in quantum dot orbitals in Fock-Darwin basis set
! with quantum numbers n (quasi-Landau level) and m (angular mom.)
! The definition of "Landau level" is different in quantum Hall litt
! and quantum dot litt. For quantum dots we use the notation n_fd
! (for Fock-Darwin), for projected composite fermions we use n_cf
! just for convenience.
! For Fock-Darwin states, zex is more than an exponential parameter.
! It serves as a multiplicatif optimization factor for the spring constant.
      use all_tools_mod
      use coefs_mod
      use basis1_mod
      use const_mod
      use basis3_mod
      implicit real*8(a-h,o-z)

      common /compferm/ emagv,nv,idot

      call alloc ('n_fd', n_fd, nbasis)
      call alloc ('m_fd', m_fd, nbasis)
      call alloc ('n_cf', n_cf, nbasis)

! read the total number of "quasi-Landau" levels:
      read(5,*) nlandau
! next nlandau lines to read represent: n #m m1 m2 m3 ....
! for instance if we want landau levels with in the first LL; m=0,2,3, and
! in the second LL; m=1,2 then the input file should be:
!      03 0 2 3
!      12 1 2
      icount=0
      ncfmax=0
      do 10 in=1,nlandau
        read(5,*) n,nm,(m_fd(i),i=icount+1,icount+nm)
        do 5 i=icount+1,icount+nm
          n_fd(i)=n
          n_cf(i)=n
          if(m_fd(i).lt.0) n_cf(i)=n-m_fd(i)
          if(n_cf(i).gt.ncfmax) ncfmax=n_cf(i)
    5   enddo
        icount=icount+nm
   10 enddo
      if(icount.ne.nbasis) stop 'nbasis doesnt match the basis set'
!JT      if(ncfmax.gt.MBASIS) stop 'ncfmax.gt.MBASIS. this is a problem in cbasis_fns.f'
!      if(ncfmax.gt.6 .and. idot.eq.3)
      if(ncfmax.gt.6) write(6,'(''WARNING: landau levels 7 and 8 can cause numerical problems in projected cfs'')') ! idot not defined at this point

! read orbital coefficients
!      write(6,'(/,(12a10))') (n_fd(j),m_fd(j),j=1,nbasis)
      write(6,'(''orbital coefficients'')')
      do 20 iorb=1,norb
        read(5,*) (coef(j,iorb,1),j=1,nbasis)
        write(6,'(12f10.6)') (coef(j,iorb,1),j=1,nbasis)
   20 enddo

      write(6,'(''screening constants'')')
      read(5,*) (zex(i,1),i=1,nbasis)
      write(6,'(12f10.6)') (zex(i,1),i=1,nbasis)
      do 30 i=1,nbasis
! zex only used for idot=0, and should not be < 0. To be safe set it to 1.
! idot not defined at this point
        if(zex(i,1).le.0.d0) then
          write(6,'(''WARNING: exponent zex set to 1'')')
          zex(i,1)=1
        endif
   30 enddo


      return
      end

!-----------------------------------------------------------------------

      subroutine read_orb_dot_gauss
! Written by A.D.Guclu, Apr 2006.
! Edited by Gokhan Oztarhan, Feb 2022.
! Reads in quantum dot orbitals in gaussian basis set
! the witdh of gaussians is given by zex*we

      use control_mod
      use coefs_mod
      use const_mod
      use contrl_per_mod
      use optimo_mod
      use forcepar_mod
      use orbpar_mod
      use periodic_1d_mod
      use pseudo_mod !GO
      use files_tools_mod !GO
      use dets_mod !GO
      implicit real*8(a-h,o-z)
      
      integer orb_file_unit !GO
      

      write(6,'(/,''Reading floating gaussian orbitals for dots'')')
      call alloc ('oparm', oparm, notype, nbasis, nwf)
      do it=1,notype
        read(5,*)  (oparm(it,ib,1),ib=1,nbasis)
        if(ibasis.eq.4) then
          if(it.eq.1) write(6,'(''Floating gaussian x-positions:'')')
          if(it.eq.2) write(6,'(''Floating gaussian y-positions:'')')
          if(it.eq.3) write(6,'(''Floating gaussian widths:'')')
        elseif(ibasis.eq.5) then
          if(it.eq.1) write(6,'(''Floating gaussian radial positions:'')')
          if(it.eq.2) write(6,'(''Floating gaussian angular positions:'')')
          if(it.eq.3) write(6,'(''Floating gaussian radial widths:'')')
          if(it.eq.4) write(6,'(''Floating gaussian angular widths:'')')
        elseif(ibasis.eq.6) then
          if(it.eq.1) write(6,'(''Floating gaussian x-positions:'')')
          if(it.eq.2) write(6,'(''Floating gaussian y-positions:'')')
          if(it.eq.3) write(6,'(''Floating gaussian x-widths:'')')
          if(it.eq.4) write(6,'(''Floating gaussian y-widths:'')')
        elseif(ibasis.eq.7) then
!       Periodic Gaussians: do modulo math
          if(it.eq.1) then
            write(6,'(''Floating gaussian x-positions:'')')
            do ib=1,nbasis
              oparm(it,ib,1) = modulo(oparm(it,ib,1), alattice)
              if (oparm(it,ib,1).ge.(alattice/2.)) oparm(it,ib,1) = oparm(it,ib,1) - alattice
            enddo
          endif
          if(it.eq.2) write(6,'(''Floating gaussian y-positions:'')')
!     TO DO:  We should add in a check here to make sure that gaussians
!          aren't too much wider than alattice (i.e., cell size)
          if(it.eq.3) then
            write(6,'(''Floating gaussian x-widths:'')')
          endif
          if(it.eq.4) write(6,'(''Floating gaussian y-widths:'')')
        else
          write(6,'(''ibasis must be 4, 5, 6, or 7 in read_orb_dot_gauss'')')
          stop 'ibasis must be 4, 5, 6, or 7 in read_orb_dot_gauss'
        endif
        write(6,'(1000f12.6)') (oparm(it,ib,1),ib=1,nbasis)
      enddo

      do ib=1,nbasis
        if(oparm(3,ib,1).le.0.d0) then
          write(6,'(''WARNING: exponent oparm(3,ib,1) set to 1'')')
          oparm(3,ib,1)=1
        endif
        if(ibasis.eq.5 .or. ibasis.eq.6 .or. ibasis.eq.7) then
          if(oparm(4,ib,1).le.0.d0) then
            write(6,'(''oparm(4,ib,1) must be  > 0'')')
            stop 'oparm(4,ib,1) must be  > 0'
          endif
        endif
      enddo

! read orbital coefficients
!      write(6,'(/,(12a10))') (n_fd(j),m_fd(j),j=1,nbasis)
!      if(norb.le.100) then
     
      write(6,*)
      
      ! For backward compatibility with the older inputs, if 'orb_dot_coef' file exists, 
      ! program reads the orbital coefficients from the file. Otherwise, basis set=orbitals.
      if ( (nloc.eq.-6 .or. nloc.eq.-7) .and. file_exist('orb_dot_coef') ) then !GO
      
        write(6,'(''orbital coefficients - rows: norb, columns: nbasis'')')
      
        orb_file_unit = -1 ! -1 is here to get a unique file unit
        call open_file_or_die('orb_dot_coef', orb_file_unit)
        
        do iorb=1,norb 
            read(orb_file_unit,*) (coef(j,iorb,1),j=1,nbasis)
            write(6,'(10000f12.6)') (coef(j,iorb,1),j=1,nbasis)
        enddo
   
        close(orb_file_unit)
        
      else
      
        if(norb.ne.nbasis) stop &
     &  'norb must be equal to nbasis in read_orb_dot_gauss'
     
        write(6,'(''orbital coefficients'')')
        write(6,'(''Assuming basis set=orbitals for 2D-gaussian orbitals'')')
        do 40 iorb=1,norb
          do 30 j=1,nbasis
            if(iorb.eq.j) then
              coef(j,iorb,1)=1
            else
              coef(j,iorb,1)=0
            endif
   30     enddo
   40   enddo

      endif

      return
      end
