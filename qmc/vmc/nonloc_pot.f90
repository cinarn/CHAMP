      subroutine nonloc_pot(x,rshift,rvec_en,r_en,vpsp,psid,pe)
! Written by Claudia Filippi, modified by Cyrus Umrigar
! Calculates the local and nonlocal components of the pseudopotential
! pe_en(loc) is computed in distances and pe_en(nonloc) here in nonloc_pot if nloc !=0 and iperiodic!=0.

      use control_mod
      use deriv_orb_mod
      use eloc_mod
      use periodic_jastrow_mod  !WAS
      use atom_mod
      use dets_mod
      use const_mod
      use pseudo_mod
      use contrl_per_mod
      use slater_mod, only: detu, detd
      implicit real*8(a-h,o-z)

      dimension x(3,*),rshift(3,nelec,ncent),rvec_en(3,nelec,ncent),r_en(nelec,ncent)
!    &,detu(*),detd(*),slmui(nupdn_square,*),slmdi(nupdn_square,*)

! Calculate local and nonlocal pseudopotentials for all electrons, nuclei and l-components
! and store in vps
      do 20 i=1,nelec
        if(nloc.eq.1) then
          call getvps_fahy(r_en,i)
         elseif(nloc.ge.2 .and. nloc.le.5) then
          call getvps_champ(r_en,i)
         elseif(nloc.eq.6) then
          call getvps_gauss(r_en,i)
         elseif(nloc.ge.7) then
          stop 'nloc >= 7'
        endif
   20 continue

      if(ipr.ge.4) write(6,'(''r_en='',9f9.5)') ((r_en(iel,ic),iel=1,nelec),ic=1,ncent)

! local component
      if(iperiodic.eq.0) then
        if(ipr.ge.4) write(6,'(''nonloc_pot: pe before local psp'',f9.5)') pe
        do 30 ic=1,ncent
          do 30 i=1,nelec
            if(ipr.ge.4) write(6,'(''i,ic,vps,pe'',2i3,9f9.5)')i,ic,vps(i,ic,lpotp1(iwctype(ic))),pe+vps(i,ic,lpotp1(iwctype(ic)))
   30       pe=pe+vps(i,ic,lpotp1(iwctype(ic)))
        if(ipr.ge.4) write(6,'(''nonloc_pot: pe after local psp'',f12.5)') pe
      endif

!     Total ee and local en potential
      eloc_pot_loc = pe                      !JT
      call object_modified_by_index (eloc_pot_loc_index)  !JT

! non-local component (division by the Jastrow already in nonloc)
      call nonloc(x,rshift,rvec_en,r_en,vpsp)
      pe=pe+vpsp/psid
      if(ipr.ge.1) write(6,'(''nonloc_pot: vpsp, psid, pe after nonlocal psp'',9f12.5)') vpsp, psid, pe
      if(ipr.ge.4) write(6,'(''nonloc_pot: pe,vpsp/psid,vpsp,psid,detu(1),detd(1),r_en(1,1)='',2f9.4,9d12.4)') &
     &pe,vpsp/psid,vpsp,psid,detu(1),detd(1),r_en(1,1)

      return
      end
