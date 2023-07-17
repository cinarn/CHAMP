      subroutine finwrt
! Written by Cyrus Umrigar
! routine to print out final results
      use all_tools_mod
      use constants_mod
      use control_mod
      use montecarlo_mod
      use main_menu_mod
      use mpi_mod
      use atom_mod
      use contrl_mod
      use const_mod
      use dim_mod
      use contr2_mod
      use gradhess_mod
      use contrl_opt2_mod
      use forcepar_mod
!     use doefp_mod
      use pseudo_mod
      use contrl_per_mod
      use contr3_mod
      use forcest_mod
      use denupdn_mod
      use stepv_mod
      use pairden_mod
      use fourier_mod
      use forcewt_mod
      use est2cm_mod
      use estsig_mod
      use estcum_mod
      use estsum_mod
      use zigzag_mod
      implicit real*8(a-h,o-z)
      character*80 fmt

      common /dot/ w0,we,bext,emag,emaglz,emagsz,glande,p1,p2,p3,p4,rring
      common /compferm/ emagv,nv,idot


!MS Jellium sphere
      common /jel_sph1/ dn_background,rs_jel,radius_b ! RM

      dimension zzfin(nzzvars), zzerr(nzzvars)

      err(x,x2,i)=dsqrt(abs(x2/wcum(i)-(x/wcum(i))**2)/iblk)
      err1(x,x2)=dsqrt(dabs(x2/passes-(x/passes)**2)/passes)
!     err1s(x,x2,i)=dsqrt(dabs(x2/wcum1(i)-(x/wcum1(i))**2)/passes)
      err1s(x,x2,i)=dsqrt(dabs(x2/wcum(i)-(x/wcum(i))**2)/passes)

      if(index(mode,'mpi').ne.0) then
        call finwrt_mpi
        call grad_hess_jas_mpi
      endif

!     write(6,'(''ecum1,ecum(1)'',9d16.8)') ecum1,ecum(1)
!     write(6,'(''wcum1,wcum(1)'',9d16.8)') wcum1,wcum(1)

!     if(idtask.ne.0) return

      if(index(mode,'mpi').eq.0) then
        passes=dfloat(iblk)*dfloat(nstep)
       else
        passes=dfloat(iblk)*dfloat(nstep)*dfloat(nproc)
      endif
      rtpass=dsqrt(passes)

!     if(index(mode,'mov1').eq.0) then
!       accfin=acccum/passes
!      else
!       accfin=acccum/(passes*nelec)
!     endif

      efin=ecum(1)/passes
      energy(1)=efin
      pefin=pecum/passes
      peifin=peicum/passes
      tpbfin=tpbcum/passes
      tjffin=tjfcum/passes
      r1fin=r1cum/passes
      r2fin=r2cum/passes
      r3fin=r3cum/passes
      r4fin=r4cum/passes
      if(izigzag.gt.0) then
        zzfin(:)=zzcum(:)/passes
      endif
      accfin=acccum/passes

! In all-electron move algorithm, eerr1 differs from sigma in that eerr1 contains
! p*new+q*old, so eerr1 is a bit smaller than sigma.  sigma is a property
! of the wavefn. only, whereas eerr1 depends on how quickly
! one evolves the system.  In the calculation of T_corr, if one
! uses Tcorr=(eerr/eerr1)^2 then Tcorr=1 when nstep=1, whereas if one
! uses Tcorr=(eerr/sigma)^2 then Tcorr will be a bit < 1 when nstep=1.
! However, it makes sense to use the latter definition because
! p*new+q*old does reduce Tcorr and that is precisely what is being
! reflected when we get Tcorr < 1.
      eerr1=err1(ecum1,ecm21)
      eer1s=err1(ecum1s(1),ecm21s(1))
      eerr=err(ecum(1),ecm2(1),1)
      peerr=err(pecum,pecm2,1)
      peierr=err(peicum,peicm2,1)
      tpberr=err(tpbcum,tpbcm2,1)
      tjferr=err(tjfcum,tjfcm2,1)
      r1err=err(r1cum,r1cm2,1)
      r2err=err(r2cum,r2cm2,1)
      r3err=err(r3cum,r3cm2,1)
      r4err=err(r4cum,r4cm2,1)
      if(izigzag.gt.0) then
       do iz=1,nzzvars
        zzerr(iz)=err(zzcum(iz),zzcm2(iz),1)
       enddo
      endif
!     tcsq=eerr/eerr1
      tcsq=eerr/eer1s
      call alloc ('eloc_tc', eloc_tc, nforce)
      eloc_tc (1) = tcsq**2 !JT
!      sigma=eer1s*rtpass

      call object_modified ('eerr')  !JT
!      call object_modified ('sigma') !JT

! separate "magnetic energy" for quantum dots:
      if(ndim.eq.2) then
        temp=0.25d0*bext*bext/(we*we)
        tmfin=(pefin-peifin-emag)*temp
        tmerr=(peerr+peierr)*temp
        pefin=pefin-tmfin-emag
!       peerr=peerr+tmerr                          ! is this correct?
! note that temp is always smaller than 1
        peerr=peerr*(1-temp)+peierr*temp
      endif

! save energy, energy_sigma and energy_err for optimization
      energy(1)=efin
      energy_sigma(1)=sigma
      energy_err(1)=eerr

!     write(6,*) 'before grad_hess_jas_fin'
      if(igradhess.ge.1) call grad_hess_jas_fin(passes,efin)
!     write(6,*) 'after grad_hess_jas_fin'

      trysum=0.d0 !JT
      sucsum=0.d0 !JT
      do 90 i=1,NRAD
        trysum=trysum+try(i)
   90   sucsum=sucsum+suc(i)

      if(print_radial_probability .and. iperiodic.eq.0 .and. ncent.eq.1 .and. ipr.ge.-4) then
        if(ndim.eq.3) write(6,'(/,'' r  4*pi*r^2*rho 4*pi*r^2*rhoup 4*pi*r^2*rhodn'')')
        if(ndim.eq.2) write(6,'(/,'' r  2*pi*r*rho 2*pi*r*rhoup 2*pi*r*rhodn'')')
        delr=one/delri
        term=one/(passes*delr)
        do 100 i=1,NRAD
! 100     write(6,'(f5.3,3f10.6)') delr*(i-half),rprob(i)*term,rprobup(i)*term,rprobdn(i)*term
  100     write(6,'(f8.4,3f10.6)') delr*(i-half),rprob(i)*term,rprobup(i)*term,rprobdn(i)*term
      endif

!      if(nefp.gt.0) call writebas(passes,ecum1)

      write(fmt,'(''(/,a12,2x,a'',i3,'')'')') len_trim(title)
      write(6,fmt) mode,title
      if(index(mode,'mpi').eq.0) then
        write(6,'(a,f12.0,a,i6,a,i6,a)') 'Final results after ',passes,' passes (nstep = ',nstep,', nblk = ',iblk,')'
       else
        write(6,'(a,f12.0,a,i6,a,i6,a,i6,a)') 'Final results after ',passes,' passes (nstep = ',nstep,', nproc = ',nproc,', nblk = ' &
     &  ,iblk,')'
      endif
     ! write(6,'(''physical variable'',t20,''average'',t34,''rms error'' &
     !&,t47,''rms er*rt(pass)'',t65,''sigma'',t86,''Tcor'')')  !JT

      write(6,'(a17,1x,a14,1x,"  ",1x,a11,1x,a9,1x,a11,1x,a9,1x,"  ",1x,a9,1x,a8)') & !GO
     & 'physical variable','average','rms error','rms','er*rt(pass)','sigma','','Tcor'

      if(nforce.eq.1) then
        write(6,'("total E =        ",1x,f14.7,1x,"+-",1x,f11.7,1x,f9.5,1x,f11.5,1x,f9.5,1x,"+-",1x,f9.5,1x,f8.2)') & !GO
     &  efin,eerr,eerr*rtpass,eerr1*rtpass,sigma,error_sigma,tcsq*tcsq

!MS Jellium sphere
! eb_slf = Self-energy of the background charge
! ebz    = Interaction between background and central charges
        if(nloc.eq.-3) then
          eb_slf=3.d0/5.d0*(dn_background)**(5.d0/3.d0)/rs_jel
          ebz   =3.d0/2.d0*(dn_background)**(2.d0/3.d0)*znuc(1)/rs_jel
          write(6,'(''background self E   ='',f15.7)') eb_slf
          write(6,'(''z-background    E   ='',f15.7)') ebz
          write(6,'(''total E+eb_self+ebz ='',f15.7)') efin+eb_slf+ebz
        endif

       else
        ifr=1
!       write(6,'(''total E'',i3,'' ='',t17,f12.7,'' +-'',f11.7,3f9.5,t82,f8.2)') !JT
!    &  ifr,efin,eerr,eerr*rtpass,eerr1*rtpass,sigma,tcsq*tcsq
        write(6,'(''total E'',i3,'' ='',f16.7,1x,'' +-'',1x,f11.7,1x,f9.5,1x,f9.5,1x,f9.5,1x,'' +-'',1x,f9.5,1x,f8.2)') &
     &  ifr,efin,eerr,eerr*rtpass,eerr1*rtpass,sigma,error_sigma,tcsq*tcsq
      endif

      wcum(1)=passes
      call alloc ('eloc_tc', eloc_tc, nforce)
      do 110 ifr=2,nforce
        efin=ecum(ifr)/wcum(ifr)
        eerr=err(ecum(ifr),ecm2(ifr),ifr)
        eer1s=err1s(ecum1s(ifr),ecm21s(ifr),ifr)
        tcsq=eerr/eer1s
        eloc_tc (ifr) = tcsq**2 !JT
        sigma=eer1s*rtpass
        ffin=efin-ecum(1)/passes
        if(deltot(ifr).ne.0.d0) then
!         ffin=ffin/deltot(ifr)
!         ferr=err(fcum(ifr),fcm2(ifr),1)/abs(deltot(ifr))
          ferr=err(fcum(ifr),fcm2(ifr),1)
         else
          ferr=err(fcum(ifr),fcm2(ifr),1)
        endif
! save energy, force, energy_sigma, energy_err and force_err for optimization
! force and force_err are really the energy difference and the error in the energy difference.
        energy(ifr)=efin
        energy_err(ifr)=eerr
        energy_sigma(ifr)=sigma
        force(ifr)=ffin
        force_err(ifr)=ferr
        write(6,'(''total E'',i3,'' ='',f16.7,'' +-'',f11.7,f9.5,f18.5,t82,f8.2)') & !JT
       &ifr,efin,eerr,eerr*rtpass,sigma,tcsq*tcsq
  110   write(6,'(''E_diff'',i3,'' ='',t13,f16.7,'' +-'',f11.7,f9.5)') &
     &  ifr,ffin,ferr,ferr*rtpass
      write(6,'("potential E =    ",1x,f14.7,1x,"+-",1x,f11.7,1x,f9.5)') pefin,peerr,peerr*rtpass !GO
      write(6,'("interaction E =  ",1x,f14.7,1x,"+-",1x,f11.7,1x,f9.5)') peifin,peierr,peierr*rtpass !GO
      write(6,'("jf kinetic E =   ",1x,f14.7,1x,"+-",1x,f11.7,1x,f9.5)') tjffin,tjferr,tjferr*rtpass !GO
      write(6,'("pb kinetic E =   ",1x,f14.7,1x,"+-",1x,f11.7,1x,f9.5)') tpbfin,tpberr,tpberr*rtpass !GO
!JT         write(6,'(a,f12.7,a,f11.7)') 'Potential energy    =',pefin,' +-',peerr
!JT         write(6,'(a,f12.7,a,f11.7)') 'Kinetic energy (JF) =',tjffin,' +-',tjferr
!JT         write(6,'(a,f12.7,a,f11.7)') 'Kinetic energy (PB) =',tpbfin,' +-',tpberr

      call object_modified ('eloc_tc')

      if(ndim.eq.2) then
        write(6,'("radial mag. E =  ",1x,f14.7,1x,"+-",1x,f11.7,1x,f9.5)') tmfin,tmerr,tmerr*rtpass !GO
        write(6,'("orbital mag. E = ",1x,f14.7)') emaglz+emagv !GO
        write(6,'("Zeeman E =       ",1x,f14.7)') emagsz !GO
      endif

      if(iperiodic.eq.0.and.ncent.eq.1) then
        write(6,'(''<r> ='',t8,f21.7,'' +-'',f11.7,f9.5)') r1fin,r1err,r1err*rtpass
        write(6,'(''<r2> ='',t8,f21.7,'' +-'',f11.7,f9.5)') r2fin,r2err,r2err*rtpass
        write(6,'(''<r3> ='',t8,f21.4,'' +-'',f11.5,f9.2)') r3fin,r3err,r3err*rtpass
        write(6,'(''<r4> ='',t8,f21.2,'' +-'',f11.3,f9.1)') r4fin,r4err,r4err*rtpass
      endif

      if(izigzag.ge.1) then
        call print_zigzag_vars(zzfin,zzerr,rtpass)
!       write(6,'(''<ZigZag Amp> ='',t17,f12.7,'' +-'',f11.7,f9.5)') zzfin(3),zzerr(3),zzerr(3)*rtpass
!       write(6,'(''<|ZigZag Amp|> ='',t17,f12.7,'' +-'',f11.7,f9.5)') zzfin(1),zzerr(1),zzerr(1)*rtpass
!       write(6,'(''<ZigZag Amp^2> ='',t17,f12.7,'' +-'',f11.7,f9.5)') zzfin(2),zzerr(2),zzerr(2)*rtpass
!       write(6,'(''<ZigZag Amp (red)>='',t22,f12.7,'' +-'',f11.7,f9.5)') zzfin(6),zzerr(6),zzerr(6)*rtpass
!       write(6,'(''<|ZigZag Amp| (red)>='',t22,f12.7,'' +-'',f11.7,f9.5)') zzfin(4),zzerr(4),zzerr(4)*rtpass
!       write(6,'(''<ZigZag Amp^2 (red)>='',t22,f12.7,'' +-'',f11.7,f9.5)') zzfin(5),zzerr(5),zzerr(5)*rtpass
!       write(6,'(''<ZigZag rand Amp>='',t22,f12.7,'' +-'',f11.7,f9.5)') zzfin(9),zzerr(9),zzerr(9)*rtpass
!       write(6,'(''<|ZigZag rand Amp|>='',t22,f12.7,'' +-'',f11.7,f9.5)') zzfin(7),zzerr(7),zzerr(7)*rtpass
!       write(6,'(''<ZigZag rand Amp^2>='',t22,f12.7,'' +-'',f11.7,f9.5)') zzfin(8),zzerr(8),zzerr(8)*rtpass
      endif
      if(print_radial_probability .and. index(mode,'mov1').ne.0.and.iperiodic.eq.0.and.ncent.eq.1) then
        write(6,'(''acceptance          ='',t17,2f12.7)') accfin,sucsum/trysum
       else
        !write(6,'(''acceptance          ='',t17,2f12.7)') accfin
        write(6,'("acceptance       ",1x,f14.7)') accfin !GO
      endif

      if(ifixe.ne.0 .or. ifourier.ne.0 .or. izigzag.ne.0) call den2dwrt(passes,r1fin)

      call routines_write_final
      call reinit_routines_write_block
      call reinit_routines_write_final

      return
      end
