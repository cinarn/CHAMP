      subroutine acuest_dmc_movall
! Written by Cyrus Umrigar and Claudia Filippi
      use constants_mod
      use control_mod
      use atom_mod
      use dets_mod
      use contrl_mod
      use const_mod
      use dim_mod
      use forcepar_mod
      use pseudo_mod
      use contrl_per_mod
      use delocc_mod
      use force_dmc_mod
      use iterat_mod
      use jacobsave_mod
      use forcest_dmc_mod
      use denupdn_mod
      use stepv_mod
      use config_dmc_mod
      use branch_mod
      use estsum_dmc_mod
      use estcum_dmc_mod
      use div_v_dmc_mod
      use contrldmc_mod
      use estcm2_mod
      use stats_mod
      use age_mod
      use pairden_mod
      use fourier_mod
      use pop_control_mod, only : ffn
      use zigzag_mod
      implicit real*8(a-h,o-z)

! routine to accumulate estimators for energy etc.

      common /dot/ w0,we,bext,emag,emaglz,emagsz,glande,p1,p2,p3,p4,rring
!     common /compferm/ emagv,nv,idot

      dimension zznow(nzzvars)

! statement function for error calculation
      rn_eff(w,w2)=w**2/w2
      error(x,x2,w,w2)=dsqrt(max((x2/w-(x/w)**2)/(rn_eff(w,w2)-1),0.d0))
      errg(x,x2,i)=error(x,x2,wgcum(i),wgcm2(i))

! wt   = weight of configurations
! xsum = sum of values of x from dmc
! xnow = average of values of x from dmc
! xcum = accumulated sums of xnow
! xcm2 = accumulated sums of xnow**2
! xave = current average value of x
! xerr = current error of x

      iblk=iblk+1
      npass=iblk*nstep

!     wnow=wsum/nstep
!     wfnow=wfsum/nstep
      enow=esum/wsum
      efnow=efsum/wfsum
      ei1now=wfsum/wdsum
      ei2now=wgsum(1)/wgdsum
      rinow=risum/wgsum(1)
      r1now=r1sum/wgsum(1)
      r2now=r2sum/wgsum(1)
      r3now=r3sum/wgsum(1)
      r4now=r4sum/wgsum(1)
      if(izigzag.gt.0) then
       zznow(:)=zzsum(:)/wgsum(1)
      endif

      wcm2=wcm2+wsum**2
      wfcm2=wfcm2+wfsum**2
      ecm2=ecm2+esum*enow
      efcm2=efcm2+efsum*efnow
      ei1cm2=ei1cm2+ei1now**2
      ei2cm2=ei2cm2+ei2now**2
      r1cm2=r1cm2+r1sum*r1now
      r2cm2=r2cm2+r2sum*r2now
      r3cm2=r3cm2+r3sum*r3now
      r4cm2=r4cm2+r4sum*r4now
      ricm2=ricm2+risum*rinow
      if(izigzag.gt.0) then
       zzcm2(:)=zzcm2(:)+zzsum(:)*zznow(:)
      endif

      wcum=wcum+wsum
      wfcum=wfcum+wfsum
      wdcum=wdcum+wdsum
      wgdcum=wgdcum+wgdsum
      ecum=ecum+esum
      efcum=efcum+efsum
      ei1cum=ei1cum+ei1now
      ei2cum=ei2cum+ei2now
      r1cum=r1cum+r1sum
      r2cum=r2cum+r2sum
      r3cum=r3cum+r3sum
      r4cum=r4cum+r4sum
      ricum=ricum+risum
      if(izigzag.gt.0) then
       zzcum(:)=zzcum(:)+zzsum(:)
      endif

      do 15 ifr=1,nforce

!       wgnow=wgsum(ifr)/nstep
        egnow=egsum(ifr)/wgsum(ifr)
        penow=pesum(ifr)/wgsum(ifr)
        peinow=peisum(ifr)/wgsum(ifr)
        tpbnow=tpbsum(ifr)/wgsum(ifr)
        tjfnow=tjfsum(ifr)/wgsum(ifr)

        wgcm2(ifr)=wgcm2(ifr)+wgsum(ifr)**2
        egcm2(ifr)=egcm2(ifr)+egsum(ifr)*egnow
        pecm2(ifr)=pecm2(ifr)+pesum(ifr)*penow
        peicm2(ifr)=peicm2(ifr)+peisum(ifr)*peinow
        tpbcm2(ifr)=tpbcm2(ifr)+tpbsum(ifr)*tpbnow
        tjfcm2(ifr)=tjfcm2(ifr)+tjfsum(ifr)*tjfnow

        wgcum(ifr)=wgcum(ifr)+wgsum(ifr)
        egcum(ifr)=egcum(ifr)+egsum(ifr)
        pecum(ifr)=pecum(ifr)+pesum(ifr)
        peicum(ifr)=peicum(ifr)+peisum(ifr)
        tpbcum(ifr)=tpbcum(ifr)+tpbsum(ifr)
        tjfcum(ifr)=tjfcum(ifr)+tjfsum(ifr)

        if(iblk.eq.1) then
          egerr=0
          peerr=0
          peierr=0
          tpberr=0
          tjferr=0
         else
          egerr=errg(egcum(ifr),egcm2(ifr),ifr)
          peerr=errg(pecum(ifr),pecm2(ifr),ifr)
          peierr=errg(peicum(ifr),peicm2(ifr),ifr)
          tpberr=errg(tpbcum(ifr),tpbcm2(ifr),ifr)
          tjferr=errg(tjfcum(ifr),tjfcm2(ifr),ifr)
        endif

        egave=egcum(ifr)/wgcum(ifr)
        peave=pecum(ifr)/wgcum(ifr)
        peiave=peicum(ifr)/wgcum(ifr)
        tpbave=tpbcum(ifr)/wgcum(ifr)
        tjfave=tjfcum(ifr)/wgcum(ifr)

        if(ifr.gt.1) then
          fgcum(ifr)=fgcum(ifr)+wgsum(1)*(egnow-egsum(1)/wgsum(1))
          fgcm2(ifr)=fgcm2(ifr)+wgsum(1)*(egnow-egsum(1)/wgsum(1))**2
          fgave=egcum(1)/wgcum(1)-egcum(ifr)/wgcum(ifr)
          if(iblk.eq.1) then
            fgerr=0
           else
            fgerr=errg(fgcum(ifr),fgcm2(ifr),1)
          endif
          if(deltot(ifr).ne.0.d0) then
            fgave=fgave/abs(deltot(ifr))
            fgerr=fgerr/abs(deltot(ifr))
          endif
        endif

!       if(iblk.eq.1) write(6,*) ecum,ecm2,dsqrt(ecm2),wcum,ecm2/wcum,(ecum/wcum)**2

! write out header first time

        if(iblk.eq.1.and.ifr.eq.1) then
          if(ibasis.eq.3) then
            write(6,'(t5,''egnow'',t15,''egave'',t21,''(egerr)'' ,t32,''peave'',t38,''(peerr)'',t49,''tpbave'',t55,''(tpberr)'',t66 &
     &  ,''tjfave'',t72,''(tjferr)'',t83,''emave'',t89,''(emave)'',t100,''fgave'',t106,''(fgerr)'', &
     &  t118,''npass'',t128,''wgsum'',t138,''ioldest'')')
           else
            write(6,'(t5,''egnow'',t15,''egave'',t21,''(egerr)'' ,t32,''peave'',t38,''(peerr)'',t49,''tpbave'',t55,''(tpberr)'',t66 &
     &  ,''tjfave'',t72,''(tjferr)'',t83,''fgave'',t89,''(fgerr)'',t101,''npass'',t111,''wgsum'',t121,''ioldest'')')
          endif
        endif

! write out current values of averages etc.

        if(ndim.eq.2) then
          iegerr=nint(10000000*egerr)
          ipeerr=nint(10000000*peerr)
          ipeierr=nint(10000000*peierr)
          itpber=nint(10000000*tpberr)
          itjfer=nint(10000000*tjferr)
          ifgerr=nint(10000000*fgerr)
         else
          iegerr=nint(100000*egerr)
          ipeerr=nint(100000*peerr)
          itpber=nint(100000*tpberr)
          itjfer=nint(100000*tjferr)
          ifgerr=nint(100000*fgerr)
        endif

! magnetic energy for quantum dots...
! right definition of the potential energy does not include magnetic energy.
        if(ndim.eq.2) then
!         emave=0.125*bext*bext*r2cum/wgcum(ifr)+emaglz+emagsz+emagv
          temp=0.25d0*bext*bext/(we*we)
          emave=(peave-peiave-emag)*temp+emag
          emerr=(peerr+peierr)*temp
          iemerr=nint(10000000*emerr)
          peave=peave-emave
!         ipeerr=ipeerr+iemerr
          ipeerr=nint(10000000*(peerr*(1-temp)+temp*peierr))
        endif

        if(ifr.eq.1) then
          if(ndim.eq.2) then
            write(6,'(f14.7,5(f14.7,''('',i7,'')''),17x,3i10)') egnow,egave,iegerr,peave,ipeerr,tpbave,itpber,tjfave,itjfer, &
     &      emave,iemerr,npass,nint(wgsum(ifr)),ioldest
           else
            write(6,'(f10.5,4(f10.5,''('',i5,'')''),17x,3i10)') egnow,egave,iegerr,peave,ipeerr,tpbave,itpber,tjfave,itjfer, &
     &      npass,nint(wgsum(ifr)),ioldest
          endif
         else
          if(ndim.eq.2) then
            write(6,'(f14.7,5(f14.7,''('',i7,'')''),17x,3i10)') egnow,egave,iegerr,peave,ipeerr,tpbave,itpber,tjfave,itjfer, &
     &      emave,iemerr,nint(wgsum(ifr))
           else
            write(6,'(f10.5,5(f10.5,''('',i5,'')''),10x,3i10)') egnow,egave,iegerr,peave,ipeerr,tpbave,itpber,tjfave,itjfer, &
     &      fgave,ifgerr,nint(wgsum(ifr))
          endif
        endif
   15 continue

! zero out xsum variables

      wsum=zero
      wfsum=zero
      wdsum=zero
      wgdsum=zero
      esum=zero
      efsum=zero
      ei1sum=zero
      ei2sum=zero
      r1sum=zero
      r2sum=zero
      r3sum=zero
      r4sum=zero
      risum=zero
      if(izigzag.gt.0) then
       zzsum(:)=zero
      endif

      do 20 ifr=1,nforce
        egsum(ifr)=zero
        wgsum(ifr)=zero
        pesum(ifr)=zero
        peisum(ifr)=zero
        tpbsum(ifr)=zero
   20   tjfsum(ifr)=zero

      call systemflush(6)

      return

      entry acues1_dmc_movall
! statistical fluctuations without blocking

      if(ipr.gt.-2) then
         write(11,'(i8,f11.8,f15.8,f13.8,i5)') ipass,ffn,wsum1(1),esum1(1)/wsum1(1),nwalk
      end if

      wcum1=wcum1+wsum1(1)
      wfcum1=wfcum1+wfsum1
      ecum1=ecum1+esum1(1)
      efcum1=efcum1+efsum1
      ei3cum=ei3cum+wfsum1/wdsum1

      wcm21=wcm21+wsum1(1)**2
      wfcm21=wfcm21+wfsum1**2
      ecm21=ecm21+esum1(1)**2/wsum1(1)
      efcm21=efcm21+efsum1**2/wfsum1
      ei3cm2=ei3cm2+(wfsum1/wdsum1)**2

      wfsum1=zero
      wdsum1=zero
      efsum1=zero

      do 22 ifr=1,nforce
        wgcum1(ifr)=wgcum1(ifr)+wgsum1(ifr)
        egcum1(ifr)=egcum1(ifr)+egsum1(ifr)
        wgcm21(ifr)=wgcm21(ifr)+wgsum1(ifr)**2
        if(wgsum1(ifr).ne.0.d0) then
          egcm21(ifr)=egcm21(ifr)+egsum1(ifr)**2/wgsum1(ifr)
         else
          egcm21(ifr)=0
        endif
        wsum1(ifr)=zero
        wgsum1(ifr)=zero
        esum1(ifr)=zero
   22   egsum1(ifr)=zero

      return

      entry zeres0_dmc_movall
! Initialize various quantities at beginning of run
! the initial values of energy psi etc. are calculated here
      ipass=0

! set quadrature points

!     if(nloc.gt.0) call gesqua(nquad,xq,yq,zq,wq)
      if(nloc.gt.0) call rotqua

      eigv=one
      eest=etrial
      nwalk=nconf
      wdsumo=nconf_global
      wgdsumo=nconf_global
      fprod=one
      do 70 i=0,nfprod
        wtgen(i)=nconf_global
   70   ff(i)=one

      do 80 iw=1,nconf
        current_walker=iw !TA
        call object_modified_by_index (current_walker_index) !TA
        wt(iw)=one
!       if(istrech.eq.0) then
!         do 71 ifr=2,nforce
!           do 71 ie=1,nelec
!             do 71 k=1,ndim
!  71           xoldw(k,ie,iw,ifr)=xoldw(k,ie,iw,1)
!       endif
        do 72 ifr=1,nforce
          if(nforce.gt.1) then
            call strech(xoldw(1,1,iw,1),xoldw(1,1,iw,ifr),ajacob,ifr,1)
           else
            ajacob=one
          endif
          ajacold(iw,ifr)=ajacob
          call hpsi(xoldw(1,1,iw,ifr),psidow(iw,ifr),psijow(iw,ifr),voldw(1,1,iw,ifr),div_vow(1,iw),d2ow(iw,ifr),peow(iw,ifr), &
     &    peiow(iw,ifr),eoldw(iw,ifr),denergy,ifr)
          pwt(iw,ifr)=one
          do 72 ip=0,nwprod-1
            wthist(iw,ip,ifr)=one
   72   continue
        if(psidow(iw,1).lt.zero) then
          do 76 ifr=1,nforce
            psidow(iw,ifr)=-psidow(iw,ifr)
            if(nup.gt.1) then
               do 74 k=1,ndim
               temp=voldw(k,1,iw,ifr)
               voldw(k,1,iw,ifr)=voldw(k,2,iw,ifr)
               voldw(k,2,iw,ifr)=temp
               temp=div_vow(1,iw)
               div_vow(1,iw)=div_vow(2,iw)
               div_vow(2,iw)=temp
               temp=xoldw(k,1,iw,ifr)
               xoldw(k,1,iw,ifr)=xoldw(k,2,iw,ifr)
   74          xoldw(k,2,iw,ifr)=temp
             else if(ndn.gt.1) then
               do 75 k=1,ndim
               temp=voldw(k,nup+1,iw,ifr)
               voldw(k,nup+1,iw,ifr)=voldw(k,nup+2,iw,ifr)
               voldw(k,nup+2,iw,ifr)=temp
               temp=div_vow(nup+1,iw)
               div_vow(nup+1,iw)=div_vow(nup+2,iw)
               div_vow(nup+2,iw)=temp
               temp=xoldw(k,nup+1,iw,ifr)
               xoldw(k,nup+1,iw,ifr)=xoldw(k,nup+2,iw,ifr)
   75          xoldw(k,nup+2,iw,ifr)=temp
             else
               write(6,'(5x,''negative psi for boson wave function'')')
!              stop
            endif
   76     continue
        endif
   80 continue

      entry zerest_dmc_movall
! entry point to zero out all averages etc. after equilibration runs

      iblk=0

! zero out estimators

      wcum1=zero
      wfcum1=zero
      wcum=zero
      wfcum=zero
      wdcum=zero
      wgdcum=zero
      ecum1=zero
      efcum1=zero
      ecum=zero
      efcum=zero
      ei1cum=zero
      ei2cum=zero
      ei3cum=zero
      r1cum=zero
      r2cum=zero
      r3cum=zero
      r4cum=zero
      ricum=zero
      if(izigzag.gt.0) then
       zzcum(:)=zero
      endif

      wcm21=zero
      wfcm21=zero
      wcm2=zero
      wfcm2=zero
      wdcm2=zero
      wgdcm2=zero
      ecm21=zero
      efcm21=zero
      ecm2=zero
      efcm2=zero
      ei1cm2=zero
      ei2cm2=zero
      ei3cm2=zero
      r1cm2=zero
      r2cm2=zero
      r3cm2=zero
      r4cm2=zero
      ricm2=zero
      if(izigzag.gt.0) then
       zzcm2(:)=zero
      endif

      wfsum1=zero
      wsum=zero
      wfsum=zero
      wdsum=zero
      wgdsum=zero
      efsum1=zero
      esum=zero
      efsum=zero
      ei1sum=zero
      ei2sum=zero
      ei3sum=zero
      r1sum=zero
      r2sum=zero
      r3sum=zero
      r4sum=zero
      risum=zero
      if(izigzag.gt.0) then
       zzsum(:)=zero
      endif

      call alloc ('fgcum', fgcum, nforce)
      call alloc ('fgcm2', fgcm2, nforce)
      call alloc ('wgcm2', wgcm2, nforce)
      call alloc ('wgcm21', wgcm21, nforce)
      call alloc ('egcm2', egcm2, nforce)
      call alloc ('egcm21', egcm21, nforce)
      call alloc ('pecm2', pecm2, nforce)
      call alloc ('tpbcm2', tpbcm2, nforce)
      call alloc ('tjfcm2', tjfcm2, nforce)
      call alloc ('peicm2', peicm2, nforce)
      call alloc ('wgcum', wgcum, nforce)
      call alloc ('wgcum1', wgcum1, nforce)
      call alloc ('egcum', egcum, nforce)
      call alloc ('egcum1', egcum1, nforce)
      call alloc ('pecum', pecum, nforce)
      call alloc ('peicum', peicum, nforce)
      call alloc ('tpbcum', tpbcum, nforce)
      call alloc ('tjfcum', tjfcum, nforce)
      call alloc ('wgsum', wgsum, nforce)
      call alloc ('wsum1', wsum1, nforce)
      call alloc ('wgsum1', wgsum1, nforce)
      call alloc ('esum1', esum1, nforce)
      call alloc ('egsum', egsum, nforce)
      call alloc ('egsum1', egsum1, nforce)
      call alloc ('pesum', pesum, nforce)
      call alloc ('peisum', peisum, nforce)
      call alloc ('tpbsum', tpbsum, nforce)
      call alloc ('tjfsum', tjfsum, nforce)
      do 85 ifr=1,nforce
        wgcum1(ifr)=zero
        wgcum(ifr)=zero
        egcum1(ifr)=zero
        egcum(ifr)=zero
        wgcm21(ifr)=zero
        wgcm2(ifr)=zero
        egcm21(ifr)=zero
        egcm2(ifr)=zero
        wsum1(ifr)=zero
        wgsum1(ifr)=zero
        wgsum(ifr)=zero
        esum1(ifr)=zero
        egsum1(ifr)=zero
        egsum(ifr)=zero
        pecum(ifr)=zero
        peicum(ifr)=zero
        tpbcum(ifr)=zero
        tjfcum(ifr)=zero
        pecm2(ifr)=zero
        peicm2(ifr)=zero
        tpbcm2(ifr)=zero
        tjfcm2(ifr)=zero
        pesum(ifr)=zero
        peisum(ifr)=zero
        tpbsum(ifr)=zero
        tjfsum(ifr)=zero
        fgcum(ifr)=zero
   85   fgcm2(ifr)=zero

      nbrnch=0

      call alloc ('taueff', taueff, nforce)
! **Warning** taueff temporarily set low.  Not any more
      if(try_int.eq.0) then
        if(idmc.eq.1 .or. idmc.eq.2) then
          taueff(1)=tau/(one+(znuc(iwctype(1))**2*tau)/10)
        elseif(idmc.eq.3) then
          taueff(1)=tau
        else
          taueff(1)=0
        endif
        write(6,'(''taueff set equal to'',f9.5)') taueff(1)
        do 86 ifr=2,nforce
          if(itau_eff.ge.1) then
            taueff(ifr)=taueff(1)
           else
            taueff(ifr)=tau
          endif
   86   continue
       else
        taueff(1)=tau*dfus2ac/dfus2unf(1)
        write(6,'(''various possibilities for mult tau are:'',3f9.5)') acc/try_int,dfus2ac/dfus2unf(1),dr2ac/dr2un
        write(6,'(''taueff set equal to tau*'',f12.5,'' ='',f9.5)') dfus2ac/dfus2unf(1),taueff(1)
        if(itausec.eq.1) then
          do 87 ifr=2,nforce
            if(itau_eff.ge.1) then
              taueff(ifr)=taueff(1)*dfus2unf(ifr)/dfus2unf(1)
             else
              taueff(ifr)=tau*dfus2unf(ifr)/dfus2unf(1)
            endif
   87     continue
          write(6,'(''secondary taueff set equal to '',20f9.5)') (taueff(ifr),ifr=2,nforce)
         else
          do 88 ifr=2,nforce
   88       taueff(ifr)=taueff(1)
        endif
      endif
      dr2ac=zero
      dr2un=zero
      dfus2ac=zero
      call alloc ('dfus2unf', dfus2unf, nforce)
      do 89 ifr=1,nforce
   89   dfus2unf(ifr)=zero
      tautot=zero
      try_int=0
      acc=0
      acc_int=0
      nodecr=0

! Zero out estimators for charge density of atom.
      do 90 i=1,NRAD
        rprobup(i)=zero
        rprobdn(i)=zero
   90   rprob(i)=zero

! Zero out estimators for pair densities:
      if (ifixe.ne.0) then
      if (.not. allocated(den2d_t)) allocate (den2d_t(-NAX:NAX,-NAX:NAX))
      if (.not. allocated(den2d_u)) allocate (den2d_u(-NAX:NAX,-NAX:NAX))
      if (.not. allocated(den2d_d)) allocate (den2d_d(-NAX:NAX,-NAX:NAX))
      if (.not. allocated(pot_ee2d_t)) allocate (pot_ee2d_t(-NAX:NAX,-NAX:NAX))
      if (.not. allocated(pot_ee2d_u)) allocate (pot_ee2d_u(-NAX:NAX,-NAX:NAX))
      if (.not. allocated(pot_ee2d_d)) allocate (pot_ee2d_d(-NAX:NAX,-NAX:NAX))
      if (.not. allocated(xx0probdt)) allocate (xx0probdt(0:NAX,-NAX:NAX,-NAX:NAX))
      if (.not. allocated(xx0probdu)) allocate (xx0probdu(0:NAX,-NAX:NAX,-NAX:NAX))
      if (.not. allocated(xx0probdd)) allocate (xx0probdd(0:NAX,-NAX:NAX,-NAX:NAX))
      if (.not. allocated(xx0probut)) allocate (xx0probut(0:NAX,-NAX:NAX,-NAX:NAX))
      if (.not. allocated(xx0probuu)) allocate (xx0probuu(0:NAX,-NAX:NAX,-NAX:NAX))
      if (.not. allocated(xx0probud)) allocate (xx0probud(0:NAX,-NAX:NAX,-NAX:NAX))
      do 100 i2=-NAX,NAX
        do 100 i3=-NAX,NAX
          den2d_t(i2,i3)=0
          den2d_u(i2,i3)=0
          den2d_d(i2,i3)=0
          pot_ee2d_t(i2,i3)=0
          pot_ee2d_u(i2,i3)=0
          pot_ee2d_d(i2,i3)=0
          do 100 i1=0,NAX
            xx0probdt(i1,i2,i3)=0
            xx0probdu(i1,i2,i3)=0
            xx0probdd(i1,i2,i3)=0
            xx0probut(i1,i2,i3)=0
            xx0probuu(i1,i2,i3)=0
  100       xx0probud(i1,i2,i3)=0
      endif
      if (ifourier.ne.0) then
      if (.not. allocated(fourierrk_t)) allocate(fourierrk_t(-NAX:NAX,0:NAK1))
      if (.not. allocated(fourierrk_u)) allocate(fourierrk_u(-NAX:NAX,0:NAK1))
      if (.not. allocated(fourierrk_d)) allocate(fourierrk_d(-NAX:NAX,0:NAK1))
      if (.not. allocated(fourierkk_t)) allocate(fourierkk_t(-NAK2:NAK2,-NAK2:NAK2))
      if (.not. allocated(fourierkk_u)) allocate(fourierkk_u(-NAK2:NAK2,-NAK2:NAK2))
      if (.not. allocated(fourierkk_d)) allocate(fourierkk_d(-NAK2:NAK2,-NAK2:NAK2))
      do 110 i1=-NAX,NAX
        do 110 i2=0,NAK1
          fourierrk_t(i1,i2)=0
          fourierrk_u(i1,i2)=0
  110     fourierrk_d(i1,i2)=0
      do 120 i1=-NAK2,NAK2
        do 120 i2=-NAK2,NAK2
          fourierkk_t(i1,i2)=0
          fourierkk_u(i1,i2)=0
  120     fourierkk_d(i1,i2)=0
      endif
      if (izigzag.ne.0) then
        call alloc_range ('zzpairden_t', zzpairden_t, -NAX, NAX, -NAX, NAX)
        call alloc_range ('zzpairdenij_t', zzpairdenij_t, -NAX, NAX, 0, (nelec-1))
        call alloc_range ('zzcorr', zzcorr, 0, NAX)
        call alloc_range ('zzcorr1', zzcorr1, 0, NAX)
        call alloc_range ('zzcorr2', zzcorr2, 0, NAX)
        call alloc_range ('zzcorrij', zzcorrij, 0, (nelec-1))
        call alloc_range ('zzcorrij1', zzcorrij1, 0, (nelec-1))
        call alloc_range ('zzcorrij2', zzcorrij2, 0, (nelec-1))
        call alloc_range ('yycorr', yycorr, 0, NAX)
        call alloc_range ('yycorr1', yycorr1, 0, NAX)
        call alloc_range ('yycorr2', yycorr2, 0, NAX)
        call alloc_range ('yycorrij', yycorrij, 0, (nelec-1))
        call alloc_range ('yycorrij1', yycorrij1, 0, (nelec-1))
        call alloc_range ('yycorrij2', yycorrij2, 0, (nelec-1))
        call alloc_range ('znncorr', znncorr, 0, NAX)
        call alloc_range ('zn2ncorr', zn2ncorr, 0, NAX)
        zzpairden_t(:,:) = 0
        zzpairdenij_t(:,:) = 0
        zzcorr(:) = 0
        zzcorr1(:) = 0
        zzcorr2(:) = 0
        zzcorrij(:) = 0
        zzcorrij1(:) = 0
        zzcorrij2(:) = 0
        yycorr(:) = 0
        yycorr1(:) = 0
        yycorr2(:) = 0
        yycorrij(:) = 0
        yycorrij1(:) = 0
        yycorrij2(:) = 0
        znncorr(:) = 0
        zn2ncorr(:) = 0
      endif



      return
      end
