      subroutine corbitals_loc_ana(iel,rvec_en,r_en,corb,cdorb,cddorb)
c orbitals_loc_ana adapted to complex orbitals by A.D.Guclu, Feb2004
c Calculate localized orbitals and derivatives for all or 1 electrons

      use orbitals_mod, only: orb_tot_nb
      use coefs_mod
      use const_mod
      implicit real*8(a-h,o-z)
!JT      include 'vmc.h'
!JT      include 'force.h'
!JT      include 'numbas.h'

      complex*16 corb,cdorb,cddorb
      complex*16 cphin,cdphin,cd2phin

      common /dim/ ndim
!JT      common /const/ pi,hb,etrial,delta,deltai,fbias,nelec,imetro,ipr
      common /cphifun/ cphin(MBASIS,MELEC),cdphin(3,MBASIS,MELEC,MELEC)
     &,cd2phin(MBASIS,MELEC)
!JT      common /coefs/ coef(MBASIS,MORB,MWF),nbasis,norb
      common /compferm/ emagv,nv,idot
      common /numbas/ exp_h_bas(MCTYPE),r0_bas(MCTYPE)
     &,rwf(MRWF_PTS,MRWF,MCTYPE,MWF),d2rwf(MRWF_PTS,MRWF,MCTYPE,MWF)
     &,numr,nrbas(MCTYPE),igrid(MCTYPE),nr(MCTYPE),iwrwf(MBASIS_CTYPE,MCTYPE)
      common /wfsec/ iwftype(MFORCE),iwf,nwftype

      dimension rvec_en(3,MELEC,MCENT),r_en(MELEC,MCENT)
     &,corb(nelec,orb_tot_nb),cdorb(3,nelec,nelec,orb_tot_nb),cddorb(nelec,orb_tot_nb)

c Decide whether we are computing all or one electron
      if(iel.eq.0) then
        nelec1=1
        nelec2=nelec
       else
        nelec1=iel
        nelec2=iel
      endif

c get basis functions
c nej controls whether if we have correlated basis set
      if(numr.eq.1) then
	call cbasis_fns_num(iel,rvec_en,r_en)
	nej=1

      else if(idot.eq.3) then
        if(iel.ne.0) stop '1 electron move not possible with projected comp.ferm.'
        call cbasis_fns_cf(rvec_en)
        nej=nelec
      else
        call cbasis_fns_fd(iel,rvec_en,r_en)
        nej=1
      endif

c      do 5 ib=1,nbasis
c        write(6,'(''ib,cphin='',i3,(30f9.5))') ib,(cphin(i,ib),i=1,nelec)
c        write(6,'(''ib,cdphin1='',i3,(30f9.5))') ib,(cdphin(1,i,ib),i=1,nelec)
c        write(6,'(''ib,cdphin2='',i3,(30f9.5))') ib,(cdphin(1,i,ib),i=1,nelec)
c 5      write(6,'(''ib,cd2phin='',i3,(30f9.5))') ib,(cd2phin(i,ib),i=1,nelec)

      do 25 iorb=1,norb
        do 25 ie=nelec1,nelec2
          corb(ie,iorb)=dcmplx(0,0)
            do 10 idim=1,ndim
              do 10 je=1,nej
   10         cdorb(idim,ie,je,iorb)=dcmplx(0,0)
          cddorb(ie,iorb)=dcmplx(0,0)

          do 25 m=1,nbasis
            corb(ie,iorb)=corb(ie,iorb)+coef(m,iorb,iwf)*cphin(m,ie)
            do 15 idim=1,ndim
              do 15 je=1,nej
   15           cdorb(idim,ie,je,iorb)=cdorb(idim,ie,je,iorb)+coef(m,iorb,iwf)*cdphin(idim,m,ie,je)
   25       cddorb(ie,iorb)=cddorb(ie,iorb)+coef(m,iorb,iwf)*cd2phin(m,ie)

      return
      end
c-----------------------------------------------------------------------

