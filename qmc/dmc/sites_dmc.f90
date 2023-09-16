!-----------------------------------------------------------------------
      subroutine sites_dmc
!-----------------------------------------------------------------------
!     Not used at present
!     Description: random configurations for DMC on the model of the
!     routine 'sites' used in VMC
!     Created: J. Toulouse, 18 Mar 2007
!-----------------------------------------------------------------------
      use atom_mod
      use contrl_mod
      use const_mod
      use dim_mod
      use config_dmc_mod
      use site_pref_mod
      implicit real*8(a-h,o-z)

      dimension nsite(ncent)

!     sites
      ! distribute electrons randomly among the centers
      l = nelec !NC
      do i=1, ncent
         nsite(i) = 0
      end do
      do while (l .gt. 0)
         k = 1+(1+isign(1,l-nup-1))/2
         i = int(rannyu(0)*ncent + 0.5d0)
         if (rannyu(0) .le. site_prob(k,i)) then
            nsite(i) = nsite(i) + 1
            if (nsite(i) .le. znuc(iwctype(i))) then
               l = l - 1
            else
               nsite(i) = nsite(i) - 1
            end if
         end if
      end do


! loop over spins and centers. If odd number of electrons on all
! atoms then the up-spins have an additional electron.

      l=0
      do 10 ispin=1,2
        do 10 i=1,ncent
          ju=(nsite(i)+2-ispin)/2
          if(znuc(iwctype(i)).eq.0.d0) stop 'znuc should not be 0 in sites'
          do 10 j=1,ju
            l=l+1
            if(l.gt.nelec) return
            if(j.eq.1) then
              sitsca=1/znuc(iwctype(i))
             elseif(j.le.5) then
              sitsca=2/(znuc(iwctype(i))-2)
             else
              sitsca=3/(znuc(iwctype(i))-10)
            endif

!           sample position from exponentials around center
            do 10 iconf=1,nconf_global
             do 10 k=1,ndim
             site=-dlog(rannyu(0))
             site=sign(site,(rannyu(0)-0.5d0))
   10        xoldw(k,l,iconf,1)=sitsca*site+cent(k,i)

      write(6,'(i4,a,i3,a)') nconf_global,' configurations for',l,' electrons have been generated by sites.'
      write(6,'()')


      open(9,file='mc_configs_sites',status='unknown')
!write(fmt,'(a1,i2,a21)')'(',ndim*nelec,'f14.8,i3,d12.4,f12.5)'

      do iconf=1,nconf_global
       write(9,'(1000f12.8)')  ((xoldw(k,j,iconf,1),k=1,ndim),j=1,nelec)
      enddo
      close(9)

      if(l.lt.nelec) stop 'bad input to sites'
      return
      end
