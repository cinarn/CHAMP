      subroutine bffden2d(p,q,xold,xnew)

! Body-Fixed Frame Density utilizing the M-fold Orientational Order Parameter.
! Uses strict Laboratory Frame coordinates (no Center of Mass translation).
! Written by Gokhan Oztarhan, Sep 2026.

      use dets_mod
      use const_mod
      use dim_mod
      use pairden_mod
      implicit real*8(a-h,o-z)

      dimension xold(3,nelec),xnew(3,nelec)
      
      real*8 :: p, q
      real*8 :: rad_o(nelec), rad_n(nelec)
      integer :: idx_o(nelec), idx_n(nelec)
      integer :: ref_o(M_bff), ref_n(M_bff)
      real*8 :: C_o, S_o, C_n, S_n
      real*8 :: Theta_body_o, Theta_body_n, theta_j
      real*8 :: xbody_o(nelec), ybody_o(nelec)
      real*8 :: xbody_n(nelec), ybody_n(nelec)
      integer :: n_skip_bff, i, j, k, itemp, iring, m_shell
      real*8 :: mag_o, mag_n
      
      ! Local variables for Inter-Ring Phase Correlation
      real*8 :: Cin_o_bff, Sin_o_bff, Cout_o_bff, Sout_o_bff
      real*8 :: Cin_n_bff, Sin_n_bff, Cout_n_bff, Sout_n_bff
      real*8 :: Thetain_o_bff, Thetaout_o_bff, dTheta_o_bff
      real*8 :: Thetain_n_bff, Thetaout_n_bff, dTheta_n_bff
      integer :: Min_bff, Mout_bff, ibin_o_bff, ibin_n_bff
      integer :: n_skip_in_bff, n_skip_out_bff
      real*8 :: range_half_bff
      
      integer :: ix1roto, ix2roto, ix1rotn, ix2rotn

      ! 1. Lab frame radii (No Center of Mass translation)
      do i = 1, nelec
         rad_o(i) = dsqrt(xold(1,i)**2 + xold(2,i)**2)
         idx_o(i) = i

         rad_n(i) = dsqrt(xnew(1,i)**2 + xnew(2,i)**2)
         idx_n(i) = i
      end do

      ! 2. Radial Sorting (Bubble Sort)
      do i = 1, nelec - 1
         do j = i + 1, nelec
            if (rad_o(idx_o(i)) .gt. rad_o(idx_o(j))) then
               itemp = idx_o(i)
               idx_o(i) = idx_o(j)
               idx_o(j) = itemp
            end if
            if (rad_n(idx_n(i)) .gt. rad_n(idx_n(j))) then
               itemp = idx_n(i)
               idx_n(i) = idx_n(j)
               idx_n(j) = itemp
            end if
         end do
      end do
      
      ! 3. Accumulate Order Parameters for ALL valid shells in conf_bff
      n_skip_bff = 0
      do iring = 1, nrings_bff
         m_shell = conf_bff(iring)
         ! For a singleton this gives |Q_1|=1, which carries no ordering information.
         if (m_shell .gt. 0) then
            C_o = 0.d0; S_o = 0.d0
            C_n = 0.d0; S_n = 0.d0
            do k = 1, m_shell
               j = idx_o(n_skip_bff + k)
               theta_j = datan2(xold(2,j), xold(1,j))
               C_o = C_o + dcos(m_shell * theta_j)
               S_o = S_o + dsin(m_shell * theta_j)

               j = idx_n(n_skip_bff + k)
               theta_j = datan2(xnew(2,j), xnew(1,j))
               C_n = C_n + dcos(m_shell * theta_j)
               S_n = S_n + dsin(m_shell * theta_j)
            end do
            mag_o = dsqrt(C_o**2 + S_o**2) / m_shell
            mag_n = dsqrt(C_n**2 + S_n**2) / m_shell
            psi_M_acc_bff(iring) = psi_M_acc_bff(iring) + q * mag_o + p * mag_n
         end if
         n_skip_bff = n_skip_bff + m_shell
      end do

      ! ---------------------------------------------------------
      ! Inter-Ring Phase Correlation (Executes if 2 or more rings)
      ! ---------------------------------------------------------
      if (nrings_bff .ge. 2) then
         ! Automatically target the two outermost rings
         Mout_bff = conf_bff(nrings_bff)
         Min_bff = conf_bff(nrings_bff - 1)
         
         n_skip_in_bff = 0
         do i = 1, nrings_bff - 2
            n_skip_in_bff = n_skip_in_bff + conf_bff(i)
         end do
         n_skip_out_bff = n_skip_in_bff + Min_bff

         Cin_o_bff = 0.d0; Sin_o_bff = 0.d0
         Cout_o_bff = 0.d0; Sout_o_bff = 0.d0
         Cin_n_bff = 0.d0; Sin_n_bff = 0.d0
         Cout_n_bff = 0.d0; Sout_n_bff = 0.d0

         ! Inner ring
         do k = 1, Min_bff
            j = idx_o(n_skip_in_bff + k)
            theta_j = datan2(xold(2,j), xold(1,j))
            Cin_o_bff = Cin_o_bff + dcos(Min_bff * theta_j)
            Sin_o_bff = Sin_o_bff + dsin(Min_bff * theta_j)

            j = idx_n(n_skip_in_bff + k)
            theta_j = datan2(xnew(2,j), xnew(1,j))
            Cin_n_bff = Cin_n_bff + dcos(Min_bff * theta_j)
            Sin_n_bff = Sin_n_bff + dsin(Min_bff * theta_j)
         end do
         Thetain_o_bff = (1.d0 / Min_bff) * datan2(Sin_o_bff, Cin_o_bff)
         Thetain_n_bff = (1.d0 / Min_bff) * datan2(Sin_n_bff, Cin_n_bff)

         ! Outer ring
         do k = 1, Mout_bff
            j = idx_o(n_skip_out_bff + k)
            theta_j = datan2(xold(2,j), xold(1,j))
            Cout_o_bff = Cout_o_bff + dcos(Mout_bff * theta_j)
            Sout_o_bff = Sout_o_bff + dsin(Mout_bff * theta_j)

            j = idx_n(n_skip_out_bff + k)
            theta_j = datan2(xnew(2,j), xnew(1,j))
            Cout_n_bff = Cout_n_bff + dcos(Mout_bff * theta_j)
            Sout_n_bff = Sout_n_bff + dsin(Mout_bff * theta_j)
         end do
         Thetaout_o_bff = (1.d0 / Mout_bff) * datan2(Sout_o_bff, Cout_o_bff)
         Thetaout_n_bff = (1.d0 / Mout_bff) * datan2(Sout_n_bff, Cout_n_bff)

         ! Phase difference
         dTheta_o_bff = Thetaout_o_bff - Thetain_o_bff
         dTheta_n_bff = Thetaout_n_bff - Thetain_n_bff

         ! Fold by the common harmonic of these two rings, not the selected BFF ring.
         dTheta_o_bff = modulo(dTheta_o_bff + pi/dble(M_ir_bff), 2.d0*pi/dble(M_ir_bff)) - pi/dble(M_ir_bff)
         dTheta_n_bff = modulo(dTheta_n_bff + pi/dble(M_ir_bff), 2.d0*pi/dble(M_ir_bff)) - pi/dble(M_ir_bff)

         ! Set binning range to the modulo bounds.
         range_half_bff = pi / dble(M_ir_bff)

         ibin_o_bff = int((dTheta_o_bff + range_half_bff) / (2.d0 * range_half_bff) * NIRBINS_bff) + 1
         if (ibin_o_bff .lt. 1) ibin_o_bff = 1
         if (ibin_o_bff .gt. NIRBINS_bff) ibin_o_bff = NIRBINS_bff
         irphase_bff(ibin_o_bff) = irphase_bff(ibin_o_bff) + q

         ibin_n_bff = int((dTheta_n_bff + range_half_bff) / (2.d0 * range_half_bff) * NIRBINS_bff) + 1
         if (ibin_n_bff .lt. 1) ibin_n_bff = 1
         if (ibin_n_bff .gt. NIRBINS_bff) ibin_n_bff = NIRBINS_bff
         irphase_bff(ibin_n_bff) = irphase_bff(ibin_n_bff) + p
      end if

      ! ---------------------------------------------------------
      ! M-fold Target Angle & Coordinate Transformation
      ! ---------------------------------------------------------
      ! Partition only the selected main shell; main-shell correlation is unchanged.
      call bff_subring_indices(xold,rad_o,idx_o,ref_o)
      call bff_subring_indices(xnew,rad_n,idx_n,ref_n)

! Retain the original main-shell |Q| accumulators and phase correlations above.
! Only the reference angle used for the density rotation changes.
      call bff_geometric_angle(xold,ref_o,Theta_body_o)
      call bff_geometric_angle(xnew,ref_n,Theta_body_n)

      ! Transform configurations into the Body-Fixed Frame using Lab-Frame origin
      do i = 1, nelec
         xbody_o(i) = xold(1,i) * dcos(-Theta_body_o) - xold(2,i) * dsin(-Theta_body_o)
         ybody_o(i) = xold(1,i) * dsin(-Theta_body_o) + xold(2,i) * dcos(-Theta_body_o)

         xbody_n(i) = xnew(1,i) * dcos(-Theta_body_n) - xnew(2,i) * dsin(-Theta_body_n)
         ybody_n(i) = xnew(1,i) * dsin(-Theta_body_n) + xnew(2,i) * dcos(-Theta_body_n)
      end do

      ! 5. UNCONDITIONAL Binning in the Body-Fixed Frame
      do ie2 = 1, nelec
         ! Old configuration binning
         ix1roto = nint(delxi(1) * xbody_o(ie2))
         ix2roto = nint(delxi(2) * ybody_o(ie2))

         if (ix1roto .ge. -NAX .and. ix1roto .le. NAX .and. ix2roto .ge. -NAX .and. ix2roto .le. NAX) then
            bffden2d_t(ix1roto,ix2roto) = bffden2d_t(ix1roto,ix2roto) + q
            if (ie2 .le. nup) then
               bffden2d_u(ix1roto,ix2roto) = bffden2d_u(ix1roto,ix2roto) + q
            else
               bffden2d_d(ix1roto,ix2roto) = bffden2d_d(ix1roto,ix2roto) + q
            endif
         end if

         ! New configuration binning
         ix1rotn = nint(delxi(1) * xbody_n(ie2))
         ix2rotn = nint(delxi(2) * ybody_n(ie2))

         if (ix1rotn .ge. -NAX .and. ix1rotn .le. NAX .and. ix2rotn .ge. -NAX .and. ix2rotn .le. NAX) then
            bffden2d_t(ix1rotn,ix2rotn) = bffden2d_t(ix1rotn,ix2rotn) + p
            if (ie2 .le. nup) then
               bffden2d_u(ix1rotn,ix2rotn) = bffden2d_u(ix1rotn,ix2rotn) + p
            else
               bffden2d_d(ix1rotn,ix2rotn) = bffden2d_d(ix1rotn,ix2rotn) + p
            endif
         end if
      end do

      return
      end
      subroutine prepare_bff_input
! Input lists need storage before reading the optional namelist, even if disabled.
      use const_mod, only: nelec
      use pairden_mod
      implicit none
      if(allocated(conf_bff)) deallocate(conf_bff)
      if(allocated(sub_conf_bff)) deallocate(sub_conf_bff)
      if(allocated(psi_M_acc_bff)) deallocate(psi_M_acc_bff)
      allocate(conf_bff(nelec),sub_conf_bff(nelec),psi_M_acc_bff(nelec))
      conf_bff=0
      sub_conf_bff=0
      psi_M_acc_bff=0.d0
      sub_id_bff=0
      nrings_bff=0
      nsubrings_bff=0
      M_bff=0
      M_ir_bff=0
      end subroutine prepare_bff_input


      subroutine initialize_bff
! Validate the fixed configuration and allocate maps and assignment caches once.
      use const_mod, only: nelec, pi
      use pairden_mod
      implicit none
      integer :: iring, isub, first, population, ia, ib, irem, m, s, nrow, ncol

      M_bff=0
      M_ir_bff=0
      nrings_bff=0
      nsubrings_bff=0
      if(sub_id_bff.lt.0) stop 'sub_id_bff must be zero or a positive sub-ring index'
      if(sub_id_bff.eq.0) return
      nrings_bff=count(conf_bff.gt.0)
      if(nrings_bff.eq.0) stop 'conf_bff must contain the main shell populations'
      if(any(conf_bff(1:nrings_bff).le.0).or.any(conf_bff(nrings_bff+1:).ne.0)) &
     &  stop 'conf_bff must contain positive populations followed by zeros'
      if(sum(conf_bff).ne.nelec) stop 'conf_bff populations must sum to nelec'
      if(all(sub_conf_bff.eq.0)) sub_conf_bff=conf_bff
      nsubrings_bff=count(sub_conf_bff.gt.0)
      if(any(sub_conf_bff(1:nsubrings_bff).le.0).or.any(sub_conf_bff(nsubrings_bff+1:).ne.0)) &
     &  stop 'sub_conf_bff must contain positive populations followed by zeros'
      if(sum(sub_conf_bff).ne.nelec) stop 'sub_conf_bff populations must sum to nelec'
      if(sub_id_bff.gt.nsubrings_bff) stop 'sub_id_bff exceeds the number of sub-rings'

      if(allocated(first_sub_bff)) deallocate(first_sub_bff,nsub_bff,sub_shell_bff)
      allocate(first_sub_bff(nrings_bff),nsub_bff(nrings_bff),sub_shell_bff(nsubrings_bff))
      isub=1
      do iring=1,nrings_bff
        first=isub
        first_sub_bff(iring)=first
        population=0
        do while(population.lt.conf_bff(iring))
          if(isub.gt.nsubrings_bff) stop 'sub_conf_bff does not cover every main shell'
          if(sub_conf_bff(isub).ne.sub_conf_bff(first)) &
     &      stop 'sub-rings within one main shell must have equal populations'
          population=population+sub_conf_bff(isub)
          if(population.gt.conf_bff(iring)) stop 'sub_conf_bff crosses a main shell boundary'
          sub_shell_bff(isub)=iring
          isub=isub+1
        enddo
        nsub_bff(iring)=isub-first
      enddo
      M_bff=sub_conf_bff(sub_id_bff)
      if(allocated(bff_angle_order)) deallocate(bff_angle_order,bff_angle_theta,bff_vertex_cos,bff_vertex_sin)
      allocate(bff_angle_order(M_bff),bff_angle_theta(M_bff),bff_vertex_cos(M_bff),bff_vertex_sin(M_bff))
      do s=1,M_bff
        bff_vertex_cos(s)=dcos(2.d0*pi*dble(s-1)/dble(M_bff))
        bff_vertex_sin(s)=dsin(2.d0*pi*dble(s-1)/dble(M_bff))
      enddo

      if(allocated(bff_score)) then
        deallocate(bff_slot_cos,bff_slot_sin,bff_score,bff_unit_x,bff_unit_y,bff_mean_radius, &
     &    bff_u,bff_v,bff_minv,bff_members,bff_group_order,bff_available,bff_assignment, &
     &    bff_rows,bff_p,bff_way,bff_used,bff_visited)
      endif
      nrow=max(1,maxval(sub_conf_bff)-1)
      ncol=max(1,maxval(conf_bff)-1)
      allocate(bff_slot_cos(nrow,nrings_bff),bff_slot_sin(nrow,nrings_bff),bff_score(nrow,ncol), &
     &  bff_unit_x(nelec),bff_unit_y(nelec),bff_mean_radius(nelec), &
     &  bff_u(0:nrow),bff_v(0:ncol),bff_minv(ncol), &
     &  bff_members(nelec),bff_group_order(nelec),bff_available(ncol),bff_assignment(nrow), &
     &  bff_rows(nrow),bff_p(0:ncol),bff_way(ncol),bff_used(nelec),bff_visited(0:ncol))

! Pair opposite slots explicitly so reflection permutes identical cached scores.
      bff_slot_cos=0.d0
      bff_slot_sin=0.d0
      do iring=1,nrings_bff
        m=sub_conf_bff(first_sub_bff(iring))
        do s=1,(m-1)/2
          bff_slot_cos(s,iring)=dcos(2.d0*pi*dble(s)/dble(m))
          bff_slot_sin(s,iring)=dsin(2.d0*pi*dble(s)/dble(m))
          bff_slot_cos(m-s,iring)=bff_slot_cos(s,iring)
          bff_slot_sin(m-s,iring)=-bff_slot_sin(s,iring)
        enddo
        if(mod(m,2).eq.0) then
          bff_slot_cos(m/2,iring)=-1.d0
          bff_slot_sin(m/2,iring)=0.d0
        endif
      enddo

! Inter-ring correlation always refers to the two outermost MAIN shells.
      if(nrings_bff.ge.2) then
        ia=conf_bff(nrings_bff-1)
        ib=conf_bff(nrings_bff)
        do while(ib.ne.0)
          irem=mod(ia,ib)
          ia=ib
          ib=irem
        enddo
        M_ir_bff=(conf_bff(nrings_bff-1)/ia)*conf_bff(nrings_bff)
      endif
      write(6,'(a,i5,a,i5,a,i5)') 'BFF reference sub-ring=',sub_id_bff,', derived M=',M_bff, &
     &  ', main shell=',sub_shell_bff(sub_id_bff)
      write(6,'(a,*(i5))') 'BFF main shells: ',conf_bff(1:nrings_bff)
      write(6,'(a,*(i5))') 'BFF sub-rings:   ',sub_conf_bff(1:nsubrings_bff)
      end subroutine initialize_bff


      subroutine bff_subring_indices(x,rad,idx,iref)
! Exact distinct-partner assignment for each innermost-unused anchor.
! The result is a sequential group construction, not a whole-shell optimization.
      use const_mod, only: nelec
      use pairden_mod
      implicit none
      real*8, intent(in) :: x(3,nelec), rad(nelec)
      integer, intent(in) :: idx(nelec)
      integer, intent(out) :: iref(M_bff)
      real*8 :: cos_delta, sin_delta
      integer :: iring, first, last, ngroup, igroup, offset, slot, col, nfree
      integer :: k, j, anchor, partner, itemp

      iring=sub_shell_bff(sub_id_bff)
      first=sum(conf_bff(1:iring-1))+1
      last=first+conf_bff(iring)-1
      ngroup=nsub_bff(iring)
      if(ngroup.eq.1) then
        iref=idx(first:last)
        return
      endif
      bff_used=.false.
      if(M_bff.gt.1) then
        do k=first,last
          j=idx(k)
          if(rad(j).gt.0.d0) then
            bff_unit_x(j)=x(1,j)/rad(j)
            bff_unit_y(j)=x(2,j)/rad(j)
          else
! The origin has no azimuth; retain the base atan2(0,0)=0 convention.
            bff_unit_x(j)=1.d0
            bff_unit_y(j)=0.d0
          endif
        enddo
      endif
      do igroup=1,ngroup
        do k=first,last
          anchor=idx(k)
          if(.not.bff_used(anchor)) exit
        enddo
        offset=(igroup-1)*M_bff
        bff_members(offset+1)=anchor
        bff_used(anchor)=.true.
        bff_mean_radius(igroup)=rad(anchor)
        if(M_bff.gt.1) then
          nfree=0
          do k=first,last
            j=idx(k)
            if(bff_used(j)) cycle
            nfree=nfree+1
            bff_available(nfree)=j
          enddo
          if(nfree.eq.M_bff-1) then
! Final membership is forced; slot permutations cannot change its M-fold angle.
            do slot=1,M_bff-1
              bff_assignment(slot)=slot
            enddo
          else
            do col=1,nfree
              j=bff_available(col)
              cos_delta=bff_unit_x(j)*bff_unit_x(anchor)+bff_unit_y(j)*bff_unit_y(anchor)
              sin_delta=bff_unit_y(j)*bff_unit_x(anchor)-bff_unit_x(j)*bff_unit_y(anchor)
              do slot=1,M_bff-1
                bff_score(slot,col)=cos_delta*bff_slot_cos(slot,iring)+sin_delta*bff_slot_sin(slot,iring)
              enddo
            enddo
            call bff_assign_partners(M_bff-1,nfree)
          endif
          do slot=1,M_bff-1
            partner=bff_available(bff_assignment(slot))
            bff_members(offset+slot+1)=partner
            bff_used(partner)=.true.
            bff_mean_radius(igroup)=bff_mean_radius(igroup)+rad(partner)
          enddo
        endif
        bff_mean_radius(igroup)=bff_mean_radius(igroup)/dble(M_bff)
        bff_group_order(igroup)=igroup
      enddo
      do igroup=1,ngroup-1
        do k=igroup+1,ngroup
          if(bff_mean_radius(bff_group_order(igroup)).gt.bff_mean_radius(bff_group_order(k))) then
            itemp=bff_group_order(igroup)
            bff_group_order(igroup)=bff_group_order(k)
            bff_group_order(k)=itemp
          endif
        enddo
      enddo
      igroup=bff_group_order(sub_id_bff-first_sub_bff(iring)+1)
      offset=(igroup-1)*M_bff
      iref=bff_members(offset+1:offset+M_bff)
      end subroutine bff_subring_indices


      subroutine bff_assign_partners(nrow,ncol)
! Rectangular Hungarian assignment, minimizing cost=-cosine score.
! Canonical row ordering removes dependence of tie resolution on slot numbering.
      use pairden_mod
      implicit none
      integer, intent(in) :: nrow,ncol
      integer :: i, j, k, row, j0, j1, i0, itemp
      real*8 :: delta, cur
      logical :: before

      if(nrow.eq.1) then
! One slot is already an exact maximum-score assignment.
        bff_assignment(1)=maxloc(bff_score(1,1:ncol),dim=1)
        return
      endif
      do i=1,nrow
        bff_rows(i)=i
      enddo
      do i=2,nrow
        row=bff_rows(i)
        k=i-1
        do while(k.ge.1)
          before=.false.
          do j=1,ncol
            if(bff_score(row,j).gt.bff_score(bff_rows(k),j)) then
              before=.true.
              exit
            elseif(bff_score(row,j).lt.bff_score(bff_rows(k),j)) then
              exit
            endif
          enddo
          if(.not.before) exit
          bff_rows(k+1)=bff_rows(k)
          k=k-1
        enddo
        bff_rows(k+1)=row
      enddo
      bff_u=0.d0
      bff_v=0.d0
      bff_p=0
      do i=1,nrow
        bff_p(0)=i
        j0=0
        bff_minv(1:ncol)=huge(1.d0)
        bff_visited=.false.
        do
          bff_visited(j0)=.true.
          i0=bff_p(j0)
          delta=huge(1.d0)
          do j=1,ncol
            if(bff_visited(j)) cycle
            cur=-bff_score(bff_rows(i0),j)-bff_u(i0)-bff_v(j)
            if(cur.lt.bff_minv(j)) then
              bff_minv(j)=cur
              bff_way(j)=j0
            endif
            if(bff_minv(j).lt.delta) then
              delta=bff_minv(j)
              j1=j
            endif
          enddo
          do j=0,ncol
            if(bff_visited(j)) then
              bff_u(bff_p(j))=bff_u(bff_p(j))+delta
              bff_v(j)=bff_v(j)-delta
            elseif(j.gt.0) then
              bff_minv(j)=bff_minv(j)-delta
            endif
          enddo
          j0=j1
          if(bff_p(j0).eq.0) exit
        enddo
        do
          j1=bff_way(j0)
          bff_p(j0)=bff_p(j1)
          j0=j1
          if(j0.eq.0) exit
        enddo
      enddo
      do j=1,ncol
        if(bff_p(j).gt.0) bff_assignment(bff_rows(bff_p(j)))=j
      enddo
      end subroutine bff_assign_partners

      subroutine bff_geometric_angle(x,iref,theta_body)
! Rotation-only fit to distinct regular-polygon vertices in circular order.
! Cyclic relabeling changes the angle by an equivalent 2*pi/M_bff.
      use const_mod, only: nelec, pi
      use pairden_mod, only: M_bff, bff_angle_order, bff_angle_theta, bff_vertex_cos, bff_vertex_sin
      implicit none
      real*8, intent(in) :: x(3,nelec)
      integer, intent(in) :: iref(M_bff)
      real*8, intent(out) :: theta_body
      real*8 :: angle, c, s, period
      integer :: i, j, itemp, ie

      do i=1,M_bff
        ie=iref(i)
        bff_angle_theta(i)=modulo(datan2(x(2,ie),x(1,ie)),2.d0*pi)
        bff_angle_order(i)=i
      enddo
      do i=2,M_bff
        itemp=bff_angle_order(i)
        angle=bff_angle_theta(itemp)
        j=i-1
        do while(j.ge.1)
          if(bff_angle_theta(bff_angle_order(j)).le.angle) exit
          bff_angle_order(j+1)=bff_angle_order(j)
          j=j-1
        enddo
        bff_angle_order(j+1)=itemp
      enddo
      c=0.d0
      s=0.d0
      do i=1,M_bff
        ie=iref(bff_angle_order(i))
        c=c+x(1,ie)*bff_vertex_cos(i)+x(2,ie)*bff_vertex_sin(i)
        s=s+x(2,ie)*bff_vertex_cos(i)-x(1,ie)*bff_vertex_sin(i)
      enddo
      period=2.d0*pi/dble(M_bff)
      theta_body=modulo(datan2(s,c)+period/2.d0,period)-period/2.d0
      end subroutine bff_geometric_angle
