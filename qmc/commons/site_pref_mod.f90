module site_pref_mod !NC

    implicit none
    save

    double precision, allocatable  :: site_prob(:,:)
    logical :: site_prob_given = .false.

contains

    subroutine sample_site_occupations(ncent, nelec, nup, capacity, weights, nsite, nsite_spin)
      ! Without nsite_spin preserve the legacy total-count sampler.
      ! Preserve the old spin-weight schedule: down first, then up.
      use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
      implicit none
      integer, intent(in) :: ncent, nelec, nup
      double precision, intent(in) :: capacity(ncent), weights(2,ncent)
      integer, intent(out) :: nsite(ncent)
      integer, intent(out), optional :: nsite_spin(ncent,2)
      double precision, external :: rannyu
      double precision :: mass(ncent), scale, total, threshold, cumulative, u
      integer :: remaining, spin, i, selected
      integer :: slots(ncent), available(3), rest(3), left(2)
      logical :: supported(3,ncent)

      if (ncent < 1 .or. nelec < 0 .or. nup < 0 .or. nup > nelec) &
        error stop 'site sampling: invalid electron/site counts'
      if (any(.not. ieee_is_finite(capacity))) &
        error stop 'site sampling: invalid capacities'
      if (any(capacity < 0.d0)) error stop 'site sampling: invalid capacities'
      if (any(.not. ieee_is_finite(weights))) &
        error stop 'site sampling: invalid weights'
      if (any(weights < 0.d0)) error stop 'site sampling: invalid weights'
      nsite = 0
      if (present(nsite_spin)) then
        nsite_spin = 0
        left = [nup, nelec-nup]
        slots = int(min(capacity, dble(nelec)))
        supported(1,:) = weights(1,:) > 0.d0
        supported(2,:) = weights(2,:) > 0.d0
        supported(3,:) = supported(1,:) .or. supported(2,:)
        do i = 1, 3
          available(i) = sum(slots, mask=supported(i,:))
        enddo
      endif
      do remaining = nelec, 1, -1
        spin = 1
        if (remaining > nup) spin = 2
        if (present(nsite_spin)) left(spin) = left(spin) - 1
        mass = 0.d0
        do i = 1, ncent
          ! Test the proposed count, preserving fractional-capacity semantics.
          if (dble(nsite(i)) + 1.d0 > capacity(i) .or. weights(spin,i) <= 0.d0) cycle
          if (present(nsite_spin)) then
            ! Reserving capacity for up, down and their union is necessary
            ! and sufficient for two spin groups with shared site capacities.
            rest = available - merge(1, 0, supported(:,i))
            if (any(rest(1:2) < left) .or. rest(3) < sum(left)) cycle
          endif
          mass(i) = weights(spin,i)
        enddo
        scale = maxval(mass)
        if (scale <= 0.d0) error stop 'site sampling: no available positive-weight site'
        ! Normalize by the largest weight before summation to avoid overflow.
        mass = mass / scale
        total = sum(mass)
        u = rannyu(0)
        threshold = u * total
        cumulative = 0.d0
        selected = 0
        do i = 1, ncent
          if (mass(i) <= 0.d0) cycle
          selected = i
          cumulative = cumulative + mass(i)
          if (threshold < cumulative) exit
        enddo
        ! selected remains the last positive site if rounding reaches total.
        nsite(selected) = nsite(selected) + 1
        if (present(nsite_spin)) then
          nsite_spin(selected,spin) = nsite_spin(selected,spin) + 1
          available = available - merge(1, 0, supported(:,selected))
        endif
      enddo
    end subroutine sample_site_occupations

    subroutine initial_site_occupations(ncent, nelec, nup, nloc, capacity, nsite, nsite_spin)
      integer, intent(in) :: ncent, nelec, nup, nloc
      double precision, intent(in) :: capacity(ncent)
      integer, intent(out) :: nsite(ncent), nsite_spin(ncent,2)
      if (nloc == -8 .and. site_prob_given) then
        call sample_site_occupations(ncent, nelec, nup, capacity, site_prob, nsite, nsite_spin)
      else
        call sample_site_occupations(ncent, nelec, nup, capacity, site_prob, nsite)
        nsite_spin(:,1) = (nsite + 1)/2
        nsite_spin(:,2) = nsite/2
      endif
    end subroutine initial_site_occupations

  end module site_pref_mod