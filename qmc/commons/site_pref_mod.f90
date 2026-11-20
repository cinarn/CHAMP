module site_pref_mod !NC

    implicit none
    save

    double precision, allocatable  :: site_prob(:,:)

contains

    subroutine sample_site_occupations(ncent, nelec, nup, capacity, weights, nsite)
      ! Direct categorical proposals, equivalent to the old rejection sampler
      ! for finite weights in [0,1]. Keep nsite as total counts, not spin counts.
      ! Preserve the old spin-weight schedule: down first, then up.
      use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
      implicit none
      integer, intent(in) :: ncent, nelec, nup
      double precision, intent(in) :: capacity(ncent), weights(2,ncent)
      integer, intent(out) :: nsite(ncent)
      double precision, external :: rannyu
      double precision :: mass(ncent), scale, total, threshold, cumulative, u
      integer :: remaining, spin, i, selected

      if (ncent < 1 .or. nelec < 0 .or. nup < 0 .or. nup > nelec) &
        error stop 'site sampling: invalid electron/site counts'
      if (any(.not. ieee_is_finite(capacity))) &
        error stop 'site sampling: invalid capacities'
      if (any(capacity < 0.d0)) error stop 'site sampling: invalid capacities'
      if (any(.not. ieee_is_finite(weights))) &
        error stop 'site sampling: invalid weights'
      if (any(weights < 0.d0)) error stop 'site sampling: invalid weights'
      nsite = 0
      do remaining = nelec, 1, -1
        spin = 1
        if (remaining > nup) spin = 2
        mass = 0.d0
        do i = 1, ncent
          ! Test the proposed count, preserving fractional-capacity semantics.
          if (dble(nsite(i)) + 1.d0 <= capacity(i)) mass(i) = weights(spin,i)
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
      enddo
    end subroutine sample_site_occupations

  end module site_pref_mod