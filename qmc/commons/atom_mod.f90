module atom_mod

 use constants_mod

 integer               :: nctype,ncent
 integer,  allocatable :: iwctype(:),ncent_ctype(:)
 real(dp), allocatable :: cent(:,:)
 real(dp), allocatable :: znuc(:)
 real(dp)              :: pecent, znuc_tot

end module atom_mod
