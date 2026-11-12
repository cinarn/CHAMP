module contrl_mod

 implicit none
 save

 integer               :: nstep,nblk,nblkeq,nconf,nconf_global,nconf_new,isite,idump,irstar
 integer :: irot_order=3, irot_every=0
 integer, dimension(4) :: irand_seed

end module contrl_mod
