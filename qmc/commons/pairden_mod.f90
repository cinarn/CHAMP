module pairden_mod

 implicit none
 save
 
 double precision, allocatable :: xx0probut(:,:,:),xx0probuu(:,:,:),xx0probud(:,:,:),xx0probdt(:,:,:)
 double precision, allocatable :: xx0probdu(:,:,:),xx0probdd(:,:,:),den2d_t(:,:),den2d_d(:,:),den2d_u(:,:)
 double precision, allocatable :: pot_ee2d_t(:,:), pot_ee2d_u(:,:), pot_ee2d_d(:,:)       
 double precision, allocatable :: dos(:)
 double precision, allocatable :: thetafix(:) !GO
 integer, allocatable          :: imeshfix1(:), imeshfix2(:) !GO
 double precision              :: dos_dele
 double precision              :: delxi(3),xmax,xfix(3)
 double precision              :: pair_hits_u, pair_hits_d  ! Exact trigger trackers
 integer ithetafix !GO
 integer ifixe
 integer, parameter :: NAX = 50
 
! Main radial shells and mean-radius-ranked angular sub-rings.
! Only conf_bff, sub_conf_bff and sub_id_bff are input parameters.
 integer :: sub_id_bff, M_bff, nrings_bff, nsubrings_bff
 integer, allocatable :: conf_bff(:), sub_conf_bff(:)
 integer, allocatable :: first_sub_bff(:), nsub_bff(:), sub_shell_bff(:)
 double precision, allocatable :: bffden2d_t(:,:), bffden2d_u(:,:), bffden2d_d(:,:)
 double precision, allocatable :: psi_M_acc_bff(:)

! Fixed slot angles and sequential assignment workspace, allocated once per setup.
 double precision, allocatable :: bff_slot_cos(:,:), bff_slot_sin(:,:), bff_score(:,:)
 double precision, allocatable :: bff_unit_x(:), bff_unit_y(:), bff_mean_radius(:)
 double precision, allocatable :: bff_u(:), bff_v(:), bff_minv(:)
 integer, allocatable :: bff_members(:), bff_group_order(:), bff_available(:), bff_assignment(:)
 integer, allocatable :: bff_rows(:), bff_p(:), bff_way(:)
 logical, allocatable :: bff_used(:), bff_visited(:)
 
! Variables for Inter-Ring Phase Correlation
 integer :: M_ir_bff ! Common harmonic of the two outermost shells.
 integer, parameter :: NIRBINS_bff = 100
 double precision, allocatable :: irphase_bff(:)

end module pairden_mod
