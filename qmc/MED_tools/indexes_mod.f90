module indexes_mod

! indexes of objects
  integer     :: xold_index
  integer     :: xoldw_index
  integer     :: denergy_index
  integer     :: det_ex_unq_up_index
  integer     :: det_ex_unq_dn_index
  integer     :: gvalue_index
  integer     :: g_index
  integer     :: d2d2a_index
  integer     :: d1d2a_index
  integer     :: d1d2b_index
  integer     :: d2d2b_index
  integer     :: vj_index
  integer     :: sum_lap_lnj_index
  integer     :: orb_index
  integer     :: dorb_index
  integer     :: ddorb_index
  integer     :: deti_det_index
  integer     :: psi_det_index
  integer     :: psi_jas_index
  integer     :: psijo_index
  integer     :: eloc_index
  integer     :: eloc_vmc_index
  integer     :: eloc_dmc_index
  integer     :: slmui_index
  integer     :: slmdi_index
  integer     :: detu_index
  integer     :: detd_index
  integer     :: electron_index
  integer     :: vpsp_ex_index
  integer     :: dvpsp_exp_index !fp
  integer     :: orbe_index
  integer     :: detn_index
  integer     :: slmin_index !fp
  integer     :: eloc_bav_index
  integer     :: eloc_av_index
  integer     :: eloc_av_err_index
  integer     :: eloc_pot_nloc_index
  integer     :: vold_index
  integer     :: voldw_index
  integer     :: div_vo_index
  integer     :: div_vow_index
  integer     :: dpsi_csf_index
  integer     :: dpsi_jas_index
  integer     :: dpsi_orb_index
  integer     :: dpsi_geo_index
  integer     :: dpsi_csf_av_index
  integer     :: dpsi_jas_av_index
  integer     :: dpsi_orb_av_index
  integer     :: deloc_csf_index
  integer     :: deloc_jas_index
  integer     :: deloc_orb_index
  integer     :: deloc_geo_index
  integer     :: deloc_av_index
  integer     :: eloc_pot_index
  integer     :: eloc_pot_loc_index
  integer     :: eloc_pot_nloc_ex_index
  integer     :: eloc_pot_nloc_exp_index !fp
  integer     :: jas_pairs_nb_index
  integer     :: d2psi_jas_index
  integer     :: d2eloc_jas_index
  integer     :: opt_nwt_nb_index
  integer     :: opt_lin_nb_index
  integer     :: opt_ptb_nb_index
  integer     :: grad_nwt_index
  integer     :: grad_index
  integer     :: grad_ptb_index
  integer     :: det1_det_index
  integer     :: intra_sp_histo_av_index
  integer     :: intra_sp_histo_av_err_index
  integer     :: intra_sp_zv1_av_index
  integer     :: intra_sp_zv1_av_err_index
  integer     :: intra_sp_zv2_av_index
  integer     :: intra_sp_zv2_av_err_index
  integer     :: intra_sp_zv3_av_index
  integer     :: intra_sp_zv3_av_err_index
  integer     :: intra_sp_zv4_av_index
  integer     :: intra_sp_zv4_av_err_index
  integer     :: intra_sp_zv5_av_index
  integer     :: intra_sp_zv5_av_err_index
  integer     :: intra_sp_zvzb1_av_index
  integer     :: intra_sp_zvzb1_av_err_index
  integer     :: intra_sp_zvzb3_av_index
  integer     :: intra_sp_zvzb3_av_err_index
  integer     :: intra_sp_zvzb4_av_index
  integer     :: intra_sp_zvzb4_av_err_index
  integer     :: intra_sp_zv1zb3_av_index
  integer     :: intra_sp_zv1zb3_av_err_index
  integer     :: intra_sp_zvzb5_av_index
  integer     :: intra_sp_zvzb5_av_err_index
  integer     :: hess_sor_index
  integer     :: hess_lzr_index
  integer     :: hess_uf_index
  integer     :: hess_tu_index
  integer     :: hess_lin_index
  integer     :: nparmcsf_index
  integer     :: nparmj_index
  integer     :: param_orb_nb_index
  integer     :: delta_jas_nwt_index
  integer     :: delta_jas_lin_index
  integer     :: delta_jas_ptb_index
  integer     :: delta_csf_nwt_index
  integer     :: delta_csf_lin_index
  integer     :: delta_csf_ptb_index
  integer     :: delta_coef_ex_nwt_index
  integer     :: delta_coef_ex_lin_index
  integer     :: delta_coef_ex_ptb_index
  integer     :: walker_weights_index
  integer     :: wt_index
  integer     :: fprod_index
  integer     :: eold_index
  integer     :: eoldw_index
  integer     :: dist_ee_min_index
  integer     :: dist_ee_max_index
  integer     :: nwalk_index
  integer     :: phin_index
  integer     :: dphin_index
  integer     :: d2phin_index
  integer     :: r_en_index
  integer     :: rvec_en_index
  integer     :: param_exp_nb_index
  integer     :: ddet_dexp_unq_up_index
  integer     :: ddet_dexp_unq_dn_index
  integer     :: dpsi_exp_index
  integer     :: dpsi_lnexp_index
  integer     :: deloc_exp_index
  integer     :: deloc_lnexp_index
  integer     :: gradient_variance_index
  integer     :: hessian_variance_lm_index
  integer     :: hessian_variance_lmcov_index
  integer     :: hessian_variance_lin_index
  integer     :: hessian_variance_index
  integer     :: pe_ee_index
  integer     :: pe_en_index
  integer     :: dphin_dz_index
  integer     :: grd_dphin_dz_index
  integer     :: lap_dphin_dz_index
  integer     :: dphin_norm_dz_index
  integer     :: grd_dphin_norm_dz_index
  integer     :: lap_dphin_norm_dz_index
  integer     :: dphin_ortho_dz_index
  integer     :: grd_dphin_ortho_dz_index
  integer     :: lap_dphin_ortho_dz_index
  integer     :: coef_index
  integer     :: coef_orb_on_norm_basis_index
  integer     :: coef_orb_on_ortho_basis_index
  integer     :: psid_ex_in_x_index
  integer     :: add_diag_mult_exp_index
  integer     :: is_param_type_geo_index
  integer     :: dpsi_av_index
  integer     :: dpsi_var_index
  integer     :: deloc_var_index

! indexes of nodes
  integer     :: eloc_bld_index
  integer     :: dpsi_bld_index
  integer     :: dpsi_nwt_bld_index
  integer     :: dpsi_lin_bld_index
  integer     :: dpsi_ptb_bld_index
  integer     :: deloc_bld_index
  integer     :: deloc_nwt_bld_index
  integer     :: deloc_lin_bld_index
  integer     :: deloc_ptb_bld_index
  integer     :: d2psi_bld_index
  integer     :: d2psi_nwt_bld_index
  integer     :: d2eloc_nwt_bld_index
  integer     :: eloc_pot_ex_bld_index
  integer     :: eloc_pot_exp_bld_index
  integer     :: gradient_bld_index
  integer     :: intra_sp_bld_index
  integer     :: hess_nwt_energy_bld_index
  integer     :: hess_nwt_bld_index
  integer     :: param_nb_bld_index
  integer     :: delta_jas_bld_index
  integer     :: delta_csf_bld_index
  integer     :: delta_coef_ex_bld_index
  integer     :: walker_weights_bld_index
  integer     :: walker_weights_sum_bld_index
  integer     :: wgcum_index
  integer     :: eloc_wlk_bld_index
  integer     :: eloc_wlk_test_bld_index
  integer     :: eloc_test2_bld_index
  integer     :: eloc_test3_bld_index
  integer     :: eloc_test4_bld_index
  integer     :: eloc_wlk_test2_bld_index
  integer     :: coord_elec_bld_index
  integer     :: coord_elec_wlk_bld_index
  integer     :: grd_psi_over_psi_wlk_bld_index
  integer     :: grd_psi_over_psi_sq_wlk_bld_index
  integer     :: div_grd_psi_over_psi_wlk_bld_index
  integer     :: hessian_variance_bld_index
  integer     :: dorb_dexp_bld_index
  integer     :: grd_dorb_dexp_bld_index
  integer     :: lap_dorb_dexp_bld_index
  integer     :: gradient_energy_bld_index
  !! pjasen
  integer    :: param_pjasen_nb_index
  integer    :: grad_dpsi_pjasen_index
  integer    :: lap_dpsi_pjasen_index
  integer    :: deloc_pjasen_index
  integer    :: dpsi_pjasen_index
  integer    :: d2psi_pjasen_index
!! WAS
  !! pjasee
  integer    :: param_pjasee_nb_index
  integer    :: grad_dpsi_pjasee_index
  integer    :: lap_dpsi_pjasee_index
  integer    :: deloc_pjasee_index
  integer    :: dpsi_pjasee_index
  integer    :: d2psi_pjasee_index
  !! pjas
  integer    :: param_pjas_nb_index
  integer    :: grad_dpsi_pjas_index
  integer    :: lap_dpsi_pjas_index
  integer    :: deloc_pjas_index
  integer    :: dpsi_pjas_index
  integer    :: d2psi_pjas_index
  integer    :: vd_index
  integer    :: star_en_index
  integer    :: cos_star_en_index
  integer    :: sin_star_en_index
  integer    :: grad_star_en_index
  integer    :: grad_cos_star_en_index
  integer    :: grad_sin_star_en_index
  integer    :: psid_pjas_index
  integer    :: star_sum_en_index
  integer    :: cos_star_ee_index
  integer    :: grad_cos_star_ee_index

  integer    :: sin_star_ee_index
  integer    :: grad_sin_star_ee_index
  integer    :: sin_star_sum_ee_index

  integer    :: star_ee_index
  integer    :: grad_star_ee_index

  integer    :: dvpsp_pjas_index
  integer    :: deloc_pot_nloc_pjas_index
!!!
  integer     :: param_geo_nb_index
!! solid orbitals
!  integer    :: ngvec_orb_index
  integer     :: current_walker_index
  integer     :: current_walker_weight_index
  integer     :: total_iterations_block_nb_index
  integer     :: total_iterations_nb_index
  integer     :: walker_weights_sum_block_index
  integer     :: walker_weights_sum_index
  integer     :: walker_weights_sq_sum_block_index
  integer     :: walker_weights_sq_sum_index

end module indexes_mod
