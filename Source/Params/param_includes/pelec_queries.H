
// This file is automatically created by parse_pelec_params.py.  To update
// or add runtime parameters, please edit _cpp_parameters and then run
// mk_params.sh

pp.query("state_interp_order", state_interp_order);
pp.query("lin_limit_state_interp", lin_limit_state_interp);
pp.query("state_nghost", state_nghost);
pp.query("do_reflux", do_reflux);
pp.query("do_avg_down", do_avg_down);
pp.query("init_pltfile", init_pltfile);
pp.query("init_pltfile_massfrac_tol", init_pltfile_massfrac_tol);
pp.query("init_pltfile_coarse_levels", init_pltfile_coarse_levels);
pp.query("dump_old", dump_old);
pp.query("difmag", difmag);
pp.query("small_pres", small_pres);
pp.query("do_hydro", do_hydro);
pp.query("do_mol", do_mol);
pp.query("use_hybrid_weno", use_hybrid_weno);
pp.query("weno_scheme", weno_scheme);
pp.query("bc_nscbc", bc_nscbc);
pp.query("bc_nscbc_sigma", bc_nscbc_sigma);
pp.query("bc_nscbc_relax_u", bc_nscbc_relax_u);
pp.query("bc_nscbc_relax_t", bc_nscbc_relax_t);
pp.query("bc_nscbc_beta", bc_nscbc_beta);
pp.query("bc_nscbc_beta_s", bc_nscbc_beta_s);
pp.query("bc_nscbc_order", bc_nscbc_order);
pp.query("bc_nscbc_pin_farfield", bc_nscbc_pin_farfield);
pp.query("bc_nscbc_extrap_temperature", bc_nscbc_extrap_temperature);
pp.query("bc_nscbc_extrap_material", bc_nscbc_extrap_material);
pp.query("bc_nscbc_backflow_material", bc_nscbc_backflow_material);
pp.query("bc_nscbc_nri", bc_nscbc_nri);
pp.query("bc_nscbc_ndnr", bc_nscbc_ndnr);
pp.query("add_ext_src", add_ext_src);
{
  amrex::Vector<amrex::Real> tmp(AMREX_SPACEDIM, 0.0);
  pp.queryarr("external_forcing", tmp, 0, AMREX_SPACEDIM);
  for (int i = 0; i < tmp.size(); i++) {
    external_forcing[i] = tmp[i];
  }
}
pp.query("add_forcing_src", add_forcing_src);
pp.query("forcing_u0", forcing_u0);
pp.query("forcing_v0", forcing_v0);
pp.query("forcing_w0", forcing_w0);
pp.query("forcing_force", forcing_force);
pp.query("ppm_type", ppm_type);
pp.query("ppm_trace_sources", ppm_trace_sources);
pp.query("plm_iorder", plm_iorder);
pp.query("mol_iorder", mol_iorder);
pp.query("use_laxf_flux", use_laxf_flux);
pp.query("use_flattening", use_flattening);
pp.query("dual_energy_update_E_from_e", dual_energy_update_E_from_e);
pp.query("dual_energy_eta2", dual_energy_eta2);
pp.query("allow_negative_energy", allow_negative_energy);
pp.query("allow_small_energy", allow_small_energy);
pp.query("xl_ext_bc_type", xl_ext_bc_type);
pp.query("xr_ext_bc_type", xr_ext_bc_type);
pp.query("yl_ext_bc_type", yl_ext_bc_type);
pp.query("yr_ext_bc_type", yr_ext_bc_type);
pp.query("zl_ext_bc_type", zl_ext_bc_type);
pp.query("zr_ext_bc_type", zr_ext_bc_type);
pp.query("do_rf", do_rf);
pp.query("rf_omega", rf_omega);
pp.query("rf_axis", rf_axis);
pp.query("rf_axis_x", rf_axis_x);
pp.query("rf_axis_y", rf_axis_y);
pp.query("rf_axis_z", rf_axis_z);
pp.query("rf_rad", rf_rad);
pp.query("diffuse_temp", diffuse_temp);
pp.query("diffuse_enth", diffuse_enth);
pp.query("diffuse_spec", diffuse_spec);
pp.query("diffuse_vel", diffuse_vel);
pp.query("transport_harmonic_mean", transport_harmonic_mean);
pp.query("do_isothermal_walls", do_isothermal_walls);
{
  amrex::Vector<amrex::Real> tmp(AMREX_SPACEDIM, -1.0);
  pp.queryarr("domlo_isothermal_temp", tmp, 0, AMREX_SPACEDIM);
  for (int i = 0; i < tmp.size(); i++) {
    domlo_isothermal_temp[i] = tmp[i];
  }
}
{
  amrex::Vector<amrex::Real> tmp(AMREX_SPACEDIM, -1.0);
  pp.queryarr("domhi_isothermal_temp", tmp, 0, AMREX_SPACEDIM);
  for (int i = 0; i < tmp.size(); i++) {
    domhi_isothermal_temp[i] = tmp[i];
  }
}
pp.query("do_les", do_les);
pp.query("use_explicit_filter", use_explicit_filter);
pp.query("Cs", Cs);
pp.query("Cw", Cw);
pp.query("CI", CI);
pp.query("PrT", PrT);
pp.query("eb_boundary_T", eb_boundary_T);
pp.query("eb_isothermal", eb_isothermal);
pp.query("eb_noslip", eb_noslip);
pp.query("redistribution_type", redistribution_type);
pp.query("eb_clean_massfrac", eb_clean_massfrac);
pp.query("eb_clean_massfrac_threshold", eb_clean_massfrac_threshold);
pp.query("eb_srd_max_order", eb_srd_max_order);
pp.query("eb_weights_type", eb_weights_type);
pp.query("eb_zero_body_state", eb_zero_body_state);
pp.query("eb_problem_state", eb_problem_state);
pp.query("do_mms", do_mms);
pp.query("masa_solution_name", masa_solution_name);
pp.query("fixed_dt", fixed_dt);
pp.query("initial_dt", initial_dt);
pp.query("dt_cutoff", dt_cutoff);
pp.query("max_dt", max_dt);
pp.query("cfl", cfl);
pp.query("init_shrink", init_shrink);
pp.query("change_max", change_max);
pp.query("sdc_iters", sdc_iters);
pp.query("mol_iters", mol_iters);
pp.query("do_react", do_react);
pp.query("chem_integrator", chem_integrator);
pp.query("bndry_func_thread_safe", bndry_func_thread_safe);
pp.query("print_energy_diagnostics", print_energy_diagnostics);
pp.query("sum_interval", sum_interval);
pp.query("track_extrema", track_extrema);
pp.query("extrema_spec_name", extrema_spec_name);
pp.query("sum_per", sum_per);
pp.query("hard_cfl_limit", hard_cfl_limit);
pp.query("job_name", job_name);
pp.query("flame_trac_name", flame_trac_name);
pp.query("fuel_name", fuel_name);
