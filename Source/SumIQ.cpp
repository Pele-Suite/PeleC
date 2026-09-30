#include <iomanip>

#include "PeleC.H"

void
PeleC::sum_integrated_quantities()
{
  BL_PROFILE("PeleC::sum_integrated_quantities()");

  if (verbose <= 0) {
    return;
  }

  nscbc_report_diagnostics();

  bool local_flag = true;

  int finest_level = parent->finestLevel();
  amrex::Real time = state[State_Type].curTime();
  amrex::Real mass = 0.0;
  amrex::Real mom[3] = {0.0};
  amrex::Real rho_e = 0.0;
  amrex::Real rho_K = 0.0;
  amrex::Real rho_E = 0.0;
  amrex::Real fuel_prod = 0;
  amrex::Real temp = 0;

  for (int lev = 0; lev <= finest_level; lev++) {
    PeleC& pc_lev = getLevel(lev);

    mass += pc_lev.volWgtSum("density", time, local_flag);
    mom[0] += pc_lev.volWgtSum("xmom", time, local_flag);
    mom[1] += pc_lev.volWgtSum("ymom", time, local_flag);
    mom[2] += pc_lev.volWgtSum("zmom", time, local_flag);
    rho_e += pc_lev.volWgtSum("rho_e", time, local_flag);
    rho_K += pc_lev.volWgtSum("kineng", time, local_flag);
    rho_E += pc_lev.volWgtSum("rho_E", time, local_flag);

    if (!fuel_name.empty()) {
      fuel_prod += pc_lev.volWgtSum("rho_omega_" + fuel_name, time, local_flag);
    }

    temp += pc_lev.volWgtSum("Temp", time, local_flag);
  }

  if (verbose > 0) {
    const int nfoo = 9;
    amrex::Real foo[nfoo] = {mass,  mom[0], mom[1],    mom[2], rho_e,
                             rho_K, rho_E,  fuel_prod, temp};
    amrex::ParallelDescriptor::ReduceRealSum(
      foo, nfoo, amrex::ParallelDescriptor::IOProcessorNumber());

    if (amrex::ParallelDescriptor::IOProcessor()) {
      int i = 0;
      mass = foo[i++];
      mom[0] = foo[i++];
      mom[1] = foo[i++];
      mom[2] = foo[i++];
      rho_e = foo[i++];
      rho_K = foo[i++];
      rho_E = foo[i++];
      fuel_prod = foo[i++];
      temp = foo[i++];

      amrex::Print() << '\n';
      amrex::Print() << "TIME = " << time << " MASS        = " << mass << '\n';
      amrex::Print() << "TIME = " << time << " XMOM        = " << mom[0]
                     << '\n';
      amrex::Print() << "TIME = " << time << " YMOM        = " << mom[1]
                     << '\n';
      amrex::Print() << "TIME = " << time << " ZMOM        = " << mom[2]
                     << '\n';
      amrex::Print() << "TIME = " << time << " RHO*e       = " << rho_e << '\n';
      amrex::Print() << "TIME = " << time << " RHO*K       = " << rho_K << '\n';
      amrex::Print() << "TIME = " << time << " RHO*E       = " << rho_E << '\n';
      amrex::Print() << "TIME = " << time << " FUEL PROD   = " << fuel_prod
                     << '\n';

      const int log_index = find_datalog_index("datlog");
      if (log_index >= 0) {
        std::ostream& data_log1 = parent->DataLog(log_index);
        if (data_log1.good()) {
          const int datwidth = 24;
          if (time == 0.0) {
            data_log1 << std::setw(datwidth) << "          time";
            data_log1 << std::setw(datwidth) << "          mass";
            data_log1 << std::setw(datwidth) << "          xmom";
            data_log1 << std::setw(datwidth) << "          ymom";
            data_log1 << std::setw(datwidth) << "          zmom";
            data_log1 << std::setw(datwidth) << "         rho_K";
            data_log1 << std::setw(datwidth) << "         rho_e";
            data_log1 << std::setw(datwidth) << "         rho_E";
            data_log1 << std::setw(datwidth) << "     fuel_prod";
            data_log1 << std::setw(datwidth) << "          temp";
            data_log1 << std::endl;
          }

          // Write the quantities at this time
          const int datprecision = 16;
          data_log1 << std::setw(datwidth) << time;
          data_log1 << std::setw(datwidth) << std::setprecision(datprecision)
                    << mass;
          data_log1 << std::setw(datwidth) << std::setprecision(datprecision)
                    << mom[0];
          data_log1 << std::setw(datwidth) << std::setprecision(datprecision)
                    << mom[1];
          data_log1 << std::setw(datwidth) << std::setprecision(datprecision)
                    << mom[2];
          data_log1 << std::setw(datwidth) << std::setprecision(datprecision)
                    << rho_K;
          data_log1 << std::setw(datwidth) << std::setprecision(datprecision)
                    << rho_e;
          data_log1 << std::setw(datwidth) << std::setprecision(datprecision)
                    << rho_E;
          data_log1 << std::setw(datwidth) << std::setprecision(datprecision)
                    << fuel_prod;
          data_log1 << std::setw(datwidth) << std::setprecision(datprecision)
                    << temp;
          data_log1 << std::endl;
        }
      }
    }
  }
}

void
PeleC::monitor_extrema()
{
  BL_PROFILE("PeleC::monitor_extrema()");

  if (verbose <= 0) {
    return;
  }

  const bool local_flag = true;
  const int finest_level = parent->finestLevel();
  const amrex::Real time = state[State_Type].curTime();
  amrex::Vector<std::string> extrema_vars = {
    "density", "x_velocity", "y_velocity", "z_velocity", "eint_e",
    "Temp",    "pressure",   "massfrac",   "sumYminus1"};

  int nspec_extrema = 2;
  bool use_all_spec = false;
  if (extrema_spec_name == "ALL") {
    use_all_spec = true;
    nspec_extrema += NUM_SPECIES;
    extrema_vars.insert(
      extrema_vars.end(), PeleC::spec_names.begin(), PeleC::spec_names.end());
  } else {
    if (!fuel_name.empty()) {
      nspec_extrema++;
      extrema_vars.push_back(fuel_name);
    }
    if (!flame_trac_name.empty()) {
      nspec_extrema++;
      extrema_vars.push_back(flame_trac_name);
    }
    if (!extrema_spec_name.empty()) {
      nspec_extrema++;
      extrema_vars.push_back(extrema_spec_name);
    }
  }

  const auto nextrema = static_cast<int>(extrema_vars.size());
  constexpr amrex::Real neg_huge = std::numeric_limits<amrex::Real>::lowest();
  constexpr amrex::Real huge = std::numeric_limits<amrex::Real>::max();
  amrex::Vector<amrex::Real> minima(nextrema, huge), maxima(nextrema, neg_huge);

  for (int lev = 0; lev <= finest_level; lev++) {
    PeleC& pc_lev = getLevel(lev);
    for (int ii = 0; ii < nextrema - nspec_extrema; ++ii) {
      maxima[ii] = amrex::max<amrex::Real>(
        maxima[ii], pc_lev.maxDerive(extrema_vars[ii], time, local_flag));
      minima[ii] = amrex::min<amrex::Real>(
        minima[ii], pc_lev.minDerive(extrema_vars[ii], time, local_flag));
    }

    // Handle species separately
    auto mf = pc_lev.derive("massfrac", time, 0);
    BL_ASSERT(!(mf == nullptr));
    amrex::MultiFab sumY(pc_lev.grids, pc_lev.dmap, 1, 0);
    sumY.setVal(0.0);
    const auto idx_massfrac = nextrema - nspec_extrema;

    for (int ispec = 0; ispec < NUM_SPECIES; ispec++) {
      amrex::Real maxval = mf->max(ispec, 0, local_flag);
      amrex::Real minval = mf->min(ispec, 0, local_flag);
      amrex::MultiFab::Add(sumY, *mf, ispec, 0, 1, 0);

      // "massfrac" gets the extrema across all mass fractions
      maxima[idx_massfrac] =
        amrex::max<amrex::Real>(maxima[idx_massfrac], maxval);
      minima[idx_massfrac] =
        amrex::min<amrex::Real>(minima[idx_massfrac], minval);

      // Get values for any individual species if relevant
      if (use_all_spec) {
        const auto idx_spec = idx_massfrac + ispec + 2;
        maxima[idx_spec] = amrex::max<amrex::Real>(maxima[idx_spec], maxval);
        minima[idx_spec] = amrex::min<amrex::Real>(minima[idx_spec], minval);
      } else {
        for (int iext = 0; iext < nspec_extrema - 2; iext++) {
          const auto idx_spec = idx_massfrac + iext + 2;
          if (extrema_vars[idx_spec] == PeleC::spec_names[ispec]) {
            maxima[idx_spec] =
              amrex::max<amrex::Real>(maxima[idx_spec], maxval);
            minima[idx_spec] =
              amrex::min<amrex::Real>(minima[idx_spec], minval);
          }
        }
      }
    }
    // sum of mass fractions
    maxima[idx_massfrac + 1] = amrex::max<amrex::Real>(
      maxima[idx_massfrac + 1], sumY.max(0, 0, local_flag) - 1.0);
    minima[idx_massfrac + 1] = amrex::min<amrex::Real>(
      minima[idx_massfrac + 1], sumY.min(0, 0, local_flag) - 1.0);
  }

  if (verbose > 0) {
    amrex::ParallelDescriptor::ReduceRealMax(
      maxima.data(), nextrema, amrex::ParallelDescriptor::IOProcessorNumber());

    amrex::ParallelDescriptor::ReduceRealMin(
      minima.data(), nextrema, amrex::ParallelDescriptor::IOProcessorNumber());

    if (amrex::ParallelDescriptor::IOProcessor()) {

      amrex::Print() << std::endl;
      for (int ii = 0; ii < nextrema; ++ii) {
        const int datwidth = 15;
        const int datwidth_txt = 10;
        const int datprecision = 8;
        amrex::Print() << "TIME = " << time << " " << std::left
                       << std::setw(datwidth_txt) << extrema_vars[ii]
                       << "  MIN = " << std::setw(datwidth)
                       << std::setprecision(datprecision) << minima[ii]
                       << "  MAX = " << std::setw(datwidth)
                       << std::setprecision(datprecision) << maxima[ii]
                       << std::endl;
      }

      const int log_index = find_datalog_index("extremalog");
      if (log_index >= 0) {
        std::ostream& data_log1 = parent->DataLog(log_index);
        if (data_log1.good()) {
          const int datwidth = 18;
          if (time == 0.0) {
            data_log1 << std::setw(datwidth) << "          time";
            for (int ii = 0; ii < nextrema; ++ii) {
              data_log1 << std::setw(datwidth - 4) << extrema_vars[ii]
                        << "-min";
              data_log1 << std::setw(datwidth - 4) << extrema_vars[ii]
                        << "-max";
            }
            data_log1 << std::endl;
          }

          // Write the quantities at this time
          const int datprecision = 10;
          data_log1 << std::setw(datwidth) << time;
          for (int ii = 0; ii < nextrema; ++ii) {
            data_log1 << std::setw(datwidth) << std::setprecision(datprecision)
                      << minima[ii];
            data_log1 << std::setw(datwidth) << std::setprecision(datprecision)
                      << maxima[ii];
          }
          data_log1 << std::endl;
        }
      }
    }
  }
}
