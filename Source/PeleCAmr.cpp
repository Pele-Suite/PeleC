#include "PeleCAmr.H"

#ifdef PELE_USE_SPRAY
#include "SprayParticles.H"
#endif

void
PeleCAmr::writePlotFile()
{
  if (!Plot_Files_Output()) {
    return;
  }
  BL_PROFILE("PeleCAmr::writePlotFile()");

  if (first_plotfile) {
    first_plotfile = false;
    amr_level[0]->setPlotVariables();
  }

  // Don't continue if we have no variables to plot.
  if (statePlotVars().empty() && derivePlotVars().empty()) {
    return;
  }

  const std::string& pltfile =
    amrex::Concatenate(plot_file_root, level_steps[0], file_name_digits);

  if (verbose > 0) {
    amrex::Print() << "PLOTFILE: file = " << pltfile << '\n';
  }

  if ((record_run_info != 0) && amrex::ParallelDescriptor::IOProcessor()) {
    runlog << "PLOTFILE: file = " << pltfile << '\n';
  }

  bool write_hdf5_plots = false;
  std::string hdf5_compression{"None@0"};
#ifdef AMREX_USE_HDF5
  amrex::ParmParse pp("pelec");
  pp.query("write_hdf5_plots", write_hdf5_plots);
#ifdef AMREX_USE_HDF5_ZFP
  pp.query("hdf5_compression", hdf5_compression);
#endif
#endif
  writePlotFileDoit(pltfile, true, write_hdf5_plots, hdf5_compression);
}

void
PeleCAmr::writeSmallPlotFile()
{
  if (!Plot_Files_Output()) {
    return;
  }

  BL_PROFILE("PeleCAmr::writeSmallPlotFile()");

  if (first_smallplotfile) {
    first_smallplotfile = false;
    amr_level[0]->setSmallPlotVariables();
  }

  // Don't continue if we have no variables to plot.
  if (stateSmallPlotVars().empty() && deriveSmallPlotVars().empty()) {
    return;
  }

  const std::string& pltfile =
    amrex::Concatenate(small_plot_file_root, level_steps[0], file_name_digits);

  if (verbose > 0) {
    amrex::Print() << "SMALL PLOTFILE: file = " << pltfile << '\n';
  }

  if ((record_run_info != 0) && amrex::ParallelDescriptor::IOProcessor()) {
    runlog << "SMALL PLOTFILE: file = " << pltfile << '\n';
  }

  writePlotFileDoit(pltfile, false);
}

void
PeleCAmr::constructPlotMF(
  const bool regular,
  amrex::Vector<std::unique_ptr<amrex::MultiFab>>& plotMFs,
  amrex::Vector<std::string>& plt_var_names)
{

  const auto& desc_lst = amrex::AmrLevel::get_desc_lst();
  amrex::Vector<std::pair<int, int>> plot_var_map;
  for (int typ = 0; typ < desc_lst.size(); typ++) {
    for (int comp = 0; comp < desc_lst[typ].nComp(); comp++) {
      const std::string& name = desc_lst[typ].name(comp);
      const bool is_plot_var = regular ? amrex::Amr::isStatePlotVar(name)
                                       : amrex::Amr::isStateSmallPlotVar(name);
      if (
        is_plot_var &&
        desc_lst[typ].getType() == amrex::IndexType::TheCellType()) {
        plot_var_map.push_back(std::pair<int, int>(typ, comp));
      }
    }
  }

  int num_derive = 0;
  auto& derive_lst = amrex::AmrLevel::get_derive_lst();
  std::list<std::string> derive_names;
  const std::list<amrex::DeriveRec>& dlist = derive_lst.dlist();
  for (const auto& it : dlist) {
    const bool is_plot_var = regular
                               ? amrex::Amr::isDerivePlotVar(it.name())
                               : amrex::Amr::isDeriveSmallPlotVar(it.name());
    if (is_plot_var) {
      derive_names.push_back(it.name());
      num_derive += it.numDerive();
    }
  }
#ifdef PELE_USE_SPRAY
  // Add spray derive variables
  if (regular && (SprayParticleContainer::NumDeriveVars() > 0)) {
    num_derive += SprayParticleContainer::NumDeriveVars();
  }
#endif

  // Decide to plot vfrac
  bool plot_vfrac =
    ebInDomain() && (amrex::EB2::TopIndexSpaceIfPresent() != nullptr);
  amrex::ParmParse pp("pelec");
  pp.query("plot_vfrac", plot_vfrac);

  const auto n_data_items =
    plot_var_map.size() + num_derive + static_cast<int>(plot_vfrac);

  const int nGrow = 0;
  const amrex::Real cur_time =
    (amr_level[0]->get_state_data(State_Type)).curTime();

  const int nlevels = finestLevel() + 1;
  for (int lev = 0; lev < nlevels; ++lev) {

    plotMFs[lev] = std::make_unique<amrex::MultiFab>(
      boxArray(lev), DistributionMap(lev), n_data_items, nGrow, amrex::MFInfo(),
      amr_level[lev]->Factory());

    // Cull data from state variables -- use no ghost cells.
    int cnt = 0;
    for (auto& pv : plot_var_map) {
      int typ = pv.first;
      int comp = pv.second;
      amrex::MultiFab& this_dat = amr_level[lev]->get_new_data(typ);
      amrex::MultiFab::Copy(*plotMFs[lev], this_dat, comp, cnt, 1, nGrow);
      cnt++;
    }

    // Cull data from derived variables.
    if ((!derive_names.empty())) {
      for (const auto& derive_name : derive_names) {
        const amrex::DeriveRec* rec = derive_lst.get(derive_name);
        int ncomp = rec->numDerive();

        auto derive_dat = amr_level[lev]->derive(derive_name, cur_time, nGrow);
        amrex::MultiFab::Copy(*plotMFs[lev], *derive_dat, 0, cnt, ncomp, nGrow);
        cnt += ncomp;
      }
    }

#ifdef PELE_USE_SPRAY
    if (SprayParticleContainer::NumDeriveVars() > 0 && regular) {
      const int num_spray_derive = SprayParticleContainer::NumDeriveVars();
      PeleC::setupVirtualParticles(lev, finestLevel());
      plotMFs[lev]->setVal(0., cnt, num_spray_derive);
      // Compute derived spray variables for active particles
      PeleC::SprayPC->computeDerivedVars(*plotMFs[lev], lev, cnt);
      if (lev < finestLevel()) {
        amrex::MultiFab tmp_plt(
          boxArray(lev), DistributionMap(lev), num_spray_derive, 0,
          amrex::MFInfo(), amr_level[lev]->Factory());
        tmp_plt.setVal(0.);
        // Compute derived spray variables for virtual particles under refined
        // regions
        PeleC::VirtPC->computeDerivedVars(tmp_plt, lev, 0);
        amrex::MultiFab::Add(
          *plotMFs[lev], tmp_plt, 0, cnt, num_spray_derive, 0);
      }
      PeleC::removeVirtualParticles(lev);
      cnt += num_spray_derive;
    }
#endif

    if (plot_vfrac) {
      const auto& ebfactory = dynamic_cast<amrex::EBFArrayBoxFactory const&>(
        amr_level[lev]->Factory());
      amrex::MultiFab::Copy(
        *plotMFs[lev], ebfactory.getVolFrac(), 0, cnt, 1, nGrow);
    }
  }

  for (auto& pv : plot_var_map) {
    int typ = pv.first;
    int comp = pv.second;
    plt_var_names.push_back(desc_lst[typ].name(comp));
  }

  for (const auto& derive_name : derive_names) {
    const amrex::DeriveRec* rec = derive_lst.get(derive_name);
    for (int i = 0; i < rec->numDerive(); i++) {
      plt_var_names.push_back(rec->variableName(i));
    }
  }

#ifdef PELE_USE_SPRAY
  if (SprayParticleContainer::NumDeriveVars() > 0 && regular) {
    for (const auto& spray_derive_name :
         SprayParticleContainer::DeriveVarNames()) {
      plt_var_names.push_back(spray_derive_name);
    }
  }
#endif

  if (plot_vfrac) {
    plt_var_names.push_back("vfrac");
  }

  AMREX_ASSERT(n_data_items == plt_var_names.size());
}

void
PeleCAmr::writePlotFileDoit(
  const std::string& pltfile,
  const bool regular,
  const bool write_hdf5_plots,
  const std::string& hdf5_compression)
{
  auto dPlotFileTime0 = amrex::second();

  const int nlevels = finestLevel() + 1;
  amrex::Vector<std::string> plt_var_names;
  amrex::Vector<std::unique_ptr<amrex::MultiFab>> plotMFs(nlevels);
  constructPlotMF(regular, plotMFs, plt_var_names);

  amrex::Vector<const amrex::MultiFab*> plotMFs_constvec;
  plotMFs_constvec.reserve(nlevels);
  for (int lev = 0; lev < nlevels; ++lev) {
    plotMFs_constvec.push_back(
      static_cast<const amrex::MultiFab*>(plotMFs[lev].get()));
  }

  const amrex::Real cur_time =
    (amr_level[0]->get_state_data(State_Type)).curTime();
  amrex::Vector<int> istep(nlevels);
  for (int lev = 0; lev < nlevels; ++lev) {
    istep[lev] = levelSteps(lev);
  }

#ifdef AMREX_USE_HDF5
  if (write_hdf5_plots) {
    amrex::WriteMultiLevelPlotfileHDF5SingleDset(
      pltfile, nlevels, plotMFs_constvec, plt_var_names, Geom(), cur_time,
      istep, refRatio(), hdf5_compression);
  } else {
#endif
    (void)hdf5_compression; // Avoid unused warning
    amrex::WriteMultiLevelPlotfile(
      pltfile, nlevels, plotMFs_constvec, plt_var_names, Geom(), cur_time,
      istep, refRatio());
#ifdef AMREX_USE_HDF5
  }
#endif

  amrex::VisMF::IO_Buffer io_buffer(amrex::VisMF::GetIOBufferSize());
  std::ofstream HeaderFile;
  if (!write_hdf5_plots) {
    HeaderFile.rdbuf()->pubsetbuf(io_buffer.dataPtr(), io_buffer.size());
    if (amrex::ParallelDescriptor::IOProcessor()) {
      // Only the IOProcessor() writes to the header file.
      std::string HeaderFileName(pltfile + "/Header");
      HeaderFile.open(HeaderFileName.c_str(), std::ios::app | std::ios::binary);
      if (!HeaderFile.good()) {
        amrex::FileOpenFailed(HeaderFileName);
      }
    }

    if (amrex::EB2::TopIndexSpaceIfPresent() != nullptr) {
      for (int lev = 0; lev < nlevels; ++lev) {
        HeaderFile << "1.0e-6\n";
      }
    }
  }

  if (regular) {
    if (!write_hdf5_plots) {
      for (int lev = 0; lev < nlevels; ++lev) {
        amr_level[lev]->writePlotFilePost(pltfile, HeaderFile);
      }
    }

    last_plotfile = level_steps[0];
  } else {
    last_smallplotfile = level_steps[0];
  }

  if ((amrex::ParallelDescriptor::IOProcessor()) && (HeaderFile.is_open())) {
    HeaderFile.close();
  }

#ifdef PELE_USE_SPRAY
  if (PeleC::SprayPC != nullptr) {
    if (regular) {
      for (int lev = 0; lev < nlevels; ++lev) {
        PeleC::SprayPC->SprayParticleIO(lev, false, pltfile);
      }
    } else {
      amrex::Abort("HDF5 particle writing incomplete");
    }
  }
#endif

  if (verbose > 0) {
    const int IOProc = amrex::ParallelDescriptor::IOProcessorNumber();
    auto dPlotFileTime = amrex::second() - dPlotFileTime0;
    amrex::ParallelDescriptor::ReduceRealMax(dPlotFileTime, IOProc);
    if (regular) {
      amrex::Print() << "Write plotfile time = " << dPlotFileTime << "  seconds"
                     << "\n\n";
    } else {
      amrex::Print() << "Write small plotfile time = " << dPlotFileTime
                     << "  seconds"
                     << "\n\n";
    }
  }
}

#ifdef AMREX_USE_ASCENT
void
PeleCAmr::doInSituViz(const int step)
{
  BL_PROFILE("PeleCAmr::doInSituViz()");

  // Output only on given frequency
  if (!(step % pele_ascent.plot_int == 0)) {
    return;
  }

  if (step == 0) {
    amr_level[0]->setPlotVariables();
  }

  if (statePlotVars().empty() && derivePlotVars().empty()) {
    return;
  }

  auto dPlotFileTime0 = amrex::second();

  const int nlevels = finestLevel() + 1;
  amrex::Vector<std::string> plt_var_names;
  amrex::Vector<std::unique_ptr<amrex::MultiFab>> plotMFs(nlevels);
  constructPlotMF(true, plotMFs, plt_var_names);

  amrex::Vector<const amrex::MultiFab*> plotMFs_constvec;
  plotMFs_constvec.reserve(nlevels);
  for (int lev = 0; lev < nlevels; ++lev) {
    plotMFs_constvec.push_back(
      static_cast<const amrex::MultiFab*>(plotMFs[lev].get()));
  }

  const amrex::Real cur_time =
    (amr_level[0]->get_state_data(State_Type)).curTime();
  amrex::Vector<int> istep(nlevels);
  for (int lev = 0; lev < nlevels; ++lev) {
    istep[lev] = levelSteps(lev);
  }

  conduit::Node bp_mesh;
  amrex::MultiLevelToBlueprint(
    nlevels, plotMFs_constvec, plt_var_names, Geom(), cur_time, istep,
    refRatio(), bp_mesh);

  ascent::Ascent ascent;
  conduit::Node open_opts;

#ifdef AMREX_USE_MPI
  open_opts["mpi_comm"] =
    MPI_Comm_c2f(amrex::ParallelDescriptor::Communicator());
#endif
  ascent.open(open_opts);
  conduit::Node verify_info;
  if (!conduit::blueprint::mesh::verify(bp_mesh, verify_info)) {
    ASCENT_INFO("Error: Mesh Blueprint Verify Failed!");
    verify_info.print();
  }

  conduit::Node actions;
  ascent.publish(bp_mesh);

  ascent.execute(actions);

  ascent.close();

  if (verbose > 0) {
    const int IOProc = amrex::ParallelDescriptor::IOProcessorNumber();
    auto dPlotFileTime = amrex::second() - dPlotFileTime0;
    amrex::ParallelDescriptor::ReduceRealMax(dPlotFileTime, IOProc);
    amrex::Print() << "Ascent write time = " << dPlotFileTime << "  seconds"
                   << std::endl;
  }
}
#endif
