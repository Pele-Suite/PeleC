#include "Geometry.H"

namespace pele::pelec {

void
FlatPlate::build(
  const amrex::Geometry& /*geom*/, const int /*max_coarsening_level*/)
{
  amrex::Print() << "flat plate  geometry not currently supported. \n";
  amrex::Abort();
}

void
Ramp::build(const amrex::Geometry& geom, const int max_coarsening_level)
{
  amrex::Print() << "ramp geometry\n";
  int upDir;
  int indepVar;
  amrex::Real startPt;
  amrex::Real slope;
  amrex::ParmParse pp("eb2");
  pp.get("up_dir", upDir);
  pp.get("indep_var", indepVar);
  pp.get("start_pt", startPt);
  pp.get("ramp_slope", slope);

  amrex::RealArray normal;
  normal.fill(0.0);
  normal[upDir] = 1.0;
  normal[indepVar] = -slope;

  amrex::RealArray point;
  point.fill(0.0);
  point[upDir] = -slope * startPt;

  amrex::EB2::PlaneIF ramp(point, normal);
  auto gshop = amrex::EB2::makeShop(ramp);
  amrex::EB2::Build(gshop, geom, max_coarsening_level, max_coarsening_level);
}

void
Combustor::build(const amrex::Geometry& geom, const int max_coarsening_level)
{
  amrex::ParmParse pp("combustor");

  amrex::Real fwl;
  pp.get("far_wall_loc", fwl);

  amrex::EB2::PlaneIF farwall(
    {AMREX_D_DECL(fwl, 0., 0.)}, {AMREX_D_DECL(1., 0., 0.)});

  amrex::Vector<amrex::Real> pl1pt;
  amrex::Vector<amrex::Real> pl2pt;
  amrex::Vector<amrex::Real> pl2nm;
  amrex::Vector<amrex::Real> pl3pt;
  pp.getarr("ramp_plane1_point", pl1pt);
  pp.getarr("ramp_plane2_point", pl2pt);
  pp.getarr("ramp_plane2_normal", pl2nm);
  pp.getarr("ramp_plane3_point", pl3pt);

  amrex::EB2::PlaneIF r0(
    {AMREX_D_DECL(pl1pt[0], pl1pt[1], 0.)}, {AMREX_D_DECL(0., -1., 0.)});
  amrex::EB2::PlaneIF r1(
    {AMREX_D_DECL(pl2pt[0], pl2pt[1], 0.)},
    {AMREX_D_DECL(pl2nm[0], pl2nm[1], 0.)});
  amrex::EB2::PlaneIF r2(
    {AMREX_D_DECL(pl3pt[0], pl3pt[1], 0.)}, {AMREX_D_DECL(1., 0., 0.)});
  auto ramp = amrex::EB2::makeIntersection(r0, r1, r2);

  amrex::Vector<amrex::Real> pipelo;
  amrex::Vector<amrex::Real> pipehi;
  pp.getarr("pipe_lo", pipelo);
  pp.getarr("pipe_hi", pipehi);

  amrex::EB2::BoxIF pipe(
    {AMREX_D_DECL(pipelo[0], pipelo[1], -1.)},
    {AMREX_D_DECL(pipehi[0], pipehi[1], 1.)}, false);

  // where does plane 1 and plane 2 intersect?
  amrex::Real k2 = std::abs(pl2nm[0] / pl2nm[1]);
  amrex::Real secty = pl2pt[1] + k2 * (pl3pt[0] - pl2pt[0]);
  // How much do we cut?
  amrex::Real dx = geom.CellSize(0);
  amrex::Real dycut =
    4. * (1. + max_coarsening_level) * amrex::min<amrex::Real>(dx, k2 * dx);
  amrex::EB2::BoxIF flat_corner(
    {AMREX_D_DECL(pl3pt[0], 0., -1.)}, {AMREX_D_DECL(1.e10, secty + dycut, 1.)},
    false);

  auto polys = amrex::EB2::makeUnion(farwall, ramp, pipe, flat_corner);

  amrex::Real lenx = amrex::DefaultGeometry().ProbLength(0);
  amrex::Real leny = amrex::DefaultGeometry().ProbLength(1);
  auto pr = amrex::EB2::translate(
    amrex::EB2::lathe(polys), {AMREX_D_DECL(
                                static_cast<amrex::Real>(lenx * 0.5),
                                static_cast<amrex::Real>(leny * 0.5), 0.)});

  auto gshop = amrex::EB2::makeShop(pr);
  amrex::EB2::Build(gshop, geom, max_coarsening_level, max_coarsening_level);
}

void
ICEPistonBowl::build(
  const amrex::Geometry& geom, const int max_coarsening_level)
{
  amrex::RealArray center({AMREX_D_DECL(0.04 - 0.0125 - 0.02, 0.0, 0.0)});

  const amrex::Real radius = 0.02;

  const bool has_fluid_inside = false;
  amrex::EB2::SphereIF sf(radius, center, has_fluid_inside);

  amrex::EB2::CylinderIF cf1(
    0.04, 0.09, 0, {AMREX_D_DECL(0.045, 0.0, 0.0)}, true);
  amrex::EB2::CylinderIF cf2(
    0.04, 0.125, 0, {AMREX_D_DECL(-0.0125 - 0.02, 0.0, 0.0)}, false);
  amrex::EB2::CylinderIF cf3(
    0.03, 0.125, 0, {AMREX_D_DECL(-0.0125 - 0.02, 0.0, 0.0)}, false);
  auto pipe = amrex::EB2::makeDifference(cf2, cf3);
  amrex::EB2::CylinderIF cf4(
    0.03, 0.10, 0, {AMREX_D_DECL(-0.0125 - 0.02, 0.0, 0.0)}, false);

  amrex::RealArray center2({AMREX_D_DECL(0.0, 0.0, 0.0)});

  const amrex::Real radius2 = 0.09;

  const bool has_fluid_inside2 = true;
  amrex::EB2::SphereIF sf2(radius2, center2, has_fluid_inside2);

  auto polys = amrex::EB2::makeUnion(cf1, pipe, cf4, sf, sf2);

  auto gshop = amrex::EB2::makeShop(polys);
  amrex::EB2::Build(gshop, geom, max_coarsening_level, max_coarsening_level, 4);
}
void
ExtrudedTriangles::build(
  const amrex::Geometry& geom, const int max_coarsening_level)
{
  // setting some constants
  // the polygon is triangle
  // we can only do a maximum of 5 triangles (change if needed)
  const int npts_in_tri = 3;
  const int max_tri = 5;

  // number of user defined triangles
  int num_tri;

  amrex::ParmParse pp("extruded_triangles");
  amrex::Vector<amrex::Array<amrex::Real, AMREX_SPACEDIM>> alltri(
    npts_in_tri * max_tri);

  // initialize all triangles with some dummy values
  // that fall outside of the domain
  const amrex::Real* problo;
  const amrex::Real* probhi;
  amrex::Real maxlen;

  problo = geom.ProbLo();
  probhi = geom.ProbHi();

  AMREX_D_TERM(
    maxlen = geom.ProbLength(0);
    , maxlen = amrex::max<amrex::Real>(maxlen, geom.ProbLength(1));
    , maxlen = amrex::max<amrex::Real>(maxlen, geom.ProbLength(2)););

  // setting all triangles to be waaay outside the domain initially
  for (int itri = 0; itri < max_tri; itri++) {
    AMREX_D_TERM(
      alltri[npts_in_tri * itri + 0][0] = problo[0] + 100.0 * maxlen;
      , alltri[npts_in_tri * itri + 0][1] = problo[1] - 100.0 * maxlen;
      , alltri[npts_in_tri * itri + 0][2] = 0.0;);

    AMREX_D_TERM(
      alltri[npts_in_tri * itri + 1][0] = probhi[0] + 100.0 * maxlen;
      , alltri[npts_in_tri * itri + 1][1] = problo[1] - 100.0 * maxlen;
      , alltri[npts_in_tri * itri + 1][2] = 0.0;);

    AMREX_D_TERM(
      alltri[npts_in_tri * itri + 2][0] = probhi[0] + 100.0 * maxlen;
      , alltri[npts_in_tri * itri + 2][1] = problo[1] + 100.0 * maxlen;
      , alltri[npts_in_tri * itri + 2][2] = 0.0;);
  }

  // get user defined number of triangles
  pp.get("num_tri", num_tri);

  for (int itri = 0; itri < num_tri; itri++) {
    amrex::Array<amrex::Real, AMREX_SPACEDIM> point{
      AMREX_D_DECL(0.0, 0.0, 0.0)};

    for (int ipt = 0; ipt < npts_in_tri; ipt++) {
      std::string pointstr =
        "tri_" + convertIntGG(itri) + "_point_" + convertIntGG(ipt);
      amrex::Vector<amrex::Real> vecpt;
      pp.getarr(pointstr, vecpt, 0, AMREX_SPACEDIM);
      for (int dir = 0; dir < AMREX_SPACEDIM; dir++) {
        point[dir] = vecpt[dir];
      }
      alltri[npts_in_tri * itri + ipt] = point;
    }
  }

  // intersection of the 3 planes in a triangle for all triangles
  amrex::Vector<std::unique_ptr<amrex::EB2::IntersectionIF<
    amrex::EB2::PlaneIF, amrex::EB2::PlaneIF, amrex::EB2::PlaneIF>>>
    impfunc_triangles(max_tri);

  for (int itri = 0; itri < max_tri; itri++) {
    // make sure points are in anti clockwise direction to set the inside of
    // the triangle as solid phase correctly
    amrex::Array<amrex::Real, AMREX_SPACEDIM> norm0;
    amrex::Array<amrex::Real, AMREX_SPACEDIM> norm1;
    amrex::Array<amrex::Real, AMREX_SPACEDIM> norm2;

    amrex::Array<amrex::Real, AMREX_SPACEDIM> point0;
    amrex::Array<amrex::Real, AMREX_SPACEDIM> point1;
    amrex::Array<amrex::Real, AMREX_SPACEDIM> point2;

    point0 = alltri[npts_in_tri * itri + 0];
    point1 = alltri[npts_in_tri * itri + 1];
    point2 = alltri[npts_in_tri * itri + 2];

    AMREX_D_TERM(
      norm0[0] = -(point1[1] - point0[1]);, norm0[1] = (point1[0] - point0[0]);
      , norm0[2] = 0.0;);

    AMREX_D_TERM(
      norm1[0] = -(point2[1] - point1[1]);, norm1[1] = (point2[0] - point1[0]);
      , norm1[2] = 0.0;);

    AMREX_D_TERM(
      norm2[0] = -(point0[1] - point2[1]);, norm2[1] = (point0[0] - point2[0]);
      , norm2[2] = 0.0;);

    // normalize so that magnitude is 1
    amrex::Real norm = sqrt(norm0[0] * norm0[0] + norm0[1] * norm0[1]);
    norm0[0] = norm0[0] / norm;
    norm0[1] = norm0[1] / norm;

    // normalize so that magnitude is 1
    norm = sqrt(norm1[0] * norm1[0] + norm1[1] * norm1[1]);
    norm1[0] = norm1[0] / norm;
    norm1[1] = norm1[1] / norm;

    // normalize so that magnitude is 1
    norm = sqrt(norm2[0] * norm2[0] + norm2[1] * norm2[1]);
    norm2[0] = norm2[0] / norm;
    norm2[1] = norm2[1] / norm;

    amrex::EB2::PlaneIF plane0(point0, norm0);
    amrex::EB2::PlaneIF plane1(point1, norm1);
    amrex::EB2::PlaneIF plane2(point2, norm2);

    impfunc_triangles[itri] = std::make_unique<amrex::EB2::IntersectionIF<
      amrex::EB2::PlaneIF, amrex::EB2::PlaneIF, amrex::EB2::PlaneIF>>(

      plane0, plane1, plane2);
  }

  auto alltri_IF = amrex::EB2::makeUnion(
    *impfunc_triangles[0], *impfunc_triangles[1], *impfunc_triangles[2],
    *impfunc_triangles[3], *impfunc_triangles[4]);

  auto gshop = amrex::EB2::makeShop(alltri_IF);
  amrex::EB2::Build(gshop, geom, max_coarsening_level, max_coarsening_level);
}
void
PolygonRevolution::build(
  const amrex::Geometry& /*geom*/, const int /*max_coarsening_level*/)
{
  amrex::Print() << "polygon_revolution  geometry not currently supported. "
                    " combustor?\n";
  amrex::Abort();
}

void
QuarterCircle::build(
  const amrex::Geometry& geom, const int max_coarsening_level)
{
  amrex::Real r_inner = 1.0;
  amrex::Real r_outer = 2.0;
  amrex::ParmParse pp("eb2");
  pp.query("r_inner", r_inner);
  pp.query("r_outer", r_outer);

  amrex::EB2::CylinderIF inner(r_inner, 10, 2, {AMREX_D_DECL(0, 0, 0)}, false);
  amrex::EB2::CylinderIF outer(r_outer, 10, 2, {AMREX_D_DECL(0, 0, 0)}, true);

  auto polys = amrex::EB2::makeUnion(inner, outer);
  auto gshop = amrex::EB2::makeShop(polys);
  amrex::EB2::Build(
    gshop, geom, max_coarsening_level, max_coarsening_level, 4, false);
}

void
RotatedCylinder::build(
  const amrex::Geometry& geom, const int max_coarsening_level)
{
  bool inside = true;
  amrex::Real radius = 0.0002;
  int direction = 0;
  amrex::Vector<amrex::Real> centervec(3, 0.0);
  amrex::Real rotation = 0;
  int rotation_axe = 0;

  amrex::ParmParse pp("eb2");
  pp.query("cylinder_has_fluid_inside", inside);
  pp.query("cylinder_radius", radius);
  pp.query("cylinder_direction", direction);
  pp.query("cylinder_rotation", rotation);
  pp.query("cylinder_rotation_axe", rotation_axe);
  pp.getarr("cylinder_center", centervec, 0, 3);
  amrex::Array<amrex::Real, AMREX_SPACEDIM> center = {
    AMREX_D_DECL(centervec[0], centervec[1], centervec[2])};

  rotation = (rotation / 180.) * constants::PI();

  amrex::EB2::CylinderIF my_cyl(radius, direction, center, inside);

  auto my_cyl_rot = amrex::EB2::rotate(my_cyl, rotation, rotation_axe);
  auto gshop = amrex::EB2::makeShop(my_cyl_rot);
  amrex::EB2::Build(
    gshop, geom, max_coarsening_level, max_coarsening_level, 4, false);
}

void
RotatedBox::build(const amrex::Geometry& geom, const int max_coarsening_level)
{

  amrex::ParmParse pp("eb2");
  amrex::RealArray lo;
  pp.get("box_lo", lo);

  amrex::RealArray hi;
  pp.get("box_hi", hi);

  bool has_fluid_inside;
  pp.get("box_has_fluid_inside", has_fluid_inside);

  amrex::Real rotation = 0;
  int rotation_axe = 0;
  pp.query("box_rotation", rotation);
  pp.query("box_rotation_axe", rotation_axe);

  rotation = (rotation / 180.) * constants::PI();

  amrex::EB2::BoxIF bf(lo, hi, has_fluid_inside);

  auto bf_rot = amrex::EB2::rotate(bf, rotation, rotation_axe);
  auto gshop = amrex::EB2::makeShop(bf_rot);
  amrex::EB2::Build(
    gshop, geom, max_coarsening_level, max_coarsening_level, 4, false);
}

void
CheckpointFile::build(
  const amrex::Geometry& geom, const int max_coarsening_level)
{
  std::string chkfile = "chk_geom";
  amrex::ParmParse pp("eb2");
  pp.query("chkfile", chkfile);

  amrex::EB2::BuildFromChkptFile(chkfile, geom, 0, max_coarsening_level);
}
} // namespace pele::pelec
