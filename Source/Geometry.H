#ifndef GEOMETRY_H
#define GEOMETRY_H
#include "Factory.H"
#include "EB.H"
#include "PeleCUtilities.H"

namespace pele::pelec {

class Geometry : public physics::Factory<Geometry>
{
public:
  static std::string base_identifier() { return "Geometry"; }
  virtual void
  build(const amrex::Geometry& geom, const int max_coarsening_level) = 0;
  ~Geometry() override = default;
};

class FlatPlate : public Geometry::Register<FlatPlate>
{
public:
  static std::string identifier() { return "flat_plate"; }

  void
  build(const amrex::Geometry& geom, const int max_coarsening_level) override;
};

class Ramp : public Geometry::Register<Ramp>
{
public:
  static std::string identifier() { return "ramp"; }

  void
  build(const amrex::Geometry& geom, const int max_coarsening_level) override;
};

class Combustor : public Geometry::Register<Combustor>
{
public:
  static std::string identifier() { return "combustor"; }

  void
  build(const amrex::Geometry& geom, const int max_coarsening_level) override;
};

class ICEPistonBowl : public Geometry::Register<ICEPistonBowl>
{
public:
  static std::string identifier() { return "ICE_PistonBowl"; }

  void
  build(const amrex::Geometry& geom, const int max_coarsening_level) override;
};

class ExtrudedTriangles : public Geometry::Register<ExtrudedTriangles>
{
public:
  static std::string identifier() { return "extruded_triangles"; }

  void
  build(const amrex::Geometry& geom, const int max_coarsening_level) override;
};

class PolygonRevolution : public Geometry::Register<PolygonRevolution>
{
public:
  static std::string identifier() { return "polygon_revolution"; }

  void
  build(const amrex::Geometry& geom, const int max_coarsening_level) override;
};

class QuarterCircle : public Geometry::Register<QuarterCircle>
{
public:
  static std::string identifier() { return "quarter-circle"; }

  void
  build(const amrex::Geometry& geom, const int max_coarsening_level) override;
};

class RotatedCylinder : public Geometry::Register<RotatedCylinder>
{
public:
  static std::string identifier() { return "rotated_cylinder"; }

  void
  build(const amrex::Geometry& geom, const int max_coarsening_level) override;
};

class RotatedBox : public Geometry::Register<RotatedBox>
{
public:
  static std::string identifier() { return "rotated_box"; }

  void
  build(const amrex::Geometry& geom, const int max_coarsening_level) override;
};

class CheckpointFile : public Geometry::Register<CheckpointFile>
{
public:
  static std::string identifier() { return "chkfile"; }

  void
  build(const amrex::Geometry& geom, const int max_coarsening_level) override;
};

} // namespace pele::pelec
#endif
