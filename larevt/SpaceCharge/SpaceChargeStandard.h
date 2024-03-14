////////////////////////////////////////////////////////////////////////
// \file SpaceChargeStandard.h
//
// \brief header of class for storing/accessing space charge distortions
//
// \author mrmooney@bnl.gov
//
////////////////////////////////////////////////////////////////////////
#ifndef SPACECHARGE_SPACECHARGESTANDARD_H
#define SPACECHARGE_SPACECHARGESTANDARD_H

// LArSoft libraries
#include "larcoreobj/SimpleTypesAndConstants/geo_vectors.h"
#include "larevt/SpaceCharge/SpaceCharge.h"

// FHiCL libraries
namespace fhicl {
  class ParameterSet;
}

// ROOT includes
#include "TF1.h"
class TGraph;

// C/C++ standard libraries
#include <stdint.h>
#include <string>
#include <vector>

// Forward declarations
namespace geo {
  class TPCID;
}

namespace spacecharge {

  class SpaceChargeStandard : public SpaceCharge {

  public:
    explicit SpaceChargeStandard(fhicl::ParameterSet const& pset);
    SpaceChargeStandard(SpaceChargeStandard const&) = delete;
    virtual ~SpaceChargeStandard() = default;

    bool Configure(fhicl::ParameterSet const& pset);
    bool Update(uint64_t ts = 0);

    bool EnableSimSpatialSCE() const override;
    bool EnableSimEfieldSCE() const override;
    bool EnableCorrSCE() const override;
    bool EnableCalSpatialSCE() const override;
    bool EnableCalEfieldSCE() const override;

    /// GetPosOffset(true point) + true point = apparent position after drift
    /// All values are in global coordinates, so method is aware of drift direction
    /// Note that this method must calculate the TPD ID, so is potentially expensive
    virtual geo::Vector_t GetPosOffsets(geo::Point_t const& point) const override;
   
    /// GetPosOffset(true point) + true point = apparent position after drift
    /// Preferred method when TPC ID is known
    /// All values are in global coordinates, so method is aware of drift direction
    virtual geo::Vector_t GetPosOffsets(geo::Point_t const& point, geo::TPCID const& tpcID) const override;
   
    /// Nom Efield vector + mag(Efield) * GetEfieldOffsets(true point) = Actual Efield at true point
    /// All values are in global coordinates
    /// Note that this method must calculate the TPD ID, so is potentially expensive
    virtual geo::Vector_t GetEfieldOffsets(geo::Point_t const& point) const override;
    
    /// Nom Efield vector + mag(Efield) * GetEfieldOffsets(true point) = Actual Efield at true point
    /// Preferred method when the TPC ID is known
    /// All values are in global coordinates, so method is aware of drift direction
    virtual geo::Vector_t GetEfieldOffsets(geo::Point_t const& point, geo::TPCID const& tpcID) const override;

    /// GetCalPosOffsets(measured point, TPC ID) + measured point = estimated true position
    /// All values are in global coorinates, so method is aware of drift direction for given TPC ID
    virtual geo::Vector_t GetCalPosOffsets(geo::Point_t const& point, geo::TPCID const& tpcID) const override;
    
    // To get corrected field at measured point, first calculate estimated true position, 
    // then get Efield offsets for that position 
    // Preferred method is to get estimated true position, then call GetEfieldOffsets
    geo::Vector_t GetCalEfieldOffsets(geo::Point_t const& point, geo::TPCID const& tpcID) const override;
  
  private:
  protected:
    std::vector<double> GetPosOffsetsParametric(double xVal, double yVal, double zVal, geo::TPCID const& id) const;
    
    // Note that this example code ignores TPCID when returning position offsets
    double GetOnePosOffsetParametric(double xVal, double yVal, double zVal, std::string axis) const;
    std::vector<double> GetEfieldOffsetsParametric(double xVal, double yVal, double zVal, geo::TPCID const& tpcid) const;
    double GetOneEfieldOffsetParametric(double xVal,
                                        double yVal,
                                        double zVal,
                                        std::string axis) const;
    double TransformX(double xVal) const;
    double TransformY(double yVal) const;
    double TransformZ(double zVal) const;
    bool IsInsideBoundaries(double xVal, double yVal, double zVal) const;

    bool fEnableSimSpatialSCE;
    bool fEnableSimEfieldSCE;
    bool fEnableCalSpatialSCE;
    bool fEnableCalEfieldSCE;
    bool fEnableCorrSCE;

    std::string fRepresentationType;
    std::string fInputFilename;

    TGraph** g1_x = new TGraph*[7];
    TGraph** g2_x = new TGraph*[7];
    TGraph** g3_x = new TGraph*[7];
    TGraph** g4_x = new TGraph*[7];
    TGraph** g5_x = new TGraph*[7];

    TGraph** g1_y = new TGraph*[7];
    TGraph** g2_y = new TGraph*[7];
    TGraph** g3_y = new TGraph*[7];
    TGraph** g4_y = new TGraph*[7];
    TGraph** g5_y = new TGraph*[7];
    TGraph** g6_y = new TGraph*[7];

    TGraph** g1_z = new TGraph*[7];
    TGraph** g2_z = new TGraph*[7];
    TGraph** g3_z = new TGraph*[7];
    TGraph** g4_z = new TGraph*[7];

    TF1* f1_x = new TF1("f1_x", "pol6");
    TF1* f2_x = new TF1("f2_x", "pol6");
    TF1* f3_x = new TF1("f3_x", "pol6");
    TF1* f4_x = new TF1("f4_x", "pol6");
    TF1* f5_x = new TF1("f5_x", "pol6");
    TF1* fFinal_x = new TF1("fFinal_x", "pol4");

    TF1* f1_y = new TF1("f1_y", "pol5");
    TF1* f2_y = new TF1("f2_y", "pol5");
    TF1* f3_y = new TF1("f3_y", "pol5");
    TF1* f4_y = new TF1("f4_y", "pol5");
    TF1* f5_y = new TF1("f5_y", "pol5");
    TF1* f6_y = new TF1("f6_y", "pol5");
    TF1* fFinal_y = new TF1("fFinal_y", "pol5");

    TF1* f1_z = new TF1("f1_z", "pol4");
    TF1* f2_z = new TF1("f2_z", "pol4");
    TF1* f3_z = new TF1("f3_z", "pol4");
    TF1* f4_z = new TF1("f4_z", "pol4");
    TF1* fFinal_z = new TF1("fFinal_z", "pol3");

    TGraph** g1_Ex = new TGraph*[7];
    TGraph** g2_Ex = new TGraph*[7];
    TGraph** g3_Ex = new TGraph*[7];
    TGraph** g4_Ex = new TGraph*[7];
    TGraph** g5_Ex = new TGraph*[7];

    TGraph** g1_Ey = new TGraph*[7];
    TGraph** g2_Ey = new TGraph*[7];
    TGraph** g3_Ey = new TGraph*[7];
    TGraph** g4_Ey = new TGraph*[7];
    TGraph** g5_Ey = new TGraph*[7];
    TGraph** g6_Ey = new TGraph*[7];

    TGraph** g1_Ez = new TGraph*[7];
    TGraph** g2_Ez = new TGraph*[7];
    TGraph** g3_Ez = new TGraph*[7];
    TGraph** g4_Ez = new TGraph*[7];

    TF1* f1_Ex = new TF1("f1_Ex", "pol6");
    TF1* f2_Ex = new TF1("f2_Ex", "pol6");
    TF1* f3_Ex = new TF1("f3_Ex", "pol6");
    TF1* f4_Ex = new TF1("f4_Ex", "pol6");
    TF1* f5_Ex = new TF1("f5_Ex", "pol6");
    TF1* fFinal_Ex = new TF1("fFinal_Ex", "pol4");

    TF1* f1_Ey = new TF1("f1_Ey", "pol5");
    TF1* f2_Ey = new TF1("f2_Ey", "pol5");
    TF1* f3_Ey = new TF1("f3_Ey", "pol5");
    TF1* f4_Ey = new TF1("f4_Ey", "pol5");
    TF1* f5_Ey = new TF1("f5_Ey", "pol5");
    TF1* f6_Ey = new TF1("f6_Ey", "pol5");
    TF1* fFinal_Ey = new TF1("fFinal_Ey", "pol5");

    TF1* f1_Ez = new TF1("f1_Ez", "pol4");
    TF1* f2_Ez = new TF1("f2_Ez", "pol4");
    TF1* f3_Ez = new TF1("f3_Ez", "pol4");
    TF1* f4_Ez = new TF1("f4_Ez", "pol4");
    TF1* fFinal_Ez = new TF1("fFinal_Ez", "pol3");

  }; // class SpaceChargeStandard
} //namespace spacecharge
#endif // SPACECHARGE_SPACECHARGESTANDARD_H
