////////////////////////////////////////////////////////////////////////
// \file SpaceCharge.h
//
// \brief pure virtual base interface for space charge distortions
//
// \author mrmooney@bnl.gov
//
////////////////////////////////////////////////////////////////////////
#ifndef SPACECHARGE_SPACECHARGE_H
#define SPACECHARGE_SPACECHARGE_H

// C/C++ standard libraries
#include "larcoreobj/SimpleTypesAndConstants/geo_vectors.h"

// Forward declarations
namespace geo {
  class TPCID;
}

namespace spacecharge {

  class SpaceCharge {
  public:
    SpaceCharge(const SpaceCharge&) = delete;
    SpaceCharge(SpaceCharge&&) = delete;
    SpaceCharge& operator=(const SpaceCharge&) = delete;
    SpaceCharge& operator=(SpaceCharge&&) = delete;
    virtual ~SpaceCharge() = default;

    virtual bool EnableSimSpatialSCE() const = 0;
    virtual bool EnableSimEfieldSCE() const = 0;
    virtual bool EnableCorrSCE() const = 0;
    virtual bool EnableCalSpatialSCE() const = 0;
    virtual bool EnableCalEfieldSCE() const = 0;

    /// GetPosOffset(true point) + true point = apparent position after drift
    /// All values are in global coordinates, so method is aware of drift direction
    /// Note that this method must calculate the TPD ID, so is potentially expensive
    virtual geo::Vector_t GetPosOffsets(geo::Point_t const& point) const = 0;
   
    /// GetPosOffset(true point) + true point = apparent position after drift
    /// Preferred method when TPC ID is known
    /// All values are in global coordinates, so method is aware of drift direction
    virtual geo::Vector_t GetPosOffsets(geo::Point_t const& point, geo::TPCID const& tpcID) const = 0;
   
    /// Nom Efield vector + mag(Efield) * GetEfieldOffsets(true point) = Actual Efield at true point
    /// All values are in global coordinates, so method is aware of drift direction
    /// Note that this method must calculate the TPD ID, so is potentially expensive
    virtual geo::Vector_t GetEfieldOffsets(geo::Point_t const& point) const = 0;
    
    /// Nom Efield vector + mag(Efield) * GetEfieldOffsets(true point) = Actual Efield at true point
    /// Preferred method when the TPC ID is known
    /// All values are in global coordinates, so method is aware of drift direction
    virtual geo::Vector_t GetEfieldOffsets(geo::Point_t const& point, geo::TPCID const& tpcID) const = 0;

    /// GetCalPosOffsets(measured point, TPC ID) + measured point = estimated true position
    /// All values are in global coorinates, so method is aware of drift direction for given TPC ID
    virtual geo::Vector_t GetCalPosOffsets(geo::Point_t const& point, geo::TPCID const& tpcID) const = 0;
    
    // To get corrected field at measured point, first calculate estimated true position, 
    // then get Efield offsets for that position
    // Preferred method is to get estimated true position, then call GetEfieldOffsets
    virtual geo::Vector_t GetCalEfieldOffsets(geo::Point_t const& point,
                                              geo::TPCID const& tpcID) const = 0;

  protected:
    SpaceCharge() = default;

  }; // class SpaceCharge
} //namespace spacecharge

#endif // SPACECHARGE_SPACECHARGE_H
