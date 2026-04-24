//==============================================================================
//
//  This file is part of GNSSTk, the ARL:UT GNSS Toolkit.
//
//  The GNSSTk is free software; you can redistribute it and/or modify
//  it under the terms of the GNU Lesser General Public License as published
//  by the Free Software Foundation; either version 3.0 of the License, or
//  any later version.
//
//  The GNSSTk is distributed in the hope that it will be useful,
//  but WITHOUT ANY WARRANTY; without even the implied warranty of
//  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//  GNU Lesser General Public License for more details.
//
//  You should have received a copy of the GNU Lesser General Public
//  License along with GNSSTk; if not, write to the Free Software Foundation,
//  Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110, USA
//
//  This software was developed by Applied Research Laboratories at the
//  University of Texas at Austin.
//  Copyright 2004-2022, The Board of Regents of The University of Texas System
//
//==============================================================================

#ifndef PPP_OBS_PREPROCESSOR_HPP
#define PPP_OBS_PREPROCESSOR_HPP

#include <map>
#include <set>
#include <string>
#include <vector>

#include "CommonTime.hpp"
#include "Rinex3ObsData.hpp"
#include "Rinex3ObsHeader.hpp"
#include "RinexSatID.hpp"
#include "SatID.hpp"

#include "PPPConstants.hpp"

namespace gnsstk
{
namespace ppp
{

/// Observation types needed for PPP processing
struct PPPRequiredObs
{
    bool hasP1;     ///< Frequency-1 pseudorange available
    bool hasP2;     ///< Frequency-2 pseudorange available
    bool hasL1;     ///< Frequency-1 carrier phase available
    bool hasL2;     ///< Frequency-2 carrier phase available
    double P1;      ///< L1 pseudorange (m)
    double P2;      ///< L2 pseudorange (m)
    double L1;      ///< L1 carrier phase (m)
    double L2;      ///< L2 carrier phase (m)
    double lambda1; ///< L1 wavelength (m)
    double lambda2; ///< L2 wavelength (m)
};

/// Cycle slip detection flags for a satellite
struct PPPCycleSlipFlags
{
    bool mwSlip;    ///< Melbourne-Wubbena slip detected
    bool gfSlip;    ///< Geometry-free slip detected
    bool tdSlip;    ///< Time-difference slip detected
    bool anySlip() const { return mwSlip || gfSlip || tdSlip; }
};

/**
 * @class PPPObsPreprocessor
 * @brief Preprocess RINEX observations for PPP.
 *
 * Responsibilities:
 *   - Extract required observations (P1, P2, L1, L2) per satellite per system
 *   - Detect cycle slips using MW, GF, and TD combinations
 *   - Detect gross outliers in pseudorange
 *   - Apply elevation mask
 *   - Flag invalid or slipped observations
 */
class PPPObsPreprocessor
{
  public:
    /// Default constructor
    PPPObsPreprocessor();

    /**
     * @brief Constructor with configuration.
     * @param mode      Processing mode (IF or UC)
     * @param elevMask  Minimum elevation angle (degrees)
     * @param systems   Systems to process
     */
    PPPObsPreprocessor(PPPMode mode, double elevMask,
                       const std::vector<SatelliteSystem>& systems);

    /// Destructor
    ~PPPObsPreprocessor() = default;

    // -------------------------------------------------------------------------
    // Configuration
    // -------------------------------------------------------------------------

    /// Set minimum elevation mask (degrees)
    void setElevationMask(double mask)
    {
        elevMask = mask;
    }

    /// Set cycle slip detection thresholds
    void setCycleSlipThresholds(double mwThresh, double gfThresh, double tdThresh)
    {
        mwThreshold = mwThresh;
        gfThreshold = gfThresh;
        tdThreshold = tdThresh;
    }

    /// Set processing mode
    void setMode(PPPMode m)
    {
        mode = m;
    }

    // -------------------------------------------------------------------------
    // Processing
    // -------------------------------------------------------------------------

    /**
     * @brief Process a RINEX observation epoch.
     * @param header   RINEX observation header (for obs type mapping)
     * @param obsData  RINEX observation data for this epoch
     * @param elevMap  Pre-computed elevation angles (sat -> elev deg)
     * @return Map of satellite -> extracted observations
     */
    std::map<SatID, PPPRequiredObs> processEpoch(
        const Rinex3ObsHeader& header,
        const Rinex3ObsData& obsData,
        const std::map<SatID, double>& elevMap);

    /**
     * @brief Check if a satellite has cycle slips this epoch.
     * @param sat Satellite ID
     * @return Cycle slip flags
     */
    PPPCycleSlipFlags getCycleSlipFlags(const SatID& sat) const;

    /// @return true if the satellite was rejected this epoch
    bool isRejected(const SatID& sat) const;

    /// @return Number of satellites with detected cycle slips this epoch
    size_t numCycleSlips() const;

    /// @return Number of satellites rejected this epoch
    size_t numRejected() const;

    /// Clear epoch-specific state (call at start of each epoch)
    void clearEpochState();

  private:
    PPPMode mode;
    double elevMask;
    std::vector<SatelliteSystem> systems;

    // Cycle slip detection thresholds
    double mwThreshold;   ///< MW slip threshold (m)
    double gfThreshold;   ///< GF slip threshold (m)
    double tdThreshold;   ///< TD slip threshold (m)

    // Epoch state
    std::map<SatID, PPPCycleSlipFlags> slipFlags;
    std::set<SatID> rejectedSats;

    // History for cycle slip detection (sat -> last epoch data)
    struct EpochHistory
    {
        CommonTime time;
        double mw;      ///< MW combination (m)
        double gf;      ///< GF combination (m)
        double L1;      ///< L1 phase (m)
        double L2;      ///< L2 phase (m)
        double P1;      ///< L1 code (m)
        double P2;      ///< L2 code (m)
    };
    std::map<SatID, EpochHistory> history;

    // -------------------------------------------------------------------------
    // Internal helpers
    // -------------------------------------------------------------------------

    /// Extract P1, P2, L1, L2 from RINEX data for a satellite
    PPPRequiredObs extractObservations(const Rinex3ObsHeader& header,
                                       const Rinex3ObsData& obsData,
                                       const SatID& sat) const;

    /// Detect cycle slips for a satellite
    PPPCycleSlipFlags detectCycleSlips(const SatID& sat, const PPPRequiredObs& obs);

    /// Compute Melbourne-Wubbena combination
    static double computeMW(const PPPRequiredObs& obs);

    /// Compute Geometry-Free combination (in meters)
    static double computeGF(const PPPRequiredObs& obs);

    /// Check if observation is a gross outlier
    bool isOutlier(const PPPRequiredObs& obs) const;

    /// @return true if system is in the active list
    bool isActiveSystem(SatelliteSystem sys) const;

    /// Get frequency-dependent wavelength for a system
    static bool getWavelengths(SatelliteSystem sys, double& lambda1, double& lambda2);
};

} // namespace ppp
} // namespace gnsstk

#endif // PPP_OBS_PREPROCESSOR_HPP
