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

#ifndef PPP_CORRECTIONS_HPP
#define PPP_CORRECTIONS_HPP

#include <map>
#include <string>

#include "AntexData.hpp"
#include "CommonTime.hpp"
#include "Position.hpp"
#include "SatID.hpp"
#include "Triple.hpp"

// Forward declarations for optional corrections
namespace gnsstk
{
class SolarSystem;
} // namespace gnsstk

namespace gnsstk
{
namespace ppp
{

/**
 * @class PPPCorrections
 * @brief Manages and applies PPP-specific measurement corrections.
 *
 * This class wraps the existing GNSSTk correction facilities into a unified
 * interface suitable for PPP processing. Corrections include:
 *   - Satellite antenna PCO/PCV (from ANTEX)
 *   - Receiver antenna PCO/PCV (from ANTEX)
 *   - Solid Earth tides
 *   - Pole tides (optional)
 *   - Ocean loading tides (optional)
 *   - Phase windup (carrier phase only)
 *
 * Usage:
 *   PPPCorrections corr;
 *   corr.loadAntex("igs14.atx");
 *   double corrected = corr.apply(obs, sat, time, rxPos, antType, isPhase);
 */
class PPPCorrections
{
  public:
    /// Default constructor
    PPPCorrections();

    /// Destructor
    ~PPPCorrections();

    // -------------------------------------------------------------------------
    // Load correction data
    // -------------------------------------------------------------------------

    /**
     * @brief Load an ANTEX file for PCO/PCV corrections.
     * @param filename Path to ANTEX file (e.g., "igs14.atx")
     * @return true if loaded successfully
     */
    bool loadAntex(const std::string& filename);

    /**
     * @brief Set receiver antenna type for PCO/PCV lookup.
     * @param antType Antenna type string from RINEX header (e.g., "TRM59800.00")
     */
    void setReceiverAntenna(const std::string& antType)
    {
        rxAntennaType = antType;
    }

    /**
     * @brief Set SolarSystem object for tide and geometry computations.
     * If not set, tide corrections that require Sun/Moon positions will be skipped.
     */
    void setSolarSystem(SolarSystem* ss)
    {
        solarSys = ss;
    }

    // -------------------------------------------------------------------------
    // Apply corrections
    // -------------------------------------------------------------------------

    /**
     * @brief Apply all available corrections to an observation.
     * @param obs      Raw observation (m)
     * @param sat      Satellite ID
     * @param t        Observation epoch
     * @param rxPos    Receiver position (ECEF)
     * @param satPos   Satellite position (ECEF)
     * @param freqStr  Frequency string for PCV lookup (e.g., "G01", "E05")
     * @param isPhase  true for carrier phase, false for pseudorange
     * @return Corrected observation (m)
     */
    double apply(double obs, const SatID& sat, const CommonTime& t,
                 const Position& rxPos, const Position& satPos,
                 const std::string& freqStr, bool isPhase) const;

    /**
     * @brief Compute the sum of all geometric/site displacement corrections.
     * @param t        Epoch
     * @param rxPos    Receiver position (ECEF)
     * @return Total site displacement vector (ECEF XYZ, meters)
     */
    Triple getSiteDisplacement(const CommonTime& t, const Position& rxPos) const;

    // -------------------------------------------------------------------------
    // Individual corrections (for diagnostics/output)
    // -------------------------------------------------------------------------

    /**
     * @brief Get satellite PCO + PCV correction (m).
     * @param sat      Satellite ID
     * @param t        Epoch
     * @param satPos   Satellite position (ECEF)
     * @param rxPos    Receiver position (ECEF)
     * @param freqStr  Frequency string
     * @return PCO+PCV correction along line-of-sight (m)
     */
    double getSatPCO_PCVCorr(const SatID& sat, const CommonTime& t,
                             const Position& satPos, const Position& rxPos,
                             const std::string& freqStr) const;

    /**
     * @brief Get receiver PCO + PCV correction (m).
     * @param t        Epoch
     * @param rxPos    Receiver position (ECEF)
     * @param satPos   Satellite position (ECEF)
     * @param freqStr  Frequency string
     * @return PCO+PCV correction along line-of-sight (m)
     */
    double getRxPCO_PCVCorr(const CommonTime& t,
                            const Position& rxPos, const Position& satPos,
                            const std::string& freqStr) const;

    /**
     * @brief Get solid Earth tide displacement at receiver (m).
     * @param t     Epoch
     * @param rxPos Receiver position (ECEF)
     * @return Displacement vector (ECEF XYZ, meters)
     */
    Triple getSolidEarthTide(const CommonTime& t, const Position& rxPos) const;

    /**
     * @brief Get phase windup correction (m).
     * @param sat    Satellite ID
     * @param t      Epoch
     * @param rxPos  Receiver position (ECEF)
     * @param satPos Satellite position (ECEF)
     * @param isLHC  true for Left Hand Circular polarization (e.g., Galileo E5)
     * @return Phase windup in meters (to be ADDED to carrier phase)
     */
    double getPhaseWindup(const SatID& sat, const CommonTime& t,
                          const Position& rxPos, const Position& satPos,
                          bool isLHC = false) const;

    /// @return true if ANTEX data is loaded
    bool hasAntex() const
    {
        return !antexData.empty();
    }

    /// @return Last error message
    std::string getLastError() const
    {
        return lastError;
    }

  private:
    // ANTEX data: key = satellite PRN string (e.g., "G01") or antenna type
    mutable std::map<std::string, AntexData> antexData;

    std::string rxAntennaType;   ///< Receiver antenna type from RINEX
    SolarSystem* solarSys;       ///< Optional SolarSystem for tides/geometry
    mutable std::string lastError; ///< Last error message (mutable for const methods)

    /// Parse ANTEX file and populate antexData
    bool parseAntexFile(const std::string& filename);

    /// Find AntexData for a satellite by PRN
    const AntexData* findSatAntex(const SatID& sat) const;

    /// Find AntexData for receiver by antenna type
    const AntexData* findRxAntex() const;

    /// Convert SatID to ANTEX key string (e.g., "G01")
    static std::string satToAntexKey(const SatID& sat);
};

} // namespace ppp
} // namespace gnsstk

#endif // PPP_CORRECTIONS_HPP
