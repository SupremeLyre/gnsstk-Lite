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

#ifndef PPP_STATE_VECTOR_HPP
#define PPP_STATE_VECTOR_HPP

#include <map>
#include <set>
#include <string>
#include <vector>

#include "Matrix.hpp"
#include "Namelist.hpp"
#include "SatID.hpp"
#include "SatelliteSystem.hpp"
#include "Vector.hpp"

#include "PPPConstants.hpp"

namespace gnsstk
{
namespace ppp
{

/// Forward declaration
class PPPStateVector;

/**
 * @class PPPStateVector
 * @brief Manages the PPP state vector layout for both IF and UC modes.
 *
 * The state vector is dynamically sized based on the number of satellites
 * and systems. It provides index mapping from physical parameters to
 * vector indices, and manages the covariance matrix and Namelist labels.
 *
 * IF Mode layout:
 *   [dX, dY, dZ][dT_gps][ISB_sys2...ISB_sysN][ZWD][N_IF_sat1...N_IF_satM]
 *
 * UC Mode layout:
 *   [dX, dY, dZ][dT_sys1...dT_sysN][ZWD][I1_sat1...I1_satM][N1_sat1...N1_satM][N2_sat1...N2_satM]
 */
class PPPStateVector
{
  public:
    /// Default constructor
    PPPStateVector();

    /**
     * @brief Constructor with mode and system list.
     * @param mode PPP processing mode (IF or UC)
     * @param systems List of satellite systems to process
     * @param refSys Reference system for clock (default GPS)
     */
    PPPStateVector(PPPMode mode, const std::vector<SatelliteSystem>& systems,
                   SatelliteSystem refSys = SatelliteSystem::GPS);

    /// Destructor
    ~PPPStateVector() = default;

    // -------------------------------------------------------------------------
    // Setup and configuration
    // -------------------------------------------------------------------------

    /**
     * @brief (Re-)initialize the state vector with given mode and systems.
     * Clears all existing satellites.
     */
    void initialize(PPPMode mode, const std::vector<SatelliteSystem>& systems,
                    SatelliteSystem refSys = SatelliteSystem::GPS);

    /**
     * @brief Set the list of active satellites.
     * This will resize the state vector and covariance, preserving existing
     * states where possible.
     */
    void setSatellites(const std::vector<SatID>& sats);

    /**
     * @brief Add a single satellite to the state vector.
     * @return true if the satellite was newly added
     */
    bool addSatellite(const SatID& sat);

    /**
     * @brief Remove a satellite from the state vector.
     * @return true if the satellite was found and removed
     */
    bool removeSatellite(const SatID& sat);

    /**
     * @brief Remove all satellites that are not in the given set.
     * @return Number of satellites removed
     */
    size_t pruneSatellites(const std::set<SatID>& activeSats);

    // -------------------------------------------------------------------------
    // Index queries
    // -------------------------------------------------------------------------

    /// @return Index of dX (always 0)
    size_t posIndex() const
    {
        return 0;
    }

    /// @return Index of clock for given system (IF: reference system first; UC: all systems)
    size_t clkIndex(SatelliteSystem sys) const;

    /// @return Index of ZWD
    size_t zwdIndex() const;

    /// @return Index of ISB for given system (IF mode only; sys must not be reference)
    size_t isbIndex(SatelliteSystem sys) const;

    /// @return Index of IF ambiguity for given satellite (IF mode only)
    size_t ambIFIndex(const SatID& sat) const;

    /// @return Index of L1 ionosphere for given satellite (UC mode only)
    size_t ionoIndex(const SatID& sat) const;

    /// @return Index of freq-1 ambiguity for given satellite (UC mode only)
    size_t amb1Index(const SatID& sat) const;

    /// @return Index of freq-2 ambiguity for given satellite (UC mode only)
    size_t amb2Index(const SatID& sat) const;

    /// @return true if the satellite is present in the state vector
    bool hasSatellite(const SatID& sat) const;

    // -------------------------------------------------------------------------
    // Dimension and mode
    // -------------------------------------------------------------------------

    /// @return Total dimension of the state vector
    size_t size() const
    {
        return numFixedStates() + satellites.size() * numSatStates();
    }

    /// @return Number of satellite systems
    size_t numSystems() const
    {
        return systems.size();
    }

    /// @return Number of satellites
    size_t numSatellites() const
    {
        return satellites.size();
    }

    /// @return Current processing mode
    PPPMode getMode() const
    {
        return mode;
    }

    /// @return Number of fixed states (position + clocks + ZWD, no ISB)
    size_t numFixedStates() const;

    /// @return Number of satellite-dependent states per satellite
    size_t numSatStates() const;

    // -------------------------------------------------------------------------
    // Data access
    // -------------------------------------------------------------------------

    /// State vector
    Vector<double> state;

    /// Covariance matrix
    Matrix<double> cov;

    /// Namelist labels for each state element
    Namelist names;

    // -------------------------------------------------------------------------
    // Helpers
    // -------------------------------------------------------------------------

    /// Build the Namelist based on current state layout
    void buildNamelist();

    /// Resize state and covariance to current dimension, preserving existing values
    void resize();

    /// @return String representation of state layout
    std::string dump() const;

  private:
    PPPMode mode;                           ///< Processing mode
    SatelliteSystem refSystem;              ///< Reference clock system
    std::vector<SatelliteSystem> systems;   ///< Active systems (reference first in IF mode)
    std::vector<SatID> satellites;          ///< Active satellites in deterministic order

    // Caches for quick lookup
    std::map<SatelliteSystem, size_t> sysIndexMap;   ///< System -> clock/ISB index
    std::map<SatID, size_t> satIndexMap;             ///< Sat -> satellite block index

    // Offsets into the state vector
    size_t offsetClocks;    ///< Start of clock/ISB block
    size_t offsetZWD;       ///< Start of ZWD
    size_t offsetSatStates; ///< Start of satellite-dependent states

    /// Rebuild index caches after satellite/system list changes
    void rebuildIndices();

    /// Sort systems: reference system first, then others in deterministic order
    void sortSystems();
};

} // namespace ppp
} // namespace gnsstk

#endif // PPP_STATE_VECTOR_HPP
