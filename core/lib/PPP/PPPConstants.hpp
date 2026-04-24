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

#ifndef PPP_CONSTANTS_HPP
#define PPP_CONSTANTS_HPP

#include <vector>
#include <string>
#include "SatelliteSystem.hpp"

namespace gnsstk
{
namespace ppp
{

/// PPP processing mode
enum class PPPMode
{
    IonosphereFree,   ///< Ionosphere-free combination (IF)
    Uncombined        ///< Uncombined (UC), estimate ionosphere per satellite
};

/// Default process noise sigmas for EKF (m or m/sqrt(s) as appropriate)
struct PPPProcessNoise
{
    /// Position random constant (effectively zero process noise)
    static constexpr double POS = 0.0;
    /// Receiver clock white noise (m per sqrt(s)) ~ 1e9 m/s * 1e-10 s = 0.1m/sqrt(s)
    static constexpr double CLK = 1.0e3;
    /// ZWD random walk (m per sqrt(s)) ~ 1-5 mm/sqrt(h) = 0.5-2.5e-3 mm/sqrt(s)
    static constexpr double ZWD = 5.0e-4;
    /// Ionosphere random walk (m per sqrt(s)) ~ 1-10 mm/sqrt(h)
    static constexpr double IONO = 1.0e-3;
    /// Ambiguity random constant (zero process noise)
    static constexpr double AMB = 0.0;
    /// ISB random walk (m per sqrt(s))
    static constexpr double ISB = 1.0e2;
};

/// Default initial standard deviations (m) for apriori state
struct PPPInitStd
{
    static constexpr double POS = 100.0;      ///< Position (m)
    static constexpr double CLK = 1.0e4;      ///< Clock (m)
    static constexpr double ZWD = 0.5;        ///< ZWD (m)
    static constexpr double IONO = 5.0;       ///< Ionosphere (m)
    static constexpr double AMB = 100.0;      ///< Ambiguity (m)
    static constexpr double ISB = 1.0e3;      ///< ISB (m)
};

/// Default measurement noise sigmas (m)
struct PPPMeasNoise
{
    static constexpr double PCODE = 0.3;      ///< Pseudocode (m)
    static constexpr double PHASE = 0.003;    ///< Carrier phase (m)
    static constexpr double IF_PC = 0.3;      ///< IF pseudorange (m)
    static constexpr double IF_LC = 0.003;    ///< IF carrier phase (m)
};

/// Supported GNSS systems for PPP (default set)
inline std::vector<SatelliteSystem> defaultPPPSystems()
{
    return {
        SatelliteSystem::GPS,
        SatelliteSystem::Galileo,
        SatelliteSystem::BeiDou,
        SatelliteSystem::Glonass,
        SatelliteSystem::QZSS,
        SatelliteSystem::IRNSS
    };
}

/// Reference system for clock estimation (usually GPS)
inline SatelliteSystem referenceSystem()
{
    return SatelliteSystem::GPS;
}

} // namespace ppp
} // namespace gnsstk

#endif // PPP_CONSTANTS_HPP
