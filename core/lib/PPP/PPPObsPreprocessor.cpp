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

#include "PPPObsPreprocessor.hpp"

#include <cmath>

#include "GNSSconstants.hpp"
#include "RinexObsID.hpp"
#include "SatID.hpp"

namespace gnsstk
{
namespace ppp
{

using namespace std;

// --------------------------------------------------------------------------------
// Construction
// --------------------------------------------------------------------------------

PPPObsPreprocessor::PPPObsPreprocessor()
    : mode(PPPMode::IonosphereFree),
      elevMask(7.5),
      mwThreshold(3.0),
      gfThreshold(0.05),
      tdThreshold(0.05)
{
}

PPPObsPreprocessor::PPPObsPreprocessor(PPPMode m, double mask,
                                       const vector<SatelliteSystem>& sysList)
    : mode(m),
      elevMask(mask),
      systems(sysList),
      mwThreshold(3.0),
      gfThreshold(0.05),
      tdThreshold(0.05)
{
}

// --------------------------------------------------------------------------------
// Epoch processing
// --------------------------------------------------------------------------------

map<SatID, PPPRequiredObs> PPPObsPreprocessor::processEpoch(
    const Rinex3ObsHeader& header,
    const Rinex3ObsData& obsData,
    const map<SatID, double>& elevMap)
{
    clearEpochState();

    map<SatID, PPPRequiredObs> result;

    // Loop over all satellites in the observation data
    for (const auto& obsPair : obsData.obs)
    {
        const SatID& sat = obsPair.first;

        // Check if system is active
        if (!isActiveSystem(sat.system))
            continue;

        // Check elevation mask
        auto elevIt = elevMap.find(sat);
        if (elevIt == elevMap.end() || elevIt->second < elevMask)
        {
            rejectedSats.insert(sat);
            continue;
        }

        // Extract observations
        PPPRequiredObs robs = extractObservations(header, obsData, sat);

        // Check if we have enough observations for the mode
        if (mode == PPPMode::IonosphereFree)
        {
            if (!robs.hasP1 || !robs.hasP2 || !robs.hasL1 || !robs.hasL2)
            {
                rejectedSats.insert(sat);
                continue;
            }
        }
        else // UC mode
        {
            if ((!robs.hasP1 && !robs.hasP2) || (!robs.hasL1 && !robs.hasL2))
            {
                rejectedSats.insert(sat);
                continue;
            }
        }

        // Detect cycle slips
        PPPCycleSlipFlags flags = detectCycleSlips(sat, robs);
        if (flags.anySlip())
        {
            slipFlags[sat] = flags;
            // Don't reject - just flag for ambiguity reset
        }

        // Check for gross outliers
        if (isOutlier(robs))
        {
            rejectedSats.insert(sat);
            continue;
        }

        result[sat] = robs;
    }

    return result;
}

// --------------------------------------------------------------------------------
// Query state
// --------------------------------------------------------------------------------

PPPCycleSlipFlags PPPObsPreprocessor::getCycleSlipFlags(const SatID& sat) const
{
    auto it = slipFlags.find(sat);
    if (it != slipFlags.end())
        return it->second;
    return PPPCycleSlipFlags();
}

bool PPPObsPreprocessor::isRejected(const SatID& sat) const
{
    return rejectedSats.find(sat) != rejectedSats.end();
}

size_t PPPObsPreprocessor::numCycleSlips() const
{
    return slipFlags.size();
}

size_t PPPObsPreprocessor::numRejected() const
{
    return rejectedSats.size();
}

void PPPObsPreprocessor::clearEpochState()
{
    slipFlags.clear();
    rejectedSats.clear();
}

// --------------------------------------------------------------------------------
// Extract observations
// --------------------------------------------------------------------------------

PPPRequiredObs PPPObsPreprocessor::extractObservations(
    const Rinex3ObsHeader& header,
    const Rinex3ObsData& obsData,
    const SatID& sat) const
{
    PPPRequiredObs robs;
    robs.hasP1 = robs.hasP2 = robs.hasL1 = robs.hasL2 = false;
    robs.P1 = robs.P2 = robs.L1 = robs.L2 = 0.0;

    // Get system character via RinexSatID
    RinexSatID rsat(sat);
    string sysChar = string(1, rsat.systemChar());
    auto typeIt = header.mapObsTypes.find(sysChar);
    if (typeIt == header.mapObsTypes.end())
        return robs;

    const vector<RinexObsID>& obsTypes = typeIt->second;

    // Find the observation vector for this satellite
    auto dataIt = obsData.obs.find(sat);
    if (dataIt == obsData.obs.end())
        return robs;

    const vector<RinexDatum>& data = dataIt->second;
    if (data.size() != obsTypes.size())
        return robs;

    // Get wavelengths for this system
    double lambda1, lambda2;
    if (!getWavelengths(sat.system, lambda1, lambda2))
        return robs;
    robs.lambda1 = lambda1;
    robs.lambda2 = lambda2;

    // Map observation types to indices
    for (size_t i = 0; i < obsTypes.size(); ++i)
    {
        const RinexObsID& rid = obsTypes[i];
        string typeStr = rid.asString();  // e.g., "C1C", "L1X"
        if (typeStr.empty())
            continue;
        char type = typeStr[0];    // C, L, D, S
        char band = typeStr[1];    // 1, 2, 5, etc
        // char code = typeStr[2];    // C, P, W, X, etc

        if (type != 'C' && type != 'L')
            continue;

        double value = data[i].data;
        if (value == 0.0)
            continue;

        // Try to map to primary frequencies
        // For GPS/GAL/QZSS: band '1' = L1/E1, band '2' = L2/E5a
        // For BDS: band '2' = B1I, band '6' = B3I
        // For GLONASS: band '1' = G1, band '2' = G2

        bool isPrimaryFreq = false;
        bool isSecondaryFreq = false;

        switch (sat.system)
        {
            case SatelliteSystem::GPS:
            case SatelliteSystem::Galileo:
            case SatelliteSystem::QZSS:
                isPrimaryFreq = (band == '1');
                isSecondaryFreq = (band == '2' || band == '5');
                break;
            case SatelliteSystem::BeiDou:
                isPrimaryFreq = (band == '2');  // B1I
                isSecondaryFreq = (band == '6' || band == '7'); // B3I/B2a
                break;
            case SatelliteSystem::Glonass:
                isPrimaryFreq = (band == '1');
                isSecondaryFreq = (band == '2');
                break;
            case SatelliteSystem::IRNSS:
                isPrimaryFreq = (band == '5');
                isSecondaryFreq = (band == '9');
                break;
            default:
                break;
        }

        if (type == 'C')
        {
            if (isPrimaryFreq && !robs.hasP1)
            {
                robs.P1 = value;
                robs.hasP1 = true;
            }
            else if (isSecondaryFreq && !robs.hasP2)
            {
                robs.P2 = value;
                robs.hasP2 = true;
            }
        }
        else if (type == 'L')
        {
            if (isPrimaryFreq && !robs.hasL1)
            {
                robs.L1 = value * lambda1;  // Convert cycles to meters
                robs.hasL1 = true;
            }
            else if (isSecondaryFreq && !robs.hasL2)
            {
                robs.L2 = value * lambda2;  // Convert cycles to meters
                robs.hasL2 = true;
            }
        }
    }

    return robs;
}

// --------------------------------------------------------------------------------
// Cycle slip detection
// --------------------------------------------------------------------------------

PPPCycleSlipFlags PPPObsPreprocessor::detectCycleSlips(const SatID& sat,
                                                        const PPPRequiredObs& obs)
{
    PPPCycleSlipFlags flags;
    flags.mwSlip = flags.gfSlip = flags.tdSlip = false;

    if (!obs.hasL1 || !obs.hasL2 || !obs.hasP1 || !obs.hasP2)
        return flags; // Can't detect slips without dual-frequency

    auto it = history.find(sat);
    if (it == history.end())
    {
        // First epoch for this satellite - initialize history, no slip
        EpochHistory eh;
        eh.time = CommonTime::BEGINNING_OF_TIME;
        eh.mw = computeMW(obs);
        eh.gf = computeGF(obs);
        eh.L1 = obs.L1;
        eh.L2 = obs.L2;
        eh.P1 = obs.P1;
        eh.P2 = obs.P2;
        history[sat] = eh;
        return flags;
    }

    EpochHistory& prev = it->second;

    // Melbourne-Wubbena slip detection
    double mw = computeMW(obs);
    if (fabs(mw - prev.mw) > mwThreshold)
        flags.mwSlip = true;

    // Geometry-free slip detection
    double gf = computeGF(obs);
    if (fabs(gf - prev.gf) > gfThreshold)
        flags.gfSlip = true;

    // Time-difference check: detect phase jump if L1 or L2 changed by > threshold
    // (simplified, no dt available here; use absolute difference)
    if (obs.hasL1 && prev.L1 != 0.0)
    {
        if (fabs(obs.L1 - prev.L1) > tdThreshold)
            flags.tdSlip = true;
    }
    if (obs.hasL2 && prev.L2 != 0.0)
    {
        if (fabs(obs.L2 - prev.L2) > tdThreshold)
            flags.tdSlip = true;
    }

    // Update history
    prev.mw = mw;
    prev.gf = gf;
    prev.L1 = obs.L1;
    prev.L2 = obs.L2;
    prev.P1 = obs.P1;
    prev.P2 = obs.P2;

    return flags;
}

// --------------------------------------------------------------------------------
// Static combination computations
// --------------------------------------------------------------------------------

double PPPObsPreprocessor::computeMW(const PPPRequiredObs& obs)
{
    if (!obs.hasL1 || !obs.hasL2 || !obs.hasP1 || !obs.hasP2)
        return 0.0;

    double f1 = C_MPS / obs.lambda1;
    double f2 = C_MPS / obs.lambda2;

    // Wide-lane phase (m)
    double L_w = (f1 * obs.L1 - f2 * obs.L2) / (f1 - f2);
    // Narrow-lane code (m)
    double P_n = (f1 * obs.P1 + f2 * obs.P2) / (f1 + f2);

    return L_w - P_n;
}

double PPPObsPreprocessor::computeGF(const PPPRequiredObs& obs)
{
    if (!obs.hasL1 || !obs.hasL2)
        return 0.0;

    // Geometry-free combination in meters
    return obs.L1 - obs.L2;
}

// --------------------------------------------------------------------------------
// Outlier detection
// --------------------------------------------------------------------------------

bool PPPObsPreprocessor::isOutlier(const PPPRequiredObs& obs) const
{
    // Check for unreasonably large pseudorange values
    if (obs.hasP1 && (obs.P1 < 1.0e6 || obs.P1 > 5.0e7))
        return true;
    if (obs.hasP2 && (obs.P2 < 1.0e6 || obs.P2 > 5.0e7))
        return true;

    // Check for large code-phase divergence (indicates bad data)
    if (obs.hasP1 && obs.hasL1)
    {
        double diff = fabs(obs.P1 - obs.L1);
        if (diff > 100.0) // 100m divergence
            return true;
    }

    return false;
}

// --------------------------------------------------------------------------------
// System helpers
// --------------------------------------------------------------------------------

bool PPPObsPreprocessor::isActiveSystem(SatelliteSystem sys) const
{
    for (const auto& s : systems)
    {
        if (s == sys)
            return true;
    }
    return false;
}

bool PPPObsPreprocessor::getWavelengths(SatelliteSystem sys, double& lambda1, double& lambda2)
{
    switch (sys)
    {
        case SatelliteSystem::GPS:
        case SatelliteSystem::QZSS:
            lambda1 = C_MPS / (L1_MULT_GPS * OSC_FREQ_GPS);
            lambda2 = C_MPS / (L2_MULT_GPS * OSC_FREQ_GPS);
            return true;
        case SatelliteSystem::Galileo:
            lambda1 = C_MPS / (154.0 * OSC_FREQ_GPS);  // E1
            lambda2 = C_MPS / (115.0 * OSC_FREQ_GPS);  // E5a
            return true;
        case SatelliteSystem::Glonass:
            lambda1 = C_MPS / (160.0 * OSC_FREQ_GPS);  // G1 center freq
            lambda2 = C_MPS / (125.0 * OSC_FREQ_GPS);  // G2 center freq
            return true;
        case SatelliteSystem::BeiDou:
            lambda1 = C_MPS / (152.6 * OSC_FREQ_GPS);  // B1I ~1561.098 MHz
            lambda2 = C_MPS / (127.0 * OSC_FREQ_GPS);  // B3I ~1268.520 MHz
            return true;
        case SatelliteSystem::IRNSS:
            lambda1 = C_MPS / (115.0 * OSC_FREQ_GPS);  // L5
            lambda2 = C_MPS / (115.0 * OSC_FREQ_GPS);  // S-band (same as placeholder)
            return true;
        default:
            lambda1 = lambda2 = 0.0;
            return false;
    }
}

} // namespace ppp
} // namespace gnsstk
