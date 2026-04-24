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

#include "PPPStateVector.hpp"

#include <algorithm>
#include <iomanip>
#include <sstream>

#include "StringUtils.hpp"

namespace gnsstk
{
namespace ppp
{

using namespace std;

// --------------------------------------------------------------------------------
// Construction / Destruction
// --------------------------------------------------------------------------------

PPPStateVector::PPPStateVector()
    : mode(PPPMode::IonosphereFree),
      refSystem(SatelliteSystem::GPS),
      offsetClocks(0),
      offsetZWD(0),
      offsetSatStates(0)
{
}

PPPStateVector::PPPStateVector(PPPMode m, const vector<SatelliteSystem>& sysList,
                               SatelliteSystem refSys)
    : mode(m),
      refSystem(refSys),
      offsetClocks(0),
      offsetZWD(0),
      offsetSatStates(0)
{
    initialize(m, sysList, refSys);
}

// --------------------------------------------------------------------------------
// Initialization
// --------------------------------------------------------------------------------

void PPPStateVector::initialize(PPPMode m, const vector<SatelliteSystem>& sysList,
                                SatelliteSystem refSys)
{
    mode = m;
    refSystem = refSys;
    systems = sysList;
    satellites.clear();
    sysIndexMap.clear();
    satIndexMap.clear();

    sortSystems();
    rebuildIndices();

    state = Vector<double>(size(), 0.0);
    cov = Matrix<double>(size(), size(), 0.0);
    buildNamelist();
}

// --------------------------------------------------------------------------------
// Satellite management
// --------------------------------------------------------------------------------

void PPPStateVector::setSatellites(const vector<SatID>& sats)
{
    // Preserve existing states for satellites that remain
    Vector<double> oldState = state;
    Matrix<double> oldCov = cov;
    vector<SatID> oldSats = satellites;

    satellites = sats;
    // Sort to ensure deterministic order
    sort(satellites.begin(), satellites.end());

    rebuildIndices();
    resize();

    // Copy back old values where possible
    for (size_t i = 0; i < oldSats.size(); ++i)
    {
        auto it = satIndexMap.find(oldSats[i]);
        if (it != satIndexMap.end())
        {
            size_t newIdx = it->second;
            size_t oldIdx = i;
            size_t nSatState = numSatStates();

            // Copy satellite-dependent states
            for (size_t j = 0; j < nSatState; ++j)
            {
                state(offsetSatStates + newIdx * nSatState + j) =
                    oldState(offsetSatStates + oldIdx * nSatState + j);

                for (size_t k = 0; k < nSatState; ++k)
                {
                    cov(offsetSatStates + newIdx * nSatState + j,
                        offsetSatStates + newIdx * nSatState + k) =
                        oldCov(offsetSatStates + oldIdx * nSatState + j,
                               offsetSatStates + oldIdx * nSatState + k);
                }
            }
        }
    }
}

bool PPPStateVector::addSatellite(const SatID& sat)
{
    if (hasSatellite(sat))
        return false;

    vector<SatID> newSats = satellites;
    newSats.push_back(sat);
    setSatellites(newSats);
    return true;
}

bool PPPStateVector::removeSatellite(const SatID& sat)
{
    auto it = find(satellites.begin(), satellites.end(), sat);
    if (it == satellites.end())
        return false;

    vector<SatID> newSats = satellites;
    newSats.erase(it);
    setSatellites(newSats);
    return true;
}

size_t PPPStateVector::pruneSatellites(const set<SatID>& activeSats)
{
    size_t removed = 0;
    vector<SatID> newSats;
    for (const auto& sat : satellites)
    {
        if (activeSats.find(sat) != activeSats.end())
            newSats.push_back(sat);
        else
            ++removed;
    }
    if (removed > 0)
        setSatellites(newSats);
    return removed;
}

// --------------------------------------------------------------------------------
// Index queries
// --------------------------------------------------------------------------------

size_t PPPStateVector::clkIndex(SatelliteSystem sys) const
{
    if (mode == PPPMode::IonosphereFree)
    {
        if (sys == refSystem)
            return offsetClocks; // Reference clock
        auto it = sysIndexMap.find(sys);
        if (it != sysIndexMap.end())
            return offsetClocks + it->second; // ISB offset
        return string::npos;
    }
    else // Uncombined
    {
        auto it = sysIndexMap.find(sys);
        if (it != sysIndexMap.end())
            return offsetClocks + it->second;
        return string::npos;
    }
}

size_t PPPStateVector::zwdIndex() const
{
    return offsetZWD;
}

size_t PPPStateVector::isbIndex(SatelliteSystem sys) const
{
    if (mode != PPPMode::IonosphereFree)
        return string::npos;
    if (sys == refSystem)
        return string::npos; // Not an ISB
    auto it = sysIndexMap.find(sys);
    if (it != sysIndexMap.end())
        return offsetClocks + it->second;
    return string::npos;
}

size_t PPPStateVector::ambIFIndex(const SatID& sat) const
{
    if (mode != PPPMode::IonosphereFree)
        return string::npos;
    auto it = satIndexMap.find(sat);
    if (it != satIndexMap.end())
        return offsetSatStates + it->second; // 1 state per sat in IF mode
    return string::npos;
}

size_t PPPStateVector::ionoIndex(const SatID& sat) const
{
    if (mode != PPPMode::Uncombined)
        return string::npos;
    auto it = satIndexMap.find(sat);
    if (it != satIndexMap.end())
        return offsetSatStates + it->second * 3; // I1 is first of 3 states
    return string::npos;
}

size_t PPPStateVector::amb1Index(const SatID& sat) const
{
    if (mode != PPPMode::Uncombined)
        return string::npos;
    auto it = satIndexMap.find(sat);
    if (it != satIndexMap.end())
        return offsetSatStates + it->second * 3 + 1; // N1 is second
    return string::npos;
}

size_t PPPStateVector::amb2Index(const SatID& sat) const
{
    if (mode != PPPMode::Uncombined)
        return string::npos;
    auto it = satIndexMap.find(sat);
    if (it != satIndexMap.end())
        return offsetSatStates + it->second * 3 + 2; // N2 is third
    return string::npos;
}

bool PPPStateVector::hasSatellite(const SatID& sat) const
{
    return satIndexMap.find(sat) != satIndexMap.end();
}

// --------------------------------------------------------------------------------
// Dimension queries
// --------------------------------------------------------------------------------

size_t PPPStateVector::numFixedStates() const
{
    // Position (3) + Clock block + ZWD (1)
    if (mode == PPPMode::IonosphereFree)
        return 3 + systems.size() + 1; // dT_ref + ISB for each non-ref system
    else
        return 3 + systems.size() + 1; // dT for each system
}

size_t PPPStateVector::numSatStates() const
{
    if (mode == PPPMode::IonosphereFree)
        return 1; // N_IF per satellite
    else
        return 3; // I1, N1, N2 per satellite
}

// --------------------------------------------------------------------------------
// Internal helpers
// --------------------------------------------------------------------------------

void PPPStateVector::sortSystems()
{
    // Put reference system first, then sort the rest
    vector<SatelliteSystem> sorted;
    sorted.push_back(refSystem);
    for (const auto& sys : systems)
    {
        if (sys != refSystem)
            sorted.push_back(sys);
    }
    systems = sorted;
}

void PPPStateVector::rebuildIndices()
{
    sysIndexMap.clear();
    satIndexMap.clear();

    // Build system index map
    if (mode == PPPMode::IonosphereFree)
    {
        // In IF mode: index 0 is reference clock, rest are ISBs
        size_t idx = 0;
        for (const auto& sys : systems)
        {
            if (sys == refSystem)
                sysIndexMap[sys] = 0;
            else
                sysIndexMap[sys] = ++idx;
        }
    }
    else
    {
        // In UC mode: each system has its own clock
        for (size_t i = 0; i < systems.size(); ++i)
            sysIndexMap[systems[i]] = i;
    }

    // Build satellite index map
    for (size_t i = 0; i < satellites.size(); ++i)
        satIndexMap[satellites[i]] = i;

    // Compute offsets
    offsetClocks = 3; // After position
    if (mode == PPPMode::IonosphereFree)
        offsetZWD = offsetClocks + systems.size(); // dT_ref + ISBs
    else
        offsetZWD = offsetClocks + systems.size(); // dT per system
    offsetSatStates = offsetZWD + 1;
}

void PPPStateVector::resize()
{
    size_t dim = size();
    state.resize(dim);
    cov.resize(dim, dim);
}

void PPPStateVector::buildNamelist()
{
    names.clear();

    // Position
    names += "dX";
    names += "dY";
    names += "dZ";

    // Clocks / ISBs
    if (mode == PPPMode::IonosphereFree)
    {
        names += "dT_" + StringUtils::asString(refSystem);
        for (const auto& sys : systems)
        {
            if (sys != refSystem)
                names += "ISB_" + StringUtils::asString(sys);
        }
    }
    else
    {
        for (const auto& sys : systems)
            names += "dT_" + StringUtils::asString(sys);
    }

    // ZWD
    names += "ZWD";

    // Satellite-dependent states
    for (const auto& sat : satellites)
    {
        ostringstream satOss;
        satOss << convertSatelliteSystemToString(sat.system) << sat.id;
        string satStr = satOss.str();
        if (mode == PPPMode::IonosphereFree)
        {
            names += "N_IF_" + satStr;
        }
        else
        {
            names += "I1_" + satStr;
            names += "N1_" + satStr;
            names += "N2_" + satStr;
        }
    }
}

// --------------------------------------------------------------------------------
// String representation
// --------------------------------------------------------------------------------

string PPPStateVector::dump() const
{
    ostringstream oss;
    oss << "PPPStateVector: mode=" << (mode == PPPMode::IonosphereFree ? "IF" : "UC")
        << ", dim=" << size() << ", systems=" << systems.size() << ", sats=" << satellites.size() << endl;
    oss << "  Fixed states: pos(3) + clocks/ISBs(" << (mode == PPPMode::IonosphereFree ? systems.size() : systems.size())
        << ") + ZWD(1) = " << numFixedStates() << endl;
    oss << "  Sat states per sat: " << numSatStates() << ", total sat states: "
        << satellites.size() * numSatStates() << endl;
    oss << "  Layout: [dX,dY,dZ][clocks/ISBs][ZWD][sat-dependent]" << endl;
    oss << "  Offsets: clocks=" << offsetClocks << ", zwd=" << offsetZWD
        << ", satStates=" << offsetSatStates << endl;

    for (size_t i = 0; i < names.size(); ++i)
    {
        oss << "    [" << setw(3) << i << "] " << left << setw(15) << names.getName(i)
            << " = " << fixed << setprecision(4) << state(i) << endl;
    }
    return oss.str();
}

} // namespace ppp
} // namespace gnsstk
