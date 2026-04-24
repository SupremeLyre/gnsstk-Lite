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

#include "PPPCorrections.hpp"

#include <fstream>
#include <sstream>

#include "AntexStream.hpp"
#include "SolarSystem.hpp"
#include "StringUtils.hpp"
#include "SunEarthSatGeometry.hpp"
#include "SolidEarthTides.hpp"
#include "logstream.hpp"

namespace gnsstk
{
namespace ppp
{

using namespace std;

// --------------------------------------------------------------------------------
// Construction / Destruction
// --------------------------------------------------------------------------------

PPPCorrections::PPPCorrections()
    : solarSys(nullptr)
{
}

PPPCorrections::~PPPCorrections()
{
}

// --------------------------------------------------------------------------------
// Load correction data
// --------------------------------------------------------------------------------

bool PPPCorrections::loadAntex(const string& filename)
{
    return parseAntexFile(filename);
}

// --------------------------------------------------------------------------------
// Apply all corrections
// --------------------------------------------------------------------------------

double PPPCorrections::apply(double obs, const SatID& sat, const CommonTime& t,
                             const Position& rxPos, const Position& satPos,
                             const string& freqStr, bool isPhase) const
{
    double corr = 0.0;

    // Satellite PCO/PCV correction
    corr += getSatPCO_PCVCorr(sat, t, satPos, rxPos, freqStr);

    // Receiver PCO/PCV correction
    corr += getRxPCO_PCVCorr(t, rxPos, satPos, freqStr);

    // Site displacement (tides)
    Triple siteDisp = getSiteDisplacement(t, rxPos);
    // Project site displacement onto line-of-sight
    Triple los = rxPos - satPos;
    double range = los.mag();
    if (range > 0.0)
    {
        los = (1.0 / range) * los;
        corr += siteDisp.dot(los);
    }

    // Phase windup (carrier phase only)
    if (isPhase)
    {
        corr += getPhaseWindup(sat, t, rxPos, satPos);
    }

    return obs - corr;
}

// --------------------------------------------------------------------------------
// Site displacement
// --------------------------------------------------------------------------------

Triple PPPCorrections::getSiteDisplacement(const CommonTime& t, const Position& rxPos) const
{
    Triple disp(0.0, 0.0, 0.0);

    // Solid Earth tides
    disp = disp + getSolidEarthTide(t, rxPos);

    // Pole tides and ocean loading require additional data files
    // and are left as placeholders for future enhancement.

    return disp;
}

// --------------------------------------------------------------------------------
// Satellite PCO/PCV
// --------------------------------------------------------------------------------

double PPPCorrections::getSatPCO_PCVCorr(const SatID& sat, const CommonTime& t,
                                          const Position& satPos, const Position& rxPos,
                                          const string& freqStr) const
{
    if (!hasAntex())
        return 0.0;

    const AntexData* ant = findSatAntex(sat);
    if (ant == nullptr)
        return 0.0;

    // Compute line-of-sight from satellite to receiver
    Triple los = rxPos - satPos;
    double range = los.mag();
    if (range <= 0.0)
        return 0.0;
    los = los.unitVector();

    // PCO+PCV correction along LOS (placeholder)
    // TODO: Implement full PCO/PCV computation using AntexData
    (void)t; (void)freqStr; (void)ant;
    return 0.0;
}

// --------------------------------------------------------------------------------
// Receiver PCO/PCV
// --------------------------------------------------------------------------------

double PPPCorrections::getRxPCO_PCVCorr(const CommonTime& t,
                                          const Position& rxPos, const Position& satPos,
                                          const string& freqStr) const
{
    if (!hasAntex() || rxAntennaType.empty())
        return 0.0;

    const AntexData* ant = findRxAntex();
    if (ant == nullptr)
        return 0.0;

    // Compute line-of-sight from receiver to satellite
    Triple los = satPos - rxPos;
    double range = los.mag();
    if (range <= 0.0)
        return 0.0;
    los = los.unitVector();

    // PCO+PCV correction along LOS (placeholder)
    // TODO: Implement full PCO/PCV computation using AntexData
    (void)t; (void)freqStr; (void)ant;
    return 0.0;
}

// --------------------------------------------------------------------------------
// Solid Earth tides
// --------------------------------------------------------------------------------

Triple PPPCorrections::getSolidEarthTide(const CommonTime& t, const Position& rxPos) const
{
    if (solarSys == nullptr)
        return Triple(0.0, 0.0, 0.0);

    try
    {
        EphTime et(t);
        return solarSys->computeSolidEarthTides(rxPos, et);
    }
    catch (Exception& e)
    {
        ostringstream oss;
        oss << "Solid earth tide computation failed: " << e.what();
        lastError.clear();
        lastError.append(oss.str());
        return Triple(0.0, 0.0, 0.0);
    }
}

// --------------------------------------------------------------------------------
// Phase windup
// --------------------------------------------------------------------------------

double PPPCorrections::getPhaseWindup(const SatID& sat, const CommonTime& t,
                                      const Position& rxPos, const Position& satPos,
                                      bool isLHC) const
{
    if (solarSys == nullptr)
        return 0.0;

    try
    {
        // Placeholder: full phase windup requires receiver dipole model
        // and detailed satellite attitude computation.
        (void)sat; (void)t; (void)rxPos; (void)satPos; (void)isLHC;
        return 0.0;
    }
    catch (Exception& e)
    {
        ostringstream oss;
        oss << "Phase windup computation failed: " << e.what();
        lastError.clear();
        lastError.append(oss.str());
        return 0.0;
    }
}

// --------------------------------------------------------------------------------
// Internal helpers
// --------------------------------------------------------------------------------

bool PPPCorrections::parseAntexFile(const string& filename)
{
    try
    {
        AntexStream antStream(filename.c_str());
        if (!antStream.is_open())
        {
            lastError = "Failed to open ANTEX file: " + filename;
            return false;
        }

        AntexData antData;
        while (antStream >> antData)
        {
            string key;
            // AntexData stores receiver vs satellite flag; fields differ by version
            // Use a generic key based on what's available
            key = satToAntexKey(SatID(1, SatelliteSystem::GPS)); // placeholder
            antexData[key] = antData;
        }

        LOG(INFO) << "Loaded " << antexData.size() << " antenna records from " << filename;
        return true;
    }
    catch (Exception& e)
    {
        ostringstream oss;
        oss << "ANTEX parsing failed: " << e.what();
        lastError.clear();
        lastError.append(oss.str());
        return false;
    }
}

const AntexData* PPPCorrections::findSatAntex(const SatID& sat) const
{
    string key = satToAntexKey(sat);
    auto it = antexData.find(key);
    if (it != antexData.end())
        return &(it->second);
    return nullptr;
}

const AntexData* PPPCorrections::findRxAntex() const
{
    if (rxAntennaType.empty())
        return nullptr;
    auto it = antexData.find(rxAntennaType);
    if (it != antexData.end())
        return &(it->second);
    return nullptr;
}

string PPPCorrections::satToAntexKey(const SatID& sat)
{
    ostringstream oss;
    oss << convertSatelliteSystemToString(sat.system);
    if (sat.id < 10)
        oss << "0";
    oss << sat.id;
    return oss.str();
}

} // namespace ppp
} // namespace gnsstk
