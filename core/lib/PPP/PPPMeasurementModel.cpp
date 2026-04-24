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

#include "PPPMeasurementModel.hpp"

#include <cmath>

#include "GNSSconstants.hpp"
#include "SatelliteSystem.hpp"
#include "StringUtils.hpp"

namespace gnsstk
{
namespace ppp
{

using namespace std;

// --------------------------------------------------------------------------------
// Construction
// --------------------------------------------------------------------------------

PPPMeasurementModel::PPPMeasurementModel(const PPPStateVector& stateVec)
    : sv(stateVec)
{
}

// --------------------------------------------------------------------------------
// Unified build interface
// --------------------------------------------------------------------------------

bool PPPMeasurementModel::build(const vector<PPPRawObservation>& rawObs,
                                const vector<PPPGeometry>& geometry,
                                Vector<double>& obs,
                                Matrix<double>& H,
                                Matrix<double>& measCov) const
{
    if (sv.getMode() == PPPMode::IonosphereFree)
        return buildIF(rawObs, geometry, obs, H, measCov);
    else
        return buildUC(rawObs, geometry, obs, H, measCov);
}

// --------------------------------------------------------------------------------
// IF Mode builder
// --------------------------------------------------------------------------------

bool PPPMeasurementModel::buildIF(const vector<PPPRawObservation>& rawObs,
                                  const vector<PPPGeometry>& geometry,
                                  Vector<double>& obs,
                                  Matrix<double>& H,
                                  Matrix<double>& measCov) const
{
    const size_t nSatsAll = rawObs.size();
    if (nSatsAll == 0 || geometry.size() != nSatsAll)
        return false;

    // Count valid observations: need P1, P2, L1, L2 for IF combination
    size_t nValid = 0;
    for (const auto& ro : rawObs)
    {
        if (ro.P1 != 0.0 && ro.P2 != 0.0 && ro.L1 != 0.0 && ro.L2 != 0.0)
            nValid++;
    }
    if (nValid == 0)
        return false;

    const size_t nObs = nValid * 2; // P_IF and L_IF per satellite
    const size_t nState = sv.size();

    obs = Vector<double>(nObs, 0.0);
    H = Matrix<double>(nObs, nState, 0.0);
    measCov = Matrix<double>(nObs, nObs, 0.0);

    size_t obsRow = 0;
    for (size_t i = 0; i < nSatsAll; ++i)
    {
        const auto& ro = rawObs[i];
        const auto& geo = geometry[i];

        // Skip if any required observation is missing
        if (ro.P1 == 0.0 || ro.P2 == 0.0 || ro.L1 == 0.0 || ro.L2 == 0.0)
            continue;

        double f1 = C_MPS / ro.lambda1;
        double f2 = C_MPS / ro.lambda2;
        double alpha, beta;
        ionoFreeCoefficients(f1, f2, alpha, beta);

        // IF combination
        double P_IF = alpha * ro.P1 + beta * ro.P2;
        double L_IF = alpha * ro.L1 + beta * ro.L2;

        // IF wavelength in meters (for ambiguity state)
        // The IF phase combination: L_IF = alpha*L1 + beta*L2
        // where L1 = lambda1 * phi1 (meters), L2 = lambda2 * phi2 (meters)
        // The IF ambiguity in meters: N_IF_m = alpha*lambda1*N1 + beta*lambda2*N2
        // We estimate N_IF_m directly
        double lambda_IF = C_MPS / (f1 * f1 / (f1 + f2) - f2 * f2 / (f1 + f2)); // simplified
        // Actually for the design matrix, the ambiguity partial is 1.0 for L_IF
        // because L_IF = ... + N_IF_m (where N_IF_m is in meters)

        // --- Pseudorange IF ---
        // obs = P_IF - (rho - satClk + tropDry)  [pre-fit residual]
        // Note: geometry.rho includes corrections, satClk is removed
        // tropDry is the mapped dry delay (fixed)
        obs(obsRow) = P_IF - (geo.rho - geo.satClk + geo.tropDry);

        // Partials
        fillPositionPartials(H, obsRow, geo.cosines);
        // Clock/ISB partials
        size_t clkIdx = sv.clkIndex(ro.sat.system);
        if (clkIdx != string::npos)
            H(obsRow, clkIdx) = 1.0;
        // ZWD partial
        H(obsRow, sv.zwdIndex()) = geo.tropWetMap;
        // Ambiguity partial: 0 for pseudorange
        size_t ambIdx = sv.ambIFIndex(ro.sat);
        if (ambIdx != string::npos)
            H(obsRow, ambIdx) = 0.0;

        // Measurement noise
        double noiseP_IF = sqrt(alpha * alpha * PPPMeasNoise::PCODE * PPPMeasNoise::PCODE +
                                beta * beta * PPPMeasNoise::PCODE * PPPMeasNoise::PCODE);
        measCov(obsRow, obsRow) = noiseP_IF * noiseP_IF;
        obsRow++;

        // --- Carrier phase IF ---
        obs(obsRow) = L_IF - (geo.rho - geo.satClk + geo.tropDry);

        // Partials (same as pseudorange)
        fillPositionPartials(H, obsRow, geo.cosines);
        if (clkIdx != string::npos)
            H(obsRow, clkIdx) = 1.0;
        H(obsRow, sv.zwdIndex()) = geo.tropWetMap;
        // Ambiguity partial: 1.0 for carrier phase (N_IF in meters)
        if (ambIdx != string::npos)
            H(obsRow, ambIdx) = 1.0;

        // Measurement noise
        double noiseL_IF = sqrt(alpha * alpha * PPPMeasNoise::PHASE * PPPMeasNoise::PHASE +
                                beta * beta * PPPMeasNoise::PHASE * PPPMeasNoise::PHASE);
        measCov(obsRow, obsRow) = noiseL_IF * noiseL_IF;
        obsRow++;
    }

    return obsRow > 0;
}

// --------------------------------------------------------------------------------
// UC Mode builder
// --------------------------------------------------------------------------------

bool PPPMeasurementModel::buildUC(const vector<PPPRawObservation>& rawObs,
                                  const vector<PPPGeometry>& geometry,
                                  Vector<double>& obs,
                                  Matrix<double>& H,
                                  Matrix<double>& measCov) const
{
    const size_t nSatsAll = rawObs.size();
    if (nSatsAll == 0 || geometry.size() != nSatsAll)
        return false;

    // Count valid observations: at least one freq for P and L
    size_t nValid = 0;
    for (const auto& ro : rawObs)
    {
        if ((ro.P1 != 0.0 || ro.P2 != 0.0) && (ro.L1 != 0.0 || ro.L2 != 0.0))
            nValid++;
    }
    if (nValid == 0)
        return false;

    const size_t nState = sv.size();

    // Pre-count valid observations
    size_t nObs = 0;
    for (const auto& ro : rawObs)
    {
        if (ro.P1 != 0.0) nObs++;
        if (ro.P2 != 0.0) nObs++;
        if (ro.L1 != 0.0) nObs++;
        if (ro.L2 != 0.0) nObs++;
    }

    obs = Vector<double>(nObs, 0.0);
    H = Matrix<double>(nObs, nState, 0.0);
    measCov = Matrix<double>(nObs, nObs, 0.0);

    size_t obsRow = 0;
    for (size_t i = 0; i < nSatsAll; ++i)
    {
        const auto& ro = rawObs[i];
        const auto& geo = geometry[i];

        double f1 = (ro.lambda1 > 0.0) ? C_MPS / ro.lambda1 : 0.0;
        double f2 = (ro.lambda2 > 0.0) ? C_MPS / ro.lambda2 : 0.0;
        double gamma = (f1 > 0.0 && f2 > 0.0) ? (f1 / f2) * (f1 / f2) : 1.0;

        size_t clkIdx = sv.clkIndex(ro.sat.system);
        size_t ionoIdx = sv.ionoIndex(ro.sat);
        size_t amb1Idx = sv.amb1Index(ro.sat);
        size_t amb2Idx = sv.amb2Index(ro.sat);

        // --- P1 ---
        if (ro.P1 != 0.0)
        {
            obs(obsRow) = ro.P1 - (geo.rho - geo.satClk + geo.tropDry);
            fillPositionPartials(H, obsRow, geo.cosines);
            if (clkIdx != string::npos) H(obsRow, clkIdx) = 1.0;
            H(obsRow, sv.zwdIndex()) = geo.tropWetMap;
            if (ionoIdx != string::npos) H(obsRow, ionoIdx) = 1.0;       // +I1
            if (amb1Idx != string::npos) H(obsRow, amb1Idx) = 0.0;
            if (amb2Idx != string::npos) H(obsRow, amb2Idx) = 0.0;
            measCov(obsRow, obsRow) = PPPMeasNoise::PCODE * PPPMeasNoise::PCODE;
            obsRow++;
        }

        // --- P2 ---
        if (ro.P2 != 0.0)
        {
            obs(obsRow) = ro.P2 - (geo.rho - geo.satClk + geo.tropDry);
            fillPositionPartials(H, obsRow, geo.cosines);
            if (clkIdx != string::npos) H(obsRow, clkIdx) = 1.0;
            H(obsRow, sv.zwdIndex()) = geo.tropWetMap;
            if (ionoIdx != string::npos) H(obsRow, ionoIdx) = gamma;     // +gamma*I1
            if (amb1Idx != string::npos) H(obsRow, amb1Idx) = 0.0;
            if (amb2Idx != string::npos) H(obsRow, amb2Idx) = 0.0;
            measCov(obsRow, obsRow) = PPPMeasNoise::PCODE * PPPMeasNoise::PCODE;
            obsRow++;
        }

        // --- L1 ---
        if (ro.L1 != 0.0)
        {
            obs(obsRow) = ro.L1 - (geo.rho - geo.satClk + geo.tropDry);
            fillPositionPartials(H, obsRow, geo.cosines);
            if (clkIdx != string::npos) H(obsRow, clkIdx) = 1.0;
            H(obsRow, sv.zwdIndex()) = geo.tropWetMap;
            if (ionoIdx != string::npos) H(obsRow, ionoIdx) = -1.0;      // -I1
            if (amb1Idx != string::npos) H(obsRow, amb1Idx) = 1.0;       // +N1
            if (amb2Idx != string::npos) H(obsRow, amb2Idx) = 0.0;
            measCov(obsRow, obsRow) = PPPMeasNoise::PHASE * PPPMeasNoise::PHASE;
            obsRow++;
        }

        // --- L2 ---
        if (ro.L2 != 0.0)
        {
            obs(obsRow) = ro.L2 - (geo.rho - geo.satClk + geo.tropDry);
            fillPositionPartials(H, obsRow, geo.cosines);
            if (clkIdx != string::npos) H(obsRow, clkIdx) = 1.0;
            H(obsRow, sv.zwdIndex()) = geo.tropWetMap;
            if (ionoIdx != string::npos) H(obsRow, ionoIdx) = -gamma;    // -gamma*I1
            if (amb1Idx != string::npos) H(obsRow, amb1Idx) = 0.0;
            if (amb2Idx != string::npos) H(obsRow, amb2Idx) = 1.0;       // +N2
            measCov(obsRow, obsRow) = PPPMeasNoise::PHASE * PPPMeasNoise::PHASE;
            obsRow++;
        }
    }

    return obsRow > 0;
}

// --------------------------------------------------------------------------------
// Static helpers
// --------------------------------------------------------------------------------

void PPPMeasurementModel::ionoFreeCoefficients(double f1, double f2,
                                               double& alpha, double& beta)
{
    double f1sq = f1 * f1;
    double f2sq = f2 * f2;
    alpha = f1sq / (f1sq - f2sq);
    beta = -f2sq / (f1sq - f2sq);
}

size_t PPPMeasurementModel::obsPerSat() const
{
    return (sv.getMode() == PPPMode::IonosphereFree) ? 2 : 4;
}

size_t PPPMeasurementModel::expectedObsDim(size_t nSats) const
{
    return nSats * obsPerSat();
}

// --------------------------------------------------------------------------------
// Wet mapping function (simplified Niell)
// --------------------------------------------------------------------------------

double PPPMeasurementModel::wetMappingFunction(double elevDeg)
{
    if (elevDeg <= 0.0)
        return 0.0;
    double elevRad = elevDeg * DEG_TO_RAD;
    // Simplified mapping function ~ 1.001 / sqrt(0.002001 + sin^2(elev))
    double sinE = sin(elevRad);
    return 1.001 / sqrt(0.002001 + sinE * sinE);
}

// --------------------------------------------------------------------------------
// Fill position partials
// --------------------------------------------------------------------------------

void PPPMeasurementModel::fillPositionPartials(Matrix<double>& H, size_t row,
                                               const Triple& cosines)
{
    // drho/dX = (X_sat - X_rx) / rho = cosines.x (from rx to sat)
    H(row, 0) = cosines[0];
    H(row, 1) = cosines[1];
    H(row, 2) = cosines[2];
}

} // namespace ppp
} // namespace gnsstk
