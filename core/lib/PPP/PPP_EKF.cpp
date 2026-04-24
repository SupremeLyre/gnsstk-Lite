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

#include "PPP_EKF.hpp"

#include <cmath>
#include <algorithm>

#include "GNSSconstants.hpp"
#include "SatelliteSystem.hpp"

namespace gnsstk
{
namespace ppp
{

using namespace std;

// --------------------------------------------------------------------------------
// Construction
// --------------------------------------------------------------------------------

PPP_EKF::PPP_EKF(PPPMode m, const vector<SatelliteSystem>& systems,
                 const Position& refPos, SatelliteSystem refSys)
    : mode(m),
      nominalPos(refPos),
      refSystem(refSys),
      epochDT(30.0),
      stateResized(false)
{
    sv.initialize(m, systems, refSys);
    // Set Namelist in KalmanFilter base class
    Reset(sv.names);
}

// --------------------------------------------------------------------------------
// Data input
// --------------------------------------------------------------------------------

void PPP_EKF::setEpochData(double dt, const vector<PPPRawObservation>& obs,
                           const vector<PPPGeometry>& geo)
{
    epochDT = dt;
    currentObs = obs;
    currentGeo = geo;
}

// --------------------------------------------------------------------------------
// Result queries
// --------------------------------------------------------------------------------

Position PPP_EKF::getPosition() const
{
    Position pos(nominalPos);
    if (sv.size() >= 3)
    {
        pos.setECEF(pos.X() + sv.state(sv.posIndex()),
                    pos.Y() + sv.state(sv.posIndex() + 1),
                    pos.Z() + sv.state(sv.posIndex() + 2));
    }
    return pos;
}

Matrix<double> PPP_EKF::getPositionCov() const
{
    Matrix<double> cov3(3, 3);
    for (size_t i = 0; i < 3; ++i)
        for (size_t j = 0; j < 3; ++j)
            cov3(i, j) = sv.cov(i, j);
    return cov3;
}

double PPP_EKF::getZWD() const
{
    return sv.state(sv.zwdIndex());
}

double PPP_EKF::getClock(SatelliteSystem sys) const
{
    if (mode == PPPMode::IonosphereFree)
    {
        if (sys == refSystem)
            return sv.state(sv.clkIndex(refSystem));
        // For non-ref system, clock = ref_clock + ISB
        size_t isbIdx = sv.isbIndex(sys);
        if (isbIdx != string::npos)
            return sv.state(sv.clkIndex(refSystem)) + sv.state(isbIdx);
        return 0.0;
    }
    else
    {
        size_t idx = sv.clkIndex(sys);
        if (idx != string::npos)
            return sv.state(idx);
        return 0.0;
    }
}

double PPP_EKF::getISB(SatelliteSystem sys) const
{
    if (mode != PPPMode::IonosphereFree)
        return 0.0;
    size_t idx = sv.isbIndex(sys);
    if (idx != string::npos)
        return sv.state(idx);
    return 0.0;
}

// --------------------------------------------------------------------------------
// defineInitial - Initialize filter with apriori
// --------------------------------------------------------------------------------

int PPP_EKF::defineInitial(double& T0, Vector<double>& X, Matrix<double>& Cov)
{
    T0 = 0.0;

    // State vector initially all zeros (deviations from nominal)
    X = Vector<double>(sv.size(), 0.0);

    // Apriori covariance: diagonal with default standard deviations
    Cov = Matrix<double>(sv.size(), sv.size(), 0.0);

    // Position
    Cov(0, 0) = PPPInitStd::POS * PPPInitStd::POS;
    Cov(1, 1) = PPPInitStd::POS * PPPInitStd::POS;
    Cov(2, 2) = PPPInitStd::POS * PPPInitStd::POS;

    // Clocks / ISBs
    size_t nClk = (mode == PPPMode::IonosphereFree) ? sv.numSystems() : sv.numSystems();
    for (size_t i = 0; i < nClk; ++i)
    {
        size_t idx = sv.clkIndex(sv.getMode() == PPPMode::IonosphereFree ? refSystem : refSystem);
        // Set all clock/ISB variances
        if (mode == PPPMode::IonosphereFree)
        {
            Cov(3, 3) = PPPInitStd::CLK * PPPInitStd::CLK;          // Ref clock
            for (size_t j = 1; j < sv.numSystems(); ++j)
            {
                size_t isbIdx = 3 + j;
                if (isbIdx < sv.size())
                    Cov(isbIdx, isbIdx) = PPPInitStd::ISB * PPPInitStd::ISB;
            }
        }
        else
        {
            for (size_t j = 0; j < sv.numSystems(); ++j)
            {
                size_t clkIdx = sv.clkIndex(refSystem) + j; // Simple offset
                if (clkIdx < sv.size())
                    Cov(clkIdx, clkIdx) = PPPInitStd::CLK * PPPInitStd::CLK;
            }
        }
        break; // Only need to do this once
    }

    // ZWD
    size_t zwdIdx = sv.zwdIndex();
    if (zwdIdx < sv.size())
        Cov(zwdIdx, zwdIdx) = PPPInitStd::ZWD * PPPInitStd::ZWD;

    // Satellite-dependent states (ambiguities, ionosphere)
    size_t nSatState = sv.numSatStates();
    for (size_t i = 0; i < sv.numSatellites(); ++i)
    {
        for (size_t j = 0; j < nSatState; ++j)
        {
            size_t idx = sv.zwdIndex() + 1 + i * nSatState + j;
            if (idx >= sv.size())
                break;
            if (mode == PPPMode::IonosphereFree)
                Cov(idx, idx) = PPPInitStd::AMB * PPPInitStd::AMB;
            else
            {
                // UC: I1, N1, N2
                if (j == 0)
                    Cov(idx, idx) = PPPInitStd::IONO * PPPInitStd::IONO;
                else
                    Cov(idx, idx) = PPPInitStd::AMB * PPPInitStd::AMB;
            }
        }
    }

    return 1; // Return state and covariance (not information form)
}

// --------------------------------------------------------------------------------
// defineMeasurements - Build observation equations
// --------------------------------------------------------------------------------

KalmanFilter::KalmanReturn PPP_EKF::defineMeasurements(double& T, const Vector<double>& X,
                                                        const Matrix<double>& C, bool useFlag)
{
    // Update time
    T = time + epochDT;

    // Sync KalmanFilter state/cov to PPPStateVector for measurement building
    sv.state = State;
    sv.cov = Cov;

    // Check for satellite changes before measurement update
    if (satellitesChanged())
    {
        updateStateForNewSatellites();
    }

    // Build measurement equations
    PPPMeasurementModel mm(sv);
    Vector<double> obsVec;
    Matrix<double> H;
    Matrix<double> measCov;

    bool ok = mm.build(currentObs, currentGeo, obsVec, H, measCov);
    if (!ok)
    {
        // No valid observations this epoch
        return KalmanReturn::SkipThisEpoch;
    }

    // Update KalmanFilter members for measurement update
    Data = obsVec;
    Partials = H;
    MCov = measCov;

    // Zero the state for extended Kalman filter (deviations from nominal)
    State = Vector<double>(State.size(), 0.0);

    return KalmanReturn::Process;
}

// --------------------------------------------------------------------------------
// defineTimestep - State transition and process noise
// --------------------------------------------------------------------------------

void PPP_EKF::defineTimestep(double T, double DT, const Vector<double>& StateIn,
                             const Matrix<double>& CovIn, bool useFlag)
{
    // Phi = I (identity) for PPP: all states are random constants or random walks
    const size_t n = sv.size();
    PhiInv = Matrix<double>(n, n, 0.0);
    for (size_t i = 0; i < n; ++i)
        PhiInv(i, i) = 1.0;

    // Process noise: build G and Rw
    Matrix<double> Gmat;
    Matrix<double> Rwmat;
    buildProcessNoise(DT, Gmat, Rwmat);

    G = Gmat;
    Rw = Rwmat;

    // Control vector (zero for unbiased process noise)
    Control = Vector<double>(Gmat.cols(), 0.0);
}

// --------------------------------------------------------------------------------
// defineInterim - Handle satellite changes and position update
// --------------------------------------------------------------------------------

int PPP_EKF::defineInterim(int which, double Time)
{
    if (which == 1)
    {
        // Before measurement update
        // Check satellite changes (handled in defineMeasurements)
        if (stateResized)
        {
            // After satellite change, the SRI was reshaped.
            // Need to invert to get valid State and Cov.
            doInversions = true;
            stateResized = false;
        }
    }
    else if (which == 2)
    {
        // Between MU and TU
        // Nothing special needed
    }
    else if (which == 3)
    {
        // After TU
        // Update nominal position with the estimated correction
        updateNominalPosition();
        // Reset position deviations to zero (extended KF)
        State(0) = 0.0;
        State(1) = 0.0;
        State(2) = 0.0;
        // Sync back to sv
        sv.state = State;
        sv.cov = Cov;
    }
    return 0;
}

// --------------------------------------------------------------------------------
// Internal helpers
// --------------------------------------------------------------------------------

bool PPP_EKF::satellitesChanged() const
{
    set<SatID> currentSet;
    for (const auto& ro : currentObs)
        currentSet.insert(ro.sat);

    set<SatID> prevSet(previousSats.begin(), previousSats.end());

    return currentSet != prevSet;
}

void PPP_EKF::updateStateForNewSatellites()
{
    // Extract current satellite list from observations
    vector<SatID> newSats;
    for (const auto& ro : currentObs)
        newSats.push_back(ro.sat);
    sort(newSats.begin(), newSats.end());

    // Get old dimension and states
    size_t oldDim = sv.size();
    Vector<double> oldState = sv.state;
    Matrix<double> oldCov = sv.cov;

    // Update state vector with new satellite list
    sv.setSatellites(newSats);

    // Reshape the SRI in KalmanFilter base class
    // We need to extend/shrink the SRI to match new Namelist
    srif.reshape(sv.names);

    previousSats = newSats;
    stateResized = true;
}

void PPP_EKF::updateNominalPosition()
{
    // Add estimated position deviation to nominal position
    if (sv.size() >= 3)
    {
        double dx = State(0);
        double dy = State(1);
        double dz = State(2);
        nominalPos.setECEF(nominalPos.X() + dx,
                           nominalPos.Y() + dy,
                           nominalPos.Z() + dz);
    }
}

void PPP_EKF::buildProcessNoise(double DT, Matrix<double>& Gout, Matrix<double>& Rwout)
{
    // Number of noise states: clocks, ZWD, ionosphere (position and ambiguity have no noise)
    size_t nNoise = countNoiseStates();
    size_t nState = sv.size();

    // G is N x nNoise (maps noise to state)
    Gout = Matrix<double>(nState, nNoise, 0.0);
    // Rw is nNoise x nNoise (square root of inverse process noise covariance)
    Rwout = Matrix<double>(nNoise, nNoise, 0.0);

    size_t noiseIdx = 0;

    // Clock process noise (white noise)
    for (size_t i = 0; i < sv.numSystems(); ++i)
    {
        size_t clkIdx = sv.clkIndex(refSystem) + i;
        if (clkIdx < sv.size())
        {
            Gout(clkIdx, noiseIdx) = 1.0;
            // White noise variance = sigma^2 * dt
            double var = PPPProcessNoise::CLK * PPPProcessNoise::CLK * DT;
            Rwout(noiseIdx, noiseIdx) = 1.0 / sqrt(var);
            noiseIdx++;
        }
    }

    // ZWD random walk
    size_t zwdIdx = sv.zwdIndex();
    if (zwdIdx < sv.size())
    {
        Gout(zwdIdx, noiseIdx) = 1.0;
        double var = PPPProcessNoise::ZWD * PPPProcessNoise::ZWD * DT;
        Rwout(noiseIdx, noiseIdx) = 1.0 / sqrt(var);
        noiseIdx++;
    }

    // Ionosphere random walk (UC mode only)
    if (mode == PPPMode::Uncombined)
    {
        for (size_t i = 0; i < previousSats.size(); ++i)
        {
            size_t ionoIdx = sv.ionoIndex(previousSats[i]);
            if (ionoIdx < sv.size())
            {
                Gout(ionoIdx, noiseIdx) = 1.0;
                double var = PPPProcessNoise::IONO * PPPProcessNoise::IONO * DT;
                Rwout(noiseIdx, noiseIdx) = 1.0 / sqrt(var);
                noiseIdx++;
            }
        }
    }
}

size_t PPP_EKF::countNoiseStates() const
{
    size_t count = sv.numSystems(); // Clock noise for each system
    count += 1;                     // ZWD
    if (mode == PPPMode::Uncombined)
        count += sv.numSatellites(); // Ionosphere per satellite
    return count;
}

} // namespace ppp
} // namespace gnsstk
