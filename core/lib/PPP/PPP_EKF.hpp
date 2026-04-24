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

#ifndef PPP_EKF_HPP
#define PPP_EKF_HPP

#include <set>
#include <vector>

#include "KalmanFilter.hpp"
#include "Position.hpp"

#include "PPPConstants.hpp"
#include "PPPMeasurementModel.hpp"
#include "PPPStateVector.hpp"

namespace gnsstk
{
namespace ppp
{

/**
 * @class PPP_EKF
 * @brief Extended Kalman Filter for PPP processing.
 *
 * Inherits from gnsstk::KalmanFilter and implements the PPP-specific
 * measurement model, process noise, and state management.
 *
 * Usage:
 *   PPP_EKF ppp(PPPMode::IonosphereFree, systems, refPos);
 *   ppp.initializeFilter();
 *   for each epoch:
 *     ppp.setEpochData(dt, rawObs, geometry);
 *     ppp.ForwardFilter(ppp.time + dt, dt);
 *     Position pos = ppp.getPosition();
 */
class PPP_EKF : public KalmanFilter
{
  public:
    /**
     * @brief Constructor.
     * @param mode      PPP processing mode (IF or UC)
     * @param systems   List of satellite systems to process
     * @param refPos    Reference/nominal receiver position (ECEF, m)
     * @param refSys    Reference system for clock (default GPS)
     */
    PPP_EKF(PPPMode mode, const std::vector<SatelliteSystem>& systems,
            const Position& refPos,
            SatelliteSystem refSys = SatelliteSystem::GPS);

    /// Destructor
    ~PPP_EKF() override = default;

    // -------------------------------------------------------------------------
    // Data input (call before ForwardFilter for each epoch)
    // -------------------------------------------------------------------------

    /**
     * @brief Set the observation data for the upcoming epoch.
     * @param dt   Time interval from previous epoch (s)
     * @param obs  Raw observations (one per satellite)
     * @param geo  Pre-computed geometry (one per satellite)
     */
    void setEpochData(double dt, const std::vector<PPPRawObservation>& obs,
                      const std::vector<PPPGeometry>& geo);

    // -------------------------------------------------------------------------
    // Result queries (valid after measurement update)
    // -------------------------------------------------------------------------

    /// @return Current estimated position (ECEF, m)
    Position getPosition() const;

    /// @return Current position covariance (3x3, m^2)
    Matrix<double> getPositionCov() const;

    /// @return Reference to the PPP state vector (includes state layout)
    const PPPStateVector& getStateVector() const
    {
        return sv;
    }

    /// @return Current ZWD estimate (m)
    double getZWD() const;

    /// @return Clock bias for given system (m)
    double getClock(SatelliteSystem sys) const;

    /// @return ISB for given system relative to reference (m), IF mode only
    double getISB(SatelliteSystem sys) const;

    /// @return Current filter time (seconds since start)
    double getFilterTime() const
    {
        return time;
    }

    /// @return Number of processed epochs
    int getEpochCount() const
    {
        return NMU;
    }

    // -------------------------------------------------------------------------
    // KalmanFilter virtual overrides
    // -------------------------------------------------------------------------

    /** Initialize the filter with apriori state and covariance. */
    int defineInitial(double& T0, Vector<double>& X, Matrix<double>& Cov) override;

    /** Build measurement equations for the current epoch. */
    KalmanReturn defineMeasurements(double& T, const Vector<double>& X,
                                    const Matrix<double>& C, bool useFlag) override;

    /** Build state transition and process noise matrices. */
    void defineTimestep(double T, double DT, const Vector<double>& State,
                        const Matrix<double>& Cov, bool useFlag) override;

    /**
     * Handle satellite additions/removals and position update.
     * @param which 1=before MU, 2=between MU and TU, 3=after TU
     * @return >=0 to continue, -1 to skip epoch
     */
    int defineInterim(int which, double Time) override;

  private:
    // PPP configuration
    PPPStateVector sv;                    ///< State vector layout manager
    PPPMode mode;                         ///< Processing mode
    Position nominalPos;                  ///< Nominal/reference position (ECEF)
    SatelliteSystem refSystem;            ///< Reference clock system

    // Epoch data (set before ForwardFilter)
    double epochDT;                       ///< Epoch interval (s)
    std::vector<PPPRawObservation> currentObs;  ///< Current epoch observations
    std::vector<PPPGeometry> currentGeo;        ///< Current epoch geometry

    // Satellite tracking
    std::vector<SatID> previousSats;      ///< Satellites from previous epoch
    bool stateResized;                    ///< True if state was resized this epoch

    // -------------------------------------------------------------------------
    // Internal helpers
    // -------------------------------------------------------------------------

    /// Check if satellite list changed from previous epoch
    bool satellitesChanged() const;

    /// Update state vector and SRI when satellites change
    void updateStateForNewSatellites();

    /// Update nominal position after measurement update (position correction)
    void updateNominalPosition();

    /// Build process noise matrices G and Rw
    void buildProcessNoise(double DT, Matrix<double>& Gout, Matrix<double>& Rwout);

    /// @return Number of states with process noise
    size_t countNoiseStates() const;
};

} // namespace ppp
} // namespace gnsstk

#endif // PPP_EKF_HPP
