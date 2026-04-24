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

#ifndef PPP_MEASUREMENT_MODEL_HPP
#define PPP_MEASUREMENT_MODEL_HPP

#include <vector>

#include "Matrix.hpp"
#include "Position.hpp"
#include "SatID.hpp"
#include "Triple.hpp"
#include "Vector.hpp"

#include "PPPConstants.hpp"
#include "PPPStateVector.hpp"

namespace gnsstk
{
namespace ppp
{

/// Geometry information for a satellite, pre-computed from PreciseRange
struct PPPGeometry
{
    SatID sat;          ///< Satellite ID
    double rho;         ///< Geometric range (m)
    double elev;        ///< Elevation angle (deg)
    double azim;        ///< Azimuth angle (deg)
    Triple cosines;     ///< Direction cosines (unit vector from sat to rx)
    double satClk;      ///< Satellite clock bias (m)
    double tropDry;     ///< Dry troposphere delay (m), mapped to elevation
    double tropWetMap;  ///< Wet mapping function value
};

/// Raw observation for a satellite, extracted from RINEX
struct PPPRawObservation
{
    SatID sat;          ///< Satellite ID
    double P1;          ///< Frequency-1 pseudorange (m), 0 if invalid
    double P2;          ///< Frequency-2 pseudorange (m), 0 if invalid
    double L1;          ///< Frequency-1 carrier phase (m), 0 if invalid
    double L2;          ///< Frequency-2 carrier phase (m), 0 if invalid
    double lambda1;     ///< Frequency-1 wavelength (m)
    double lambda2;     ///< Frequency-2 wavelength (m)
};

/**
 * @class PPPMeasurementModel
 * @brief Builds the PPP observation equations for IF or UC mode.
 *
 * Given raw observations and pre-computed geometry, this class constructs:
 *   - Observation vector (pre-fit residuals)
 *   - Design matrix (partials w.r.t. state parameters)
 *   - Measurement covariance matrix
 *
 * IF mode: 2 obs per satellite (P_IF, L_IF)
 * UC mode: 4 obs per satellite (P1, P2, L1, L2)
 */
class PPPMeasurementModel
{
  public:
    /**
     * @brief Constructor.
     * @param sv Reference to the PPP state vector (defines state layout)
     */
    explicit PPPMeasurementModel(const PPPStateVector& sv);

    /// Destructor
    ~PPPMeasurementModel() = default;

    // -------------------------------------------------------------------------
    // Unified interface
    // -------------------------------------------------------------------------

    /**
     * @brief Build observation equations according to current mode.
     * @param rawObs   Raw observations (one per satellite)
     * @param geometry Pre-computed geometry (one per satellite, same order)
     * @param obs      Output: observation (pre-fit residual) vector
     * @param H        Output: design matrix (partials)
     * @param measCov  Output: measurement covariance matrix
     * @return true if at least one valid observation was built
     */
    bool build(const std::vector<PPPRawObservation>& rawObs,
               const std::vector<PPPGeometry>& geometry,
               Vector<double>& obs,
               Matrix<double>& H,
               Matrix<double>& measCov) const;

    // -------------------------------------------------------------------------
    // Mode-specific builders
    // -------------------------------------------------------------------------

    /**
     * @brief Build IF-mode observation equations.
     * Requires P1, P2, L1, L2 for each satellite. Invalid observations are skipped.
     */
    bool buildIF(const std::vector<PPPRawObservation>& rawObs,
                 const std::vector<PPPGeometry>& geometry,
                 Vector<double>& obs,
                 Matrix<double>& H,
                 Matrix<double>& measCov) const;

    /**
     * @brief Build UC-mode observation equations.
     * Requires P1, P2, L1, L2 for each satellite. Invalid observations are skipped.
     */
    bool buildUC(const std::vector<PPPRawObservation>& rawObs,
                 const std::vector<PPPGeometry>& geometry,
                 Vector<double>& obs,
                 Matrix<double>& H,
                 Matrix<double>& measCov) const;

    // -------------------------------------------------------------------------
    // Helpers
    // -------------------------------------------------------------------------

    /**
     * @brief Compute IF combination coefficients.
     * @param f1 Frequency 1 (Hz)
     * @param f2 Frequency 2 (Hz)
     * @param alpha Output: alpha coefficient for f1
     * @param beta  Output: beta coefficient for f2
     */
    static void ionoFreeCoefficients(double f1, double f2,
                                     double& alpha, double& beta);

    /// @return Number of observations per satellite for current mode
    size_t obsPerSat() const;

    /// @return Expected observation dimension for N satellites
    size_t expectedObsDim(size_t nSats) const;

  private:
    const PPPStateVector& sv;   ///< Reference state vector

    /// Simple wet mapping function (approximate Niell)
    static double wetMappingFunction(double elevDeg);

    /// Fill the position partials (direction cosines) into a row of H
    static void fillPositionPartials(Matrix<double>& H, size_t row,
                                     const Triple& cosines);
};

} // namespace ppp
} // namespace gnsstk

#endif // PPP_MEASUREMENT_MODEL_HPP
