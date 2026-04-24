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

//==============================================================================
//
//  PPPSolve - Precise Point Positioning using GNSSTk
//
//==============================================================================

/** \page apps
 * - \subpage PPPSolve - Compute PPP position solution from RINEX
 * \page PPPSolve
 * \tableofcontents
 *
 * \section PPPSolve_name NAME
 * PPPSolve - Compute Precise Point Positioning solution from RINEX
 *
 * \section PPPSolve_synopsis SYNOPSIS
 * \b PPPSolve [\argarg{OPTION}] ...
 *
 * \section PPPSolve_description DESCRIPTION
 * The application reads one or more RINEX observation files, plus SP3
 * ephemeris and clock files, and computes a precise point positioning
 * solution using extended Kalman filtering.
 *
 * \dictionary
 * \dicterm{\--obs \argarg{FN}}
 * \dicdef{RINEX observation file name(s) [repeatable] ()}
 * \dicterm{\--eph \argarg{FN}}
 * \dicdef{SP3 ephemeris file name(s) [repeatable] ()}
 * \dicterm{\--clk \argarg{FN}}
 * \dicdef{RINEX clock file name(s) [repeatable] ()}
 * \dicterm{\--nav \argarg{FN}}
 * \dicdef{RINEX navigation file name(s) [repeatable] ()}
 * \dicterm{\--ant \argarg{FN}}
 * \dicdef{ANTEX antenna file name (e.g., igs14.atx) ()}
 * \dicterm{\--met \argarg{FN}}
 * \dicdef{RINEX meteorological file name(s) [repeatable] ()}
 * \dicterm{\--mode \argarg{M}}
 * \dicdef{Processing mode: IF (ionosphere-free) or UC (uncombined) (IF)}
 * \dicterm{\--systems \argarg{S}}
 * \dicdef{GNSS systems: G,R,E,C,J,I (default: G,R,E,C)}
 * \dicterm{\--elev \argarg{DEG}}
 * \dicdef{Minimum elevation angle (deg) (7.5)}
 * \dicterm{\--start \argarg{T[:F]}}
 * \dicdef{Start processing at this epoch ([Beginning of dataset])}
 * \dicterm{\--stop \argarg{T[:F]}}
 * \dicdef{Stop processing at this epoch ([End of dataset])}
 * \dicterm{\--decimate \argarg{DT}}
 * \dicdef{Decimate data to time interval dt (0: no decimation) (0.00)}
 * \dicterm{\--out \argarg{FN}}
 * \dicdef{Output PPP solution file name (ppp.out)}
 * \dicterm{\--log \argarg{FN}}
 * \dicdef{Output log file name (ppp.log)}
 * \dicterm{\--ref \argarg{P[:F]}}
 * \dicdef{Known position for residual computation ()}
 * \dicterm{\--verbose}
 * \dicdef{Print extended output (don't)}
 * \dicterm{\--help}
 * \dicdef{Print this syntax page and quit (don't)}
 * \enddictionary
 */

#include <algorithm>
#include <fstream>
#include <iostream>
#include <map>
#include <string>
#include <vector>

#include "CommandLine.hpp"
#include "CommonTime.hpp"
#include "Epoch.hpp"
#include "Exception.hpp"
#include "GPSWeekSecond.hpp"
#include "NavLibrary.hpp"
#include "NavSearchOrder.hpp"
#include "NewNavInc.h"
#include "Position.hpp"
#include "Rinex3NavBase.hpp"
#include "Rinex3NavData.hpp"
#include "Rinex3NavHeader.hpp"
#include "Rinex3NavStream.hpp"
#include "Rinex3ObsData.hpp"
#include "Rinex3ObsHeader.hpp"
#include "Rinex3ObsStream.hpp"
#include "RinexMetData.hpp"
#include "RinexMetHeader.hpp"
#include "RinexMetStream.hpp"
#include "RinexObsID.hpp"
#include "RinexSatID.hpp"
#include "SP3Data.hpp"
#include "SP3Header.hpp"
#include "SP3NavDataFactory.hpp"
#include "SP3Stream.hpp"
#include "TimeString.hpp"
#include "expandtilde.hpp"
#include "logstream.hpp"
#include "singleton.hpp"

#include "MultiFormatNavDataFactory.hpp"

#include "EphemerisRange.hpp"
#include "PreciseRange.hpp"

#include "PPPConstants.hpp"
#include "PPP_EKF.hpp"
#include "PPPCorrections.hpp"
#include "PPPMeasurementModel.hpp"
#include "PPPObsPreprocessor.hpp"
#include "PPPStateVector.hpp"

using namespace std;
using namespace gnsstk;
using namespace gnsstk::StringUtils;
using namespace gnsstk::ppp;

//------------------------------------------------------------------------------------
string Version(string("1.0 4/24/26"));

//------------------------------------------------------------------------------------
// Global configuration
class PPPConfig : public Singleton<PPPConfig>
{
  public:
    PPPConfig() noexcept
    {
        SetDefaults();
    }

    int ProcessUserInput(int argc, char** argv) noexcept;
    string BuildCommandLine(void) noexcept;
    int ExtraProcessing(string& errors, string& extras) noexcept;

    void SetDefaults(void) noexcept;

    // Command line data
    CommandLine opts;
    static const string PrgmName;
    string Title;

    vector<string> InputObsFiles;
    vector<string> InputSP3Files;
    vector<string> InputClkFiles;
    vector<string> InputNavFiles;
    vector<string> InputMetFiles;
    string AntexFile;

    string Obspath, SP3path, Clkpath, Navpath, Metpath;

    string startStr, stopStr;
    CommonTime beginTime, endTime;

    double decimate;
    double elevLimit;

    string modeStr;         // IF or UC
    string systemsStr;      // e.g., "G,R,E,C"
    PPPMode mode;
    vector<SatelliteSystem> systems;
    vector<RinexSatID> exclSat;

    string OutputFile;
    string LogFile;
    string refPosStr;

    bool help, verbose;
    int debug;

    // Output
    ofstream logstrm;

    // Stores
    NavLibrary navLib;
    shared_ptr<NavDataFactory> ndfp;
    list<RinexMetData> MetStore;

    // Reference position
    Position knownPos;

    // Time formats
    static const string calfmt, gpsfmt, longfmt;

    string PrgmDesc, cmdlineUsage, cmdlineErrors, cmdlineExtras;
    vector<string> cmdlineUnrecognized;
};

const string PPPConfig::PrgmName = string("PPPSolve");
const string PPPConfig::calfmt = string("%04Y/%02m/%02d %02H:%02M:%02S");
const string PPPConfig::gpsfmt = string("%4F %10.3g");
const string PPPConfig::longfmt = calfmt + " = %4F %w %10.3g %P";

//------------------------------------------------------------------------------------
// Main processing functions
int Initialize(string& errors);
int ProcessFiles(void);
void OutputPPPResult(const PPP_EKF& ppp, const CommonTime& t, ofstream& out);

//------------------------------------------------------------------------------------
int main(int argc, char** argv)
{
#include "NewNavInit.h"
    PPPConfig& C(PPPConfig::Instance());

    try
    {
        int iret;
        clock_t totaltime(clock());
        Epoch wallclkbeg;
        wallclkbeg.setLocalTime();

        C.Title = C.PrgmName + ", part of the GNSS Toolkit, Ver " + Version + ", Run " +
                  printTime(wallclkbeg, C.calfmt);

        for (;;)
        {
            if ((iret = C.ProcessUserInput(argc, argv)) != 0)
                break;

            string errs;
            if ((iret = Initialize(errs)) != 0)
            {
                LOG(ERROR) << "------- Input is not valid: ----------\n" << errs << "------- end errors -----------";
                break;
            }

            int nfiles = ProcessFiles();
            if (nfiles < 0)
                break;
            LOG(VERBOSE) << "Successfully processed " << nfiles << " RINEX observation file" << (nfiles > 1 ? "s." : ".");

            break;
        }

        if (iret == 0)
        {
            totaltime = clock() - totaltime;
            Epoch wallclkend;
            wallclkend.setLocalTime();
            ostringstream oss;
            oss << C.PrgmName << " timing: processing " << fixed << setprecision(3)
                << double(totaltime) / double(CLOCKS_PER_SEC) << " sec, wallclock: " << setprecision(0)
                << (wallclkend - wallclkbeg) << " sec.\n";
            LOGstrm << oss.str();
            cout << oss.str();
        }

        return iret;
    }
    catch (Exception& e)
    {
        cerr << "Exception: " << e.what() << endl;
    }
    catch (...)
    {
        cerr << "Unknown exception. Abort." << endl;
    }
    return 1;
}

//------------------------------------------------------------------------------------
int Initialize(string& errors)
{
    try
    {
        PPPConfig& C(PPPConfig::Instance());
        bool isValid(true);
        size_t nfile;
        ostringstream ossE;

        // Set up navigation store
        C.ndfp = make_shared<MultiFormatNavDataFactory>();
        auto* mfndfp = dynamic_cast<MultiFormatNavDataFactory*>(C.ndfp.get());
        C.navLib.addFactory(C.ndfp);
        C.navLib.setTypeFilter({NavMessageType::Ephemeris, NavMessageType::Clock});

        // Add paths to filenames
        include_path(C.Obspath, C.InputObsFiles);
        include_path(C.SP3path, C.InputSP3Files);
        include_path(C.Clkpath, C.InputClkFiles);
        include_path(C.Navpath, C.InputNavFiles);
        include_path(C.Metpath, C.InputMetFiles);

        expand_filename(C.InputSP3Files);
        expand_filename(C.InputClkFiles);
        expand_filename(C.InputNavFiles);
        expand_filename(C.InputMetFiles);

        // Check obs files
        if (C.InputObsFiles.empty())
        {
            ossE << "Error : no RINEX observation files specified.\n";
            isValid = false;
        }

        // Load SP3 files
        if (!C.InputSP3Files.empty())
        {
            try
            {
                for (nfile = 0; nfile < C.InputSP3Files.size(); nfile++)
                {
                    LOG(VERBOSE) << "Load SP3 file " << C.InputSP3Files[nfile];
                    C.ndfp->addDataSource(C.InputSP3Files[nfile]);
                }
            }
            catch (Exception& e)
            {
                ossE << "Error: failed to read SP3 files: " << e.getText(0) << endl;
                isValid = false;
            }
        }
        else
        {
            ossE << "Error: no SP3 ephemeris files specified.\n";
            isValid = false;
        }

        // Load clock files
        if (!C.InputClkFiles.empty())
        {
            try
            {
                for (nfile = 0; nfile < C.InputClkFiles.size(); nfile++)
                {
                    LOG(VERBOSE) << "Load Clock file " << C.InputClkFiles[nfile];
                    C.ndfp->addDataSource(C.InputClkFiles[nfile]);
                }
            }
            catch (Exception& e)
            {
                ossE << "Error: failed to read clock files: " << e.getText(0) << endl;
                isValid = false;
            }
        }

        // Load nav files (optional)
        if (!C.InputNavFiles.empty())
        {
            try
            {
                for (nfile = 0; nfile < C.InputNavFiles.size(); nfile++)
                {
                    C.ndfp->addDataSource(C.InputNavFiles[nfile]);
                }
            }
            catch (Exception& e)
            {
                ossE << "Warning: failed to read nav files: " << e.getText(0) << endl;
            }
        }

        // Configure SP3 store
        auto sp3fact = mfndfp->getFactory<SP3NavDataFactory>();
        if (sp3fact && sp3fact->size() > 0)
        {
            sp3fact->setClockLinearInterp();
            sp3fact->rejectPredPositions(true);
            sp3fact->rejectPredClocks(true);
            sp3fact->setPositionInterpOrder(10);
        }

        errors = ossE.str();
        if (!isValid)
            return -5;
        return 0;
    }
    catch (Exception& e)
    {
        GNSSTK_RETHROW(e);
    }
}

//------------------------------------------------------------------------------------
int ProcessFiles(void)
{
    try
    {
        PPPConfig& C(PPPConfig::Instance());
        int nfiles = 0;

        // Open output file
        ofstream outstrm;
        if (!C.OutputFile.empty())
        {
            outstrm.open(C.OutputFile.c_str(), ios::out);
            if (!outstrm.is_open())
            {
                LOG(WARNING) << "Warning: could not open output file " << C.OutputFile;
            }
            else
            {
                // Write header
                outstrm << "# PPPSolve output\n";
                outstrm << "# week sow x y z sx sy sz dt zwd nsv\n";
            }
        }

        // Initialize PPP EKF with nominal position
        Position nominalPos(0.0, 0.0, 0.0);
        // Parse mode and systems if not already done
        if (C.modeStr == "UC" || C.modeStr == "uc" || C.modeStr == "Uncombined")
            C.mode = PPPMode::Uncombined;
        else
            C.mode = PPPMode::IonosphereFree;

        if (!C.systemsStr.empty())
        {
            C.systems.clear();
            vector<string> sysChars = split(C.systemsStr, ',');
            for (const auto& sc : sysChars)
            {
                if (sc == "G" || sc == "GPS")
                    C.systems.push_back(SatelliteSystem::GPS);
                else if (sc == "R" || sc == "GLO")
                    C.systems.push_back(SatelliteSystem::Glonass);
                else if (sc == "E" || sc == "GAL")
                    C.systems.push_back(SatelliteSystem::Galileo);
                else if (sc == "C" || sc == "BDS")
                    C.systems.push_back(SatelliteSystem::BeiDou);
                else if (sc == "J" || sc == "QZSS")
                    C.systems.push_back(SatelliteSystem::QZSS);
                else if (sc == "I" || sc == "IRNSS")
                    C.systems.push_back(SatelliteSystem::IRNSS);
            }
        }

        PPP_EKF ppp(C.mode, C.systems, nominalPos);
        ppp.initializeFilter();

        // Initialize corrections
        PPPCorrections corrections;
        if (!C.AntexFile.empty())
            corrections.loadAntex(C.AntexFile);

        // Process each observation file
        for (size_t nfile = 0; nfile < C.InputObsFiles.size(); nfile++)
        {
            Rinex3ObsStream istrm(C.InputObsFiles[nfile].c_str(), ios::in);
            if (!istrm.is_open())
            {
                LOG(WARNING) << "Warning: could not open file " << C.InputObsFiles[nfile];
                continue;
            }

            Rinex3ObsHeader Rhead;
            Rinex3ObsData Rdata;

            try
            {
                istrm >> Rhead;
            }
            catch (Exception& e)
            {
                LOG(WARNING) << "Warning: Failed to read header from " << C.InputObsFiles[nfile];
                continue;
            }

            if (C.verbose)
            {
                LOG(VERBOSE) << "Processing " << C.InputObsFiles[nfile];
                Rhead.dump(LOGstrm);
            }

            // Preprocessor
            PPPObsPreprocessor preprocessor(C.mode, C.elevLimit, C.systems);

            // Process epochs
            while (istrm >> Rdata)
            {
                if (!istrm.good() || istrm.eof())
                    break;

                if (Rdata.epochFlag > 1 || Rdata.obs.empty())
                    continue;

                // Time limits
                if (Rdata.time < C.beginTime)
                    continue;
                if (Rdata.time > C.endTime)
                    break;

                // Decimation
                if (C.decimate > 0.0)
                {
                    double sow = GPSWeekSecond(Rdata.time).sow;
                    if (fmod(sow, C.decimate) > 0.25)
                        continue;
                }

                // Compute satellite positions and elevations using PreciseRange
                map<SatID, double> elevMap;
                map<SatID, PPPGeometry> geoMap;

                for (const auto& obsPair : Rdata.obs)
                {
                    const SatID& sat = obsPair.first;

                    try
                    {
                        CorrectedEphemerisRange CER;
                        CER.ComputeAtReceiveTime(Rdata.time, nominalPos, sat, C.navLib, NavSearchOrder::Nearest);
                        elevMap[sat] = CER.elevation;

                        PPPGeometry geo;
                        geo.sat = sat;
                        geo.rho = CER.rawrange;
                        geo.elev = CER.elevation;
                        geo.azim = CER.azimuth;
                        geo.cosines = CER.cosines;
                        geo.satClk = CER.svclkbias;
                        geo.tropDry = 0.0;  // Will be computed by trop model
                        geo.tropWetMap = 1.0 / sin(max(CER.elevation, 1.0) * DEG_TO_RAD);
                        geoMap[sat] = geo;
                    }
                    catch (Exception& e)
                    {
                        LOG(DEBUG) << "No ephemeris for " << sat << " at " << printTime(Rdata.time, C.gpsfmt);
                    }
                }

                // Preprocess observations
                auto preprocessed = preprocessor.processEpoch(Rhead, Rdata, elevMap);
                if (preprocessed.empty())
                    continue;

                // Convert to PPP raw observations and geometry
                vector<PPPRawObservation> rawObs;
                vector<PPPGeometry> geometry;

                for (const auto& pp : preprocessed)
                {
                    const SatID& sat = pp.first;
                    const PPPRequiredObs& ro = pp.second;

                    PPPRawObservation raw;
                    raw.sat = sat;
                    raw.P1 = ro.hasP1 ? ro.P1 : 0.0;
                    raw.P2 = ro.hasP2 ? ro.P2 : 0.0;
                    raw.L1 = ro.hasL1 ? ro.L1 : 0.0;
                    raw.L2 = ro.hasL2 ? ro.L2 : 0.0;
                    raw.lambda1 = ro.lambda1;
                    raw.lambda2 = ro.lambda2;
                    rawObs.push_back(raw);

                    auto geoIt = geoMap.find(sat);
                    if (geoIt != geoMap.end())
                        geometry.push_back(geoIt->second);
                }

                // Run PPP EKF for this epoch
                double dt = 30.0;  // Nominal epoch interval
                ppp.setEpochData(dt, rawObs, geometry);
                ppp.ForwardFilter(ppp.getFilterTime() + dt, dt);

                // Output results
                if (outstrm.is_open())
                    OutputPPPResult(ppp, Rdata.time, outstrm);

                // Update nominal position
                nominalPos = ppp.getPosition();
            }

            istrm.close();
            nfiles++;
        }

        outstrm.close();
        return nfiles;
    }
    catch (Exception& e)
    {
        GNSSTK_RETHROW(e);
    }
}

//------------------------------------------------------------------------------------
void OutputPPPResult(const PPP_EKF& ppp, const CommonTime& t, ofstream& out)
{
    Position pos = ppp.getPosition();
    Matrix<double> cov = ppp.getPositionCov();

    GPSWeekSecond gws(t);

    // Basic epoch/position/clock/ZWD output
    out << fixed << setprecision(4)
        << gws.week << " " << gws.sow << " "
        << pos.X() << " " << pos.Y() << " " << pos.Z() << " "
        << sqrt(cov(0, 0)) << " " << sqrt(cov(1, 1)) << " " << sqrt(cov(2, 2)) << " "
        << ppp.getClock(SatelliteSystem::GPS) / C_MPS << " "  // GPS clock in seconds
        << ppp.getZWD() << " "
        << ppp.getEpochCount();

    // Additional system clocks (ISB output for IF mode)
    PPPConfig& C(PPPConfig::Instance());
    if (ppp.getStateVector().getMode() == PPPMode::IonosphereFree)
    {
        for (const auto& sys : C.systems)
        {
            if (sys == SatelliteSystem::GPS)
                continue;
            out << " " << ppp.getISB(sys) / C_MPS;
        }
    }
    else
    {
        // UC mode: output each system clock
        for (const auto& sys : C.systems)
        {
            if (sys == SatelliteSystem::GPS)
                continue;
            out << " " << ppp.getClock(sys) / C_MPS;
        }
    }

    // Satellite-dependent states (ambiguities / ionosphere)
    const PPPStateVector& sv = ppp.getStateVector();
    if (sv.getMode() == PPPMode::IonosphereFree)
    {
        // Output number of ambiguity states and their values
        out << " " << sv.numSatellites();
        // Note: full ambiguity dump can be very long; just output count here
    }
    else
    {
        // UC mode: output ionosphere count
        out << " " << sv.numSatellites();
    }

    out << "\n";
}

//------------------------------------------------------------------------------------
// Command line processing
//------------------------------------------------------------------------------------

int PPPConfig::ProcessUserInput(int argc, char** argv) noexcept
{
    opts.DefineUsageString(PrgmName + " [options]");
    PrgmDesc = BuildCommandLine();

    int iret = opts.ProcessCommandLine(argc, argv, PrgmDesc, cmdlineUsage, cmdlineErrors, cmdlineUnrecognized);

    if (iret == -2 || iret == -3)
        return iret;

    if (opts.hasHelp())
    {
        LOG(INFO) << Title;
        LOG(INFO) << cmdlineUsage;
        return 1;
    }

    iret = ExtraProcessing(cmdlineErrors, cmdlineExtras);

    if (!cmdlineErrors.empty())
    {
        LOG(INFO) << "Errors found on command line:\n " << cmdlineErrors << "\nEnd of command line errors.";
        return 1;
    }

    return 0;
}

string PPPConfig::BuildCommandLine(void) noexcept
{
    string PrgmDesc = " Program " + PrgmName +
                      " reads one or more RINEX (v.2+) observation files, plus\n"
                      " SP3 ephemeris and optional clock files, and computes a\n"
                      " Precise Point Positioning solution using extended Kalman filtering.\n";

    opts.Add(0, "obs", "fn", true, true, &InputObsFiles, "# Required input:", "RINEX observation file name(s)");
    opts.Add(0, "eph", "fn", true, true, &InputSP3Files, "", "SP3 ephemeris file name(s)");
    opts.Add(0, "clk", "fn", true, false, &InputClkFiles, "", "RINEX clock file name(s)");
    opts.Add(0, "nav", "fn", true, false, &InputNavFiles, "", "RINEX nav file name(s)");
    opts.Add(0, "ant", "fn", false, false, &AntexFile, "", "ANTEX antenna file (e.g., igs14.atx)");
    opts.Add(0, "met", "fn", true, false, &InputMetFiles, "", "RINEX meteorological file name(s)");

    opts.Add(0, "mode", "M", false, false, &modeStr,
             "# Processing options:", "Processing mode: IF (ionosphere-free) or UC (uncombined)");
    opts.Add(0, "systems", "S", false, false, &systemsStr, "", "GNSS systems to use (e.g., G,R,E,C)");
    opts.Add(0, "elev", "deg", false, false, &elevLimit, "", "Minimum elevation angle (deg)");
    opts.Add(0, "start", "t[:f]", false, false, &startStr, "", "Start processing at this epoch");
    opts.Add(0, "stop", "t[:f]", false, false, &stopStr, "", "Stop processing at this epoch");
    opts.Add(0, "decimate", "dt", false, false, &decimate, "", "Decimate data to interval dt (s)");

    opts.Add(0, "out", "fn", false, false, &OutputFile, "# Output:", "Output PPP solution file name");
    opts.Add(0, "log", "fn", false, false, &LogFile, "", "Output log file name");
    opts.Add(0, "ref", "p[:f]", false, false, &refPosStr, "", "Known position for residual computation");

    opts.Add(0, "verbose", "", false, false, &verbose, "# Other:", "Print extended output");
    opts.Add(0, "help", "", false, false, &help, "", "Print this syntax page and quit");

    return PrgmDesc;
}

int PPPConfig::ExtraProcessing(string& errors, string& extras) noexcept
{
    // Reference position
    if (!refPosStr.empty())
    {
        vector<string> fld = split(refPosStr, ',');
        if (fld.size() == 3)
        {
            try
            {
                double X = asDouble(fld[0]);
                double Y = asDouble(fld[1]);
                double Z = asDouble(fld[2]);
                knownPos.setECEF(X, Y, Z);
            }
            catch (Exception& e)
            {
                errors += "Error: invalid position in --ref arg\n";
            }
        }
    }

    // Start/stop times
    if (!startStr.empty() && startStr != "[Beginning of dataset]")
    {
        try
        {
            Epoch ep;
            ep.scanf(startStr, "%Y,%m,%d,%H,%M,%S");
            beginTime = static_cast<CommonTime>(ep);
        }
        catch (Exception&)
        {
            errors += "Error: invalid start time\n";
        }
    }
    else
    {
        beginTime = CommonTime::BEGINNING_OF_TIME;
    }

    if (!stopStr.empty() && stopStr != "[End of dataset]")
    {
        try
        {
            Epoch ep;
            ep.scanf(stopStr, "%Y,%m,%d,%H,%M,%S");
            endTime = static_cast<CommonTime>(ep);
        }
        catch (Exception&)
        {
            errors += "Error: invalid stop time\n";
        }
    }
    else
    {
        endTime = CommonTime::END_OF_TIME;
    }

    // Open log file
    logstrm.open(LogFile.c_str(), ios::out);
    if (logstrm.is_open())
    {
        pLOGstrm = &logstrm;
        LOG(INFO) << Title;
    }

    verbose = (LOGlevel >= VERBOSE);
    debug = (LOGlevel - DEBUG);

    return 0;
}

void PPPConfig::SetDefaults(void) noexcept
{
    modeStr = "IF";
    systemsStr = "G,R,E,C";
    mode = PPPMode::IonosphereFree;
    systems = {SatelliteSystem::GPS, SatelliteSystem::Galileo,
               SatelliteSystem::BeiDou, SatelliteSystem::Glonass};
    elevLimit = 7.5;
    decimate = 0.0;
    OutputFile = "ppp.out";
    LogFile = "ppp.log";
    verbose = false;
    debug = -1;
    help = false;
}
