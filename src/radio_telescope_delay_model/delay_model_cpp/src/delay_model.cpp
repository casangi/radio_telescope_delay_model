// C++ port of the CALC11 ALMA driver (almacalc.f) -- see delay_model.hpp.
//
// The port reproduces almacalc.f statement for statement: it fills the same
// COMMON blocks (SITCM, /ephem/, EOPCM, WOBCM, reading CONVDS from CMATH
// after aINITL has set it), performs the identical scalar arithmetic in the
// identical order, and calls the identical Fortran entry points (aINITL once,
// aDRIVR per time and antenna). Bit-identical agreement with the reference
// Fortran driver is asserted in the test suite.

#include "delay_model.hpp"

#include <algorithm>
#include <cstring>
#include <mutex>
#include <stdexcept>

namespace {

constexpr int kMaxStat = 100;  // Max_Stat in cmxst11.i / param11.i

// ---------------------------------------------------------------------------
// gfortran COMMON block layouts. Each struct mirrors the declared COMMON
// order of its include file exactly (Fortran column-major arrays appear here
// with reversed index order); prefix structs are used where only leading
// members are accessed, which is well-defined because a COMMON block is a
// single storage sequence starting at the symbol.
// ---------------------------------------------------------------------------

// COMMON / SITCM / of cmxst11.i.
extern "C" struct {
    double cfrad[kMaxStat];
    double plat[kMaxStat][3];
    double plon[kMaxStat][3];
    double sitaxo[kMaxStat];
    double sitoam[kMaxStat][11];
    double sitoph[kMaxStat][11];
    double sitxyz[kMaxStat][3];
    double snrm[kMaxStat][3];
    double sitzen[kMaxStat];
    double tcrot[kMaxStat][3][3];
    double xlat[kMaxStat];
    double xlon[kMaxStat];
    double sithoa[kMaxStat][2][11];
    double sithop[kMaxStat][2][11];
    double height[kMaxStat];
    double rtrot[kMaxStat][3][3];
    double glat[kMaxStat];
    double dbtilt[kMaxStat][2];
    double rotilt[kMaxStat][3][3];
    double optl6[kMaxStat][6];
    int zero_site;
    short ktype[kMaxStat];
    short nlast[2];
    short lnsite[kMaxStat][4];
    short numsit;
} sitcm_;

// COMMON / CMATH / (set by aINITL).
extern "C" struct {
    double pi, twopi, halfpi, convd, convds, convhs, secday;
} cmath_;

// Leading members of COMMON / EOPCM / as used by almacalc.f.
extern "C" struct {
    double atmutc[3];
    double roteph[20][2];
    double a1utc[3];
    double a1diff[3];
} eopcm_;

// Leading members of COMMON / WOBCM / of cmwob11.i (through DWOBY).
extern "C" struct {
    double dwobp[2][2];
    double rwobx[3][3];
    double rwoby[3][3];
    double wobif[3];
    double wobx;
    double woby;
    double xywob[20][2];
    double dwobx;
    double dwoby;
} wobcm_;

// COMMON /ephem/ of d_input.i: a single CHARACTER*128 (raw bytes).
extern "C" char ephem_[128];

// COMMON /OBSRVN/ of cobsn11.i: the per-observation geocenter-station marker
// (set by SITG from Zero_site; reset here every call because SITG never
// clears it when no geocenter station is present).
extern "C" struct {
    int nzero;
} obsrvn_;

// Fortran entry points of the (unported) CALC11 core.
extern "C" void ainitl_();
extern "C" void adrivr_(const int* iref, const int* iremot, const double* utc,
                        const double* xjd, const double* ut1,
                        const double* dut1at, const double* cfsite,
                        const double* sitlon, const double* sitlat,
                        const double* right_asc, const double* declination,
                        const double* sur_tp, const double* sur_pr,
                        const double* sur_hm, const int* ss_obj,
                        double* delay_vac, double* dry_atm1, double* dry_atm2,
                        double* wet_atm1, double* wet_atm2);

// The CALC11 core keeps its state in COMMON blocks: one computation at a time.
std::mutex core_mutex;

}  // namespace

namespace rtdm::delay_model {

void calc_delay_model(
    const double reference_position[3], std::size_t n_antenna,
    const double* antenna_x, const double* antenna_y, const double* antenna_z,
    const double* temperature_celsius, const double* pressure_hpa,
    const double* humidity_fraction, std::size_t n_time, const double* mjd_utc,
    const double* right_ascension, const double* declination,
    const double* polar_motion_x_arcsec, const double* polar_motion_y_arcsec,
    const double* ut1_minus_utc_seconds, double leap_seconds,
    const double* axis_offset_metres, const std::string& ephemeris_path,
    double* geometric_delay, double* dry_delay, double* wet_delay,
    const int* axis_type, const double* ocean_vertical_amplitude,
    const double* ocean_vertical_phase,
    const double* ocean_horizontal_amplitude,
    const double* ocean_horizontal_phase,
    const double* ocean_pole_tide_coefficients) {
    if (n_antenna + 1 > static_cast<std::size_t>(kMaxStat))
        throw std::invalid_argument(
            "CALC11 supports at most 99 antennas plus the array reference.");
    if (n_time < 2)
        throw std::invalid_argument(
            "CALC11 needs at least two times (EOP rates are formed from the "
            "first and last time).");
    if (ephemeris_path.size() > 128)
        throw std::invalid_argument("ephemeris_path exceeds 128 characters.");
    if (axis_type != nullptr)
        for (std::size_t i = 0; i < n_antenna; ++i)
            if (axis_type[i] < 1 || axis_type[i] > 5)
                throw std::invalid_argument(
                    "axis_type entries must be CALC codes 1..5.");

    std::lock_guard<std::mutex> lock(core_mutex);

    // --- Fill the CALC site arrays (almacalc.f lines "Fill calc site...",
    //     extended with the per-station data difxcalc11's dinit.f loads) ---
    sitcm_.numsit = static_cast<short>(n_antenna + 1);
    sitcm_.sitxyz[0][0] = reference_position[0];
    sitcm_.sitxyz[0][1] = reference_position[1];
    sitcm_.sitxyz[0][2] = reference_position[2];
    sitcm_.sitaxo[0] = 0.0;
    sitcm_.ktype[0] = 3;
    // A (0,0,0) reference selects difxcalc11's geocenter mode: station 1 is
    // the geocenter (Zero_site = 1), and CALC's geocenter special-casing
    // (dummy site geometry, no station-1 atmosphere/tides) applies -- the
    // convention of difxcalc11's per-antenna .im model. Any real reference
    // position keeps the normal almacalc behaviour.
    const bool geocenter_reference = reference_position[0] == 0.0 &&
                                     reference_position[1] == 0.0 &&
                                     reference_position[2] == 0.0;
    sitcm_.zero_site = geocenter_reference ? 1 : 0;
    obsrvn_.nzero = 0;  // SITG re-marks the geocenter station per observation
    for (std::size_t n = 0; n < n_antenna + 1; ++n) {
        const std::size_t a = n - 1;  // antenna index when n >= 1
        if (n >= 1) {
            sitcm_.sitxyz[n][0] = antenna_x[a];
            sitcm_.sitxyz[n][1] = antenna_y[a];
            sitcm_.sitxyz[n][2] = antenna_z[a];
            sitcm_.sitaxo[n] = axis_offset_metres[a];
            sitcm_.ktype[n] =
                axis_type ? static_cast<short>(axis_type[a]) : short{3};
        }
        for (int k = 0; k < 11; ++k) {
            sitcm_.sitoam[n][k] =
                (n >= 1 && ocean_vertical_amplitude)
                    ? ocean_vertical_amplitude[a * 11 + k] : 0.0;
            sitcm_.sitoph[n][k] =
                (n >= 1 && ocean_vertical_phase)
                    ? ocean_vertical_phase[a * 11 + k] : 0.0;
            for (int l = 0; l < 2; ++l) {
                sitcm_.sithoa[n][l][k] =
                    (n >= 1 && ocean_horizontal_amplitude)
                        ? ocean_horizontal_amplitude[(a * 2 + l) * 11 + k] : 0.0;
                sitcm_.sithop[n][l][k] =
                    (n >= 1 && ocean_horizontal_phase)
                        ? ocean_horizontal_phase[(a * 2 + l) * 11 + k] : 0.0;
            }
        }
        for (int j = 0; j < 6; ++j)
            sitcm_.optl6[n][j] =
                (n >= 1 && ocean_pole_tide_coefficients)
                    ? ocean_pole_tide_coefficients[a * 6 + j] : 0.0;
        // Axis tilts are not applied by this CALC 11 version (caxom's ROTAXIS
        // use is disabled upstream); identity matrices mirror difxcalc11's
        // dSITI result for zero tilts.
        sitcm_.dbtilt[n][0] = 0.0;
        sitcm_.dbtilt[n][1] = 0.0;
        for (int i = 0; i < 3; ++i)
            for (int j = 0; j < 3; ++j)
                sitcm_.rotilt[n][j][i] = (i == j) ? 1.0 : 0.0;
    }

    // --- Name of the ephemeris file (CHARACTER*128, blank padded) ---
    std::memset(ephem_, ' ', sizeof(ephem_));
    std::memcpy(ephem_, ephemeris_path.data(), ephemeris_path.size());

    // --- Initialize the model modules (CALL aINITL) ---
    ainitl_();

    // --- EOP rate setup, exactly as in almacalc.f ---
    const double fjd = mjd_utc[0] + 2400000.5;
    const double deltime = (mjd_utc[n_time - 1] - mjd_utc[0]) * 86400.0;
    const double dut1at =
        1.0 + ((ut1_minus_utc_seconds[n_time - 1] - ut1_minus_utc_seconds[0]) /
               deltime);
    const double dwobxl =
        (polar_motion_x_arcsec[n_time - 1] - polar_motion_x_arcsec[0]) * 1.0e3 /
        deltime;
    const double dwobyl =
        (polar_motion_y_arcsec[n_time - 1] - polar_motion_y_arcsec[0]) * 1.0e3 /
        deltime;
    wobcm_.dwobx = dwobxl * cmath_.convds * 1.0e-3;  // radians/sec
    wobcm_.dwoby = dwobyl * cmath_.convds * 1.0e-3;  // radians/sec
    eopcm_.atmutc[0] = fjd;
    eopcm_.atmutc[1] = leap_seconds;
    eopcm_.atmutc[2] = 0.0;

    // --- Compute the delays ---
    const int iref = 1;  // array reference position (Fortran index)
    double cfsite[2][3];
    double sitlon[2];
    double sitlat[2];
    double sur_tp[2];
    double sur_pr[2];
    double sur_hm[2];
    for (std::size_t iscan = 0; iscan < n_time; ++iscan) {
        const auto ijd = static_cast<long long>(mjd_utc[iscan]);
        const double yjd = static_cast<double>(ijd);  // MJD at midnight
        const double xjd = yjd + 2400000.5;           // JD at midnight
        const double utc = mjd_utc[iscan] - yjd;      // UTC, fraction of day
        const double ut1 = utc * 86400.0 + ut1_minus_utc_seconds[iscan];
        const double wobxl = polar_motion_x_arcsec[iscan] * 1.0e3;  // mas
        const double wobyl = polar_motion_y_arcsec[iscan] * 1.0e3;  // mas
        wobcm_.wobx = wobxl * cmath_.convds * 1.0e-3;  // radians
        wobcm_.woby = wobyl * cmath_.convds * 1.0e-3;  // radians
        const double right_asc = right_ascension[iscan];
        const double decl = declination[iscan];
        sur_tp[0] = -99.0;
        sur_pr[0] = -99.0;
        sur_hm[0] = -99.0;
        const int ss_obj = 0;  // solar-system objects not supported here

        for (std::size_t iant = 0; iant < n_antenna; ++iant) {
            const int iremot = static_cast<int>(iant) + 2;  // Fortran index
            cfsite[0][0] = sitcm_.sitxyz[0][0];
            cfsite[0][1] = sitcm_.sitxyz[0][1];
            cfsite[0][2] = sitcm_.sitxyz[0][2];
            cfsite[1][0] = sitcm_.sitxyz[iant + 1][0];
            cfsite[1][1] = sitcm_.sitxyz[iant + 1][1];
            cfsite[1][2] = sitcm_.sitxyz[iant + 1][2];
            sitlon[0] = sitcm_.xlon[0];
            sitlat[0] = sitcm_.xlat[0];
            sitlon[1] = sitcm_.xlon[iant + 1];
            sitlat[1] = sitcm_.xlat[iant + 1];
            sur_tp[1] = temperature_celsius[iant];
            sur_pr[1] = pressure_hpa[iant];
            sur_hm[1] = humidity_fraction[iant];

            double delay_vac = 0.0;
            double dry_atm1 = 0.0;
            double dry_atm2 = 0.0;
            double wet_atm1 = 0.0;
            double wet_atm2 = 0.0;
            adrivr_(&iref, &iremot, &utc, &xjd, &ut1, &dut1at, &cfsite[0][0],
                    sitlon, sitlat, &right_asc, &decl, sur_tp, sur_pr, sur_hm,
                    &ss_obj, &delay_vac, &dry_atm1, &dry_atm2, &wet_atm1,
                    &wet_atm2);

            geometric_delay[iscan * n_antenna + iant] = delay_vac;
            dry_delay[iscan * n_antenna + iant] = dry_atm2;  // remote station
            wet_delay[iscan * n_antenna + iant] = wet_atm2;  // remote station
        }
    }
}

}  // namespace rtdm::delay_model
