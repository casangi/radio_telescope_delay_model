// C++ driver over the vendored difxcalc11 pipeline -- see difxcalc_core.hpp.
//
// This mirrors difxcalc11's dmain.f verbatim minus the .im output stage:
//   C_mode = 'difx  '; set file names; CALL dSTART; CALL dINITL;
//   per scan: CALL dSCAN; per 2-minute interval: CALL dDRIVR;
// then reads the per-epoch sample arrays from COMMON /OUT_C/ (c2poly.i)
// instead of fitting polynomials into a .im file. Every struct below
// mirrors the declared COMMON order of its include file exactly; offsets
// are static_asserted against the Fortran sequence-storage layout.

#include "difxcalc_core.hpp"

#include <cstddef>
#include <cstring>
#include <mutex>
#include <stdexcept>

namespace {

// --- COMMON /Contrl/ of d_input.i ---------------------------------------
extern "C" struct Contrl {
    double d_interval;
    int icalc, ijob, i_out, im_out, verbose, epoch2m;
    char calc_file_name[128], im_file_name[128], jobname[128],
        calc_out_file[128];
    char base_mode[10], l_time[10], atmdr[10], atmwt[10];
    int numjobs;
    char near_far[10];
    char overwrite[4];
    char uvw[6];
} contrl_;
static_assert(offsetof(Contrl, epoch2m) == 28);
static_assert(offsetof(Contrl, calc_file_name) == 32);
static_assert(offsetof(Contrl, base_mode) == 544);
static_assert(offsetof(Contrl, numjobs) == 584);
static_assert(offsetof(Contrl, uvw) == 602);

// --- Leading members of COMMON /Calc_input/ of d_input.i ------------------
extern "C" struct CalcInput {
    double xleap_sec;
    int jobid, numscans, numsrc, numepochs, numphcntr, pointingsrc,
        phcntrnum;
    int phcntr[500];
    int scannum, intrvls2min;
} calc_input_;
static_assert(offsetof(CalcInput, phcntr) == 36);
static_assert(offsetof(CalcInput, intrvls2min) == 2040);

// --- COMMON /mode/ of cuser11.i: CHARACTER*6 C_mode ----------------------
extern "C" char mode_[6];

// --- COMMON /RTDMPATHS/ of rtdm_paths.i ----------------------------------
extern "C" struct {
    char ocean[128], tilt[128], optl[128], de421[128], leap[128];
} rtdmpaths_;

// --- Leading members of COMMON /OUT_C/ of c2poly.i ------------------------
// Fortran: X(Max_Epochs, Nstation1, Nstation2, Max_Source) with
// Max_Epochs=6, Nstation1=1, Nstation2=254, Max_Source=1001 (the Atm*
// arrays carry a leading IB=1 dimension, which does not change the layout).
constexpr int kEp = rtdm::difxcalc::kMaxEpochsPerInterval;
constexpr int kSta = 254;
constexpr int kSrc = 1001;
using OutArray = double[kSrc][kSta][1][kEp];
extern "C" struct OutC {
    OutArray delay, rate, atmdryd, atmdryr, atmwetd, atmwetr;
    OutArray ubase, vbase, wbase;
    OutArray el, az;
    OutArray stax, stay, staz, staxt, stayt, stazt;
    int iymdhms[6][kEp];  // Iymdhms_f(Max_Epochs, 6)
    int numsite, numbaseline, numphcenter;
} out_c_;
static_assert(sizeof(OutArray) == 8ull * kEp * kSta * kSrc);
static_assert(offsetof(OutC, iymdhms) == 17ull * sizeof(OutArray));

// difxcalc11 entry points (gfortran name mangling).
extern "C" void dstart_(const int* num_scans, const int* kjob);
extern "C" void dinitl_(const int* kjob);
extern "C" void dscan_(const int* iscan, const int* kjob);
extern "C" void ddrivr_(const int* iscan, const int* interval);

std::mutex core_mutex;  // COMMON-block state: one computation at a time

void set_blank_padded(char* field, std::size_t size, const std::string& value,
                      const char* name) {
    if (value.size() > size)
        throw std::invalid_argument(std::string(name) + " exceeds " +
                                    std::to_string(size) + " characters: " +
                                    value);
    std::memset(field, ' ', size);
    std::memcpy(field, value.data(), value.size());
}

}  // namespace

namespace rtdm::difxcalc {

std::size_t run_difxcalc11(
    const std::string& calc_file, const std::string& ocean_loading_file,
    const std::string& tilt_file, const std::string& ocean_pole_tide_file,
    const std::string& de421_file, const std::string& leap_second_file,
    std::size_t n_antenna, std::size_t max_samples, int* ymdhms, double* delay,
    double* dry, double* wet, double* u, double* v, double* w) {
    if (n_antenna == 0 || n_antenna > kSta)
        throw std::invalid_argument("n_antenna out of range for difxcalc11.");

    std::lock_guard<std::mutex> lock(core_mutex);

    // The defaults difxcalc's command-line parser (GETCL in dstrt.f) would
    // set; we bypass GETCL, so replicate them exactly. Base_mode, UVW and
    // the Atm switches select model branches in dSCAN/dDRIVR.
    std::memcpy(mode_, "difx  ", 6);
    set_blank_padded(contrl_.calc_file_name, 128, calc_file, "calc_file");
    set_blank_padded(contrl_.im_file_name, 128, "", "im_file");
    set_blank_padded(contrl_.jobname, 128, "", "jobname");
    set_blank_padded(contrl_.calc_out_file, 128, "", "calc_out_file");
    set_blank_padded(contrl_.base_mode, 10, "geocenter", "base_mode");
    set_blank_padded(contrl_.l_time, 10, "dont-solve", "l_time");
    set_blank_padded(contrl_.atmdr, 10, "Add-dry", "atmdr");
    set_blank_padded(contrl_.atmwt, 10, "Add-wet", "atmwt");
    set_blank_padded(contrl_.near_far, 10, "Far-field", "near_far");
    set_blank_padded(contrl_.overwrite, 4, "no", "overwrite");
    set_blank_padded(contrl_.uvw, 6, "exact", "uvw");
    contrl_.d_interval = 24.0;
    contrl_.epoch2m = static_cast<int>(120.0001 / contrl_.d_interval) + 1;
    contrl_.verbose = 0;
    contrl_.i_out = 0;
    contrl_.im_out = 0;
    contrl_.numjobs = 0;
    set_blank_padded(rtdmpaths_.ocean, 128, ocean_loading_file, "ocean file");
    set_blank_padded(rtdmpaths_.tilt, 128, tilt_file, "tilt file");
    set_blank_padded(rtdmpaths_.optl, 128, ocean_pole_tide_file, "optl file");
    set_blank_padded(rtdmpaths_.de421, 128, de421_file, "de421 file");
    set_blank_padded(rtdmpaths_.leap, 128, leap_second_file, "leap file");

    int num_scans = 0;
    const int kjob = 1;
    dstart_(&num_scans, &kjob);
    dinitl_(&kjob);

    std::size_t sample = 0;
    for (int scan = 1; scan <= num_scans; ++scan) {
        dscan_(&scan, &kjob);
        const int intervals = calc_input_.intrvls2min;
        for (int interval = 1; interval <= intervals; ++interval) {
            ddrivr_(&scan, &interval);
            const int epochs = contrl_.epoch2m;
            if (epochs < 1 || epochs > kEp)
                throw std::runtime_error("difxcalc11: unexpected epoch2m.");
            for (int e = 0; e < epochs; ++e) {
                if (sample >= max_samples)
                    throw std::runtime_error(
                        "difxcalc11: sample buffer too small.");
                for (int c = 0; c < 6; ++c)
                    ymdhms[sample * 6 + c] = out_c_.iymdhms[c][e];
                for (std::size_t a = 0; a < n_antenna; ++a) {
                    // Fortran indices: Delay_f(e+1, 1, a+1, 1). Values are
                    // raw model samples: delays in CALC-sign seconds (d_out
                    // applies the -1e6 DiFX/.im conversion at write time),
                    // atmospheres in seconds, U/V/W in metres.
                    const std::size_t out = sample * n_antenna + a;
                    delay[out] = out_c_.delay[0][a][0][e];
                    dry[out] = out_c_.atmdryd[0][a][0][e];
                    wet[out] = out_c_.atmwetd[0][a][0][e];
                    u[out] = out_c_.ubase[0][a][0][e];
                    v[out] = out_c_.vbase[0][a][0][e];
                    w[out] = out_c_.wbase[0][a][0][e];
                }
                ++sample;
            }
        }
    }
    return sample;
}

}  // namespace rtdm::difxcalc
