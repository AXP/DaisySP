#include <functional>

#include "daisysp.h"
#if defined(__arm__)
#include "daisy_seed.h"
#include "daisy_pod.h"
#include "util/scopedirqblocker.h"
using namespace daisy;
#define DSY_SRAM1_BSS __attribute__((section(".sram1_bss")))
#elif defined(_WIN32)
#include <Windows.h>
#define DSY_SDRAM_BSS
#define DSY_SRAM1_BSS
#include <cstdint>

class ScopedIrqBlocker
{
public:
    ScopedIrqBlocker() { (void)0; }
    ~ScopedIrqBlocker() { (void)0; }
};
uint32_t dsy_tim_get_tick()
{
    LARGE_INTEGER q;
    QueryPerformanceCounter(&q);
    return q.LowPart;
}
uint32_t dsy_tim_get_ticks_per_us()
{
    LARGE_INTEGER q;
    QueryPerformanceFrequency(&q);
    uint32_t us = (uint32_t)((q.QuadPart + 500000) / 1000000);
    return us;
}
#include "../libdaisy/src/hid/logger.h"

class DaisySeed
    : public daisy::Logger<daisy::LoggerDestination::LOGGER_SEMIHOST>
{
public:
    void Configure() {}
    void Init() {}
    void SetLed(bool value) { printf("SetLed -> %d\n", value); }
};
#endif
#include "test_util.h"


using namespace daisysp;
using namespace daisy;

#define TARGET_SEED 0
#define TARGET_POD 1

#define TEST_TARGET TARGET_SEED   
static constexpr float ERROR_THRESH_DB = -138.0f;

#if TEST_TARGET == TARGET_SEED
static DaisySeed  seed;
#elif TEST_TARGET == TARGET_POD
static DaisyPod  pod;
static DaisySeed& seed = pod.seed;
#endif

static constexpr int BUF_SIZE = FFT<>::FFT_MAXSIZE;

static float DSY_SDRAM_BSS data_in[BUF_SIZE+1];
static float DSY_SDRAM_BSS data_tmp[BUF_SIZE+1];
static float DSY_SDRAM_BSS data_tmp2[BUF_SIZE+1];
static float DSY_SDRAM_BSS data_out[BUF_SIZE+1];


template <typename dut_type>
static uint32_t apply(dut_type& dut, float* __restrict pSrc, float* __restrict pDst, int repeats)
{
    uint32_t dt_min(UINT32_MAX);
    for (int i = 0; i < repeats; i++)
    {
        // disable interrupts for the duration of measurements
        ScopedIrqBlocker block;

        /* make a local copy to avoid data corruption by in-place FFT */
        memcpy(data_tmp, pSrc, dut.GetSize() * sizeof(data_tmp[0]));

        memset(pDst, 0, dut.GetSize() * sizeof(pDst[0]));
        
        const uint32_t t0 = dsy_tim_get_tick();

        dut.ForwardTransform(data_tmp, data_tmp2);
        dut.InverseTransform(data_tmp2, pDst);

        const uint32_t dt = dsy_tim_get_tick() - t0;

        dt_min = DSY_MIN(dt, dt_min);
    }

    return dt_min;
}

static FFT<4096, false> DUT_DYN;
static FFT_Generic<4096, false> GEN_DYN;

template <size_t fft_size>
static bool verify_fft_single()
{
    assert(fft_size <= DSY_COUNTOF(data_in));
    assert(fft_size <= DSY_COUNTOF(data_out));

    static FFT<fft_size, true> DUT_FIX;
    static FFT_Generic<fft_size, true> GEN_FIX;
    if (false == DUT_DYN.SetSize(fft_size))
    {
        seed.PrintLine("Error setting DUT_DYN size = %u", fft_size);
        return false;
    }
    if (false == GEN_DYN.SetSize(fft_size))
    {
        seed.PrintLine("Error setting GEN_DYN size = %u", fft_size);
        return false;
    }

    // regenerate input signal every time for some random variation
    TestUtils::GenerateSignal(data_in, fft_size);
    
    
    constexpr int repeats = 800000000 / (fft_size * get_base_2(fft_size));

    const uint32_t     dt_fix = apply(DUT_FIX, data_in, data_out, repeats);
    const float       rms_fix = TestUtils::CalcMSEdB(data_in, data_out, fft_size);
    const uint32_t     dt_dyn = apply(DUT_DYN, data_in, data_out, repeats);
    const float       rms_dyn = TestUtils::CalcMSEdB(data_in, data_out, fft_size);
    const uint32_t dt_gen_fix = apply(GEN_FIX, data_in, data_out, repeats);
    const float   rms_gen_fix = TestUtils::CalcMSEdB(data_in, data_out, fft_size);
    const uint32_t dt_gen_dyn = apply(GEN_DYN, data_in, data_out, repeats);
    const float   rms_gen_dyn = TestUtils::CalcMSEdB(data_in, data_out, fft_size);

    // verification criteria
    const bool pass =  (rms_fix < ERROR_THRESH_DB) 
                    && (rms_dyn < ERROR_THRESH_DB) 
                    && (rms_gen_fix < ERROR_THRESH_DB) 
                    && (rms_gen_dyn < ERROR_THRESH_DB);

    // produce human-readable forms
    uint32_t us_ticks = dsy_tim_get_ticks_per_us();
    const float fix_time     = (float)dt_fix / us_ticks; // 200 MHz timer
    const float dyn_time     = (float)dt_dyn / us_ticks; // 200 MHz timer
    const float gen_fix_time = (float)dt_gen_fix / us_ticks; // 200 MHz timer
    const float gen_dyn_time = (float)dt_gen_dyn / us_ticks; // 200 MHz timer

    seed.PrintLine("fft_size = %4u; RMS[dB] Fix = " FLT_FMT(3) "; Dyn = " FLT_FMT(3) "; Gen_fix = " FLT_FMT(3) "; Gen_dyn = " FLT_FMT(3) "; "\
                                  "Time[us] Fix = " FLT_FMT(2) "; Dyn = " FLT_FMT(2) "; Gen_fix = " FLT_FMT(2) "; Gen_dyn = " FLT_FMT(2) "; %s", 
        fft_size, FLT_VAR3(rms_fix), FLT_VAR3(rms_dyn), FLT_VAR3(rms_gen_fix), FLT_VAR3(rms_gen_dyn),
        FLT_VAR(2, fix_time), FLT_VAR(2, dyn_time), FLT_VAR(2, gen_fix_time), FLT_VAR(2, gen_dyn_time), TestUtils::ResultStr(pass));

    return pass;
}

static constexpr size_t fft_list[] = 
{
    32, 64, 128, 256, 512, 1024, 2048, 4096
};


template<size_t i>
class verifier
{
    public:

    bool operator()()
    {
        if (fft_list[i] <= BUF_SIZE)
        {
            return verify_fft_single<fft_list[i]>();
        }
        else
        {
            return false;
        }
    }
};


static bool verify_fft()
{
    static_for<0, DSY_COUNTOF(fft_list)> loop;
    const bool pass = loop.go_bool<verifier>();

    seed.PrintLine("Done: %s", TestUtils::ResultStr(pass));

    return pass;
}

#if 0

class DFT
{
    
    public:

        DFT()
            : size_(0)
            , scale_(0)
        {

        }
        void SetSize(size_t size)
        {
            size_ = size;
            scale_ = 6.2831853071795864 / size; // let's be more accurate than TWOPI_F
        }

        void ForwardTransform(float* __restrict buf_td, float* __restrict buf_fd)
        {
            for (size_t i = 0; i < size_/2; i++)
            {
                double acc_re(0);
                double acc_im(0);
                for (size_t j = 0; j < size_; j++)
                {
                    acc_re += cos(scale_ * i * j) * buf_td[j];
                    acc_im += sin(scale_ * i * j) * buf_td[j];
                }
                buf_fd[i*2]   = (float)acc_re;
                buf_fd[i*2+1] = (float)acc_im;
            }
            double acc_re(0);
            for (size_t j = 0; j < size_; j++)
            {
                acc_re += cos(scale_ * size_/2 * j) * buf_td[j];
            }
            buf_fd[1] = acc_re; // Nyquist freq
        }

        void InverseTransform(float* __restrict buf_fd, float* __restrict buf_td)
        {

        }
    protected:
        size_t size_;
        double scale_;
};


static void test_fft()
{
    TestUtils::GenerateSignal(data_in, 32);
    for (size_t i = 0; i < 32; i++)
    {
        //data_in[i] = sinf((TWOPI_F / 8) * (0.5f + i));
    }

    DFT ref;
    ref.SetSize(32);
    ref.ForwardTransform(data_in, data_ref);

    FFT_Generic<32> fftr;
    fftr.SetSize(32);
    fftr.ForwardTransform(data_in, data_dyn);

    FFT<32, true> fix;
    fix.ForwardTransform(data_in, data_fix);

//    FFT<1, false> dyn;
  //  dyn.SetSize(32);
    //dyn.ForwardTransform(data_in, data_dyn);

    for (size_t i = 0; i < 32; i++)
    {
        seed.PrintLine("%2u; " FLT_FMT3 "; " FLT_FMT3 "; " FLT_FMT3, i, FLT_VAR3(data_ref[i]), FLT_VAR3(data_fix[i]), FLT_VAR3(data_dyn[i]));
    }
}

#endif
int main(void)
{
    
    // initialize daisy hardware
#if TEST_TARGET == TARGET_SEED
    seed.Configure();
    seed.Init();
    seed.SetLed(true);
#elif TEST_TARGET == TARGET_POD
    Color led_ko, led_ok, led_start;
    led_ko.Init(Color::RED);
    led_ok.Init(Color::GREEN);
    led_start.Init(Color::WHITE);

    pod.Init();
    pod.led1.SetColor(led_start);
    pod.led1.Update();
#endif    

    seed.StartLog(true);

    const bool res = verify_fft();
   // test_fft(); bool res = true;


#if TEST_TARGET == TARGET_SEED
    seed.SetLed(false);
#elif TEST_TARGET == TARGET_POD
    if (res)
    {
        pod.led1.SetColor(led_ok);
    }
    else
    {
        pod.led1.SetColor(led_ko);
    }
    while(1)
    {
        pod.led1.Update();
    }
#endif
}
