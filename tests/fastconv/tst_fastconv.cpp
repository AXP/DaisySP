

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
#include "daisysp.h"



using namespace daisysp;

#define TARGET_NONE 0
#define TARGET_SEED 1
#define TARGET_POD 2

#ifdef _WIN32
#define TEST_TARGET TARGET_SEED
#else
#define TEST_TARGET TARGET_POD
#endif

static constexpr float ERROR_THRESH_DB = -60.0f;


#if TEST_TARGET == TARGET_SEED
static DaisySeed  seed;
#elif TEST_TARGET == TARGET_POD
static DaisyPod  pod;
static DaisySeed& seed = pod.seed;
#endif


//    N |      FFT  |        FIR
//-------------------------------
//   32 |      662  |      1 442
//   64 |    1 342  |      5 120
//  128 |    2 348  |     19 690
//  256 |    5 592  |     77 514
//  512 |   11 842  |    307 648
// 1024 |   23 366  |  1 227 824
// 2048 |   56 050  |  4 921 126
// 4096 |  126 628  | 21 038 706 


#if 0
static constexpr size_t filter_list[] = 
{
    32,   40,  51, */ 64,  81, 102,  128,  161, 203, 
    256, /*323, 406,*/ 512,/* 645, 813,*/ 1024, /*1290,*/ 1625, 
    2048
};
#else
static constexpr size_t filter_list[] = {
  //  32,
    64,
    128,
    256,
    512,
    1024,
    2048,
   // 4096,
   // 8192,
   // 16384,
  //  32768,
  //  65536
};
#endif
#if 0
static constexpr size_t block_list[] = 
{
    1, 2, 4, 8, 16, 21, 32, 48, 64, 100, 128, 256, 512/*, 1024, 2048 */
};
#else
static constexpr size_t block_list[] = {64};

#endif

static constexpr size_t MAX_IR_LENGTH  = TestPlatform::FindMax(filter_list);
static constexpr size_t MAX_BLOCK_SIZE = FindMax(block_list);
static constexpr size_t TEST_LENGTH   = 65536 * 16;


float DSY_SDRAM_BSS data_in[TEST_LENGTH]   = {0};
float DSY_SDRAM_BSS data_out[TEST_LENGTH]  = {0};
float DSY_SDRAM_BSS data_ref[TEST_LENGTH]  = {0};
float DSY_SDRAM_BSS data_ir[MAX_IR_LENGTH] = {0};


FFTConv<MAX_IR_LENGTH>         DSY_SRAM1_BSS DUT_DYN;
FastConvDyn<64, MAX_IR_LENGTH>  FAST_DYN;


template <typename dut_type>
static uint32_t apply_filter(dut_type& DUT, float* __restrict pFilter, float* __restrict pSrc, float* __restrict pDst, size_t filter_length, size_t signal_length, size_t block_size, size_t& dut_latency)
{
    assert(nullptr != pFilter);
    assert(nullptr != pSrc);
    assert(nullptr != pDst);
    assert(filter_length > 0);
    
    // configure fast convolution impulse response
    const bool init_res = DUT.SetIR(pFilter, filter_length);
    if (false == init_res)
    {
        seed.PrintLine("DUT init result: FAIL");
        return 0;
    }

    uint32_t dt;
    if (block_size == 1)
    {
        ScopedIrqBlocker block;

        const uint32_t t0 = dsy_tim_get_tick();

        for (size_t i = 0; i < signal_length; i++)
        {
            pDst[i] = DUT.Process(pSrc[i]);  // process sample by sample
        }

        dt = dsy_tim_get_tick() - t0;
    }
    else
    {
        ScopedIrqBlocker block;

        const uint32_t t0 = dsy_tim_get_tick();

        while(signal_length >= block_size)
        {
            DUT.ProcessBlock(pSrc, pDst, block_size);    // process whole blocks
            signal_length -= block_size;
            pSrc += block_size;
            pDst += block_size;
        }
        if (signal_length > 0)
        {
            DUT.ProcessBlock(pSrc, pDst, signal_length);    // process whatever is left 
        }

        dt = dsy_tim_get_tick() - t0;
    }
    dut_latency = DUT.GetLatency();
    return dt;
}

template <size_t filter_length>
static bool verify_fastconv_single(float* filter, size_t signal_length, size_t block_size)
{
    assert(signal_length <= DSY_COUNTOF(data_out));
    assert(signal_length <= DSY_COUNTOF(data_in));

    static FIRFilter<filter_length, MAX_BLOCK_SIZE>  REF;
    static FFTConv<get_next_power2(filter_length), true> DSY_SDRAM_BSS DUT_FIX
        = {0};
    static FastConv<64, filter_length> DSY_SDRAM_BSS DUT_FAST
        = {0};

    // regenerate input signal every time for some random variation
    TestUtils::GenerateSignal(data_in, signal_length);
    memset(data_out, 0, sizeof(data_out));
    
    uint32_t per_us = dsy_tim_get_ticks_per_us();
    size_t latency_ref(0), latency_fix(0), latency_dyn(0), latency_fast(0), latency_dfast(0);

    // apply reference FIR filter
    const uint32_t ref_dt = apply_filter(REF, filter, data_in, data_ref, filter_length, signal_length, block_size, latency_ref);

    // apply fast block convolution 
    const uint32_t dyn_dt = apply_filter(DUT_DYN, filter, data_in, data_out, filter_length, signal_length, block_size, latency_dyn);
    const float rms_dyn = TestUtils::CalcMSEdB(data_ref, data_out + latency_dyn, signal_length - latency_dyn);

    // apply fast block convolution (fixed size)
    const uint32_t fix_dt = apply_filter(DUT_FIX, filter, data_in, data_out, filter_length, signal_length, block_size, latency_fix);
    const float rms_fix = TestUtils::CalcMSEdB(data_ref, data_out + latency_fix, signal_length - latency_fix);

    // apply fast partitioned convolution 
    const uint32_t fast_dt = apply_filter(DUT_FAST, filter, data_in, data_out, filter_length, signal_length, block_size, latency_fast);
    const float rms_fast = TestUtils::CalcMSEdB(data_ref, data_out + latency_fast, signal_length - latency_fast);

    // apply fast partitioned convolution (dynamic size)
    const uint32_t dfast_dt = apply_filter(FAST_DYN, filter, data_in, data_out, filter_length, signal_length, block_size, latency_dfast);
    const float    rms_dfast = TestUtils::CalcMSEdB(data_ref, data_out + latency_dfast, signal_length - latency_dfast);

    // verification criteria
    const bool pass =    (rms_dyn   < ERROR_THRESH_DB) 
                      && (rms_fix   < ERROR_THRESH_DB)
                      && (rms_fast  < ERROR_THRESH_DB)
                      && (rms_dfast < ERROR_THRESH_DB);

    // produce human-readable forms
    const float  per_scale = 1.0f / (signal_length * per_us);
    const float   ref_time = per_scale * ref_dt;
    const float   dyn_time = per_scale * dyn_dt;
    const float   fix_time = per_scale * fix_dt;
    const float  fast_time = per_scale * fast_dt;
    const float dfast_time = per_scale * dfast_dt;

    seed.PrintLine("filter=%5u; latency=%5u; block=%4u; "\
                                            "RMS[dB]; Fix=" FLT_FMT(2) "; Dyn=" FLT_FMT(2) "; Fast=" FLT_FMT(2) "; DynFast=" FLT_FMT(2) "; " \
                   "Time[us/smp]; Ref=" FLT_FMT(3) "; Fix=" FLT_FMT(3) "; Dyn=" FLT_FMT(3) "; Fast=" FLT_FMT(3) "; DynFast=" FLT_FMT(3) "; %s", 
        filter_length, latency_fix, block_size, FLT_VAR(2, rms_fix),  FLT_VAR(2, rms_dyn),  FLT_VAR(2, rms_fast),  FLT_VAR(2, rms_dfast),
                          FLT_VAR(3, ref_time), FLT_VAR(3, fix_time), FLT_VAR(3, dyn_time), FLT_VAR(3, fast_time), FLT_VAR(3, dfast_time),
            TestUtils::ResultStr(pass));

    return pass;
}


template<size_t i>
class verifier
{
    public:

    bool operator()()
    {
        bool result = true;
        for (size_t j = 0; j < DSY_COUNTOF(block_list); j++)
        {
            result &= verify_fastconv_single<filter_list[i]>(data_ir, DSY_COUNTOF(data_in) - 1, block_list[j]);
        }
        return result;
    }
};



static bool verify_fastconv()
{
    bool pass = true;
    TestUtils::GenerateSignal(data_ir, DSY_COUNTOF(data_ir));

    static_for<0, DSY_COUNTOF(filter_list)> loop;
    pass &= loop.go_bool<verifier>();

    seed.PrintLine("Done: %s", TestUtils::ResultStr(pass));

    return pass;
}
#if 0

#define FFT_LEN(_i) (32 << (_i))

static void test_fft_perf()
{
    arm_rfft_fast_instance_f32 fft[8];
    arm_fir_instance_f32 fir[8];
    FFTReal<float> fftr[8] = {
        FFTReal<float>(FFT_LEN(0)), 
        FFTReal<float>(FFT_LEN(1)),
        FFTReal<float>(FFT_LEN(2)),
        FFTReal<float>(FFT_LEN(3)),
        FFTReal<float>(FFT_LEN(4)),
        FFTReal<float>(FFT_LEN(5)),
        FFTReal<float>(FFT_LEN(6)),
        FFTReal<float>(FFT_LEN(7)),
        };

    
    for (int i = 0; i < 8; i++)
    {
        const uint32_t t0 = dsy_tim_get_tick();
        const arm_status res = arm_rfft_fast_init_f32(&fft[i], FFT_LEN(i));
        const uint32_t dt = dsy_tim_get_tick() - t0;

        arm_fir_init_f32(&fir[i], FFT_LEN(i), data_ir, data_tmp, FFT_LEN(i));
        seed.PrintLine("arm_rfft_fast_init_f32(%4d) = %u ticks = %d", FFT_LEN(i), dt, (int)res);
        //dsy_tim_delay_ms(100);
    }

    FFTRealFixLen<5> fft32;
    FFTRealFixLen<6> fft64;
    FFTRealFixLen<7> fft128;
    FFTRealFixLen<8> fft256;
    FFTRealFixLen<9> fft512;
    FFTRealFixLen<10> fft1024;
    FFTRealFixLen<11> fft2048;
    FFTRealFixLen<12> fft4096;

    while(1) 
    {
        for (int i = 0; i < 8; i++)
        {
            const uint32_t t0 = dsy_tim_get_tick();
            arm_rfft_fast_f32(&fft[i], data_in, data_out, 0);
            const uint32_t dt = dsy_tim_get_tick() - t0;
            
            const uint32_t t1 = dsy_tim_get_tick();
            arm_fir_f32(&fir[i], data_in, data_out, FFT_LEN(i));
            const uint32_t dt2 = dsy_tim_get_tick() - t1;
            
            const uint32_t t2 = dsy_tim_get_tick();
            fftr[i].do_fft(data_in, data_out);
            const uint32_t dt3 = dsy_tim_get_tick() - t2;

            seed.PrintLine("arm_rfft_fast_f32(%4d) = %6u ticks; arm_fir_f32() = %6u tocks; fftreal() = %6u pucks", FFT_LEN(i), dt, dt2, dt3);
            dsy_tim_delay_ms(100);
        }

        uint32_t ddt0 = dsy_tim_get_tick();
        fft32.do_fft(data_in, data_out);
        ddt0 = dsy_tim_get_tick() - ddt0;

        uint32_t ddt1 = dsy_tim_get_tick();
        fft64.do_fft(data_in, data_out);
        ddt1 = dsy_tim_get_tick() - ddt1;

        uint32_t ddt2 = dsy_tim_get_tick();
        fft128.do_fft(data_in, data_out);
        ddt2 = dsy_tim_get_tick() - ddt2;

        uint32_t ddt3 = dsy_tim_get_tick();
        fft256.do_fft(data_in, data_out);
        ddt3 = dsy_tim_get_tick() - ddt3;

        uint32_t ddt4 = dsy_tim_get_tick();
        fft512.do_fft(data_in, data_out);
        ddt4 = dsy_tim_get_tick() - ddt4;

        uint32_t ddt5 = dsy_tim_get_tick();
        fft1024.do_fft(data_in, data_out);
        ddt5 = dsy_tim_get_tick() - ddt5;

        uint32_t ddt6 = dsy_tim_get_tick();
        fft2048.do_fft(data_in, data_out);
        ddt6 = dsy_tim_get_tick() - ddt6;

        uint32_t ddt7 = dsy_tim_get_tick();
        fft4096.do_fft(data_in, data_out);
        ddt7 = dsy_tim_get_tick() - ddt7;

        seed.PrintLine("%u; %u; %u; %u; %u; %u; %u; %u", ddt0, ddt1, ddt2, ddt3, ddt4, ddt5, ddt6, ddt7);

    }
}

#endif
#if 0

static float ir[128];
static void test_fastconv()
{
    for (size_t i = 0; i < DSY_COUNTOF(ir); i++)    
    {
        ir[i] = (float)(1 + i);
    }

    seed.PrintLine("IR:");
    for (size_t i = 0; i < DSY_COUNTOF(ir); i++)
    {
        seed.PrintLine("%2d; " FLT_FMT3, i, FLT_VAR3(ir[i]));
    }

    FastConv<64, 128> FC;
    FC.SetIR(ir, DSY_COUNTOF(ir));
    

    for (size_t i = 0; i < 512; i++)
    {
        float in = 0.0f;
        if (i == 100)
        {
            in = 1.0f;
        }
        float out = FC.Process(in);
        seed.PrintLine("%4u; in = " FLT_FMT3 "; out = " FLT_FMT3, i, FLT_VAR3(in), FLT_VAR3(out));
    }

}
#endif


#if 0
static FFTConvSize<2048, true>  DSY_SDRAM_BSS   FCS_fix_2048;
static FFTConvSize<500, true>    DSY_SDRAM_BSS  FCS_fix_500;
static FFTConvSize<1024, false>  DSY_SDRAM_BSS  FCS_dyn_1024;
static FFTConvSize<60, false>   DSY_SDRAM_BSS   FCS_dyn_60;

static bool verify_size_helper()
{
    constexpr size_t new_size = 100;

    bool test_pass = true;

    FCS_dyn_1024.Resize();
    FCS_dyn_60.Resize();

    seed.PrintLine("FCS_fix_2048; GetLen() = %4u; GetFFTLen() = %4u", FCS_fix_2048.GetLen(), FCS_fix_2048.GetFFTLen());
    seed.PrintLine("FCS_fix_500;  GetLen() = %4u; GetFFTLen() = %4u", FCS_fix_500.GetLen(),  FCS_fix_500.GetFFTLen());
    seed.PrintLine("FCS_dyn_1024; GetLen() = %4u; GetFFTLen() = %4u", FCS_dyn_1024.GetLen(), FCS_dyn_1024.GetFFTLen());
    seed.PrintLine("FCS_dyn_60;   GetLen() = %4u; GetFFTLen() = %4u", FCS_dyn_60.GetLen(),   FCS_dyn_60.GetFFTLen());
    
    test_pass &= FCS_fix_2048.GetLen() == 2048;
    test_pass &= FCS_fix_2048.GetFFTLen() == 4096;
    test_pass &= FCS_fix_500.GetLen() == 512;
    test_pass &= FCS_fix_500.GetFFTLen() == 1024;
    test_pass &= FCS_dyn_1024.GetLen()  == 1024;
    test_pass &= FCS_dyn_1024.GetFFTLen() == 2048;
    test_pass &= FCS_dyn_60.GetLen()  == 64;
    test_pass &= FCS_dyn_60.GetFFTLen() == 128;
    seed.PrintLine("result = %s", TestUtils::ResultStr(pass));
    
    seed.PrintLine("Resize(%u)", new_size);
    FCS_fix_2048.Resize(new_size);
    FCS_fix_500.Resize(new_size);
    FCS_dyn_1024.Resize(new_size);
    FCS_dyn_60.Resize(new_size);

    seed.PrintLine("FCS_fix_2048; GetLen() = %4u; GetFFTLen() = %4u", FCS_fix_2048.GetLen(), FCS_fix_2048.GetFFTLen());
    seed.PrintLine("FCS_fix_500;  GetLen() = %4u; GetFFTLen() = %4u", FCS_fix_500.GetLen(),  FCS_fix_500.GetFFTLen());
    seed.PrintLine("FCS_dyn_1024; GetLen() = %4u; GetFFTLen() = %4u", FCS_dyn_1024.GetLen(), FCS_dyn_1024.GetFFTLen());
    seed.PrintLine("FCS_dyn_60;   GetLen() = %4u; GetFFTLen() = %4u", FCS_dyn_60.GetLen(),   FCS_dyn_60.GetFFTLen());

    test_pass &= FCS_fix_2048.GetLen() == 2048;
    test_pass &= FCS_fix_2048.GetFFTLen() == 4096;
    test_pass &= FCS_fix_500.GetLen() == 512;
    test_pass &= FCS_fix_500.GetFFTLen() == 1024;
    test_pass &= FCS_dyn_1024.GetLen()  == 128;
    test_pass &= FCS_dyn_1024.GetFFTLen() == 256;
    test_pass &= FCS_dyn_60.GetLen()  == 64;
    test_pass &= FCS_dyn_60.GetFFTLen() == 128;

    seed.PrintLine("result = %s", TestUtils::ResultStr(pass));

    return test_pass;
}
#endif

#if 0
template<size_t... ints>
void print_sequence(std::index_sequence<ints...> int_seq)
{   
    //std::ostringstream  ss;
    //ss << "The sequence of size " << int_seq.size() << ": ";
    //((ss << ints << ' '), ...);
    //seed.PrintLine("%s", ss.str());
    seed.PrintLine("Sequence size: %u", int_seq.size());
    (seed.Print("%u ", ints), ...);
    seed.PrintLine("");
    size_t total = 0;
    ((total += ints), ...);
    seed.PrintLine("Sum = %u", total);
}

#endif

template<size_t size>
void PrintMemory()
{
    seed.PrintLine("%6u;%8u;%8u;%8u;%8u;%8u;%8u;%8u;%8u;%9u;%9u",
        size,
        sizeof(FIRFilter<size, 64>),
        sizeof(FIRFilterImplGeneric<size, 64>),
        sizeof(FFT<size, true>),
        sizeof(FFT<size, false>),
        sizeof(FFT_Generic<size, true>),
        sizeof(FFT_Generic<size, false>),
        sizeof(FFTConv<size, true>),
        sizeof(FFTConv<size, false>),
        sizeof(FastConv<64, size>),
        sizeof(FastConvDyn<64, size>));
}

template<size_t i>
struct PrintMemoryLoop
{
    bool operator()()
    {
        PrintMemory<1U << i>();
        return true;
    }
};

void PrintMemoryConsumption()
{

seed.PrintLine("%6s;%8s;%8s;%8s;%8s;%8s;%8s;%8s;%8s;%9s;%9s",
        "IR",
        "FIR ARM",
        "FIR Gen",
        "FFT ARM",
        "FFT ARM",
        "FFT Gen",
        "FFT Gen",
        "FFTConv",
        "FFTConv",
        "FastConv",
        "FastConv");
seed.PrintLine("%6s;%8s;%8s;%8s;%8s;%8s;%8s;%8s;%8s;%9s;%9s",
        "size",
        "",
        "",
        "Fix",
        "Dyn",
        "Fix",
        "Dyn",
        "Fix",
        "Dyn",
        "Fix",
        "Dyn");

    static_for<6, 17> loop;
    loop.go<PrintMemoryLoop>();

}

int main(void)
{
    //sizeof(FastConv<64, 6*48000>)
    //sizeof(FastConv<64, 2048>)
    // initialize seed hardware and daisysp modules
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
    
    //test_fastconv(); return 0;

    PrintMemoryConsumption(); //return 0;
    bool res = verify_fastconv();

 //   test_fft_perf();

  //  FIRFilter<64, 64> fir;

    //fir.SetIR(data_ir, 64);
    #if 0
    for (size_t i = 0; i < DSY_COUNTOF(data_ir); i++)
    {
        data_ir[i] = (float)i;
    }

    memset(data_in, 0, sizeof(data_in));
    memset(data_out, 0, sizeof(data_out));
    data_in[10] = 1.0f;
 
    //bool res = verify_size_helper();
#if 1
//    seed.PrintLine("before SetIR");
    FC.SetIR(data_ir, 2048);  
  //  seed.PrintLine("after SetIR");
    uint32_t t0;
    {
        ScopedIrqBlocker block;
    
        t0 = dsy_tim_get_tick();
        for (size_t i = 0; i < 2500; i++)  
        { 
            data_out[i] = FC.Process(data_in[i]);
        }
        t0 = dsy_tim_get_tick() - t0;
    }


    for (size_t i = 0; i < 2500; i++)  
    { 
        seed.PrintLine("%4u; " FLT_FMT3, i, FLT_VAR3(data_out[i]));
    }

    seed.PrintLine("ticks = %u", t0/2500);
    // ticks = 795

#endif

#endif

    
#if 0
FastConv<64, 4096>::PartPlan PP;
    size_t q = PP.PLAN_FFT_PART;
    size_t w = PP.PLAN_TOTAL;
    PP.FFT_INSTANCES;

    seed.PrintLine ("PP.IR_FFT_PART = %u", q);
    seed.PrintLine ("PP.IR_TOTAL    = %u", w);

    auto a = FC.getFFT<0>();
    auto b = FC.getFFT<1>(); 
#endif    
#if 0
    for (size_t i = 0; i < PP.StageCount(); i++) 
    {
        seed.PrintLine("Stage %u: %u", i, PP.stages_[i]);
    }

    seed.PrintLine("index_lin_:");
    print_sequence(PP.index_lin_);

    seed.PrintLine("index_:");
    print_sequence(PP.index_);
#endif

#if TEST_TARGET == TARGET_SEED

    seed.SetLed(false);
    (void)res;

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
