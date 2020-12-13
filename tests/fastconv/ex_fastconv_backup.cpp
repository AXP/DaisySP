#include "daisysp.h"
#include "daisy_seed.h"
#include <stm32h750xx.h>
#include <arm_math.h>
#include "util/logger.h"
#include "util/scopedirqblocker.h"

using namespace daisysp;
using namespace daisy;

static DaisySeed  seed;
static axp::Logger DSY_SDRAM_BSS Log;
// static void AudioCallback(float *in, float *out, size_t size)
// {
//     float saw, freq, output;
//     for(size_t i = 0; i < size; i += 2)
//     {
//         freq = 2500 + (lfo.Process() * 2500);
//         saw  = osc.Process();

//         flt.SetFreq(freq);
//         output = flt.Process(saw);

//         // left out
//         out[i] = output;

//         // right out
//         out[i + 1] = output;
//     }
// }

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












static float DSY_SDRAM_BSS data_in[65536];
static float DSY_SDRAM_BSS data_out[65536];
static float DSY_SDRAM_BSS data_ref[65536];

static float DSY_SDRAM_BSS data_ir[4096];

static float DSY_SDRAM_BSS data_tmp[4096 * 1];

FastConv<4096> DSY_SDRAM_BSS DUT;

static char* str_res[] = {"FAIL", "PASS"};
#define FFT_LEN(_i) (32 << (_i))


static void generate_signal(float* buf, size_t length)
{
    
    assert(nullptr != buf);

    if (nullptr != buf)
    {
        WhiteNoise nse;
        nse.Init();
        nse.SetAmp(1.0f);
        
        for (size_t i = 0; i < length; i++)
        {
            buf[i] = nse.Process();
        }
    }
}

static bool verify_fastconv_single(float* filter, size_t filter_length, size_t signal_length)
{
    const size_t fir_block = 256;   // used only for reference FIR processing

    assert(DSY_COUNTOF(data_in) <= DSY_COUNTOF(data_out));
    assert(DSY_COUNTOF(data_in) <= DSY_COUNTOF(data_ref));
    assert(signal_length <= DSY_COUNTOF(data_in));
    assert(filter_length + fir_block - 1 <= DSY_COUNTOF(data_tmp));

    // regenerate input signal every time for some random variation
    generate_signal(data_in, signal_length);

    // configure fast convolution impulse response
    DUT.SetIR(filter, filter_length);

    // ARM FIR implementation requires filter coefficients in reverse order
    for (size_t i = 0; i < filter_length/2; i++)
    {
        const float tmp = filter[i];
        filter[i] = filter[filter_length - 1 - i];
        filter[filter_length - 1 - i] = tmp;
    }

    // configure reference FIR filter
    arm_fir_instance_f32 REF;
    arm_fir_init_f32(&REF, filter_length, filter, data_tmp, fir_block);


    uint32_t ref_dt, dt;
    // disable interrupts for the duration of measurements
    {
        ScopedIrqBlocker block;

        // apply reference FIR filter
        const uint32_t ref_t0 = dsy_tim_get_tick();
        size_t i;
        for (i = 0; i < signal_length; i += fir_block)
        {
            arm_fir_f32(&REF, &data_in[i], &data_ref[i], fir_block);    // process whole blocks
        }
        if (i < signal_length)
        {
            arm_fir_f32(&REF, &data_in[i], &data_ref[i], signal_length - i);    // process whatever is left over
        }
        ref_dt = dsy_tim_get_tick() - ref_t0;
    

        // apply fast convolution (process several times to find the minimal processing time due to IRQ variance)
        const uint32_t t0 = dsy_tim_get_tick();
        for (size_t i = 0; i < signal_length; i++)
        {
            data_out[i] = DUT.Process(data_in[i]);  // process sample by sample
        }
        dt = dsy_tim_get_tick() - t0;
    }

    // calculate error between reference and fast algorithms
    
    const size_t latency = DUT.GetLatency(); // take fast convolution latency into account
    
    double sum_error(0);
    for (size_t i = 0; i < signal_length - latency; i++)
    {
        const double error = data_out[i + latency] - data_ref[i];
        sum_error += error * error;
    }
    sum_error /= signal_length - latency;

    // obtain human readable forms
    const float rms = 10.0f * log10f((float)sum_error); // square root operation is factored in
    const bool pass = rms < -97.0f;
    const float dut_time = (float)dt / (signal_length * 200); // 200 MHz timer
    const float ref_time = (float)ref_dt / (signal_length * 200); // 200 MHz timer

    Log.PrintF("filter %4u; signal %5u; latency %4u; RMS[dB] = " FLT_FMT3 "; Time[us/smp] Ref = " FLT_FMT(2) " DUT = " FLT_FMT(2) "; %s", 
        filter_length, signal_length, latency, FLT_VAR3(rms), FLT_VAR(2, ref_time), FLT_VAR(2, dut_time), str_res[pass]);
    return pass;
}

static bool profile_fastconv_single(float* filter, size_t filter_length, size_t signal_length)
{
    assert(DSY_COUNTOF(data_in) <= DSY_COUNTOF(data_out));
    assert(signal_length <= DSY_COUNTOF(data_in));

    // for profiling only the input data doesn't matter so its generation is omitted for speed 

    DUT.SetIR(filter, filter_length);

    // find minimum over several repeats
    uint32_t dt_min = UINT32_MAX;
    uint32_t dt_max = 0;
    {
        ScopedIrqBlocker block;
        
        for (int j = 0; j < 50; j++)
        {
            DUT.Reset();
            const uint32_t t0 = dsy_tim_get_tick();
            for (size_t i = 0; i < signal_length; i++)
            {
                data_out[i] = DUT.Process(data_in[i]);
            }
            const uint32_t dt = dsy_tim_get_tick() - t0;
            dt_min = DSY_MIN(dt_min, dt);
            dt_max = DSY_MAX(dt_max, dt);
        }
    }
    bool pass = true;
    const float time = (float)dt_min / (signal_length * 200); // 200 MHz timer
    const float time2 = (float)dt_max / (signal_length * 200); // 200 MHz timer

    Log.PrintF("filter %4u; signal %5u; Time = " FLT_FMT(2) " .. " FLT_FMT(2) " us/smp, %s", filter_length, signal_length, FLT_VAR(2, time), FLT_VAR(2, time2), str_res[pass]);
    
    return pass;
}

static void verify_fastconv()
{
    static size_t list[] = 
    {
         32,      42,     55,     64,     73,
         97,     128,    171,    227,    256,
        301,     400,    512 ,   632,    707,
        940,    1024,   1250,   1662,   2048
    };

    bool pass = true;
    generate_signal(data_ir, DSY_COUNTOF(data_ir));
#if 1
    for (size_t i = 0; i < DSY_COUNTOF(list); i++)
    {
        const size_t filter_length = list[i];
//        pass &= verify_fastconv_single(data_ir, filter_length, filter_length*4);
        pass &= verify_fastconv_single(data_ir, filter_length, DSY_COUNTOF(data_in));
    }
#endif

    for (size_t filter_length = 32; filter_length <= 2048; filter_length *= 2)
    {
        pass &= profile_fastconv_single(data_ir, filter_length, DSY_COUNTOF(data_in) - 1);  // one less to avoid calculating FFT for the next chunk
    }

    Log.PrintF("Done: %s", str_res[pass]);
}

static void test_fft_perf()
{
    arm_rfft_fast_instance_f32 fft[8];

    arm_fir_instance_f32 fir[8];

    for (int i = 0; i < 8; i++)
    {
        const uint32_t t0 = dsy_tim_get_tick();
        const arm_status res = arm_rfft_fast_init_f32(&fft[i], FFT_LEN(i));
        const uint32_t dt = dsy_tim_get_tick() - t0;

        arm_fir_init_f32(&fir[i], FFT_LEN(i), data_ir, data_tmp, FFT_LEN(i));

        Log.PrintF("arm_rfft_fast_init_f32(%4d) = %u ticks = %d", FFT_LEN(i), dt, (int)res);

        dsy_tim_delay_ms(100);
    }

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
            
            Log.PrintF("arm_rfft_fast_f32(%4d) = %6u ticks; arm_fir_f32() = %6u tocks", FFT_LEN(i), dt, dt2);
            dsy_tim_delay_ms(100);
        }
    }
}


static float ir[64];
static void test_fastconv()
{
    for (size_t i = 0; i < DSY_COUNTOF(ir); i++)    
    {
        ir[i] = 1 + i;        
    }

    Log.PrintF("IR:");
    for (size_t i = 0; i < DSY_COUNTOF(ir); i++)
    {
        Log.PrintF("%2d; " FLT_FMT3, i, FLT_VAR3(ir[i]));
    }

    FastConv<128> FC;
    FC.SetIR(ir, DSY_COUNTOF(ir));
    

    for (size_t i = 0; i < 512; i++)
    {
        float in = 0.0f;
        if (i == 100)
        {
            in = 1.0f;
        }
        float out = FC.Process(in);
        Log.PrintF("%4u; in = " FLT_FMT3 "; out = " FLT_FMT3, i, FLT_VAR3(in), FLT_VAR3(out));
    }

}


int main(void)
{
    // initialize seed hardware and daisysp modules
    float sample_rate;
    seed.Configure();
    seed.Init();
    sample_rate = seed.AudioSampleRate();
    Log.Start(UsbHandle::FS_INTERNAL, true);
    Log.PrintF("samplerate = " FLT_FMT3, FLT_VAR(3, sample_rate));

    verify_fastconv();
    



}
