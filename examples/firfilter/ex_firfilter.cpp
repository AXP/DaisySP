/******************************************************
* FIR Filter module example
*
* (c) 2021 Alexander Petrov-Savchenko
* e-mail: axp@soft-amp.com
*
******************************************************/
#include "daisysp.h"
#include "daisy_seed.h"
#include "FIRGen.h"

using namespace daisysp;
using namespace daisy;

static constexpr size_t flt_size = 512;

static DaisySeed  seed;
static FIRFilter<FIRFILTER_USER_MEMORY>  flt;
static Oscillator lfo;
static WhiteNoise osc;

static float ir_front[flt_size] = {0};
static float ir_back[flt_size+1] = {0};
static float flt_state[flt_size + 1];
static float wnd[flt_size/2] = {0};
static bool ir_update_pending = false;
static int  update_count = 0;
static int  smp_count = 0;

static void InitWindow()
{
    static CWndGenAXP<float> WndGen;

    WndGen.Init(flt_size);
    for (size_t i = 0; i < flt_size/2; i++)
    {
        wnd[i] = WndGen.Gen(i);
   }
}

static void UpdateFilter(float freq)
{
    constexpr int half = flt_size / 2;
    if (false == ir_update_pending)
    {
        const float scale = 2.0f * freq;

        for (int i = 0; i < half; i++)
        {
            const float sinc_arg = PI_F * (i - half);
            const float sinc = sinf(sinc_arg * scale) / sinc_arg;
            const float filt = sinc * wnd[i];
        
            ir_back[i] = filt;
            ir_back[flt_size - i] = filt;
        }
        ir_back[half] = scale;
        ir_update_pending = true;
    }    
}

static void ApplyFilter()
{
    if (true == ir_update_pending)
    {
        for (size_t i = 0; i < flt_size; i++)
        {
            ir_front[i] = ir_back[i];
        }
        ir_update_pending = false;
    }
}

static void AudioCallback(float *in, float *out, size_t size)
{
    float saw, output;

    ApplyFilter();

    for(size_t i = 0; i < size; i += 2)
    {
        saw  = osc.Process();

        output = flt.Process(saw);

        // left out
        out[i] = output;

        // right out
        out[i + 1] = output;
    }

    smp_count += size;
}

int main(void)
{
    // initialize seed hardware and daisysp modules
    int sample_rate;
    seed.Configure();
    seed.Init();
    seed.StartLog(false);
    sample_rate = (int) seed.AudioSampleRate();

    // initialize Biquad and set parameters
    flt.SetStateBuffer(flt_state, DSY_COUNTOF(flt_state));
    flt.SetIR(ir_front, flt_size, false);

    // set parameters for sine oscillator object
    lfo.Init(500);
    lfo.SetWaveform(Oscillator::WAVE_TRI);
    lfo.SetAmp(1.0f);
    lfo.SetFreq(.1f);

    // set parameters for sine oscillator object
    //osc.Init(sample_rate);
    //osc.SetWaveform(Oscillator::WAVE_POLYBLEP_SAW);
    //osc.SetFreq(100);
    osc.Init();
    osc.SetAmp(0.5f);
    

    InitWindow();
    UpdateFilter(0.5f);

    // start callback
    seed.StartAudio(AudioCallback);

    while(1) 
    {
        const float freq = 0.15f * (1.0f + lfo.Process());
        UpdateFilter(freq);
        update_count++;
        dsy_tim_delay_ms(1);

        if (smp_count % sample_rate == 0)
        {
            const float rate = (float)update_count * sample_rate / smp_count;
            seed.PrintLine("avg_rate = " FLT_FMT3 "upd/sec", FLT_VAR3(rate));
        }
    }
}
