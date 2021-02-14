/******************************************************
* Runtime FIR filter generator
*
* (c) 2017-2021 Alexander Petrov-Savchenko
* e-mail: axp@soft-amp.com
*
******************************************************/
#include "FIRGen.h"
#include <assert.h>


#ifdef _DEBUG
#define DBG_PRINTF printf
#define DBG_PUTS   puts
#else
#define DBG_PRINTF sizeof
#define DBG_PUTS   sizeof
#endif

CFIRGen::CFIRGen()
{
}


CFIRGen::~CFIRGen()
{
}

void CFIRGen::Create(CWndGen<float>& pWnd, float fc_norm, float* out_buf, int taps, float mul)
{
    assert((taps & (taps - 1)) == 0); // power of 2
    
    assert(out_buf);

    const float scale = 2.0 * fc_norm;
    pWnd.Init(taps);

    for (int i = 0; i < taps / 2; i++)
    {
        const float sinc_arg = PI * (i - taps / 2);
        const float wind = pWnd.Gen(i);
        const float sinc = sin(sinc_arg * scale) / sinc_arg;
        const float filt = mul * sinc * wind;

        out_buf[i] += filt;
        out_buf[taps - i] += filt;
    }
    out_buf[taps / 2] += mul * scale;
}


void CFIRGen::CreateLPF(CWndGen<float>& pWnd, float fc_norm, float* out_buf, int taps)
{
    assert((taps & (taps - 1)) == 0); // power of 2
    assert(out_buf);

    if (!out_buf)
    {
        return;
    }

    memset(out_buf, 0, sizeof(out_buf[0]) * taps);
    Create(pWnd, fc_norm, out_buf, taps, 1.0);
}

void CFIRGen::CreateHPF(CWndGen<float>& pWnd, float fc_norm, float* out_buf, int taps)
{
    assert((taps & (taps - 1)) == 0); // power of 2
    assert(out_buf);

    if (!out_buf)
    {
        return;
    }

    memset(out_buf, 0, sizeof(out_buf[0]) * taps);
    out_buf[taps / 2] = 1.0;
    Create(pWnd, fc_norm, out_buf, taps, -1.0);
}


void CFIRGen::CreateBPF(CWndGen<float> & pWnd, float fl_norm, float fh_norm, float* out_buf, int taps)
{
    assert((taps & (taps - 1)) == 0); // power of 2
    assert(out_buf);
    assert(fh_norm > fl_norm);

    if (!out_buf)
    {
        return;
    }


    memset(out_buf, 0, sizeof(out_buf[0]) * taps);
    Create(pWnd, fh_norm, out_buf, taps, 1.0);
    Create(pWnd, fl_norm, out_buf, taps, -1.0);
}



void CFIRGen::FlipSpectrum(float* buf, unsigned int length)
{
    assert(buf);
    if (buf)
    {
        for (unsigned int i = 1; i < length; i += 2)
        {
            buf[i] = -buf[i];
        }
    }
}
void CFIRGen::Scale(float* buf, unsigned int length, float factor)
{
    assert(buf);
    if (buf)
    {
        for (unsigned int i = 0; i < length; i++)
        {
            buf[i] *= factor;
        }
    }
}
float CFIRGen::GetFIRgain(float* buf, unsigned int length, float f_norm)
{
    assert(buf);
    float gain = 0.0;
    if (buf)
    {
        float re = 0.0;
        float im = 0.0;
        for (unsigned int i = 0; i < length; i++)
        {
            re += buf[i] * cosf(TWOPI_F * i * f_norm);
            im += buf[i] * sinf(TWOPI_F * i * f_norm);
        }
        gain = sqrt(re * re + im * im);
    }
    return gain;
}



void CFIRGen::Test()
{
    const int cnt = 256;
    float buf1[cnt + 1];
    float buf2[cnt + 1];
    float buf3[cnt + 1];
    float buf4[cnt + 1];
    CWndGenKaiser<float> Wnd(2.0);
    CWndGenAXP<float> WndA();
    CFIRGen::CreateLPF(Wnd, 0.1, buf1, cnt);
    CFIRGen::CreateHPF(Wnd, 0.1, buf2, cnt);
    CFIRGen::CreateBPF(Wnd, 0.1, 0.3, buf3, cnt);
    CFIRGen::CreateLPF(Wnd, 0.1, buf4, cnt);
    for (int i = 0; i < cnt; i++)
    {
        DBG_PRINTF("%4d; %g; %g; %g; %g\n", i, buf1[i], buf2[i], buf3[i], buf4[i]);
    }


}
