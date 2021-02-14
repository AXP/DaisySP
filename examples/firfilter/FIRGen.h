/******************************************************
* Runtime FIR filter generator
*
* (c) 2019-2021 Alexander Petrov-Savchenko
* e-mail: axp@soft-amp.com
*
******************************************************/
#pragma once
#include "daisysp.h"
#include "WndGen.h"

// static methods only, no instantiation required or possible
class CFIRGen
{
protected:
    CFIRGen();
    virtual ~CFIRGen();
    static void Create(CWndGen<float>& pWnd, float fc_norm, float* out_buf, int taps, float mul);
public:
    static void CreateLPF(CWndGen<float>& pWnd, float fc_norm, float* out_buf, int taps);
    static void CreateHPF(CWndGen<float>& pWnd, float fc_norm, float* out_buf, int taps); 
    static void CreateBPF(CWndGen<float>& pWnd, float fl_norm, float fh_norm, float* out_buf, int taps);

    static void Scale(float* buf, unsigned int length, float factor);
    static void FlipSpectrum(float* buf, unsigned int length);
    static float GetFIRgain(float* buf, unsigned int length, float f_norm);
    static void Test();
};
