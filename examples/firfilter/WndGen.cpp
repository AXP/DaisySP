/******************************************************
* Windowing Function generator
*
* (c) 2019-2021 Alexander Petrov-Savchenko
* e-mail: axp@soft-amp.com
*
******************************************************/
#include "WndGen.h"

template <typename TYPE>
void CWndGen<TYPE>::TestAll()
{
    CWndGenBH<TYPE> WndBH;
    CWndGenLanczos<float> WndL;
    CWndGenKaiser<double> WndK;
    CWndGenKaiser<> WndK1(1.0);
    CWndGenKaiser<> WndK2(2.0);
	CWndGenAXP<> WndAXP;
	CWndGenRect<TYPE> WndRect;

    WndBH.Test();
    WndL.Test();
    WndK1.Test();
    WndK.Test();
    WndK2.Test();
	WndAXP.Test();
	WndRect.Test();
}

template <typename TYPE>
void CWndGen<TYPE>::Test(int n)
{
    Init(n);
    printf("CWndGen::Test(%d) - %s\n", n, GetName());
    for (int i = 0; i <= n; i++)
    {
        TYPE w = Gen(i);
        printf("%4d; %g\n", i, (double)w);
    }
}

template class CWndGen<double>;
template class CWndGen<float>;