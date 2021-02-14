/******************************************************
* Windowing Function generator
*
* (c) 2019-2021 Alexander Petrov-Savchenko
* e-mail: axp@soft-amp.com
*
******************************************************/
#pragma once
#include <daisysp.h>
#include <cstdio>

// interface class

template <typename TYPE = double>
class CWndGen
{
public:
    CWndGen() {}
    virtual ~CWndGen() {}

    virtual void Init(int n) = 0;
    virtual TYPE Gen(int i) = 0;
    virtual const char* GetName() = 0;

    static void TestAll();
    void Test(int n = 100);
};

// Rectangular
template <typename TYPE = double>
class CWndGenRect: public CWndGen<TYPE>
{
public:
	CWndGenRect() {}
	virtual ~CWndGenRect() {}

	virtual void Init(int n) {};
	virtual TYPE Gen(int i) { return (TYPE)1; }
	virtual const char* GetName() {
		return "Rectangular";
	}
};


// Blackman-Harris

template <typename TYPE = double>
class CWndGenBH: public CWndGen<TYPE>
{
public:
    
    CWndGenBH() : m_Scale(0) {}
    virtual ~CWndGenBH() {}

    const char* GetName() { return "Blackman-Harris"; }
    void Init(int n)
    {
        m_Scale = (TYPE)TWOPI_F / n;
    }
    TYPE Gen(int i)
    {
        const TYPE a0 = (TYPE)0.35875;
        const TYPE a1 = (TYPE)0.48829;
        const TYPE a2 = (TYPE)0.14128;
        const TYPE a3 = (TYPE)0.01168;

        const TYPE wind_arg = m_Scale * i;

        const TYPE wind = a0 - a1 * cos(wind_arg) + a2 * cos(wind_arg * 2) - a3 * cos(wind_arg * 3);

        return wind;
    }

protected:
    TYPE   m_Scale;
};


// Lanczos

template <typename TYPE = double>
class CWndGenLanczos : public CWndGen<TYPE>
{
public:

    CWndGenLanczos() : m_Scale(0) {}
    virtual ~CWndGenLanczos() {}

    const char* GetName() { return "Lanczos"; }
    void Init(int n)
    {
        m_Scale = (TYPE)(TWOPI_F) / n;
    }
    TYPE Gen(int i)
    {
        // w[i] = sinc(2*i/N - 1), 
        // sinc x = sin(Pi*x)/(Pi*x);
        const TYPE wind_arg = m_Scale * i - (TYPE)PI_F;
        if (abs(wind_arg) < (TYPE)1.0e-5)
        {
            return (TYPE)1.0;
        }
        return sin(wind_arg) / wind_arg;
    }

protected:
    TYPE   m_Scale;
};


// Kaiser

template <typename TYPE = double>
class CWndGenKaiser : public CWndGen<TYPE>
{
public:

    CWndGenKaiser(TYPE alpha = (TYPE)1.72) : m_ScaleTap(0), m_ScaleArg(0), m_Beta(0)
    {
        SetAlpha(alpha);
    }
    ~CWndGenKaiser() {}

    const char* GetName() { return m_Name; }

    void SetAlpha(TYPE alpha)
    {
        sprintf(m_Name, "Kaiser %3.2f", alpha);
        m_Beta = alpha * (TYPE)PI_F;
        m_ScaleTap = (TYPE)1.0 / BesselFunction0(m_Beta);
    }
    void Init(int n) override
    {
        m_ScaleArg = (TYPE)2.0 / n;
    }
    TYPE Gen(int i) override
    {
        const TYPE arg = m_ScaleArg * i;
        const TYPE wind = m_ScaleTap * BesselFunction0(m_Beta * sqrt(arg * ((TYPE)2.0 - arg)));
        return wind;
    }

protected:
    TYPE BesselFunction0(TYPE z)
    {
        int i = 0;
        TYPE I0 = (TYPE)1.0;
        TYPE delta = (TYPE)1.0;

        z = (TYPE)0.25 * z * z;

        do
        {
            i++;
            delta *= z / TYPE(i * i);
            I0 += delta;
        } while (delta >= (TYPE)1.0e-8 * I0);

        return I0;
    }

    TYPE   m_ScaleTap;
    TYPE   m_ScaleArg;
    TYPE   m_Beta;
    char   m_Name[32];
};

// AXP


template <typename TYPE = double>
class CWndGenAXP : public CWndGen<TYPE>
{
public:

	CWndGenAXP(TYPE c = 0) :
		m_N(2),
		m_B((TYPE)1 / ((TYPE)0.01 + c * c)),
		m_A((TYPE)1.956)
	{
		m_K = m_B + 1;
		m_P = (TYPE)(log(m_B / m_K) / log(cos(m_A / 2)));	
	}
	virtual ~CWndGenAXP() {}

	const char* GetName() { return "AXP"; }
	void Init(int n)
	{
		m_N = (TYPE)1 / n;
	}
	TYPE Gen(int i)
	{
		

#if 0 // reference algorithm, no longer maintained
		
		const int n = i - m_N / 2; 
		if (n == 0)
		{
			return (TYPE)1;
		}

		TYPE sum = 0;
		for (int k = -m_N / 2; k <= m_N / 2; k++)
		{
			sum += k * sin(DSP_2PI * k * n / m_N);
		}

		sum *= (TYPE)(-2.0) / (m_N * m_N);

		TYPE integ = cos(DSP_PI * n) / (DSP_PI * n);

		return sum / integ;

#else // fast approximation, but parametric
		const TYPE w = pow(cos(m_A * (m_N * i - (TYPE)0.5)), m_P) * m_K - m_B;
		return w;
#endif

	}

protected:
	TYPE m_N;
	TYPE m_P;
	TYPE m_K;
	TYPE m_B;
	TYPE m_A;

};


