#pragma once
#ifndef DSY_FFT_H
#define DSY_FFT_H

#include <stdint.h>
#include "dsp.h"

#include "ffft/FFTReal.h"
#include "ffft/FFTRealFixLen.h"

#include "constvar.h"
#include "shy_fft.h"

#ifdef __arm__
#include <arm_math.h>
#include <arm_common_tables.h>
#include "fft_tables.h"
#endif

//include "daisy_seed.h"

namespace daisysp
{

template<size_t max_size, bool fixed>
class ShyFFTBase
{
public:

    using ShyFFTType = stmlib::ShyFFT<float, max_size, stmlib::RotationPhasor>;

    ShyFFTBase()
        : num_pass_((size_t)get_base_2((uint32_t)max_size))
    {
        fft_.Init();
    }

    bool SetSize(size_t new_size)
    {
        //new_size = DSY_MAX(new_size, 2);    // avoid zero length
        if (new_size <= fft_.max_size)
        {
            num_pass_ = (size_t)get_base_2((uint32_t)new_size);
            return true;
        }
        return false;
    }

    size_t GetSize() const { return 1U << num_pass_; }

    /* multiply two vectors in frequency domain */
    void Multiply(const float* buf_fdA, const float* buf_fdB, float* buf_fdDest)
    {
        const size_t mid = GetSize() / 2;

        buf_fdDest[0] = buf_fdA[0] * buf_fdB[0]; // DC
        for (size_t i = 1; i < mid; i++)
        {
            const float reA = buf_fdA[i];
            const float reB = buf_fdB[i];
            const float imA = buf_fdA[mid + i]; // actually negated, but factors out anyway
            const float imB = buf_fdB[mid + i]; // actually negated, but factors out anyway
            buf_fdDest[i] = reA * reB - imA * imB; // real
            buf_fdDest[mid + i] = reA * imB + reB * imA; // imaginary
        }
        buf_fdDest[mid] = buf_fdA[mid] * buf_fdB[mid]; // Nyquist
    }

    void ForwardTransform(float* __restrict buf_td, float* __restrict buf_fd)
    {
        fft_.Direct(buf_td, buf_fd, num_pass_);
    }
    void InverseTransform(float* __restrict buf_fd, float* __restrict buf_td)
    {
        fft_.Inverse(buf_fd, buf_td, num_pass_);
        const float scale = 1.0f / GetSize();
        
         for (size_t i = 0; i < GetSize(); i++)
         {
             buf_td[i] *= scale;
         }
    }

    
protected:
    ShyFFTType fft_;
    size_t num_pass_;
    
};


template <size_t size>
class ShyFFTBase<size, true>
{
public:
    using ShyFFTType = stmlib::ShyFFT<float, size, stmlib::RotationPhasor>;

    ShyFFTBase()
    {
        fft_.Init();
    }
    bool SetSize(size_t new_size)
    {
        return (long)new_size == GetSize();
    }
    static constexpr size_t GetSize() { return size; }

    /* multiply two vectors in frequency domain */
    static void Multiply(const float* buf_fdA, const float* buf_fdB, float* buf_fdDest)
    {
        constexpr size_t mid = GetSize() / 2;

        buf_fdDest[0] = buf_fdA[0] * buf_fdB[0]; // DC
        for (size_t i = 1; i < mid; i++)
        {
            const float reA = buf_fdA[i];
            const float reB = buf_fdB[i];
            const float imA
                = buf_fdA[mid + i]; // actually negated, but factors out anyway
            const float imB
                = buf_fdB[mid + i]; // actually negated, but factors out anyway
            buf_fdDest[i] = reA * reB - imA * imB; // real
            buf_fdDest[mid + i] = reA * imB + reB * imA; // imaginary
        }
        buf_fdDest[mid] = buf_fdA[mid] * buf_fdB[mid]; // Nyquist
    }
    void ForwardTransform(float* __restrict buf_td, float* __restrict buf_fd)
    {
        fft_.Direct(buf_td, buf_fd);
    }
    void InverseTransform(float* __restrict buf_fd, float* __restrict buf_td)
    {
        fft_.Inverse(buf_fd, buf_td);
        constexpr float scale = 1.0f / GetSize();

        for (size_t i = 0; i < GetSize(); i++)
        {
            buf_td[i] *= scale;
        }
    }

protected:
    ShyFFTType fft_;
};



template<size_t size, bool fixed>
class FFT_Shy : public ShyFFTBase<size, fixed>
{
public:
    FFT_Shy()
    {
    }
    ~FFT_Shy()
    {}
    enum Consts
    {
        FFT_MAXSIZE = 65536 // arbitrary
    };

    using ShyFFTBase<size, fixed>::SetSize;
    using ShyFFTBase<size, fixed>::GetSize;
    using ShyFFTBase<size, fixed>::fft_;
    using ShyFFTBase<size, fixed>::Multiply;
    using ShyFFTBase<size, fixed>::ForwardTransform;
    using ShyFFTBase<size, fixed>::InverseTransform;

  
};



template<size_t max_size, bool fixed>
class FFTRealBase
{
  public:
    FFTRealBase()
    : fft_(max_size) 
    {
    }
    
    bool SetSize(size_t new_size)
    {
        new_size = DSY_MAX(new_size, 2);    // avoid zero length
        if(new_size <= max_size)
        {
            fft_.set_length((long)new_size);
            return (long)new_size == fft_.get_length();
        }
        return false;
    }

    size_t GetSize() const { return (size_t)fft_.get_length(); }
    
    /* multiply two vectors in frequency domain */
    void Multiply(const float *buf_fdA, const float *buf_fdB, float *buf_fdDest)
    {
        const size_t mid = fft_.get_length() / 2;

        buf_fdDest[0] = buf_fdA[0] * buf_fdB[0]; // DC
        for(size_t i = 1; i < mid; i++)
        {
            const float reA = buf_fdA[i];
            const float reB = buf_fdB[i];
            const float imA
                = buf_fdA[mid + i]; // actually negated, but factors out anyway
            const float imB
                = buf_fdB[mid + i]; // actually negated, but factors out anyway
            buf_fdDest[i]       = reA * reB - imA * imB; // real
            buf_fdDest[mid + i] = reA * imB + reB * imA; // imaginary
        }
        buf_fdDest[mid] = buf_fdA[mid] * buf_fdB[mid]; // Nyquist
    }
  protected:
    ffft::FFTReal<float> fft_;
};

    // import implementation from common template specialization
//static constexpr auto InitFFT = ARMFFTConfig<size, false>::InitFFT;


template <size_t size>
class FFTRealBase<size, true> 
{
  public:
    using FFTRealType = ffft::FFTRealFixLen<get_base_2(size)>;

    bool SetSize(size_t new_size)
    {
        return (long)new_size == GetSize();
    }
    static constexpr size_t GetSize()  { return (size_t)FFTRealType::get_length(); }

    /* multiply two vectors in frequency domain */
    void Multiply(const float *buf_fdA, const float *buf_fdB, float *buf_fdDest)
    {
        constexpr size_t mid = fft_.get_length() / 2;

        buf_fdDest[0] = buf_fdA[0] * buf_fdB[0]; // DC
        for(size_t i = 1; i < mid; i++)
        {
            const float reA = buf_fdA[i];
            const float reB = buf_fdB[i];
            const float imA
                = buf_fdA[mid + i]; // actually negated, but factors out anyway
            const float imB
                = buf_fdB[mid + i]; // actually negated, but factors out anyway
            buf_fdDest[i]       = reA * reB - imA * imB; // real
            buf_fdDest[mid + i] = reA * imB + reB * imA; // imaginary
        }
        buf_fdDest[mid] = buf_fdA[mid] * buf_fdB[mid]; // Nyquist
    
    }
  protected:
    FFTRealType fft_;
};


template<size_t size, bool fixed>
class FFT_Generic: public FFTRealBase<size, fixed>
{
public:
    FFT_Generic()
    {
    }
    ~FFT_Generic()
    {}
    enum Consts
    {
        FFT_MAXSIZE = 65536 // arbitrary
    };

    using FFTRealBase<size, fixed>::SetSize;
    using FFTRealBase<size, fixed>::GetSize;
    using FFTRealBase<size, fixed>::fft_;
    using FFTRealBase<size, fixed>::Multiply;

    void ForwardTransform(float* __restrict buf_td, float* __restrict buf_fd)
    {
        fft_.do_fft(buf_fd, buf_td);
    }
    void InverseTransform(float* __restrict buf_fd, float* __restrict buf_td)
    {
        fft_.do_ifft(buf_fd, buf_td);
        fft_.rescale(buf_td);
    }
};





#ifdef __arm__


template <size_t size, bool fixed> 
class ARMFFTConfig
{
public:

    ARMFFTConfig()
    {
    }

    ~ARMFFTConfig()
    {
    }

    operator arm_rfft_fast_instance_f32*()
    {
        return &fft_;
    }

    // constexpr variant of CMSIS ARM FFT init
    static constexpr arm_rfft_fast_instance_f32 InitFFT(uint16_t fftLen)
    {
        arm_rfft_fast_instance_f32 S = {0};

        S.Sint.fftLen = fftLen / 2;
        S.fftLenRFFT  = fftLen;

        switch(fftLen / 2)
        {
            default:
                // leave initialized to zero, will be caught externally
                break;

            case 2048u: 
                S.Sint.bitRevLength = ARMBITREVINDEXTABLE_2048_TABLE_LENGTH;
                S.Sint.pBitRevTable = (uint16_t *)armBitRevIndexTable2048;
                S.Sint.pTwiddle     = (float32_t *)twiddleCoef_2048;
                S.pTwiddleRFFT      = (float32_t *)twiddleCoef_rfft_4096;
                break;

            case 1024u:
                S.Sint.bitRevLength = ARMBITREVINDEXTABLE_1024_TABLE_LENGTH;
                S.Sint.pBitRevTable = (uint16_t *)armBitRevIndexTable1024;
                S.Sint.pTwiddle     = (float32_t *)twiddleCoef_1024;
                S.pTwiddleRFFT      = (float32_t *)twiddleCoef_rfft_2048;
                break;

            case 512u:
                S.Sint.bitRevLength = ARMBITREVINDEXTABLE_512_TABLE_LENGTH;
                S.Sint.pBitRevTable = (uint16_t *)armBitRevIndexTable512;
                S.Sint.pTwiddle     = (float32_t *)twiddleCoef_512;
                S.pTwiddleRFFT      = (float32_t *)twiddleCoef_rfft_1024;
                break;

            case 256u:
                S.Sint.bitRevLength = ARMBITREVINDEXTABLE_256_TABLE_LENGTH;
                S.Sint.pBitRevTable = (uint16_t *)armBitRevIndexTable256;
                S.Sint.pTwiddle     = (float32_t *)twiddleCoef_256;
                S.pTwiddleRFFT      = (float32_t *)twiddleCoef_rfft_512;
                break;

            case 128u:
                S.Sint.bitRevLength = ARMBITREVINDEXTABLE_128_TABLE_LENGTH;
                S.Sint.pBitRevTable = (uint16_t *)armBitRevIndexTable128;
                S.Sint.pTwiddle     = (float32_t *)twiddleCoef_128;
                S.pTwiddleRFFT      = (float32_t *)twiddleCoef_rfft_256;
                break;

            case 64u:
                S.Sint.bitRevLength = ARMBITREVINDEXTABLE_64_TABLE_LENGTH;
                S.Sint.pBitRevTable = (uint16_t *)armBitRevIndexTable64;
                S.Sint.pTwiddle     = (float32_t *)twiddleCoef_64;
                S.pTwiddleRFFT      = (float32_t *)twiddleCoef_rfft_128;
                break;

            case 32u:
                S.Sint.bitRevLength = ARMBITREVINDEXTABLE_32_TABLE_LENGTH;
                S.Sint.pBitRevTable = (uint16_t *)armBitRevIndexTable32;
                S.Sint.pTwiddle     = (float32_t *)twiddleCoef_32;
                S.pTwiddleRFFT      = (float32_t *)twiddleCoef_rfft_64;
                break;

            case 16u:
                S.Sint.bitRevLength = ARMBITREVINDEXTABLE_16_TABLE_LENGTH;
                S.Sint.pBitRevTable = (uint16_t *)armBitRevIndexTable16;
                S.Sint.pTwiddle     = (float32_t *)twiddleCoef_16;
                S.pTwiddleRFFT      = (float32_t *)twiddleCoef_rfft_32;
                break;
        }
        return S;
    }

    bool SetSize(uint16_t fftLen)
    {
        fft_ = InitFFT(fftLen);
        return fftLen == fft_.fftLenRFFT;
    }
    uint16_t GetSize() const
    {
        return fft_.fftLenRFFT;
    }    
protected:
    ARMFFTTables        lut_;
    arm_rfft_fast_instance_f32 fft_;
};

// Specialization for fixed size configuration
template <size_t size> 
class ARMFFTConfig<size, true>
{
public:
    // import implementation from common template specialization
    static constexpr auto InitFFT = ARMFFTConfig<size, false>::InitFFT;

    ARMFFTConfig()
    {
    }

    ~ARMFFTConfig()
    {
    }

    operator arm_rfft_fast_instance_f32*()
    {
        return (arm_rfft_fast_instance_f32*) &fft_;
    }
    
    static constexpr bool SetSize(uint16_t __unused fftLen)
    {
        return fftLen <= size;
    }
    static constexpr size_t GetSize()
    {
        return size;
    }

protected: // should add "inline" keyword if C++17 is enabled
    static constexpr arm_rfft_fast_instance_f32 fft_ = InitFFT(size);
    ARMFFTTables lut_;
};

// in C++14 static variables should still be defined outside of the class
template <size_t size>
constexpr arm_rfft_fast_instance_f32 ARMFFTConfig<size, true>::fft_;



template<size_t size, bool fixed>
class FFT_ARM
{
public:
    FFT_ARM()
    {}
    ~FFT_ARM()
    {}
    enum Consts
    {
        FFT_MAXSIZE = 4096
    };
    size_t GetSize() const
    {
        return (size_t)fft_.GetSize();
    }
    bool SetSize(size_t new_size)
    {       
        return fft_.SetSize(new_size);
    }
    void ForwardTransform(float* __restrict buf_td, float* __restrict buf_fd)
    {
        arm_rfft_fast_f32(fft_, buf_td, buf_fd, FFT_FORWARD);
    }
    void InverseTransform(float* __restrict buf_fd, float* __restrict buf_td)
    {
        arm_rfft_fast_f32(fft_, buf_fd, buf_td, FFT_INVERSE);
    }

    /* multiply two vectors in frequency domain */
    void Multiply(const float *buf_fdA, const float *buf_fdB, float *buf_fdDest)
    {
        buf_fdDest[0] = buf_fdA[0] * buf_fdB[0]; // DC
        buf_fdDest[1] = buf_fdA[1] * buf_fdB[1]; // Nyquist
        arm_cmplx_mult_cmplx_f32(
            (float*)&buf_fdA[2], (float*)&buf_fdB[2], &buf_fdDest[2], GetSize()/2 - 1);
    }

  protected:
    static constexpr uint8_t FFT_FORWARD = 0;
    static constexpr uint8_t FFT_INVERSE = 1;
    
    ARMFFTConfig<size, fixed> fft_;
};

//template <size_t size = 64, bool fixed = true>
//using FFT = FFT_ARM<size, fixed>;

// !!!! Force FFT

template <size_t size = 64, bool fixed = true>
//using FFT = FFT_Generic<size, fixed>;  
using FFT = FFT_Shy<size, fixed>;  




#else // if __arm__

//template <size_t size = 64, bool fixed = true>
//using FFT = FFT_Generic<size, fixed>;

template <size_t size = 64, bool fixed = true>
using FFT = FFT_Shy<size, fixed>;


#endif  // __arm__ / generic

} // namespace daisysp

#endif // DSY_FFT_H
