#pragma once
#ifndef DSY_FFTFILTER_H
#define DSY_FFTFILTER_H

#include <stdint.h>
#include "dsp.h"
#include <type_traits>
#include <cstddef>

#include <utility>
#ifdef __arm__

#include <arm_math.h>

#else

#ifdef _MSC_VER
#define __unused
#endif
static void arm_copy_f32(float* pSrc, float* pDst, size_t blockSize) 
{
    memcpy(pDst, pSrc, blockSize * sizeof(float));
}

static void arm_add_f32(float* pSrcA, float* pSrcB, float* pDst, size_t blockSize)
{
    for(size_t i = 0; i < blockSize; i++)
    {
        pDst[i] = pSrcA[i] + pSrcB[i];
    }
}

#endif // __arm__

#include "fft.h"



namespace daisysp
{

using namespace FastConvolution;

/*
(Overlap-save algorithm for linear convolution)
h = FIR_impulse_response
M = length(h)
overlap = M − 1
N = 8 × overlap    (see next section for a better choice)
step_size = N − overlap
H = DFT(h, N)
position = 0

while position + N ≤ length(x)
    yt = IDFT(DFT(x(position+(1:N))) × H)
    y(position+(1:step_size)) = yt(M : N)    (discard M−1 y-values)
    position = position + step_size
end


(Overlap-save algorithm for linear convolution)
h = FIR_impulse_response

len = length(h)
overlap = len − 1
fft_len = 2 * len   
step_size = fft_len − overlap = len + 1
H = DFT(h, fft_len)

position = 0

while position + fft_len ≤ length(x)

    X = DFT(&x[position], fft_len)
    yt = IDFT(X × H)
    for (int j = 0; j < step_size; j++)
    {
        y[position + j] = yt[len + j];
    }
    position += step_size
end



-----xxxxxx--------------

 ffff0000 - len
[        ] FFT
     xxxxx - max number of samples
N = number of input samples in convolution: 
fft_len = N + len - 1
N = fft_len - len + 1


of which:
invalid pre-ringing
MIN(len, N) - 1
invalid post-ringing:
MIN(len, N) - 1

valid output samples:
K = fft_len - 2 * (MIN(len, N) - 1)
K = fft_len - 2 * (MIN(len, fft_len - len + 1) - 1)

len = 8
fft_len = 16
N = 9
K = 16 - 2 *(7) = 2

now, if pre-ringing and post-ringing segments are forced to overlap:
consume (MIN(len, N) - 1) more samples
len = filter
N = fft_len - len + 1 + MIN(len, fft_len - len + 1) - 1
N = fft_len - len + MIN(len, fft_len - len + 1)

now number of invalid samples should be MIN(len, N) - 1

valid output samples:
K = fft_len - (MIN(len, N) - 1)

len = 8
fft_len = 16
N = 16
invalid samples: 7
valid samples: K = 16 - 8 + 1 = 9



[-/YYYY\_]


*/








template<size_t max_size, bool fixed = false>
class FFTConvSize
{
public:
    enum
    {
        MAX_SIZE_POWER2 = get_next_power2(max_size)
    };

    FFTConvSize()
        :   fft_len_(0)
        ,   len_(0)
    {
     //   Resize(MAX_SIZE_POWER2);
    }
    void Resize(size_t size = UINT32_MAX/2)
    {
        size = (size_t)get_next_power2((uint32_t)size);
        size = DSY_MIN(size, MAX_SIZE_POWER2);
        fft_len_ = size * 2;
        len_ = size;
    }
    size_t GetLen() const noexcept
    {
        return len_;
    }
    size_t GetFFTLen() const noexcept
    {
        return fft_len_;
    }
protected:
    size_t fft_len_;
    size_t len_;
};

template<size_t size> 
struct FFTConvSize<size, true>
{
public:
    FFTConvSize()
    {
    }
    static void Resize(size_t __unused size_unused)
    {        
    }
    static constexpr size_t GetLen()
    {
        return get_next_power2(size);
    }
    static constexpr size_t GetFFTLen()
    {
        return get_next_power2(size) * 2;
    }

protected:
};




template<size_t max_size, bool fixed_size = false>
class FFTConv
{
  public:

    FFTConv(int init = 0) 
       /* : buf_ir_{0}
        , buf_io_{0}
        , buf_tmp_{0}
        , buf_copy_{0}
        , pos_(0)*/
    {
    }
    ~FFTConv()
    {
    }

    enum Consts
    {
        MAX_SIZE = get_next_power2(max_size),
    };

    using FFTModuleType = FFT<MAX_SIZE * 2, fixed_size>;
    //enum Consts2
    //{
    static constexpr size_t MAX_SUPPORTED_SIZE = FFTModuleType::FFT_MAXSIZE / 2;
    //};

    // Returns the algorithm processing latency in samples
    size_t GetLatency() const noexcept {return size_.GetLen();}

    // Resets the state of the algorithm. Initialized filter kernel is preserved.
    void Reset()
    {
        pos_ = size_.GetLen(); // start in the upper half of the buffer

     //   log_.PrintF("FFTConv:GetLen = %u, GetFFTLen = %u", size_.GetLen(), size_.GetFFTLen()); 
        memset(buf_copy_, 0, sizeof(buf_copy_));
        memset(buf_io_, 0, sizeof(buf_io_));
    }

    // Configures the algorithm to use the provided impulse response (AKA filter kernel)
    //  ir - pointer to the buffer with filter coefficients (in early-first order)
    //  len - length of the filter
    bool SetIR(const float* ir, size_t len)
    {
        size_.Resize(len); // will round up to a power of 2 internally

        if (size_.GetLen() > MAX_SIZE)
        {
            return false;
        }

        if (true != fft_.SetSize(size_.GetFFTLen()))
        {
            return false;
        }

        // copy and zero-pad to FFT length
        memcpy(buf_tmp_, ir, len * sizeof(float));
        memset(&buf_tmp_[len], 0, (size_.GetFFTLen() - len) * sizeof(float));

        //daisy::Logger<>::PrintLine("FFTConv::SetIR, len = %u, size_ = %u/%u", len, size_.GetLen(), size_.GetFFTLen());

        // transform filter kernel into frequency domain
        fft_.ForwardTransform(buf_tmp_, buf_ir_);
       // for (size_t i = 0; i < size_.GetFFTLen(); i++)
        {
         //   daisy::Logger<>::PrintLine("%4u; " FLT_FMT3, i, FLT_VAR3(buf_ir_[i]));
        }
        
        // Clean up the algorithm state
        Reset();
        
        return true;
    }

    
    void ProcessBlock(float* pSrc, float *pDst, size_t num_samples)
    {
        do // this simplest implementation also happens to be the fastest
        {
            const size_t num_proc = DSY_MIN(num_samples, size_.GetFFTLen() - pos_);

            arm_copy_f32(&buf_io_[pos_], pDst, num_proc);
            arm_copy_f32(pSrc, &buf_io_[pos_], num_proc);
            pSrc += num_proc;
            pDst += num_proc;
            num_samples -= num_proc;

            pos_ += num_proc;
            HandleNewSamples();

        } while (num_samples > 0);
    }

    void ProcessBlockAcc(float* pSrc, float* pDst, size_t num_samples)
    {
        do 
        {
            const size_t num_proc
                = DSY_MIN(num_samples, size_.GetFFTLen() - pos_);

            //arm_copy_f32(&buf_io_[pos_], pDst, num_proc);
            arm_add_f32(&buf_io_[pos_], pDst, pDst, num_proc);

            arm_copy_f32(pSrc, &buf_io_[pos_], num_proc);
            pSrc += num_proc;
            pDst += num_proc;
            num_samples -= num_proc;

            pos_ += num_proc;
            HandleNewSamples();

        } while(num_samples > 0);
    }
    // Processes one sample through the filter using Overlap&Discard algorithm and returns one sample.
    //    in - input signal 
    float Process(float in)
    {
        const float out = buf_io_[pos_];
        buf_io_[pos_] = in;
        
        pos_++;
        HandleNewSamples();
        
        return out;
    }

    

protected: // Member functions

    void Convolve()
    {
        // transform the input signal to frequency domain
        fft_.ForwardTransform(buf_io_, buf_tmp_);

        // apply filter by multiplication in frequency domain
        fft_.Multiply(buf_tmp_, buf_ir_, buf_tmp_);

        // transform the output back to time domain
        fft_.InverseTransform(buf_tmp_, buf_io_);
    }

    void HandleNewSamples()
    {
        assert(pos_ <= size_.GetFFTLen()); // no overrun allowed
        if (pos_ >= size_.GetFFTLen())
        {
            // accumulated enough samples to process a frame
            pos_ = size_.GetLen();
            // restore half of the previous frame
            arm_copy_f32(buf_copy_, buf_io_, size_.GetLen());
            // store half of the current frame for future use
            arm_copy_f32(buf_io_ + size_.GetLen(), buf_copy_, size_.GetLen());

            // perform the fast convolution of a single frame. Input and output is buf_io_.
            Convolve(); 
        }
    }


protected: // Member variables
    
    FFTConvSize<MAX_SIZE, fixed_size> size_;
    FFTModuleType fft_;    

    float buf_ir_[MAX_SIZE * 2];
    float buf_io_[MAX_SIZE * 2];
    float buf_tmp_[MAX_SIZE * 2];
    float buf_copy_[MAX_SIZE];
    
    size_t pos_;
};

} // namespace daisysp

#endif // DSY_FFTFILTER_H
