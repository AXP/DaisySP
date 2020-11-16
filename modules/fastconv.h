#pragma once
#ifndef DSY_FASTCONV_H
#define DSY_FASTCONV_H

#include <stdint.h>
#include <arm_math.h>
#include "dsp.h"
#if defined(__cplusplus) && defined(__arm__)

#include "util/logger.h"
namespace daisysp
{
/** A first-order recursive low-pass filter with variable frequency response.
*/


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

#define FASTCONV_OVERLAP_DISCARD    0
#define FASTCONV_OVERLAP_ADD        1
#define FASTCONV_MODE FASTCONV_OVERLAP_ADD



template<size_t max_size>
class FastConv
{
  public:

    FastConv() 
    {
        assert(is_power2(max_size));
    }
    ~FastConv() {}

    enum FastConvConsts
    {
        FFT_FORWARD = 0,
        FFT_INVERSE = 1,
        fft_len_ = max_size,
        len_ = max_size / 2
    };

    // Returns the algorithm processing latency in samples
    size_t GetLatency() {return len_;}

    // Resets the state of the algorithm. Initialized filter kernel is preserved.
    void Reset()
    {
        pos_ = 0;

        memset(buf_copy_, 0, sizeof(buf_copy_));
        memset(buf_io_, 0, sizeof(buf_io_));
    }

    // Configures the algorithm to use the provided impulse response (AKA filter kernel)
    //  ir - pointer to the buffer with filter coefficients (in early-first order)
    //  len - length of the filter
    bool SetIR(const float* ir, size_t len)
    {
        //len_ = get_next_power2(len);

        //fft_len_ = len_ * 2;

        if (fft_len_ > max_size)
        {
            return false;
        }

        if (ARM_MATH_SUCCESS != arm_rfft_fast_init_f32(&fft_, fft_len_))
        {
            return false;
        }

        // copy and zero-pad to FFT length
        memcpy(buf_tmp_, ir, len * sizeof(float));
        memset(&buf_tmp_[len], 0, (fft_len_ - len) * sizeof(float));

        // transform filter kernel into frequency domain
        arm_rfft_fast_f32(&fft_, buf_tmp_, buf_ir_, FFT_FORWARD);

        
#if 0   
        daisy::axp::Logger::PrintF("FFT(IR):");
#if 0
        for (size_t i = 0; i < fft_len_; i += 2)
        {
            daisy::axp::Logger::PrintF("%3u; " FLT_FMT3 "; " FLT_FMT3, i, FLT_VAR3(buf_ir_[i]), FLT_VAR3(buf_ir_[i+1]));
        }
#else
        for (size_t i = 0; i < fft_len_; i ++)
        {
            daisy::axp::Logger::PrintF("%3u; " FLT_FMT3, i, FLT_VAR3(buf_ir_[i]));
        }

#endif
        daisy::axp::Logger::PrintF("len = %u, fft_len = %u", len_, fft_len_);
#endif

        // Clean up the algorithm state
        Reset();
        pos_ = len_; // TODO decide and fix
        
        return true;
    }

    void ProcessBlock(float* pSrc, float *pDst, size_t num_samples)
    {
        do
        {
            const size_t num_proc = DSY_MIN(num_samples, fft_len_ - pos_);

            arm_copy_f32(&buf_io_[pos_], pDst, num_proc);
            arm_copy_f32(pSrc, &buf_io_[pos_], num_proc);
            pos_ += num_proc;
            pSrc += num_proc;
            pDst += num_proc;
            num_samples -= num_proc;

            if (pos_ >= fft_len_)
            {
                // accumulated enough samples to process a frame
                pos_ = len_;
                // restore half of the previous frame
                arm_copy_f32(buf_copy_, buf_io_, len_);
                // store half of the current frame for future use
                arm_copy_f32(buf_io_ + len_, buf_copy_, len_);

                // perform the fast convolution of a single frame. Input and output is buf_io_.
                Convolve(); 
            }
        } while (num_samples > 0);
    }


    // Processes one sample through the filter using Overlap&Discard algorithm and returns one sample.
    //    in - input signal 
    float Process(float in)
    {
        // assumes pos is [0..len-)

        // Note: from the time domain perspective, it would be more intuitive 
        // to write and read to/from position (pos_ + len_)
        // But here we use the periodicity of FFT transform to avoid the unnecessary addition

        const float out = buf_io_[pos_];
        buf_io_[pos_] = in;
        
        if (++pos_ >= len_)
        {
            // accumulated enough samples to process a frame
            pos_ = 0;
            // restore half of the previous frame
            arm_copy_f32(buf_copy_, buf_io_ + len_, len_);
            // store half of the current frame for future use
            arm_copy_f32(buf_io_, buf_copy_, len_);

            // perform the fast convolution of a single frame. Input and output is buf_io_.
            Convolve(); 
        }
        return out;
    }

    // Process one sample using Overlap & Add algorithm (kept for reference)
    float Process_Add(float in)
    {
        float out = buf_io_[pos_] + buf_copy_[pos_];
        buf_io_[pos_] = in;

        if (++pos_ >= len_)
        {
            // accumulated enough samples to process a frame
            pos_ = 0;
            
            // preserve copy of the tail for summation in the next frame
            arm_copy_f32(&buf_io_[len_], buf_copy_, len_);
            // zero pad the input signal
            arm_fill_f32(0.0f, &buf_io_[len_], len_);
            
            // perform the fast convolution of a single frame. Input and output is buf_io_.
            Convolve();            
        }
        return out;
    }
  private:
    void Convolve()
    {
        // transform the input signal to frequency domain
        arm_rfft_fast_f32(&fft_, buf_io_, buf_tmp_, FFT_FORWARD);

        // apply filter by multiplication in frequency domain
        buf_tmp_[0] *= buf_ir_[0]; // DC
        buf_tmp_[1] *= buf_ir_[1]; // Nyquist
        arm_cmplx_mult_cmplx_f32(&buf_tmp_[2], &buf_ir_[2], &buf_tmp_[2], len_ - 1);

        // transform the output back to time domain
        arm_rfft_fast_f32(&fft_, buf_tmp_, buf_io_, FFT_INVERSE);
    }

    float buf_ir_[max_size];
    float buf_io_[max_size];
    float buf_tmp_[max_size];
    float buf_copy_[max_size/2];
    
    size_t pos_;
    //size_t len_;
    
    //size_t fft_len_;
    arm_rfft_fast_instance_f32 fft_;

    static daisy::axp::Logger log_;
};
} // namespace daisysp
#endif // ifdef __cplusplus and __arm__
#endif // DSY_FASTCONV_H
