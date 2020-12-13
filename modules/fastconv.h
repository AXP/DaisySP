#pragma once
#ifndef DSY_FASTCONV_H
#define DSY_FASTCONV_H

#include <cstdint>
#include "dsp.h"
#include <type_traits>
#include <cstddef>

#include <utility>
#include "fastconv_plan.h"
#include "firfilter.h"
#include "fftfilter.h"


namespace daisysp
{

template<size_t base, size_t max_size>
class FastConv
{
    public:
    FastConv(int init = 0)
    {
        assert(is_power2(base));
    }
    
    static size_t GetLatency()
    {
        return 0;
    }
    void Reset()
    {
        f0_.Reset();
        for_each(fn_, [] (auto& x) 
        {
            x.Reset();
        });
    }
    float Process(float in)
    {
        float out = f0_.Process(in);

        for_each(fn_, [in, &out] (auto& x) mutable
        {
            out += x.Process(in);
        });

        return out;
    }

    void ProcessBlock(float* pSrc, float *pDst, size_t num_samples)
    {
        f0_.ProcessBlock(pSrc, pDst, num_samples);
#if 0
        for (size_t i = 0; i < num_samples; i++)
        {
            for_each(fn_, [pSrc, pDst, i] (auto& x) mutable
            {
                pDst[i] += x.Process(pSrc[i]);
            });
        }
#else
        for_each(fn_, [pSrc, pDst, num_samples] (auto& x) mutable
        {
            x.ProcessBlockAcc(pSrc, pDst, num_samples);
        });
#endif        
    }
    
    bool SetIR(const float* ir, size_t len)
    {
        assert(len <= max_size);
        assert(nullptr != ir);
        if (len > max_size || nullptr == ir)
        {
            return false;
        }

        /* zero-latency FIR stage */
        const size_t fir_len = DSY_MIN(base, len);
        bool result = f0_.SetIR(ir, fir_len);
        len -= fir_len;

        /* progressive-latency FFT stages (compile-time loop) */
        for_each(fn_, [ir, &len, &result] (auto& x) mutable 
        {
            const size_t latency = x.GetLatency();
            const size_t fft_len = DSY_MIN(latency, len);
           // daisy::Logger<>::PrintLine("latency=%u, remaining=%u, stage=%u", latency, len, fft_len);
            result &= x.SetIR(&ir[latency], fft_len);
            len -= fft_len;
        });

       // daisy::Logger<>::PrintLine("remaining=%u, result=%u", len, result);
        
        /* all samples should be consumed */
        assert(0 == len);

        return result && (0 == len);
    }

    using PartPlan = PartitionPlan<base, max_size-base, FFTConv<1>::MAX_SUPPORTED_SIZE>;

    template< size_t... ks >
    using FFTConvImpl = std::tuple<FFTConv<ks, true>...>;

    template<size_t... ks>
    static FFTConvImpl<ks...> make_storage(std::index_sequence<ks...>)
    {
        return FFTConvImpl<ks...>{};
    }
    static auto make_storage()
    {
        return make_storage(PartPlan::index_);
    }
    using FFTConvCollection = decltype(make_storage());

    protected:
    FIRFilter<base, 64>    f0_;
    FFTConvCollection       fn_;
};  //class FastConv


template <size_t base, size_t max_size>
class FastConvDyn
{
  public:
    FastConvDyn(int init = 0) 
        : num_stages_(0)
    { 
        assert(is_power2(base)); 
    }

    static size_t GetLatency() { return 0; }

    void Reset()
    {
        f0_.Reset();
        for(int i = 0; i < max_num_stages_; i++)
        {
            fn_[i].Reset();
        }        
    }
    float Process(float in)
    {
        float out = f0_.Process(in);

        for(int i = 0; i < num_stages_; i++) 
        {
            out += fn_[i].Process(in);
        }
        
        return out;
    }

    void ProcessBlock(float* pSrc, float* pDst, size_t num_samples)
    {
        f0_.ProcessBlock(pSrc, pDst, num_samples);
#if 0
        for (size_t i = 0; i < num_samples; i++)
        {
            for_each(fn_, [pSrc, pDst, i] (auto& x) mutable
            {
                pDst[i] += x.Process(pSrc[i]);
            });
        }
#else
        for(int i = 0; i < num_stages_; i++)
        {
            fn_[i].ProcessBlockAcc(pSrc, pDst, num_samples);
        };
#endif
    }

    bool SetIR(const float* ir, size_t len)
    {
        assert(len <= max_size);
        assert(nullptr != ir);
        if(len > max_size || nullptr == ir)
        {
            return false;
        }

        /* number of FFT stages */
        num_stages_ = get_base_2(get_next_power2(len) / base);

        /* zero-latency FIR stage */
        const size_t fir_len = DSY_MIN(base, len);
        bool         result  = f0_.SetIR(ir, fir_len);
        len -= fir_len;

        // base + base*2^0 + base*2^1 + base*2^2 + ... + base*2^(N-1) = len
        // base + base * (1 + 2^1 + 2^2 + ... + 2^(N-1)) = len
        // base + base * (2^(N) - 1 ) = len
        // 2^(N) = len / base

        

        //daisy::Logger<>::PrintLine("SetIR, base=%u, rem_len=%u, num_stages=%d (max=%d)", base, len, num_stages_, max_num_stages_);
        /* progressive-latency FFT stages (compile-time loop) */
        size_t stage_size = base;
        for (int i = 0; i < num_stages_; i++) 
        {
            const size_t subfilter_size = DSY_MIN(stage_size, len);
            result &= fn_[i].SetIR(&ir[stage_size], subfilter_size);
            len -= subfilter_size;
            stage_size *= 2;
        }
        
        // daisy::Logger<>::PrintLine("remaining=%u, result=%u", len, result);

        /* all samples should be consumed */
        assert(0 == len);

        return result && (0 == len);
    }

    using PartPlan
        = PartitionPlan<base, max_size - base, FFTConv<1>::MAX_SUPPORTED_SIZE>;


  protected:
    FIRFilter<base, 64> f0_;
    //FFTConvCollection   fn_;
    static constexpr int max_num_stages_
        = max_size > base
              ? get_base_2(get_next_power2((uint32_t)max_size) / base)
              : 1;
    FFTConv<max_size / 2, false> fn_[max_num_stages_];
    int                          num_stages_;
}; //class FastConv



} // namespace daisysp

#endif // DSY_FASTCONV_H
