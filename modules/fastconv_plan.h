#pragma once
#ifndef DSY_FASTCONV_PLAN_H
#define DSY_FASTCONV_PLAN_H

#include <cstdint>
#include <type_traits>
#include "dsp.h"


#if defined(__cplusplus)


namespace daisysp
{
namespace FastConvolution
{
    class IndexSeqUtils
    {
      public:
        template <size_t... Args>
        static constexpr size_t CalcSum(std::index_sequence<Args...>)
        {
            return add(Args...);
        }

      protected:
        static constexpr size_t add(size_t head) { return head; }
        template <typename... Args>
        static constexpr size_t add(size_t head, Args... tail)
        {
            return head + add(tail...);
        }

        static constexpr size_t add() { return 0; }
    };


    template <size_t base, size_t size, size_t max_block>
    class PartitionPlan
    {
      public:
        template <size_t block, size_t... Ns>
        struct plan_factory;

        template <size_t block, size_t remaining, size_t... Ns>
        struct plan_factory<block, remaining, Ns...>
        {
            using type =
                typename plan_factory<(DSY_MIN(block * 2, max_block)),
                                      (DSY_MAX((int)remaining - (int)block, 0)),
                                      Ns...,
                                      block>::type;
        };

        template <size_t block, size_t... Ns>
        struct plan_factory<block, 0, 0, Ns...>
        {
            using type = std::index_sequence<Ns...>;
        };

        using plan_t = typename plan_factory<base, size, 0>::type;
        static plan_t index_;


        enum Consts
        {
            FFT_INSTANCES = index_.size(),
            PLAN_FIR_PART = base,
            //        std::enable_if<(size > base)> PLAN_FFT_PART = IndexSeqUtils::CalcSum(index_),
            //      std::enabke_if<(size <= base)> PLAN_FFT_PART = 0,
            //PLAN_FFT_PART = size > base ? IndexSeqUtils::CalcSum(index_) : 0,
            PLAN_FFT_PART = IndexSeqUtils::CalcSum(index_),
            PLAN_TOTAL    = PLAN_FFT_PART + PLAN_FIR_PART,
        };

        constexpr PartitionPlan()
        {
            static_assert(PLAN_TOTAL >= size);
            static_assert(is_power2(max_block));
            static_assert(is_power2(base));
        }
    };


} // namespace FastConvolution
} // namespace daisysp
#endif // ifdef __cplusplus
#endif // DSY_FASTCONV_PLAN_H
