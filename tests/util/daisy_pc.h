#pragma once
#ifndef __DAISY_PC_H__
#define __DAISY_PC_H__

#if defined(_WIN32)

#include <cstdint>
#include <cstdio>
#include <cfloat>
#include <windows.h>
#include "daisysp.h"
#include "hid/logger.h"

/**   @brief Quick and dirty PC simulation for running module unit tests
 *    @author Alexander Petrov-Savchenko (axp@soft-amp.com)
 *    @date December 2020
 */


/* fix discrepnacy in msvc/gcc functions */
constexpr auto isnanf = _isnanf;

/* no need for specific allocation */
#define DSY_SDRAM_BSS   /* emtpy */

namespace daisy
{

/* Stub for blocking interrupts, ignore on Windows*/
class ScopedIrqBlocker
{
public:
    ScopedIrqBlocker() { (void)0; }
    ~ScopedIrqBlocker() { (void)0; }
};

/* timestamping services */
uint32_t dsy_tim_get_tick()
{
    LARGE_INTEGER q;
    QueryPerformanceCounter(&q);
    return q.LowPart;
}
uint32_t dsy_tim_get_ticks_per_us()
{
    LARGE_INTEGER q;
    QueryPerformanceFrequency(&q);
    /* use division with rounding */
    uint32_t us = (uint32_t)((q.QuadPart + 500000) / 1000000); 
    return us;
}

/** Simple simulatiion of Daisy platform on windows 
  * TODO: consider moving to libDaisy
  * Add services as needed
  */

class DaisyPC
{
public:
    DaisyPC()
    {}
    ~DaisyPC()
    {}
    /* redirect the log to stdout */
    using Log = daisy::Logger<daisy::LOGGER_SEMIHOST>;

    /* Import logging functions */
    static constexpr auto Print = Log::Print;
    static constexpr auto PrintLine = Log::PrintLine;

    static void StartLog(bool wait_for_pc = false)
    {
    }
};

} // namespace daisysp

#endif // _WIN32

#endif //__DAISY_PC_H__
