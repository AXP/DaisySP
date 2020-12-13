#pragma once
#ifndef DSY_FFT_TABLES_H
#define DSY_FFT_TABLES_H

#include <stdint.h>
#include <arm_math.h>



#ifdef __arm__

#define ARMBITREVINDEXTABLE_16_TABLE_LENGTH ((uint16_t)20)
#define ARMBITREVINDEXTABLE_32_TABLE_LENGTH ((uint16_t)48)
#define ARMBITREVINDEXTABLE_64_TABLE_LENGTH ((uint16_t)56)
#define ARMBITREVINDEXTABLE_128_TABLE_LENGTH ((uint16_t)208)
#define ARMBITREVINDEXTABLE_256_TABLE_LENGTH ((uint16_t)440)
#define ARMBITREVINDEXTABLE_512_TABLE_LENGTH ((uint16_t)448)
#define ARMBITREVINDEXTABLE_1024_TABLE_LENGTH ((uint16_t)1800)
#define ARMBITREVINDEXTABLE_2048_TABLE_LENGTH ((uint16_t)3808)
#define ARMBITREVINDEXTABLE_4096_TABLE_LENGTH ((uint16_t)4032)

extern const uint16_t  armBitRevIndexTable16[ARMBITREVINDEXTABLE_16_TABLE_LENGTH];
extern const uint16_t  armBitRevIndexTable32[ARMBITREVINDEXTABLE_32_TABLE_LENGTH];
extern const uint16_t  armBitRevIndexTable64[ARMBITREVINDEXTABLE_64_TABLE_LENGTH];
extern const uint16_t  armBitRevIndexTable128[ARMBITREVINDEXTABLE_128_TABLE_LENGTH];
extern const uint16_t  armBitRevIndexTable256[ARMBITREVINDEXTABLE_256_TABLE_LENGTH];
extern const uint16_t  armBitRevIndexTable512[ARMBITREVINDEXTABLE_512_TABLE_LENGTH];
extern const uint16_t  armBitRevIndexTable1024[ARMBITREVINDEXTABLE_1024_TABLE_LENGTH];
extern const uint16_t  armBitRevIndexTable2048[ARMBITREVINDEXTABLE_2048_TABLE_LENGTH];
extern const uint16_t  armBitRevIndexTable4096[ARMBITREVINDEXTABLE_4096_TABLE_LENGTH];


namespace daisysp
{
class ARMFFTTables
{
  public:
    ARMFFTTables();

  protected:
    static void InitAllTables();
    static void InitBitRev(uint16_t* tbl, int bits);
    static void InitRFFT(float* tw, int size);
    static void InitCFFT(float* tw, int size);
    static bool init_done_;
};


} // namespace daisysp


#endif // if __arm__

#endif // DSY_FFT_TABLES_H
