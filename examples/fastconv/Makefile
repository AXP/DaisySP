# Project Name
TARGET = ex_fastconv


# Sources
CPP_SOURCES = ex_fastconv.cpp
C_SOURCES = $(LIBDAISY_DIR)/Drivers/CMSIS/DSP/Source/TransformFunctions/arm_rfft_fast_f32.c \
			$(LIBDAISY_DIR)/Drivers/CMSIS/DSP/Source/TransformFunctions/arm_rfft_fast_init_f32.c \
			$(LIBDAISY_DIR)/Drivers/CMSIS/DSP/Source/TransformFunctions/arm_cfft_f32.c \
			$(LIBDAISY_DIR)/Drivers/CMSIS/DSP/Source/CommonTables/arm_common_tables.c   \
			$(LIBDAISY_DIR)/Drivers/CMSIS/DSP/Source/TransformFunctions/arm_cfft_radix8_f32.c   \
			$(LIBDAISY_DIR)/Drivers/CMSIS/DSP/Source/FilteringFunctions/arm_fir_f32.c   \
			$(LIBDAISY_DIR)/Drivers/CMSIS/DSP/Source/FilteringFunctions/arm_fir_init_f32.c  \
			$(LIBDAISY_DIR)/Drivers/CMSIS/DSP/Source/ComplexMathFunctions/arm_cmplx_mult_cmplx_f32.c \
			$(LIBDAISY_DIR)/Drivers/CMSIS/DSP/Source/SupportFunctions/arm_copy_f32.c	\
			$(LIBDAISY_DIR)/Drivers/CMSIS/DSP/Source/SupportFunctions/arm_fill_f32.c

ASM_SOURCES = $(LIBDAISY_DIR)/Drivers/CMSIS/DSP/Source/TransformFunctions/arm_bitreversal2.S


C_INCLUDES = -I$(LIBDAISY_DIR)/Drivers/CMSIS/DSP/Include/ \
			-include stm32h750xx.h

C_DEFS =  \
-DCORE_CM7  \
-DSTM32H750IB \
-DARM_MATH_CM7 \

# Library Locations
LIBDAISY_DIR ?= ../../../libdaisy
DAISYSP_DIR ?= ../../../DaisySP

# Core location, and generic Makefile.
SYSTEM_FILES_DIR = $(LIBDAISY_DIR)/core
include $(SYSTEM_FILES_DIR)/Makefile

