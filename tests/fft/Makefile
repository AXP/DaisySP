# Project Name
TARGET = tst_fft


# Sources
CPP_SOURCES = tst_fft.cpp	\
			
CMSIS_DIR ?= $(LIBDAISY_DIR)/Drivers/CMSIS
C_SOURCES += $(CMSIS_DIR)/DSP/Source/BasicMathFunctions/arm_scale_f32.c		


C_INCLUDES = -I.  -I../util/


C_DEFS =  \
-DCORE_CM7  \
-DSTM32H750IB \
-DARM_MATH_CM7 \
-fverbose-asm	\
-DNDEBUG	\


# Library Locations
LIBDAISY_DIR ?= ../../../libdaisy
DAISYSP_DIR ?= ../../../DaisySP

# Core location, and generic Makefile.
SYSTEM_FILES_DIR = $(LIBDAISY_DIR)/core
include $(SYSTEM_FILES_DIR)/Makefile

