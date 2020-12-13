# Project Name
TARGET = tst_fastconv


# Sources
CPP_SOURCES = tst_fastconv.cpp	\


	



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

