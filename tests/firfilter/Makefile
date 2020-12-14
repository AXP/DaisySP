# Project Name
TARGET = tst_firfilter


# Sources
CPP_SOURCES = tst_firfilter.cpp	\

C_INCLUDES = -I./ -I../util/

C_DEFS =  \
-fverbose-asm	\
-DNDEBUG	\


# Library Locations
LIBDAISY_DIR ?= ../../../libdaisy
DAISYSP_DIR ?= ../../../DaisySP

# Core location, and generic Makefile.
SYSTEM_FILES_DIR = $(LIBDAISY_DIR)/core
include $(SYSTEM_FILES_DIR)/Makefile

