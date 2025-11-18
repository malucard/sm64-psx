# Makefile to rebuild SM64 split image

include util.mk

.RECIPEPREFIX = >

# Default target
default: all

# Preprocessor definitions
DEFINES :=

MAKEFLAGS += --jobs=$(shell nproc) --no-builtin-rules

#==============================================================================#
# Build Options                                                                #
#==============================================================================#

# VERSION - selects the version of the game to build
#   jp - builds the 1996 Japanese version
#   us - builds the 1996 North American version
#   eu - builds the 1997 PAL version
#   sh - builds the 1997 Japanese Shindou version, with rumble pak support
VERSION ?= us
$(eval $(call validate-option,VERSION,jp us eu sh))

ifeq ($(VERSION),jp)
	DEFINES += VERSION_JP=1
	VERSION_JP_US ?= true
else ifeq ($(VERSION),us)
	DEFINES += VERSION_US=1
	VERSION_JP_US ?= true
else ifeq ($(VERSION),eu)
	DEFINES += VERSION_EU=1
	VERSION_JP_US ?= false
else ifeq ($(VERSION),sh)
	DEFINES += VERSION_SH=1
	VERSION_JP_US ?= false
endif

BIG_RAM ?= 0
BENCH ?= 0
ifneq ($(BENCH),0)
	DEFINES += BENCH=1
	BIG_RAM := 1
endif
ifneq ($(BIG_RAM),0)
	DEFINES += BIG_RAM=1
endif
MARIO_HEAD ?= 0
ifneq ($(MARIO_HEAD),0)
	DEFINES += MARIO_HEAD=1
endif

DEFINES += F3D_OLD=1 NON_MATCHING=1 AVOID_UB=1 NO_AUDIO=1

# Whether to hide commands or not
VERBOSE ?= 0
ifeq ($(VERBOSE),0)
	V := @
endif

# Whether to colorize build messages
COLOR ?= 1

# display selected options unless 'make clean' or 'make distclean' is run
ifeq ($(filter clean distclean,$(MAKECMDGOALS)),)
	$(info ==== Build Options ====)
	$(info Platform:           psx)
	$(info Region:              $(VERSION))
	$(info =======================)
endif

BUILD_DIR_BASE := build

#==============================================================================#
# Universal Dependencies                                                       #
#==============================================================================#

TOOLS_DIR := tools

# This is a bit hacky, but a lot of rules implicitly depend
# on tools and assets, and we use directory globs in the makefiles

PYTHON := python3

ifeq ($(filter clean distclean print-%,$(MAKECMDGOALS)),)
	# Make sure assets exist
	NOEXTRACT ?= 0
	ifeq ($(NOEXTRACT),0)
		DUMMY != $(PYTHON) extract_assets.py $(VERSION) >&2 || echo FAIL
		ifeq ($(DUMMY),FAIL)
			$(error Failed to extract assets)
		endif
	endif

	# Make tools if out of date
	$(info Building tools...)
	DUMMY != $(MAKE) -s -C $(TOOLS_DIR) all-except-recomp >&2 || echo FAIL
	ifeq ($(DUMMY),FAIL)
		$(error Failed to build tools)
	endif
	$(info Building...)
endif

ifeq ($(SATURN),1)
	include Makefile.ss.mk
else ifeq ($(PC),1)
	include Makefile.pc.mk
else
	include Makefile.psx.mk
endif

PSXAVENC := tools/psxavenc/build/psxavenc
MKPSXISO := tools/mkpsxiso
