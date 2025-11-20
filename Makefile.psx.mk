# this file is included by the main Makefile

#==============================================================================#
# Target Executable and Sources                                                #
#==============================================================================#

# BUILD_DIR is the location where all build artifacts are placed
BUILD_DIR      := $(BUILD_DIR_BASE)/$(VERSION)_psx
EXE            := $(BUILD_DIR)/sm64.exe
ELF            := $(BUILD_DIR)/sm64.elf
#ifeq ($(BENCH),0)
LEVEL_DIRS := $(patsubst levels/%,%,$(dir $(wildcard levels/*/header.h)))
ifneq ($(BENCH),0)
	LEVEL_DIRS := $(filter-out castle_inside/ wmotr/ wdw/ ssl/,$(LEVEL_DIRS))
endif
#else
#	LEVEL_DIRS := bob/ ttm/
#endif

# Directories containing source files
SRC_DIRS := src src/engine src/game src/audio src/menu src/buffers actors levels bin bin/$(VERSION) data assets sound
SRC_DIRS += src/port src/port/gfx src/port/dma src/port/psx
SRC_DIRS += ps1-bare-metal/ps1 ps1-bare-metal/vendor
LIBC_SRC_DIRS := ps1-bare-metal/libc

BIN_DIRS := bin bin/$(VERSION)

GODDARD_SRC_DIRS := src/goddard src/goddard/dynlists

# File dependencies and variables for specific files
include Makefile.split.mk

# Source code files
LEVEL_C_FILES        := $(foreach dir,$(LEVEL_DIRS),levels/$(dir)leveldata.c levels/$(dir)geo.c levels/$(dir)script.c)
C_FILES              := $(foreach dir,$(SRC_DIRS),$(wildcard $(dir)/*.c))
LIBC_C_FILES         := $(foreach dir,$(LIBC_SRC_DIRS),$(wildcard $(dir)/*.c))
LIBC_S_FILES         := $(foreach dir,$(LIBC_SRC_DIRS),$(wildcard $(dir)/*.s))
CXX_FILES            := $(foreach dir,$(SRC_DIRS),$(wildcard $(dir)/*.cpp))
S_FILES              := $(foreach dir,$(SRC_DIRS),$(wildcard $(dir)/*.s))
SPP_FILES            := $(foreach dir,$(SRC_DIRS),$(wildcard $(dir)/*.S))
GODDARD_C_FILES      := $(foreach dir,$(GODDARD_SRC_DIRS),$(wildcard $(dir)/*.c))
RAW_MODEL_C_FILES    := $(wildcard actors/*/model.inc.c) $(foreach dir,$(LEVEL_DIRS),$(wildcard levels/$(dir)*/model.inc.c) $(wildcard levels/$(dir)areas/*/*/model.inc.c) $(wildcard levels/$(dir)areas/*/model.inc.c))
CONV_MODEL_C_FILES   := $(RAW_MODEL_C_FILES:%.c=$(BUILD_DIR)/%.processed.c)
GENERATED_C_FILES    := $(BUILD_DIR)/sfx_defs.generated.c
#$(foreach name,$(notdir $(basename $(wildcard textures/skyboxes/*.png))),$(BUILD_DIR)/bin/$(name)_skybox.c)

C_FILES := $(filter-out src/game/main.c,$(C_FILES))

# Object files
O_FILES := \
	$(C_FILES:%.c=$(BUILD_DIR)/%.o) \
	$(CXX_FILES:%.cpp=$(BUILD_DIR)/%.o) \
	$(S_FILES:%.s=$(BUILD_DIR)/%.s.o) \
	$(GENERATED_C_FILES:%.c=%.o)

LEVEL_O_FILES := $(LEVEL_C_FILES:%.c=$(BUILD_DIR)/%.o)

GODDARD_O_FILES := $(GODDARD_C_FILES:%.c=$(BUILD_DIR)/%.o)

LIBC_O_FILES := \
	$(LIBC_C_FILES:%.c=$(BUILD_DIR)/%.o) \
	$(LIBC_S_FILES:%.s=$(BUILD_DIR)/%.s.o)
	#$(LIBC_C_FILES:%.c=$(BUILD_DIR)/%.libc.o) \
	#$(LIBC_S_FILES:%.s=$(BUILD_DIR)/%.s.libc.o)

# Automatic dependency files
DEP_FILES := $(O_FILES:.o=.d) $(LEVEL_O_FILES:.o=.d) $(GODDARD_O_FILES:.o=.d) $(LIBC_O_FILES:.o=.d)

#==============================================================================#
# Compiler Options                                                             #
#==============================================================================#

INCLUDE_DIRS := include $(BUILD_DIR) $(BUILD_DIR)/include src . include/libc ps1 ps1-bare-metal ps1-bare-metal/libc

# debugging toggles:
# USE_FLOATS: revert some of the fixed point math replacements to float versions (hasn't been needed for a while!)
# NO_INLINE: disable inlining and IPA in a few critical functions
# NO_SCRATCHPAD: disable scratchpad usage
# NO_KERNEL_RAM: disable abuse of the kernel area in ram (the first 64kb)
# BIG_RAM: make use of 8MB mode, also store all assets in main ram
# SAFE_GTE: inserts nops before GTE commands
ifeq ($(SAFE),1)
	DEFINES += NO_SCRATCHPAD=1 NO_KERNEL_RAM=1 BIG_RAM=1 SAFE_GTE=1
else ifeq ($(DEV),1)
	DEFINES += BIG_RAM=1
endif

ifeq ($(SERIAL),1)
	DEFINES += SERIAL=1
endif

DEFINES += TARGET_PSX=1 ENABLE_RUMBLE=1 RUMBLE_GRAPHIC=1
C_DEFINES := $(foreach d,$(DEFINES),-D$(d))
DEF_INC_CFLAGS := $(foreach i,$(INCLUDE_DIRS),-I$(i)) $(C_DEFINES)

ifneq      ($(call find-command,mipsel-none-elf-gcc),)
	CROSS := mipsel-none-elf-
else ifneq      ($(call find-command,mipsel-unknown-elf-gcc),)
	CROSS := mipsel-unknown-elf-
else ifneq      ($(call find-command,mipsel-linux-gnu-gcc),)
	CROSS := mipsel-linux-gnu-
else ifneq      ($(call find-command,mipsel-unknown-linux-gnu-gcc),)
	CROSS := mipsel-unknown-linux-gnu-
else ifneq      ($(call find-command,mipsel-elf-gcc),)
	CROSS := mipsel-elf-
else
	$(error Unable to detect a suitable MIPS toolchain installed)
endif

EXT_LD    := $(CROSS)ld
AS        := $(CROSS)as
AR        := $(CROSS)ar
NM		  := $(CROSS)nm
OBJDUMP   := $(CROSS)objdump
OBJCOPY   := $(CROSS)objcopy

CC := $(CROSS)gcc
LD := $(CC)

TARGET_CFLAGS := -Wall -Wextra -Werror -Wno-maybe-uninitialized -Wno-error=cpp -Wno-comment -Wno-unused-const-variable

# despite -Os/-Oz being commonly recommended due to code size savings, it runs much better with -O2
TARGET_CFLAGS += -O2

# necessary ABI and environment flags
TARGET_CFLAGS += -march=r3000 -mtune=r3000 \
	-mabi=eabi -mno-abicalls -EL -fcall-used-k0 -fcall-used-k1 -freg-struct-return \
	-mfp32 -msingle-float -fsingle-precision-constant -mno-fp-exceptions -msoft-float \
	-fno-builtin -nostdinc -nostdlib -mno-mt -fno-pic -fno-PIC -fsigned-char \
	-static -mno-shared -fomit-frame-pointer -fno-stack-protector -mno-llsc \
	-ffreestanding -mno-extern-sdata -fno-common -fno-zero-initialized-in-bss \
	-include include/stdbool.h $(DEF_INC_CFLAGS) -D_LANGUAGE_C

TARGET_CFLAGS += -ffunction-sections -fdata-sections -Wl,--gc-sections

# this can hurt performance (appears to be no longer needed)
TARGET_CFLAGS += -fconserve-stack

# objectively correct optimization flags
TARGET_CFLAGS += -fstrict-aliasing -fstrict-overflow -mno-check-zero-division -ffast-math -ffp-contract=fast \
	-flimit-function-alignment -falign-functions=16:8 -fno-align-labels -fno-align-jumps -falign-loops=16:8 \
	--param l1-cache-line-size=0 --param l1-cache-size=0 --param l2-cache-size=0 \
	-fsection-anchors -Wa,--strip-local-absolute -fallow-store-data-races -favoid-store-forwarding \
	-fno-semantic-interposition

# experiment zone
TARGET_CFLAGS += -free -fira-loop-pressure -fpredictive-commoning -fsched-pressure -fsched-spec-load \
	-ftree-pre -ftree-partial-pre -ftree-loop-im -ftree-loop-distribution -floop-interchange -freorder-blocks-algorithm=simple \
	-fweb -frename-registers -fgcse-sm -fgcse-las -fgcse-after-reload \
	-fipa-pta -fipa-icf -fipa-reorder-for-locality -fipa-bit-cp -fipa-vrp \
	-favoid-store-forwarding -fno-prefetch-loop-arrays

# -O3 territory 🥶
#TARGET_CFLAGS += -fgcse-after-reload -fipa-cp-clone -ftracer \
#	-floop-interchange -floop-unroll-and-jam -fpeel-loops \
#	-fsplit-loops -fsplit-paths -funswitch-loops -fversion-loops-for-strides

ifeq ($(SAFE),1)
	TARGET_CFLAGS += -fsanitize=unreachable,bounds-strict,pointer-overflow -UNDEBUG
else ifeq ($(DEV),1)
	TARGET_CFLAGS += -fsanitize=unreachable,bounds-strict,pointer-overflow -UNDEBUG
else
	TARGET_CFLAGS += -DNDEBUG
endif

ifeq ($(DEV),1)
	TARGET_CFLAGS += -Og
endif

# C compiler options
CFLAGS := -std=gnu2x --embed-dir=$(BUILD_DIR) $(TARGET_CFLAGS) $(if $(filter 1,$(SAFE) $(DEV)),-G8,-G32)
ifeq ($(DEV),1)
	CFLAGS += -fno-lto
else
	$(shell mkdir -p "$(BUILD_DIR)/lto_incremental")
	CFLAGS += -flto -fno-fat-lto-objects -fwhole-program -flto-incremental=$(BUILD_DIR)/lto_incremental
endif
EXT_CFLAGS := -std=gnu2x $(TARGET_CFLAGS) -G0 -fno-lto

CFLAGS_FILE := $(BUILD_DIR)/cflags.txt
$(shell mkdir -p $(dir $(CFLAGS_FILE)))
ifneq ($(CFLAGS),$(file <$(CFLAGS_FILE)))
	$(file >$(CFLAGS_FILE),$(CFLAGS))
endif

ASFLAGS := -O2 -march=r3000 -mabi=eabi -msoft-float $(foreach i,$(INCLUDE_DIRS),-I$(i)) $(foreach d,$(DEFINES),--defsym $(d))
LDFLAGS := $(CFLAGS) -EL -Wl,-Map,$(BUILD_DIR)/sm64.map -Lps1-bare-metal -T$(BUILD_DIR)/executable.preprocessed.ld

CPP      := $(CC) -E -x c
CPPFLAGS := -P -Wno-trigraphs $(DEF_INC_CFLAGS)

#==============================================================================#
# Miscellaneous Tools                                                          #
#==============================================================================#

# N64 tools
#MIO0TOOL              := $(TOOLS_DIR)/mio0
#N64CKSUM              := $(TOOLS_DIR)/n64cksum
#N64GRAPHICS           := $(TOOLS_DIR)/n64graphics
#N64GRAPHICS_CI        := $(TOOLS_DIR)/n64graphics_ci
TEXTCONV              := $(TOOLS_DIR)/textconv
#AIFF_EXTRACT_CODEBOOK := $(TOOLS_DIR)/aiff_extract_codebook
#VADPCM_ENC            := $(TOOLS_DIR)/vadpcm_enc
EXTRACT_DATA_FOR_MIO  := $(TOOLS_DIR)/extract_data_for_mio
SKYCONV               := $(TOOLS_DIR)/skyconv
PRINT = printf

ifeq ($(COLOR),1)
	NO_COL  := \033[0m
	RED     := \033[0;31m
	GREEN   := \033[0;32m
	BLUE    := \033[0;34m
	YELLOW  := \033[0;33m
	BLINK   := \033[33;5m
endif

# Common build print status function
define print
	@$(PRINT) "$(GREEN)$(1) $(YELLOW)$(2)$(GREEN) -> $(BLUE)$(3)$(NO_COL)\n"
endef

#==============================================================================#
# Main Targets                                                                 #
#==============================================================================#
ISO_OUT := $(BUILD_DIR)/sm64.iso
ifeq ($(BENCH),0)
all: $(ISO_OUT)
else
all: $(EXE)
endif

clean:
>	$(RM) -r $(BUILD_DIR_BASE)

distclean: clean
>	$(PYTHON) extract_assets.py --clean
>	$(MAKE) -C $(TOOLS_DIR) clean

#$(BUILD_DIR)/asm/ipl3_font.o:         $(IPL3_RAW_FILES)
#$(BUILD_DIR)/src/game/crash_screen.o: $(CRASH_TEXTURE_C_FILES)
#$(BUILD_DIR)/lib/rsp.o:               $(BUILD_DIR)/rsp/rspboot.bin $(BUILD_DIR)/rsp/fast3d.bin $(BUILD_DIR)/rsp/audio.bin
#$(SOUND_BIN_DIR)/sound_data.o:        $(SOUND_BIN_DIR)/sound_data.ctl.inc.c $(SOUND_BIN_DIR)/sound_data.tbl.inc.c $(SOUND_BIN_DIR)/sequences.bin.inc.c $(SOUND_BIN_DIR)/bank_sets.inc.c
$(BUILD_DIR)/levels/scripts.o:        $(BUILD_DIR)/include/level_headers.h

$(BUILD_DIR)/lib/src/math/%.o: CFLAGS += -fno-builtin

ifeq ($(VERSION),eu)
	TEXT_DIRS := text/de text/us text/fr

	# EU encoded text inserted into individual segment 0x19 files,
	# and course data also duplicated in leveldata.c
	$(BUILD_DIR)/bin/eu/translation_en.o2: $(BUILD_DIR)/text/us/define_text.inc.c
	$(BUILD_DIR)/bin/eu/translation_de.o2: $(BUILD_DIR)/text/de/define_text.inc.c
	$(BUILD_DIR)/bin/eu/translation_fr.o2: $(BUILD_DIR)/text/fr/define_text.inc.c
	$(BUILD_DIR)/levels/menu/leveldata.o2: $(BUILD_DIR)/include/text_strings.h
	$(BUILD_DIR)/levels/menu/leveldata.o2: $(BUILD_DIR)/text/us/define_courses.inc.c
	$(BUILD_DIR)/levels/menu/leveldata.o2: $(BUILD_DIR)/text/de/define_courses.inc.c
	$(BUILD_DIR)/levels/menu/leveldata.o2: $(BUILD_DIR)/text/fr/define_courses.inc.c
else
	ifeq ($(VERSION),sh)
		TEXT_DIRS := text/jp
		$(BUILD_DIR)/bin/segment2.o2: $(BUILD_DIR)/text/jp/define_text.inc.c
	else
		TEXT_DIRS := text/$(VERSION)
		# non-EU encoded text inserted into segment 0x02
		$(BUILD_DIR)/bin/segment2.o2: $(BUILD_DIR)/text/$(VERSION)/define_text.inc.c
	endif
endif

ALL_DIRS := $(BUILD_DIR) $(addprefix $(BUILD_DIR)/,$(SRC_DIRS) $(GODDARD_SRC_DIRS) $(LIBC_SRC_DIRS) $(BIN_DIRS) textures $(TEXT_DIRS) $(SOUND_SAMPLE_DIRS) $(addprefix levels/,$(c)) rsp include)

# Make sure build directory exists before compiling anything
DUMMY != mkdir -p $(ALL_DIRS)

$(BUILD_DIR)/include/text_strings.h: $(BUILD_DIR)/include/text_menu_strings.h
$(BUILD_DIR)/src/menu/file_select.o: $(BUILD_DIR)/include/text_strings.h
$(BUILD_DIR)/src/menu/star_select.o: $(BUILD_DIR)/include/text_strings.h
$(BUILD_DIR)/src/game/ingame_menu.o: $(BUILD_DIR)/include/text_strings.h

#==============================================================================#
# Texture Generation                                                           #
#==============================================================================#

# Convert PNGs to a specialized format
$(BUILD_DIR)/%.fulldata: %.png tools/convert_image_psx
>	$(call print,Converting:,$<,$@)
>	$(V)mkdir -p $(dir $@)
>	$(V)./tools/convert_image_psx 4 $< $@

PNGS_WITH_BLACK_TRANSPARENCY := levels/intro/3_tm.rgba16.png

$(PNGS_WITH_BLACK_TRANSPARENCY:%.png=$(BUILD_DIR)/%.fulldata): $(BUILD_DIR)/%.fulldata: %.png tools/convert_image_psx
>	$(call print,Converting:,$<,$@)
>	$(V)mkdir -p $(dir $@)
>	$(V)./tools/convert_image_psx 4 $< $@

ALL_PNGS := $(foreach png,$(filter-out %/cake.png %/cake_eu.png %/skyboxes/%.png,$(filter %.png,$(file <.assets-local.txt))),$(wildcard $(png))) dualshock_graphic.png
ALL_FULLDATAS := $(ALL_PNGS:%.png=$(BUILD_DIR)/%.fulldata)

$(BUILD_DIR)/tex_pack: $(ALL_FULLDATAS) tools/pack_textures.py
>	@$(PRINT) "$(GREEN)Packing all images$(NO_COL)\n"
>	$(V)rm -f $(BUILD_DIR)/fulldata_list.txt
>	$(V)for fulldata in $(ALL_FULLDATAS); do \
>		echo $$fulldata >> $(BUILD_DIR)/fulldata_list.txt ;\
>	done
>	$(V)$(PYTHON) tools/pack_textures.py $@.tmp $(BUILD_DIR)/fulldata_list.txt
>	$(V)for png in $(ALL_PNGS); do \
>		hexdump -v -e '1/1 "0x%X,"' $(BUILD_DIR)/$${png%.png}.texheader > $(BUILD_DIR)/$${png%.png}.inc.c ;\
>	done
>	$(V)mv $@.tmp $@

ALL_TEX_HEADER_FILES := $(ALL_PNGS:%.png=$(BUILD_DIR)/%.inc.c)
#ifeq ($(BENCH),0)
	CAKE_INC_C := $(if $(filter eu,$(VERSION)),$(BUILD_DIR)/levels/ending/cake_eu.inc.c,$(BUILD_DIR)/levels/ending/cake.inc.c)
#endif

$(ALL_TEX_HEADER_FILES):

#==============================================================================#
# Audio Generation                                                             #
#==============================================================================#

# background music
BGM_TRACKS := 0 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 26 27 28 29 30 31 32 33 34 35 36 37

ifeq ($(wildcard .local/0.wav),)
	$(warning music not found!)
	BGM_TRACKS := $(foreach track,$(BGM_TRACKS),dummy$(track))
endif

$(BUILD_DIR)/bgm/dummy%.track.xa: empty.wav $(PSXAVENC)
>	$(V)mkdir -p $(dir $@)
>	$(V)$(PSXAVENC) -t xa -f 37800 -b 8 -c 2 -F $(*:track=) -C 0 $< $@

$(BUILD_DIR)/bgm/%.track.xa: .local/%.wav $(PSXAVENC)
>	$(V)mkdir -p $(dir $@)
>	$(V)$(PSXAVENC) -t xa -f 37800 -b 8 -c 2 -F $(*:track=) -C 0 $< $@

$(BUILD_DIR)/bgm/pack.xa: $(TOOLS_DIR)/interleave_xa.py $(foreach track,$(BGM_TRACKS),$(BUILD_DIR)/bgm/$(track).track.xa)
>	$(V)mkdir -p $(dir $@)
>	$(V)$(PYTHON) $(TOOLS_DIR)/interleave_xa.py $(BUILD_DIR)/bgm/info.dat $(BUILD_DIR)/bgm/pack.xa $(filter-out %.py,$^)

# Sound files
SOUND_BANK_FILES    := $(wildcard sound/sound_banks/*.json)
SOUND_SAMPLE_DIRS   := $(wildcard sound/samples/*)
SOUND_SAMPLE_AIFFS  := $(foreach dir,$(SOUND_SAMPLE_DIRS),$(wildcard $(dir)/*.aiff))
#SOUND_SAMPLE_TABLES := $(foreach file,$(SOUND_SAMPLE_AIFFS),$(BUILD_DIR)/$(file:.aiff=.table))
#SOUND_SAMPLE_AIFCS  := $(foreach file,$(SOUND_SAMPLE_AIFFS),$(BUILD_DIR)/$(file:.aiff=.aifc))
#SOUND_SEQUENCE_DIRS := sound/sequences sound/sequences/$(VERSION)
# all .m64 files in SOUND_SEQUENCE_DIRS, plus all .m64 files that are generated from .s files in SOUND_SEQUENCE_DIRS
#SOUND_SEQUENCE_FILES := \
#	$(foreach dir,$(SOUND_SEQUENCE_DIRS),\
#		$(wildcard $(dir)/*.m64) \
#		$(foreach file,$(wildcard $(dir)/*.s),$(BUILD_DIR)/$(file:.s=.m64)) \
#	)

$(BUILD_DIR)/%.samplebin: %.aiff $(TOOLS_DIR)/psx_sample_gen.c
>	$(V)mkdir -p $(dir $@)
>	$(V)$(TOOLS_DIR)/psx_sample_gen $< $@

SOUND_SAMPLE_BINS := $(SOUND_SAMPLE_AIFFS:%.aiff=$(BUILD_DIR)/%.samplebin)

$(BUILD_DIR)/soundtable $(BUILD_DIR)/sounddata &: $(SOUND_SAMPLE_BINS) $(SOUND_BANK_FILES) $(TOOLS_DIR)/psx_sample_pack.py
>	$(V)$(PYTHON) $(TOOLS_DIR)/psx_sample_pack.py $(VERSION) sound/sound_banks $(BUILD_DIR)/sound/samples $(BUILD_DIR)/soundtable.tmp $(BUILD_DIR)/sounddata.tmp
>	$(V)mv $(BUILD_DIR)/soundtable.tmp $(BUILD_DIR)/soundtable
>	$(V)mv $(BUILD_DIR)/sounddata.tmp $(BUILD_DIR)/sounddata

$(BUILD_DIR)/sfx_defs.generated.c: sound/sequences/00_sound_player.s $(TOOLS_DIR)/sound_player_to_c.py
>	$(call print,Compiling sound effect definitions:,$<,$@)
>	$(V)$(PYTHON) $(TOOLS_DIR)/sound_player_to_c.py $(VERSION) $< $@

#==============================================================================#
# Segment Generation                                                           #
#==============================================================================#

# Link segment file to resolve external labels
$(BUILD_DIR)/%.elf: $(BUILD_DIR)/%.o2
>	$(call print,Linking asset ELF file (at $(SEGMENT_ADDRESS)):,$<,$@)
>	$(V)$(EXT_LD) -e 0 -Tdata=$(SEGMENT_ADDRESS) -EL -no-pie -G0 -Text_files_elf.ld --unresolved-symbols=ignore-all -Map $@.map -o $@.tmp $<
>	$(V)$(OBJCOPY) -j .data $@.tmp
>	$(V)mv $@.tmp $@
# Override for leveldata.elf, which otherwise matches the above pattern
.SECONDEXPANSION:
$(BUILD_DIR)/levels/%/leveldata.elf: $(BUILD_DIR)/levels/%/leveldata.o2 $(BUILD_DIR)/bin/$$(TEXTURE_BIN).elf
>	$(call print,Linking leveldata ELF file (at $(SEGMENT_ADDRESS)):,$<,$@)
>	$(V)$(EXT_LD) -e 0 -Tdata=$(SEGMENT_ADDRESS) -EL -no-pie -G0 -Text_files_elf.ld --unresolved-symbols=ignore-all -Map $@.map --just-symbols=$(BUILD_DIR)/bin/$(TEXTURE_BIN).elf -o $@.tmp $<
>	$(V)$(OBJCOPY) -j .data $@.tmp
>	$(V)mv $@.tmp $@

.SECONDEXPANSION:
$(BUILD_DIR)/levels/%/scriptgeo.elf: $(BUILD_DIR)/levels/%/script.o2 $(BUILD_DIR)/levels/%/geo.o2 $(BUILD_DIR)/levels/%/leveldata.elf $(GROUP_SEG_FILES) $(BUILD_DIR)/ext_files_defsym_noscriptgeo.txt
>	$(call print,Linking level script & geo ELF file (at $(SEGMENT_ADDRESS)):,$<,$@)
>	$(V)$(EXT_LD) -e 0 -Tdata=$(SEGMENT_ADDRESS) -EL -no-pie -G0 -Text_files_elf.ld --unresolved-symbols=ignore-all -Map $@.map $(addprefix --just-symbols=,$(GROUP_SEG_FILES)) --just-symbols=$(@:scriptgeo.elf=leveldata.elf) `sed "s/-Wl,//g" $(BUILD_DIR)/ext_files_defsym_noscriptgeo.txt` -o $@.tmp $< $(<:script.o2=geo.o2)
>	$(V)$(OBJCOPY) -j .data $@.tmp
>	$(V)mv $@.tmp $@

.SECONDEXPANSION:
$(BUILD_DIR)/levels/intro/scriptgeo.elf: $(BUILD_DIR)/levels/intro/script.o2 $(BUILD_DIR)/levels/intro/geo.o2 $(BUILD_DIR)/levels/menu/scriptgeo.elf $(BUILD_DIR)/levels/intro/leveldata.elf $(GROUP_SEG_FILES) $(BUILD_DIR)/ext_files_defsym_plusmenu.txt
>	$(call print,Linking intro script & geo ELF file (at 0x14000000):,$<,$@)
>	$(V)$(EXT_LD) -e 0 -Tdata=0x14000000 -EL -no-pie -G0 -Text_files_elf.ld -Map $@.map $(addprefix --just-symbols=,$(GROUP_SEG_FILES)) --just-symbols=$(@:scriptgeo.elf=leveldata.elf) --just-symbols=$(BUILD_DIR)/levels/menu/scriptgeo.elf `sed "s/-Wl,//g" $(BUILD_DIR)/ext_files_defsym_plusmenu.txt` -o $@.tmp $< $(<:script.o2=geo.o2)
>	$(V)$(OBJCOPY) -j .data $@.tmp
>	$(V)mv $@.tmp $@

.SECONDEXPANSION:
$(BUILD_DIR)/levels/menu/scriptgeo.elf: $(BUILD_DIR)/levels/menu/script.o2 $(BUILD_DIR)/levels/menu/geo.o2 $(BUILD_DIR)/levels/menu/leveldata.elf $(GROUP_SEG_FILES) $(BUILD_DIR)/ext_files_defsym_noscriptgeo.txt
>	$(call print,Linking menu script & geo ELF file (at 0x14000000):,$<,$@)
>	$(V)$(EXT_LD) -e 0 -Tdata=0x14000000 -EL -no-pie -G0 -Text_files_elf.ld -Map $@.map $(addprefix --just-symbols=,$(GROUP_SEG_FILES)) --just-symbols=$(@:scriptgeo.elf=leveldata.elf) `sed "s/-Wl,//g" $(BUILD_DIR)/ext_files_defsym_noscriptgeo.txt` -o $@.tmp $< $(<:script.o2=geo.o2)
>	$(V)$(OBJCOPY) -j .data $@.tmp
>	$(V)mv $@.tmp $@

.SECONDEXPANSION:
$(BUILD_DIR)/actors/%_geo.elf: $(BUILD_DIR)/actors/%_geo.o2 $(BUILD_DIR)/actors/%.elf
>	$(call print,Linking actor geo ELF file (at $(SEGMENT_ADDRESS)):,$<,$@)
>	$(V)$(EXT_LD) -e 0 -Tdata=$(SEGMENT_ADDRESS) -EL -no-pie -G0 -Text_files_elf.ld -Map $@.map --just-symbols=$(filter-out $<,$^) -o $@.tmp $<
>	$(V)$(OBJCOPY) -j .data $@.tmp
>	$(V)mv $@.tmp $@

$(BUILD_DIR)/%.bin: $(BUILD_DIR)/%.elf
>	$(call print,Extracting compressible data from:,$<,$@)
>	$(V)$(EXTRACT_DATA_FOR_MIO) $< $@

$(BUILD_DIR)/levels/%/leveldata.bin: $(BUILD_DIR)/levels/%/leveldata.asset.elf
>	$(call print,Extracting compressible data from:,$<,$@)
>	$(V)$(EXTRACT_DATA_FOR_MIO) $< $@

# Compress binary file
$(BUILD_DIR)/%.mio0: $(BUILD_DIR)/%.bin
>	$(call print,Compressing:,$<,$@)
>	$(V)$(MIO0TOOL) $< $@

# convert binary mio0 to object file
$(BUILD_DIR)/%.mio0.o2: $(BUILD_DIR)/%.mio0
>	$(call print,Converting MIO0 to ELF:,$<,$@)
>	$(V)$(EXT_LD) -r -b binary $< -o $@

EXT_SYMS_TXT := $(BUILD_DIR)/ext_symbols.txt
$(BUILD_DIR)/%.asset.txt: $(BUILD_DIR)/%.elf
>	@$(NM) -SUP $^ > $@.tmp
>	grep -E "[0-9A-Za-z_]+ [RD] [0-9A-Za-z]+ [0-9A-Za-z]+" $@.tmp > $@
>	@rm -f $@.tmp
$(BUILD_DIR)/levels/%/leveldata.level.txt: $(BUILD_DIR)/levels/%/leveldata.elf
>	@$(NM) -SUP $^ > $@.tmp
>	grep -E "[0-9A-Za-z_]+ [RD] [0-9A-Za-z]+ [0-9A-Za-z]+" $@.tmp > $@
>	@rm -f $@.tmp
#$(EXT_SYMS_TXT): $(MIO0_FILES_NO_LEVELS:.mio0=.asset.txt) $(LEVEL_ELF_FILES:.elf=.level.txt)
#	$(NM) -SUP $(MIO0_FILES:.mio0=.mio0.o) $(LEVEL_O_FILES) > $@.tmp
COMMA = ,
DEFSYM_PREFIX := -Wl,
$(BUILD_DIR)/%.section.txt: $(BUILD_DIR)/%.elf
>	$(V)readelf -S $^ | sed -Enz "s?.*\\.data\\s+PROGBITS\\s+[0-9A-Za-z]+\\s+([0-9A-Za-z]+)\\s+([0-9A-Za-z]+).*?$^:\1:\2!_$(basename $(notdir $^))SegmentRomStart:_$(basename $(notdir $^))SegmentRomEnd ?p" > $@
$(BUILD_DIR)/%.mio0section.txt: $(BUILD_DIR)/%.elf
>	$(V)readelf -S $^ | sed -Enz "s?.*\\.data\\s+PROGBITS\\s+[0-9A-Za-z]+\\s+([0-9A-Za-z]+)\\s+([0-9A-Za-z]+).*?$^:\1:\2!_$(basename $(notdir $^))_mio0SegmentRomStart:_$(basename $(notdir $^))_mio0SegmentRomEnd ?p" > $@
$(BUILD_DIR)/%.leveldatasection.txt: $(BUILD_DIR)/%/leveldata.elf
>	$(V)readelf -S $^ | sed -Enz "s?.*\\.data\\s+PROGBITS\\s+[0-9A-Za-z]+\\s+([0-9A-Za-z]+)\\s+([0-9A-Za-z]+).*?$^:\1:\2!_$(basename $(basename $(notdir $@)))_segment_7SegmentRomStart:_$(basename $(basename $(notdir $@)))_segment_7SegmentRomEnd ?p" > $@
$(BUILD_DIR)/%.levelscriptgeosection.txt: $(BUILD_DIR)/%/scriptgeo.elf
>	$(V)readelf -S $^ | sed -Enz "s?.*\\.data\\s+PROGBITS\\s+[0-9A-Za-z]+\\s+([0-9A-Za-z]+)\\s+([0-9A-Za-z]+).*?$^:\1:\2!_$(basename $(basename $(notdir $@)))SegmentRomStart:_$(basename $(basename $(notdir $@)))SegmentRomEnd ?p" > $@
$(BUILD_DIR)/ext_files_sections_noscriptgeo.txt: $(BIN_SEG_FILES:%.elf=%.mio0section.txt) $(GROUP_SEG_FILES:%.elf=%.mio0section.txt) $(GROUP_SEG_FILES:%.elf=%_geo.section.txt) $(LEVEL_SEG_FILES:%/leveldata.elf=%.leveldatasection.txt)
>	@cat $^ > $@
$(BUILD_DIR)/ext_files_sections_plusmenu.txt: $(BUILD_DIR)/ext_files_sections_noscriptgeo.txt $(BUILD_DIR)/levels/menu.levelscriptgeosection.txt
>	@cat $^ > $@
$(BUILD_DIR)/mario_anim_data.marioanimbin: $(BUILD_DIR)/mario_anim_data.elf $(TOOLS_DIR)/compress_mario_anims
>	@data_offset_and_size=`readelf -S $< | sed -Enz "s?.*\\.data\\s+PROGBITS\\s+[0-9A-Za-z]+\\s+([0-9A-Za-z]+)\\s+([0-9A-Za-z]+).*?0x\1 0x\2?p"` ;\
>	exec $(TOOLS_DIR)/compress_mario_anims $< $@ $$data_offset_and_size
$(BUILD_DIR)/mario_anim_data.marioanimsection.txt: $(BUILD_DIR)/mario_anim_data.marioanimbin
>	@echo -n "$^:0:`printf "%x" \`du -sb $^ | cut -f 1\``!_$(basename $(basename $(notdir $^)))SegmentRomStart:_$(basename $(basename $(notdir $^)))SegmentRomEnd " > $@

HARDCODED_SEGMENTS := -Wl,--defsym=_goddardSegmentRomStart=0 -Wl,--defsym=_goddardSegmentRomEnd=0 -Wl,--defsym=_goddardSegmentStart=0 -Wl,--defsym=_scriptsSegmentRomStart=0 -Wl,--defsym=_scriptsSegmentRomEnd=0 -Wl,--defsym=_behaviorSegmentRomStart=0 -Wl,--defsym=_behaviorSegmentRomEnd=0

$(BUILD_DIR)/ext_files_defsym_noscriptgeo.txt: $(TOOLS_DIR)/makextfiles $(BUILD_DIR)/ext_files_sections_noscriptgeo.txt
>	$(V)$(TOOLS_DIR)/makextfiles $(BUILD_DIR)/ext_files_sections_noscriptgeo.txt $(BUILD_DIR)/ext_files_noscriptgeo.dat $(BUILD_DIR)/ext_files_defsym_noscriptgeo.txt.tmp
>	@rm -f $(BUILD_DIR)/ext_files_noscriptgeo.dat
>	@echo $(HARDCODED_SEGMENTS) >> $(BUILD_DIR)/ext_files_defsym_noscriptgeo.txt.tmp
>	@mv $(BUILD_DIR)/ext_files_defsym_noscriptgeo.txt.tmp $(BUILD_DIR)/ext_files_defsym_noscriptgeo.txt

$(BUILD_DIR)/ext_files_defsym_plusmenu.txt: $(TOOLS_DIR)/makextfiles $(BUILD_DIR)/ext_files_sections_plusmenu.txt
>	$(V)$(TOOLS_DIR)/makextfiles $(BUILD_DIR)/ext_files_sections_plusmenu.txt $(BUILD_DIR)/ext_files_plusmenu.dat $(BUILD_DIR)/ext_files_defsym_plusmenu.txt.tmp
>	@rm -f $(BUILD_DIR)/ext_files_plusmenu.dat
>	@echo $(HARDCODED_SEGMENTS) >> $(BUILD_DIR)/ext_files_defsym_plusmenu.txt.tmp
>	@mv $(BUILD_DIR)/ext_files_defsym_plusmenu.txt.tmp $(BUILD_DIR)/ext_files_defsym_plusmenu.txt

ifneq ($(BENCH),0)
$(BUILD_DIR)/ext_files_sections_tmp.txt: $(BUILD_DIR)/ext_files_sections_noscriptgeo.txt $(BUILD_DIR)/levels/bob.levelscriptgeosection.txt $(BUILD_DIR)/assets/mario_anim_data.marioanimsection.txt
>	@cat $^ > $@
else
$(BUILD_DIR)/ext_files_sections_tmp.txt: $(BUILD_DIR)/ext_files_sections_plusmenu.txt $(BUILD_DIR)/levels/intro.levelscriptgeosection.txt $(LEVEL_SEG_FILES:%/leveldata.elf=%.levelscriptgeosection.txt) $(BUILD_DIR)/assets/demo_data.section.txt $(BUILD_DIR)/assets/mario_anim_data.marioanimsection.txt
>	@cat $^ > $@
endif

$(BUILD_DIR)/ext_files_sections_noaudio.txt: $(BUILD_DIR)/ext_files_sections_tmp.txt $(BUILD_DIR)/tex_pack
>	@cp $< $@.tmp
>	@echo " $(BUILD_DIR)/tex_pack:0:$$(printf "%x" $$(wc -c <"$(BUILD_DIR)/tex_pack"))!_texture_data_segment:_texture_data_segment_end " >> $@.tmp
>	@mv $@.tmp $@

$(BUILD_DIR)/ext_files_sections.txt: $(BUILD_DIR)/ext_files_sections_noaudio.txt $(BUILD_DIR)/soundtable $(BUILD_DIR)/sounddata
>	@cp $< $@.tmp
>	@echo " $(BUILD_DIR)/soundtable:0:$$(printf "%x" $$(wc -c <"$(BUILD_DIR)/soundtable"))!_audio_table_segment:_audio_table_segment_end " >> $@.tmp
>	@echo " $(BUILD_DIR)/sounddata:0:$$(printf "%x" $$(wc -c <"$(BUILD_DIR)/sounddata"))!_audio_sample_segment:_audio_sample_segment_end " >> $@.tmp
>	@mv $@.tmp $@

ifneq ($(BENCH),0)
	EXT_FILES_SECTIONS_TXT := $(BUILD_DIR)/ext_files_sections_noaudio.txt
else
	EXT_FILES_SECTIONS_TXT := $(BUILD_DIR)/ext_files_sections.txt
endif

$(BUILD_DIR)/ext_files_defsym.txt $(BUILD_DIR)/ext_files.dat &: $(TOOLS_DIR)/makextfiles $(EXT_FILES_SECTIONS_TXT)
>	$(V)$(TOOLS_DIR)/makextfiles $(EXT_FILES_SECTIONS_TXT) $(BUILD_DIR)/ext_files.dat $(BUILD_DIR)/ext_files_defsym.txt.tmp
>	@echo $(HARDCODED_SEGMENTS) >> $(BUILD_DIR)/ext_files_defsym.txt.tmp
>	@mv $(BUILD_DIR)/ext_files_defsym.txt.tmp $(BUILD_DIR)/ext_files_defsym.txt

#==============================================================================#
# Generated Source Code Files                                                  #
#==============================================================================#

# Convert binary file to a comma-separated list of byte values for inclusion in C code
$(BUILD_DIR)/%.inc.c: $(BUILD_DIR)/%
>	$(call print,Converting to C:,$<,$@)
>	$(V)hexdump -v -e '1/1 "0x%X,"' $< > $@
>	$(V)echo >> $@

$(CONV_MODEL_C_FILES): $(BUILD_DIR)/%.processed.c: %.c $(TOOLS_DIR)/preprocess_graphics.py
>	$(call print,Preprocessing graphics:,$<,$@)
>	$(V)mkdir -p $(dir $@)
>	$(V)$(PYTHON) $(TOOLS_DIR)/preprocess_graphics.py $< $@ $(C_DEFINES)

$(BUILD_DIR)/%.processed.c: %.c $(CONV_MODEL_C_FILES) $(TOOLS_DIR)/preprocess_graphics.py
>	$(call print,Preprocessing graphics:,$<,$@)
>	$(V)mkdir -p $(dir $@)
>	$(V)$(PYTHON) $(TOOLS_DIR)/preprocess_graphics.py $< $@ $(C_DEFINES)

$(BUILD_DIR)/%.processed.c: $(BUILD_DIR)/%.c $(CONV_MODEL_C_FILES) $(TOOLS_DIR)/preprocess_graphics.py
>	$(call print,Preprocessing graphics:,$<,$@)
>	$(V)mkdir -p $(dir $@)
>	$(V)$(PYTHON) $(TOOLS_DIR)/preprocess_graphics.py $< $@ $(C_DEFINES)

# Generate animation data
$(BUILD_DIR)/assets/mario_anim_data.c: $(wildcard assets/anims/*.inc.c) $(TOOLS_DIR)/mario_anims_converter.py
>	@$(PRINT) "$(GREEN)Generating animation data $(NO_COL)\n"
>	$(V)$(PYTHON) $(TOOLS_DIR)/mario_anims_converter.py > $@.tmp
>	$(V)mv $@.tmp $@

# Generate demo input data
$(BUILD_DIR)/assets/demo_data.c: assets/demo_data.json $(wildcard assets/demos/*.bin)
>	@$(PRINT) "$(GREEN)Generating demo data $(NO_COL)\n"
>	$(V)$(PYTHON) $(TOOLS_DIR)/demo_data_converter.py assets/demo_data.json $(DEF_INC_CFLAGS) > $@

# Encode in-game text strings
$(BUILD_DIR)/include/text_strings.h: include/text_strings.h.in
>	$(call print,Encoding:,$<,$@)
>	$(V)$(TEXTCONV) charmap.txt $< $@
$(BUILD_DIR)/include/text_menu_strings.h: include/text_menu_strings.h.in
>	$(call print,Encoding:,$<,$@)
>	$(V)$(TEXTCONV) charmap_menu.txt $< $@
$(BUILD_DIR)/text/%/define_courses.inc.c: text/define_courses.inc.c text/%/courses.h
>	@$(PRINT) "$(GREEN)Preprocessing: $(BLUE)$@ $(NO_COL)\n"
>	$(V)$(CPP) $(CPPFLAGS) $< -o - -I text/$*/ | $(TEXTCONV) charmap.txt - $@
$(BUILD_DIR)/text/%/define_text.inc.c: text/define_text.inc.c text/%/courses.h text/%/dialogs.h
>	@$(PRINT) "$(GREEN)Preprocessing: $(BLUE)$@ $(NO_COL)\n"
>	$(V)$(CPP) $(CPPFLAGS) $< -o - -I text/$*/ | $(TEXTCONV) charmap.txt - $@

# Level headers
$(BUILD_DIR)/include/level_headers.h: levels/level_headers.h.in
>	$(call print,Preprocessing level headers:,$<,$@)
>	$(V)$(CPP) $(CPPFLAGS) -I . $< | sed -E 's|(.+)|#include "\1"|' > $@

# Rebuild files with 'GLOBAL_ASM' if the NON_MATCHING flag changes.
$(GLOBAL_ASM_O_FILES): $(GLOBAL_ASM_DEP).$(NON_MATCHING)
$(GLOBAL_ASM_DEP).$(NON_MATCHING):
>	@$(RM) $(GLOBAL_ASM_DEP).*
>	$(V)touch $@

#==============================================================================#
# Compilation Recipes                                                          #
#==============================================================================#

# Compile C/C++ code
$(BUILD_DIR)/%.o: %.cpp $(CFLAGS_FILE)
>	$(call print,Compiling:,$<,$@)
>	@$(CXX) -fsyntax-only $(CFLAGS) -MMD -MP -MT $@ -MF $(BUILD_DIR)/$*.d $<
>	$(V)$(CXX) -c $(CFLAGS) -o $@ $<

$(BUILD_DIR)/%.o: %.c $(CFLAGS_FILE) $(BUILD_DIR)/ext_files.dat
>	$(call print,Compiling:,$<,$@)
>	$(V)$(CC) -c $(CFLAGS) -iquote $(dir $@) -MMD -MP -MT $@ -MF $(BUILD_DIR)/$*.d -o $@ $<
$(BUILD_DIR)/%.o2: $(BUILD_DIR)/%.processed.c $(TOOLS_DIR)/preprocess_graphics.py $(BUILD_DIR)/tex_pack $(CAKE_INC_C) $(CFLAGS_FILE)
>	$(call print,Compiling:,$<,$@)
>	$(V)$(CC) -c $(EXT_CFLAGS) -iquote $(dir $*) -MMD -MP -MT $@ -MF $(BUILD_DIR)/$*.d -o $@ $<
$(BUILD_DIR)/%.o: $(BUILD_DIR)/%.c $(CFLAGS_FILE)
>	$(call print,Compiling:,$<,$@)
>	$(V)$(CC) -c $(CFLAGS) -MMD -MP -MT $@ -MF $(BUILD_DIR)/$*.d -o $@ $<
$(BUILD_DIR)/%.libc.o: %.c $(CFLAGS_FILE)
>	$(call print,Compiling:,$<,$@)
>	$(V)$(CC) -c $(CFLAGS) -MMD -MP -MT $@ -MF $(BUILD_DIR)/$*.d -o $@ $<

# Assemble assembly code
$(BUILD_DIR)/%.s.o: %.s $(CFLAGS_FILE)
>	$(call print,Assembling:,$<,$@)
>	$(V)$(CPP) $(CPPFLAGS) $< -o $@.s
>	$(V)$(AS) $(ASFLAGS) -MD $(BUILD_DIR)/$*.d $@.s -o $@
$(BUILD_DIR)/%.s.libc.o: %.s $(CFLAGS_FILE)
>	$(call print,Assembling:,$<,$@)
>	$(V)$(CPP) $(CPPFLAGS) $< -o $@.s
>	$(V)$(AS) $(ASFLAGS) -MD $(BUILD_DIR)/$*.d $@.s -o $@

$(BUILD_DIR)/executable.preprocessed.ld: ps1-bare-metal/executable.ld $(CFLAGS_FILE)
>	$(V)$(CPP) $(CPPFLAGS) $< > $@

ifneq ($(BENCH),0)
	LDFLAGS += -Wl,--unresolved-symbols=ignore-all
endif

$(ELF): $(SEG_FILES_WITH_GEO) $(O_FILES) $(GODDARD_O_FILES) $(LIBC_O_FILES) $(BUILD_DIR)/ext_files_defsym.txt $(BUILD_DIR)/executable.preprocessed.ld $(EXT_FILES_DAT_IF_NEEDED)
>	@$(PRINT) "$(GREEN)Linking elf: $(BLUE)$@ $(NO_COL)\n"
>	$(V)$(LD) $(CFLAGS) -L $(BUILD_DIR) -no-pie -o $@.tmp $(addprefix -Wl$(COMMA)--just-symbols=,$(SEG_FILES_WITH_GEO)) `cat $(BUILD_DIR)/ext_files_defsym.txt` $(O_FILES) $(GODDARD_O_FILES) $(LIBC_O_FILES) $(LDFLAGS)
>	$(V)$(OBJCOPY) -R.scratchpad -R.bss -R.sbss $@.tmp
>	@mv $@.tmp $@

$(EXE): $(ELF)
>	@$(PRINT) "$(GREEN)Converting to executable: $(BLUE)$@ $(NO_COL)\n"
>	$(V)$(PYTHON) ps1-bare-metal/convertExecutable.py $< $@

$(BUILD_DIR)/psx_iso.xml: psx_iso.xml
>	$(V)cp $< $@

$(ISO_OUT): $(EXE) $(BUILD_DIR)/bgm/pack.xa $(BUILD_DIR)/psx_iso.xml system.cnf
>	@$(PRINT) "$(GREEN)Making iso file: $(BLUE)$@ $(NO_COL)\n"
>	$(V)cd $(BUILD_DIR) && ../../$(MKPSXISO) -y ./psx_iso.xml

NOPS ?= nops

IP ?= 192.168.0.32:23

send: $(EXE)
>	$(V)nops /fast /exe $< $(IP)

serve: $(EXE) $(BUILD_DIR)/ext_files.dat
>	$(V)$(PYTHON) $(TOOLS_DIR)/serial_server.py $(BUILD_DIR)/ext_files.dat $(IP)

.PHONY: all clean distclean default send serve
# with no prerequisites, .SECONDARY causes no intermediate target to be removed
.SECONDARY:

-include $(DEP_FILES)

print-% : ; $(info $* is a $(flavor $*) variable set to [$($*)]) @true
