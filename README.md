# Super Mario 64 (DualShock Version)

- This is a fork of [the full decompilation of Super Mario 64 (J), (U), (E), and (SH)](https://github.com/n64decomp/sm64).
- It is heavily modified and can no longer target Nintendo 64, only PSX and PC (for debugging).
- There are still many limitations.
- For now, it can only build from the US version.

This repo does not include all assets necessary for compiling the game.
An original copy of the game is required to extract the assets.

## Features

- Cool "DUAL SHOCK™ Compatible" graphic mimicking the original "振動パック対応" (Rumble Pak Compatible) graphic
- An analog rumble signal is now produced for the DualShock's large motor, in addition to the original modulated digital signal for the small motor and for the SCPH-1150 Dual Analog Controller
- Low-precision soft float implementation specially written for PSX to reduce the performance impact of floats
- Large amounts of code have been adapted to use fixed point math, including the 16-bit integer vectors and matrices that are standard on PSX
- Simplified rewritten render graph walker
- Tessellation (up to 2x) to reduce issues with large polygons
- RSP display lists are compiled just-in-time into a custom display list format that is more compact and faster to process
- Display list preprocessor that removes commands we won't use and optimizes meshes (TODO: make it fix more things)
- Mario's animations are compressed (from 580632 to 190324 bytes) and placed in a corner of VRAM rather than being loaded from storage (we don't have the luxury of a cartridge to copy from in the middle of a frame)
- Custom profiler
- Custom texture encoder that quantizes all textures to 4 bits per pixel
- Translucent circle-texture shadows replaced with subtractive hexagonal shadows, as the PSX doesn't support arbitrary translucency
- (TODO) Camera system adapted to rotate with the right analog stick
- (TODO) Simplified rewritten Goddard subsystem

## Known issues

- Floating trees (temporary issue due caused by a math rewrite)
- Some of Mario's animations do not play, and may even crash the game
- Music cannot be generated at build time without manually obtaining the tracks
- Sound effects work but sometimes sound odd or are missing notes
- The camera cannot be controlled in many levels due to the unfinished camera control implementation
- Crashes when entering certain levels (due to insufficient memory?)
- Ending sequence crashes on load
- When reaching the bridge in the castle grounds, Mario looks up but Lakitu never comes over
- Poles do not go down when pounded
- Textures are loaded individually, causing long stutters and loading times
- Stretched textures due to PSX limitations (the graphics preprocessor could help)
- Tessellation is not good enough to fix all large polygons (the graphics preprocessor could help)
- Some textures are rendered incorrectly (RSP JIT issues?)
- Title screen is unfinished
- Pause menu doesn't work

## Building

### Linux

1. Build and install the mipsel-none-elf-gcc toolchain. For Arch users, it is available on [AUR](https://aur.archlinux.org/packages/mipsel-none-elf-gcc-git). (You can also install it on your system from https://github.com/malucard/poeng by running `make install-gcc` from there. This may take a long time.)
2. Clone the repo: `git clone https://github.com/sm64-port/sm64-port.git`, which will create a directory `sm64-port` and then **enter** it `cd sm64-port`.
3. Place a Super Mario 64 ROM called `baserom.<VERSION>.z64` into the repository's root directory for asset extraction, ~~where `VERSION` can be `us`, `jp`, or `eu`~~. (For now, only `us` is supported.)
4. (Optional) Create a folder named `.local` in the root of the repo and place every track of the soundtrack in it as a .wav file, numbered from 0 to 37 (0.wav, 1.wav, etc).
5. Run `make` to build. To build the benchmark version without music, run `make BENCH=1`.
The disc image will be located at `build/<VERSION>_psx/sm64.<VERSION>.iso`. The benchmark version will not generate an iso, only an elf and an exe, and it will require a PSX with 8MB of RAM (an emulator or a debug unit).

### Windows (untested)

1. Install and update MSYS2, following all the directions listed on https://www.msys2.org/.
2. From the start menu, launch MSYS2 MinGW and install required packages depending on your machine (do **NOT** launch "MSYS2 MSYS"):
  * 64-bit: Launch "MSYS2 MinGW 64-bit" and install: `pacman -S git make python3 mingw-w64-x86_64-gcc mingw-w64-x86_64-meson`
  * 32-bit (will also work on 64-bit machines): Launch "MSYS2 MinGW 32-bit" and install: `pacman -S git make python3 mingw-w64-i686-gcc mingw-w64-i686-meson`
  * Do **NOT** by mistake install the packages called simply `gcc` and `meson`.
3. Build and install the mipsel-none-elf-gcc toolchain.
4. The MSYS2 terminal has a _current working directory_ that initially is `C:\msys64\home\<username>` (home directory). At the prompt, you will see the current working directory in yellow. `~` is an alias for the home directory. You can change the current working directory to `My Documents` by entering `cd /c/Users/<username>/Documents`.
5. Clone the repo: `git clone https://github.com/malucard/sm64-psx.git`, which will create a directory `sm64-psx` and then **enter** it `cd sm64-psx`.
6. Place a *Super Mario 64* ROM called `baserom.<VERSION>.z64` into the repository's root directory for asset extraction, ~~where `VERSION` can be `us`, `jp`, or `eu`~~. (For now, only `us` is supported.)
7. (Optional) Create a folder named `.local` in the root of the repo and place every track of the soundtrack in it as a .wav file, numbered from 0 to 37 (0.wav, 1.wav, etc).
8. Run `make` to build. To build the benchmark version, run `make BENCH=1`.
The disc image will be located at `build/<VERSION>_psx/sm64.<VERSION>.iso`. The benchmark version will not generate an iso, only an elf and an exe, and it will require a PSX with 8MB of RAM (an emulator or a debug unit).

#### Troubleshooting

1. If you get `make: gcc: no suitable C and C++ compiler found`, `make: gcc: command not found`, `make: gcc: No such file or directory` although the packages did successfully install, you probably launched the wrong MSYS2. Read the instructions again. The terminal prompt should contain "MINGW32" or "MINGW64" in purple text, and **NOT** "MSYS".
2. If you get `Failed to open baserom.us.z64!` you failed to place the baserom in the repository. You can write `ls` to list the files in the current working directory. If you are in the `sm64-psx` directory, make sure you see it here.
3. If you get `make: *** No targets specified and no makefile found. Stop.`, you are not in the correct directory. Make sure the yellow text in the terminal ends with `sm64-psx`. Use `cd <dir>` to enter the correct directory. If you write `ls` you should see all the project files, including `Makefile` if everything is correct.
4. If you get any error, be sure MSYS2 packages are up to date by executing `pacman -Syu` and `pacman -Su`. If the MSYS2 window closes immediately after opening it, restart your computer.
5. Check if mipsel gcc is working by executing `mipsel-none-elf-gcc -v`. If it doesn't work, you either opened the wrong MSYS start menu entry or installed the incorrect gcc package.
6. When switching between building on other platforms, run `make -C tools clean` first to allow for the tools to recompile on the new platform. This also helps when switching between shells like WSL and MSYS2.

## Project Structure

	sm64
	├── actors: object behaviors, geo layout, and display lists
	├── assets: animation and demo data
	│   ├── anims: animation data
	│   └── demos: demo data
	├── bin: C files for ordering display lists and textures
	├── build: output directory
	├── data: behavior scripts, misc. data
	├── doxygen: documentation infrastructure
	├── enhancements: example source modifications
	├── include: header files
	├── levels: level scripts, geo layout, and display lists
	├── lib: N64 SDK code
	├── sound: sequences, sound samples, and sound banks
	├── src: C source code for game
	│   ├── audio: audio code
	│   ├── buffers: stacks, heaps, and task buffers
	│   ├── engine: script processing engines and utils
	│   ├── game: behaviors and rest of game source
	│   ├── goddard: rewritten Mario intro screen
	│   ├── goddard_og: backup of original Mario intro screen
	│   ├── menu: title screen and file, act, and debug level selection menus
	│   └── port: port code, audio and video renderer
	├── text: dialog, level names, act names
	├── textures: skybox and generic texture data
	└── tools: build tools

## Contributing

Pull requests are welcome. For major changes, please open an issue first to
discuss what you would like to change.
