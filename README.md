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
- Mario's animations are compressed (from 580632 to 190324 bytes) and placed in a corner of VRAM rather than being loaded from storage (we don't have the luxury of a fast cartridge to read from in the middle of a frame)
- Custom profiler
- Custom texture encoder that quantizes all textures to 4 bits per pixel
- Translucent texture shadows replaced with subtractive hexagonal shadows, as the PSX doesn't support arbitrary translucency
- (WIP) Camera system adapted to rotate with the right analog stick
- (WIP) Simplified rewritten Goddard subsystem

## Known issues

- Some of Mario's animations do not play, and may even crash the game
- Music cannot be generated at build time without manually obtaining the tracks
- Sound effects work but sometimes sound odd or are missing notes
- The camera cannot be controlled in many levels due to the unfinished camera control implementation
- Crashes when entering certain levels (due to insufficient memory?)
- Ending sequence crashes on load
- When reaching the bridge in the castle grounds, Mario looks up but Lakitu never comes over (you need to do the Lakitu skip)
- Poles do not go down when pounded
- Textures are loaded individually, causing long stutters and loading times
- Stretched textures due to PSX limitations (the graphics preprocessor could help)
- Tessellation is not good enough to fix all large polygons (the graphics preprocessor could help)
- Some textures are rendered incorrectly (RSP JIT issues?)
- Title screen is unfinished
- Pause menu doesn't work

## Building

1. Place a *Super Mario 64* ROM called `baserom.<VERSION>.z64` into the repository's root directory for asset extraction, ~~where `VERSION` can be `us`, `jp`, or `eu`~~. (For now, only `us` is supported.)
2. (Optional) Create a folder named `.local` in the root of the repo and place every track of the soundtrack in it as a .wav file, numbered from 0 to 37 (0.wav, 1.wav, etc). Don't worry too much about this, one day this will be done automatically.
3. Now you need a properly set up environment. A Dockerfile is provided to simplify this. To use it, ensure you have [Docker](https://www.docker.com) (or a compatible alternative) installed and set up. (if on Linux, put it in [rootless mode](https://docs.docker.com/engine/security/rootless/) as well!)
	- If on Linux, I've included a convenience script for invoking a container from CLI. Run `./idc` to enter a Bash shell in the container, or `./idc <command>` to run one command in it (for example, `./idc make` or `./idc make clean`).
	- If using Visual Studio Code, you can simply install the Dev Containers extension, open the repository, and click "Reopen in Container" (either from the notification or from the Ctrl+Shift+P menu). The editor will act like it's running inside the container, including the terminal.
	- Alternatively, if you plan to do a lot of PS1 development, you can skip the container and simply have the right things installed on your system, but this is the harder option. You must be using Linux (any). You need FFMPEG's libraries, libpng, xxd, Python 3, meson, GCC or Clang, and version 15 or later of the mipsel-none-elf-gcc toolchain. To build and install mipsel-none-elf-gcc, there is a utility in [this other repo](https://github.com/malucard/poeng). Clone it and run `make install-gcc`. It will take a pretty long time.
4. Now run `make`, and when it's done, sm64.iso and sm64.cue will be in `build/us_psx/`. To build in benchmark mode, use `make BENCH=1`. The benchmark mode boots directly into a level and doesn't require a CD, but requires 8 MB of RAM and won't work on a retail console. `make clean` will remove build/, and `make distclean` will clean both build/ and tools/.

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
