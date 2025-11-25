# hydrogen_orbital

For now, it uses the time independent schrödinger equation, but in the future it should have the time dependent schrödinger one, both dirac equations and the klein-gordon.

## prerequisites
- make
- some C compiler
- something graphics related (opengl)
- emscripten (initialized on this shell session if you want to build for wasm)

## how to build
- for your pc, just make:
```bash
make -B
```
- for web:
```bash
make -B target-wasm
```

## next steps
- use nob.h or mate.h to remove the build process layer of ass
- port a subscript of raylib fully to js to remove the emscripten bloat
- remove the depencies on annoying standard libraries to make it more platform agnostic

## how to support
- don't be annoying
- please, do not scream you favorite library or cpp/rust take on this project as a pr or issue, jst lest me be happy writing dumb C code
