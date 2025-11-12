CFLAGS=-Wall -Wextra \
	   -framework CoreVideo -framework IOKit -framework Cocoa -framework GLUT -framework OpenGL \
	   -I./thirdparty -L./thirdparty ./thirdparty/libraylib.a \
	   -lm -D=PLATFORM_DESKTOP

WASM_FLAGS=--target=wasm32-wasi -nostdlib -Wl,--no-entry -Wl,--export-all

.DEFAULT_GOAL := main

main: main.c hydrogen.h
	$(CC) $(CFLAGS) -o main main.c

main-wasm: main.c hydrogen.h
	$(CC) $(CFLAGS) $(WASM_FLAGS) -o main main.c
