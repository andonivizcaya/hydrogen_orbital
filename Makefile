CFLAGS=-Wall -Wextra \
	   -framework CoreVideo -framework IOKit -framework Cocoa -framework GLUT -framework OpenGL -D=PLATFORM_DESKTOP \
	   -I./thirdparty/include

LDFLAGS=-L./thirdparty/libs ./thirdparty/libs/libraylib.a
LDLIBS=-lm

WASM_FLAGS=--target=wasm32-wasi -nostdlib -Wl,--no-entry -Wl,--export-all

.DEFAULT_GOAL := main

target-wasm: main.c hydrogen.h
	emcc -O3 -o ./wasm/game.html main.c -Os -Wall ./thirdparty/raylib-5.5_webassembly/lib/libraylib.a -I./thirdparty/include  -L./thirdparty/raylib-5.5_webassembly/lib -s USE_GLFW=3 -DPLATFORM_WEB -s ALLOW_MEMORY_GROWTH=1 -DGRAPHICS_API_OPENGL_ES2 -s FORCE_FILESYSTEM=1 --preload-file resources/fonts/Iosevka-Regular.ttc --shell-file ./wasm/shell.html

main: main.c hydrogen.h
	$(CC) $(CFLAGS) $(LDFLAGS) $(LDLIBS) -o main main.c

