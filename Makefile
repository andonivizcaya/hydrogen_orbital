SRC := src
S := $(SRC)/main.c $(SRC)/hydrogen.h

CFLAGS = -Wall -Wextra -framework CoreVideo -framework IOKit -framework Cocoa -framework GLUT \
	-framework OpenGL -D=PLATFORM_DESKTOP -I./thirdparty/include -I$(SRC)
LDFLAGS = -L./thirdparty/libs ./thirdparty/libs/libraylib.a
LDLIBS = -lm

WASM_CC ?= clang
WASM_OPT ?= -Os
WASM_RES := $(shell $(WASM_CC) --target=wasm32 -nostdlib -print-resource-dir 2>/dev/null)
WASM_INCS = -nostdlibinc -isystem "$(WASM_RES)/include"
WASM_LINKER ?=
WASM_LD := $(strip $(WASM_LINKER))

ifeq ($(WASM_LD),)
  WASM_LD := $(shell command -v wasm-ld 2>/dev/null)
endif
ifeq ($(WASM_LD),)
  WASM_LD := $(wildcard /opt/homebrew/opt/lld/bin/wasm-ld)
endif
ifeq ($(WASM_LD),)
  WASM_LD := $(wildcard /usr/local/opt/lld/bin/wasm-ld)
endif
WASM_LINK = $(if $(WASM_LD),-fuse-ld=$(WASM_LD)) -Wl,--export-table -Wl,--allow-undefined -Wl,--export=main -Wl,--export=wasm_malloc -Wl,--no-entry
WASM_RT ?=

.DEFAULT_GOAL := main

.PHONY: hydrolib-wasm wasm-resources-copy

main: $(S) Makefile
	$(CC) $(CFLAGS) $(LDFLAGS) $(LDLIBS) -o main $(SRC)/main.c

hydrolib-wasm: wasm/app.wasm

wasm-resources-copy:
	@mkdir -p wasm/resources/fonts
	@if [ -f resources/fonts/Iosevka-Regular.ttc ]; then cp resources/fonts/Iosevka-Regular.ttc wasm/resources/fonts/; fi

wasm/app.wasm: $(S) Makefile wasm/hydrolib.js wasm-resources-copy
	@if [ -z "$(WASM_LD)" ]; then \
		echo >&2 "hydrolib-wasm: wasm-ld not found (install lld / put wasm-ld on PATH)."; exit 1; fi
	$(WASM_CC) -Wall -Wextra -fno-builtin -iquote "$(SRC)/wasm_libc" -I "$(SRC)/wasm_libc" -I "$(SRC)" \
		$(WASM_INCS) --target=wasm32 -nostdlib $(WASM_OPT) -I./thirdparty/include \
		-DPLATFORM_WEB -DHYDROLIB_WEB $(SRC)/main.c $(WASM_RT) $(WASM_LINK) -o $@

clean:
	rm -f main wasm/app.wasm