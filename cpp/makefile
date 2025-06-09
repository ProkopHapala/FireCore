MAKEFLAGS += --no-print-directory
SHELL := /bin/bash

WASM_BUILD_DIR := Build-wasm
EM_CONFIG_FLAG := --em-config wasmlib/custom_cache/config


.PHONY: wasm
wasm: $(WASM_BUILD_DIR)/index.html

.PHONY: run
run: wasm
	@echo "Running local web server on port 8080, go to http://0.0.0.0:8080/MolGUIapp.html"
	cd $(WASM_BUILD_DIR)/apps/MolecularEditor && python3 -m http.server 8080



$(WASM_BUILD_DIR)/makefile: wasm-libs
	@mkdir -p $(WASM_BUILD_DIR)
	cd $(WASM_BUILD_DIR) && emcmake --em-config $(shell pwd)/wasmlib/custom_cache/config cmake .. -DWITH_SDL=ON -DWITH_FFTW=OFF -DWITH_OPENCL=OFF -DTARGET_WASM=ON -Wno-dev -DRELEASE=ON

$(WASM_BUILD_DIR)/index.html: $(WASM_BUILD_DIR)/makefile wasm-libs
	cd $(WASM_BUILD_DIR) && emmake --em-config $(shell pwd)/wasmlib/custom_cache/config make MolGUIapp -j 5



# ------------ extra libraries that have to be included or built for wasm (SIMDE, gl4es, fftw3) -------------
wasm-libs: wasmlib/custom_cache/config

wasmlib/custom_cache/config:
	@mkdir wasmlib/custom_cache/ -p
	cp /usr/share/emscripten/cache/ wasmlib/custom_cache/cache -r
	cp /usr/share/emscripten/.emscripten wasmlib/custom_cache/config
	echo -e "\n\nCACHE = '$(shell pwd)/wasmlib/custom_cache/cache'" >> wasmlib/custom_cache/config
	echo -e "\nFROZEN_CACHE = False" >> wasmlib/custom_cache/config
	embuilder build sdl2 --em-config wasmlib/custom_cache/config
	# use 'emcc --generate-config' instead?
