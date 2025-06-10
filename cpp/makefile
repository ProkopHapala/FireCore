MAKEFLAGS += --no-print-directory
SHELL := /bin/bash

WASM_BUILD_DIR := Build-wasm


.PHONY: wasm
wasm: $(WASM_BUILD_DIR)/index.html

.PHONY: run
run: wasm
	@echo "Running local web server on port 8080, go to http://0.0.0.0:8080/MolGUIapp.html"
	cd $(WASM_BUILD_DIR)/apps/MolecularEditor && python3 -m http.server 8080



$(WASM_BUILD_DIR)/makefile:
	@mkdir -p $(WASM_BUILD_DIR)
	cd $(WASM_BUILD_DIR) && emcmake cmake .. -DWITH_SDL=ON -DWITH_FFTW=OFF -DWITH_OPENCL=OFF -DTARGET_WASM=ON -Wno-dev -DRELEASE=ON

$(WASM_BUILD_DIR)/index.html: $(WASM_BUILD_DIR)/makefile
	cd $(WASM_BUILD_DIR) && emmake make MolGUIapp -j 5
