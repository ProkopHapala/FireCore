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
	cd $(WASM_BUILD_DIR) && emcmake cmake .. -DWITH_SDL=ON -DWITH_FFTW=ON -DWITH_OPENCL=OFF -DWITH_WASM=ON -Wno-dev

.PHONY: $(WASM_BUILD_DIR)/vfilesystem.js
$(WASM_BUILD_DIR)/vfilesystem.js:
	cp common_resources data -r
	/usr/share/emscripten/tools/file_packager TARGET --embed common_resources --embed data --js-output=$(WASM_BUILD_DIR)/vfilesystem.js
	rm data -r
	
$(WASM_BUILD_DIR)/index.html: $(WASM_BUILD_DIR)/makefile $(WASM_BUILD_DIR)/vfilesystem.js
	make wasm-libs
	cd $(WASM_BUILD_DIR) && emmake make MolGUIapp || echo " ----------- compilation failed, but have can use an ugly hack ---------- "

	echo ' -s FULL_ES2=1 $(shell pwd)/wasmlib/fftw3/Install/lib/libfftw3.a $(shell pwd)/wasmlib/gl4es/lib/libGL.a -sUSE_SDL=2 --pre-js $(shell pwd)/$(WASM_BUILD_DIR)/vfilesystem.js  -sFORCE_FILESYSTEM -s ALLOW_MEMORY_GROWTH=1 ' > $(WASM_BUILD_DIR)/apps/MolecularEditor/CMakeFiles/MolGUIapp.dir/linkLibs.rsp
	cd $(WASM_BUILD_DIR) && emmake make MolGUIapp	



# ------------ extra libraries that have to be included or built for wasm (SIMDE, gl4es, fftw3) -------------
wasm-libs: wasmlib/simde/simde/x86/fma.h wasmlib/fftw3/api/fftw3.h wasmlib/custom_cache/config wasmlib/gl4es/include/gl4esinit.h

wasmlib/simde/simde/x86/fma.h:
	@mkdir wasmlib -p
	cd wasmlib && git clone https://github.com/simd-everywhere/simde

wasmlib/fftw3/api/fftw3.h:
	rm -rf wasmlib/fftw3
	cd wasmlib && curl "http://www.fftw.org/fftw-3.3.10.tar.gz" -o fftw3.tar.gz
	cd wasmlib && tar -xf fftw3.tar.gz
	cd wasmlib && mv fftw-3.3.10/ fftw3/
	rm wasmlib/fftw3.tar.gz
	mkdir wasmlib/fftw3/Install
	cd wasmlib/fftw3 && emconfigure ./configure --prefix=$(shell pwd)/wasmlib/fftw3/Install --host=none --enable-silent-rules --disable-doc --disable-dependency-tracking --disable-fortran
	cd wasmlib/fftw3 && emmake make --no-print-directory -j
	cd wasmlib/fftw3 && emmake make install --no-print-directory

wasmlib/custom_cache/config:
	mkdir wasmlib/custom_cache/ -p
	cp /usr/share/emscripten/cache/ wasmlib/custom_cache/cache -r
	cp /usr/share/emscripten/.emscripten wasmlib/custom_cache/config
	echo "CACHE = '$(shell pwd)/wasmlib/custom_cache/cache'" >> wasmlib/custom_cache/config
	echo "FROZEN_CACHE = False" >> wasmlib/custom_cache/config
	# use 'emcc --generate-config' instead?

wasmlib/gl4es/include/gl4esinit.h:
	cd wasmlib && git clone https://github.com/ptitSeb/gl4es
	mkdir wasmlib/gl4es/build
	cd wasmlib/gl4es/build && emcmake cmake .. -DCMAKE_BUILD_TYPE=RelWithDebInfo -DNOX11=ON -DNOEGL=ON -DSTATICLIB=ON && make
