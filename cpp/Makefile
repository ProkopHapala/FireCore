SHELL := /bin/bash

wasm-libs: wasmlib/simde/simde/x86/fma.h wasmlib/fftw3/api/fftw3.h wasmlib/custom_cache/config


wasmlib/simde/simde/x86/fma.h:
	@mkdir wasmlib
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
