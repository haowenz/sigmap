cpp_source=sequence_batch.cc signal_batch.cc pore_model.cc cwt.cc spatial_index.cc sigmap.cc
src_dir=src
objs_dir=objs
objs+=$(patsubst %.cc,$(objs_dir)/%.o,$(cpp_source))

project_dir = $(shell pwd)
HDF5_DIR ?= ${project_dir}/extern/hdf5/build
HDF5_INCLUDE_DIR ?= ${HDF5_DIR}/include
HDF5_LIB_DIR ?= ${HDF5_DIR}/lib
HDF5_LIB ?= hdf5

cxx=g++
cxxflags=-std=c++11 -Wall -O3 -fopenmp -march=native -I${HDF5_INCLUDE_DIR}
ldflags=-L${HDF5_LIB_DIR} -Wl,-rpath=${HDF5_LIB_DIR} -l${HDF5_LIB} -lm -lz

exec=sigmap

all: hdf5 check_hdf5 dir $(exec) 
Sigmap: check_hdf5 dir $(exec) 

check_hdf5:
	@[ -f "${HDF5_INCLUDE_DIR}/H5pubconf.h" ] || { echo "HDF5 headers not found" >&2; exit 1; }
	@[ -f "${HDF5_LIB_DIR}/lib${HDF5_LIB}.so" ] || [ -f "${HDF5_LIB_DIR}/lib${HDF5_LIB}.a" ] || { echo "HDF5 library not found" >&2; exit 1; }

dir:
	mkdir -p $(objs_dir)

hdf5:
	cd extern/hdf5;\
  mkdir build;\
	./configure --enable-threadsafe --disable-hl --prefix="${HDF5_DIR}";\
	make -j;\
	make install

$(exec): $(objs)
	$(cxx) $(cxxflags) $(objs) -o $(exec) $(ldflags)
	
$(objs_dir)/%.o: $(src_dir)/%.cc
	$(cxx) $(cxxflags) -c $< -o $@ $(ldflags)

.PHONY: clean
clean:
	-rm -r $(exec) $(objs_dir)
