
# ======= Define your library paths here
PATH_MKL = /home/prokop/SW/intel/mkl
PATH_MPI = /usr/lib/x86_64-linux-gnu/openmpi
INCLUDES = -I/usr/include/mpich/
# ======= Library FLAGS
LPATHS  = -L$(PATH_MKL)/lib/intel64 -I$(PATH_MKL)/include/fftw -L$(PATH_MPI)/lib -I$(PATH_MPI)include 
LFLAGS_ = -lmkl_scalapack_lp64 -lmkl_blacs_openmpi_lp64 -lmkl_lapack95_lp64 -lmkl_intel_lp64 -lmkl_intel_thread -lmkl_core -liomp5 -lpthread -lfftw3xf_gnu -lm -Bdynamic -lmpi 
#LFLAGS = -L/home/prokop/SW/intel/mkl/lib/intel64 -lmkl_scalapack_lp64 -lmkl_blacs_openmpi_lp64 -lmkl_lapack95_lp64 -lmkl_intel_lp64 -lmkl_intel_thread -lmkl_core -liomp5 -lpthread -L/home/prokop/intel/mkl/intel64 -I/home/prokop/intel/mkl/include/fftw -lfftw3xf_gnu -lm -Bdynamic -L/usr/lib/x86_64-linux-gnu/openmpi/lib -I/usr/lib/x86_64-linux-gnu/openmpi/include -lmpi 
LFLAGS  = $(LPATHS) $(LFLAGS_) 

