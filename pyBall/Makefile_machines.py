
MKL_PATH = "/home/prokop/SW/intel"
MPI_PATH = "/usr/lib/x86_64-linux-gnu/openmpi"

def genLFLAGS( MKL_PATH=MKL_PATH, MPI_PATH=MPI_PATH ):
    path_dict = {
        "PATH_MKL"       : MKL_PATH+"/mkl/lib/intel64",
        "PATH_MKL_INCL"  : MKL_PATH+"/mkl/include",
        "PATH_FFTW_LIB"  : MKL_PATH+"/mkl/lib/intel64",
        "PATH_FFTW_INCL" : MKL_PATH+"/mkl/include/fftw",
        "PATH_MPI_LIB"   : MPI_PATH+"/lib",
        "PATH_MPI_INC"   : MPI_PATH+"/include",
    }
    #print path_dict
     
    LFLAGS   = " -lmkl_scalapack_lp64 -lmkl_blacs_openmpi_lp64 -lmkl_lapack95_lp64 -lmkl_intel_lp64 -lmkl_intel_thread -lmkl_core -liomp5 -lpthread -lfftw3xf_gnu -lm -Bdynamic -lmpi "
    #LPATHS  = "-L$(PATH_MKL) -L$(PATH_FFTW_LIB) -I$(PATH_FFTW_INCL) -L$(PATH_MPI_LIB) -I$(PATH_MPI_INC)"
    LPATHS   = " -L{PATH_MKL} -L{PATH_FFTW_LIB} -I{PATH_FFTW_INCL} -L{PATH_MPI_LIB} -I{PATH_MPI_INC} ".format(**path_dict)
    #LFLAGS   = LPATHS + LFLAGS
    return LFLAGS, LPATHS

MODLAP95   = "-I/usr/local/lib/lapack95_modules"
BLAS       = "-lblas"
MODMPI     = "-I/usr/local/mpich/include"
VISFLAGS   = ""

_FFLAGS = {
    'F90'      : " -freal-4-real-8 -ffree-form -ffree-line-length-none ",
    'F77'      : " -freal-4-real-8 ",
    'CC'       : " -freal-4-real-8 -ffree-form -ffree-line-length-none ",
    'OPT'      : " -O3 -mtune=native -ftree-vectorize ",
    'VERY_OPT' : " -Ofast -march=native -mtune=native ",
    'DEBUG'    : " -Og -g -fbounds-check -Wall -Wno-tabs",
    'DEBUGw'   : " -Og -g -fbounds-check ",
}

_LFLAGS = {
'LAPACK95' : '''-llapack -llapack95 -lg2c -L/usr/local/lib -L/usr/lib/gcc-lib/i386-redhat-linux/3.2.2''',
'USEBLAS'  : "-L${PATH_MKL} -lmkl_scalapack_lp64 -lmkl_blacs_openmpi_lp64 -lmkl_lapack95_lp64 -lmkl_intel_lp64    \
              -lmkl_intel_thread -lmkl_core -liomp5 -lpthread -L${PATH_FFTW_LIB} -I${PATH_FFTW_INCL}              \
              -lfftw3xf_gnu -lm -Bdynamic -L${PATH_MPI_LIB} -I${PATH_MPI_INC} -lmpi",
}

_PARLFLAGS = {
'SCALAPACK' : " -L/opt/intel/mkl72cluster/lib/em64t -lmkl_scalapack -lmkl_blacsF77init_intelmpi -lmkl_blacs_intelmpi \
                -lmkl_blacs_intelmpi -lmkl_blacs_intelmpi -lmkl_lapack -lmkl_def -lguide -lpthread -lg2c",
}