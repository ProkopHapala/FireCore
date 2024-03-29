'''
# === OLD BAK
PATH_MKL = /home/prokop/SW/intel/mkl
PATH_MPI = /usr/lib/x86_64-linux-gnu/openmpi
LPATHS  =  -L$(PATH_MKL)/lib/intel64 -I$(PATH_MKL)/include/fftw -L$(PATH_MPI)/lib -I$(PATH_MPI)include 
LFLAGS_ =  -lmkl_scalapack_lp64 -lmkl_blacs_openmpi_lp64 -lmkl_lapack95_lp64 -lmkl_intel_lp64 -lmkl_intel_thread -lmkl_core -liomp5 -lpthread -lfftw3xf_gnu -lm -Bdynamic -lmpi 
LFLAGS  = $(LPATHS) $(LFLAGS_) 
#LFLAGS =  -L/home/prokop/SW/intel/mkl/lib/intel64 -lmkl_scalapack_lp64 -lmkl_blacs_openmpi_lp64 -lmkl_lapack95_lp64 -lmkl_intel_lp64 -lmkl_intel_thread -lmkl_core -liomp5 -lpthread -L/home/prokop/intel/mkl/intel64 -I/home/prokop/intel/mkl/include/fftw -lfftw3xf_gnu -lm -Bdynamic -L/usr/lib/x86_64-linux-gnu/openmpi/lib -I/usr/lib/x86_64-linux-gnu/openmpi/include -lmpi 
INCLUDES = -I/usr/include/mpich/

# === NEW SIMPLYFIED
LFLAGS   = -lmkl_intel_lp64 -lmkl_intel_thread -lmkl_core -liomp5 -Bdynamic
INCLUDES = 
'''

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
    #print( path_dict )

    LFLAGS   = " -lmkl_scalapack_lp64 -lmkl_blacs_openmpi_lp64 -lmkl_lapack95_lp64 -lmkl_intel_lp64 -lmkl_intel_thread -lmkl_core -liomp5 -lpthread -lfftw3xf_gnu -lm -Bdynamic -lmpi "
    LPATHS   = " -L{PATH_MKL} -L{PATH_FFTW_LIB} -I{PATH_FFTW_INCL} -L{PATH_MPI_LIB} -I{PATH_MPI_INC} ".format(**path_dict)
    return LFLAGS, LPATHS

MODLAP95   = "-I/usr/local/lib/lapack95_modules"
BLAS       = "-lblas"
MODMPI     = "-I/usr/local/mpich/include"
VISFLAGS   = ""



'''
FFLAG_ALL = " -fPIC "
_FFLAGS = {
    'F90'      : FFLAG_ALL+" -freal-4-real-8 -ffree-form -ffree-line-length-none ",
    'F77'      : FFLAG_ALL+" -freal-4-real-8 ",
    'CC'       : FFLAG_ALL+" -freal-4-real-8 -ffree-form -ffree-line-length-none ",
    'OPT'      : FFLAG_ALL+" -O3 -mtune=native -ftree-vectorize -pg",
    'VERY_OPT' : FFLAG_ALL+" -Ofast -march=native -mtune=native",
    'DEBUGbak' : FFLAG_ALL+" -Og -g -fbounds-check -Wall -Wno-tabs",
    'DEBUG'    : FFLAG_ALL+" -Og -g -fcheck=all -fbounds-check -Wall -Wno-tabs -Wno-unused-variable -Wno-unused-label -Wno-missing-include-dirs",
    'DEBUGw'   : FFLAG_ALL+" -Og -g -fcheck=all -fbounds-check ",
}
'''

FFLAG_ALL = " -fPIC -freal-4-real-8 "
_FFLAGS = {
    'F77'      : FFLAG_ALL+" ",
    'F90'      : FFLAG_ALL+" -ffree-form -ffree-line-length-none ",
    'CC'       : FFLAG_ALL+" -ffree-form -ffree-line-length-none ",
    'OPT'      : " -O3 -mtune=native -ftree-vectorize -pg",
    'VERY_OPT' : " -Ofast -march=native -mtune=native",
    'DEBUGbak' : " -Og -g -fbounds-check -Wall -Wno-tabs",
    'DEBUG'    : " -Og -g -fcheck=all -fbounds-check -Wall -Wno-tabs -Wno-unused-variable -Wno-unused-label -Wno-missing-include-dirs",
    'DEBUGw'   : " -Og -g -fcheck=all -fbounds-check ",
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