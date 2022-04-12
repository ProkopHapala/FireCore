#ifndef  clUtils_h
#define  clUtils_h

#include <vector>
#include <CL/cl.h>

#include "OCLerrors.h"
#include "OCL_device_picker.h"

#include "datatypes.h"

//void print( cl_mem cmem ){ printf( "\n", cmem.  ); };

class OCLBuffer{
    public:
    void       * p_cpu = NULL;
    cl_mem       p_gpu    = 0;
    size_t       n        = 0;
    size_t       typesize = 0;
    bool         read_on_finish = false;
    cl_mem_flags flags    = CL_MEM_READ_WRITE;
    char       * name     = NULL;
    // following is needed only for images
    int  img_dims = 0;
    int  nImg[3]  = {0,0,0};
    cl_image_format imageFormat;

    inline size_t byteSize(){ return typesize * n; }

    inline int initOnGPU ( cl_context& context ){
        int err;
        //printf( "initOnGPU() buff_size %li | n %li typesize %li \n", byteSize(), n, typesize );
        if( (flags&CL_MEM_COPY_HOST_PTR)||(flags&CL_MEM_USE_HOST_PTR) ){
                p_gpu = clCreateBuffer(context, flags, byteSize(), p_cpu,   &err);
        }else{  p_gpu = clCreateBuffer(context, flags, byteSize(), 0,       &err); }
        //printf( "initOnGPU p_gpu: %li \n", (long)p_gpu );
        return err;
    }

    inline int initOnGPUImage( cl_context& context ){
        int err;
        //p_gpu = clCreateBuffer(context, flags, typesize * n, NULL,  &err);
        switch(img_dims){
            case 2:
                printf( " initOnGPUImage: clCreateImage2D \n" );
                p_gpu = clCreateImage2D(context, flags, &imageFormat, nImg[0],nImg[1],          0,    p_cpu, &err);   // TODO: ??? nx=nImg[0] ny=nImg[1]  ???
                break;
            case 3:
                printf( " initOnGPUImage: clCreateImage3D \n" );
                p_gpu = clCreateImage3D(context, flags, &imageFormat, nImg[0],nImg[1], nImg[2], 0, 0, p_cpu, &err);   // TODO: ??? nx=nImg[0] ny=nImg[1]  ???
                break;
        }
        printf( "initOnGPUImage img_dims: %i p_gpu: %li \n", img_dims, (long)p_gpu );
        return err;
    }

    inline int setAsArg( cl_kernel& kernel, int i   ){ return clSetKernelArg(kernel, i, sizeof(cl_mem), &p_gpu );  };


    inline int fromGPU ( cl_command_queue& commands,       void* cpu_data ){ return clEnqueueReadBuffer ( commands, p_gpu, CL_TRUE, 0, byteSize(), cpu_data, 0, NULL, NULL );  }
    inline int toGPU   ( cl_command_queue& commands, const void* cpu_data ){ return clEnqueueWriteBuffer( commands, p_gpu, CL_TRUE, 0, byteSize(), cpu_data, 0, NULL, NULL ); }
    inline int fromGPU ( cl_command_queue& commands ){ return fromGPU ( commands, p_cpu ); }
    inline int toGPU   ( cl_command_queue& commands ){ return toGPU   ( commands, p_cpu ); }
    //inline setImageParams(  );

    inline OCLBuffer(){};
    inline OCLBuffer( char* name_, size_t n_, size_t typesize_, void * p_cpu_, cl_mem_flags flags_=CL_MEM_READ_WRITE ) :n(n_),typesize(typesize_),p_cpu(p_cpu_),flags(flags_),name(name_){};
};

class OCLsystem{
    // http://stackoverflow.com/questions/20105566/advantages-of-a-program-containing-several-opencl-kernels-versus-several-program
    public:
    cl_int           err;           // error code returned from OpenCL calls
    cl_device_id     device   = 0;        // compute device id
    cl_context       context  = 0;       // compute context
    cl_command_queue commands = 0;      // compute command queue
    cl_program       program  = 0;       // compute program - TODO FIXME: There could be more than one !!!!

    std::vector<cl_kernel> kernels;
    std::vector<OCLBuffer> buffers;

    void check_programSet (){ if(program ==0){ printf("ERROR OCLsystem program  not set \n"); exit(-1); } }
    void check_contextSet (){ if(context ==0){ printf("ERROR OCLsystem context  not set \n"); exit(-1); } }
    void check_deviceSet  (){ if(device  ==0){ printf("ERROR OCLsystem device   not set \n"); exit(-1); } }
    void check_commandsSet(){ if(commands==0){ printf("ERROR OCLsystem commands not set \n"); exit(-1); } }

    int init(){
        cl_uint deviceIndex = 0;
        //parseArguments(argc, argv, &deviceIndex);
        cl_device_id devices[MAX_DEVICES];
        unsigned numDevices = getDeviceList(devices);
        if (deviceIndex >= numDevices){  printf("Invalid device index (try '--list')\n"); return -1; }
        device = devices[deviceIndex];
        char name[MAX_INFO_STRING];
        getDeviceName(device, name);
        printf("\nUsing OpenCL device: %s\n", name);
        context  = clCreateContext(0, 1, &device, NULL, NULL, &err);  OCL_checkError(err, "Creating context");
        commands = clCreateCommandQueue(context, device, 0, &err);    OCL_checkError(err, "Creating command queue");
        return err;
    }

    int newKernel( char * name ){
        check_programSet();
        int err; kernels.push_back( clCreateKernel( program, name, &err ) );  OCL_checkError(err, "clCreateKernel"); return kernels.size()-1;
    }

    int newBuffer( char* name, size_t n, size_t typesize, void * p_cpu, cl_mem_flags flags=CL_MEM_READ_WRITE ){
        check_contextSet();
        buffers.push_back( OCLBuffer( name, n, typesize, p_cpu, flags ) ); int i=buffers.size()-1; int err=buffers[i].initOnGPU(context); OCL_checkError(err, "initOnGPU"); return i;
    }

    int newBufferImage2D( char* name, size_t nx, size_t ny, size_t typesize, void * p_cpu, cl_mem_flags flags, cl_image_format imageFormat ){
        check_contextSet();
        buffers.push_back( OCLBuffer( name, nx*ny, typesize, p_cpu, flags ) );
        int i=buffers.size()-1;
        buffers[i].img_dims    = 2;
        buffers[i].nImg[0]     = nx;
        buffers[i].nImg[1]     = ny;
        buffers[i].imageFormat = imageFormat;
        int err=buffers[i].initOnGPUImage(context); OCL_checkError(err, "initOnGPUImage");
        return i;
    }

    int newBufferImage3D( char* name, size_t nx, size_t ny, size_t nz, size_t typesize, void * p_cpu, cl_mem_flags flags, cl_image_format imageFormat ){
        check_contextSet();
        buffers.push_back( OCLBuffer( name, nx*ny, typesize, p_cpu, flags ) );
        int i=buffers.size()-1;
        buffers[i].img_dims    = 3;
        buffers[i].nImg[0]     = nx;
        buffers[i].nImg[1]     = ny;
        buffers[i].nImg[2]     = nz;
        buffers[i].imageFormat = imageFormat;
        printf( "newBufferImage3D buffers[%i].img_dims %i \n", i, buffers[i].img_dims  );
        int err=buffers[i].initOnGPUImage(context); OCL_checkError(err, "initOnGPUImage" );
        return i;
    }

    int initBuffers   (){ int err = CL_SUCCESS; for(size_t i=0; i<buffers.size(); i++){  err |= buffers[i].initOnGPU ( context );     }; return err; }
    //int releaseBuffers(){ for(int i=0; i<buffers; i++){ clReleaseMemObject(buffers[i].p_gpu); } }

    char * getKernelSource(char *filename){
        FILE *file = fopen(filename, "r");
        if (!file){ fprintf(stderr, "Error: Could not open kernel source file\n"); exit(-1); }
        fseek(file, 0, SEEK_END);
        int len = ftell(file) + 1;
        rewind(file);
        char *source = (char *)calloc(sizeof(char), len);
        if (!source){ fprintf(stderr, "Error: Could not allocate memory for source string\n"); exit(-1); }
        fread(source, sizeof(char), len, file);
        fclose(file);
        return source;
    }

    // TODO : newProgram instead ?
    int buildProgram( char * fname ){
        check_deviceSet();
        char * kernelsource = getKernelSource( fname );
        // Create the comput program from the source buffer
        program = clCreateProgramWithSource(context, 1, (const char **) & kernelsource, NULL, &err);
        char tmpstr[1024];
        sprintf(tmpstr,"Creating program with %s", fname);
        OCL_checkError(err, tmpstr);
        free(kernelsource);
        //errNum = clBuildProgram(program, numDevices,deviceIDs, "-I.", NULL, NULL);
        err =      clBuildProgram(program, 0,         NULL,      "-I. -cl-std=CL2.0", NULL, NULL);
        //err = clBuildProgram(program, 0, NULL, NULL, NULL, NULL);
        if (err != CL_SUCCESS){
            printf( " ERROR in clBuildProgram %s \n", fname);
            OCL_buildProgramFailure( program, device );
            return -1;
        }
        //delete [] kernelsource; // TODO ??????
        return err;
    }

    inline int upload  (int i, const void* cpu_data ){ return buffers[i].toGPU  (commands,cpu_data); };
    inline int download(int i,       void* cpu_data ){ return buffers[i].fromGPU(commands,cpu_data); };
    //inline int upload  (int i, void* p_cpu ){ buffers[i].p_cpu=p_cpu; return buffers[i].toGPU  (commands); };
    //inline int download(int i, void* p_cpu ){ buffers[i].p_cpu=p_cpu; return buffers[i].fromGPU(commands); };

    inline int upload  (int i){ return buffers[i].toGPU  (commands); };
    inline int download(int i){ return buffers[i].fromGPU(commands); };



    inline int copy    (int from, int to, int from0, int to0, int n){ return clEnqueueCopyBuffer(commands,buffers[from].p_gpu,buffers[to].p_gpu,from0,to0,n,0,NULL,NULL); };
    inline int copyBuff(int from, int to                           ){ int n=buffers[from].n; int n_=buffers[to].n; if(n_<n)n=n_; return clEnqueueCopyBuffer(commands,buffers[from].p_gpu,buffers[to].p_gpu,0,0,n,0,NULL,NULL); };

    int download(){
        int err = CL_SUCCESS;
        for(size_t i=0; i<buffers.size(); i++ ){
            if( buffers[i].read_on_finish ){
                //printf("finish : reading buff %i \n", i);
                err |= buffers[i].fromGPU( commands );
            }
        }
        return err;
    }

    int finishRaw(){
        int err = clFinish(commands);   OCL_checkError(err, "Waiting for kernel to finish"); return err;
    }

    int finish(){
        int err;
        err = clFinish(commands);   OCL_checkError(err, "Waiting for kernel to finish");
        err |= download();
        return err;
    }

    void release_OCL(){
        clReleaseProgram(program);
        for(size_t i=0; i<kernels.size(); i++){ clReleaseKernel(kernels[i]);       }
        for(size_t i=0; i<buffers.size(); i++){ clReleaseMemObject(buffers[i].p_gpu); }
        //clReleaseKernel(kernel);
        clReleaseCommandQueue(commands);
        clReleaseContext(context);
    }
    ~OCLsystem(){ 
        release_OCL();
     }
};

#define OCL_BUFF   1
#define OCL_INT    2
#define OCL_FLOAT  3
#define OCL_LBUFF  4
#define OCL_PTR    5

#define BUFFarg(X)  OCLarg( (int)X, OCL_BUFF  )
#define INTarg(X)   OCLarg( (int)X, OCL_INT   )
#define FLOATarg(X) OCLarg( (float)X     )
#define LBUFFarg(X) OCLarg( (int)X, OCL_LBUFF )
#define PTRarg(X,n) OCLarg( (void*)X, OCL_PTR, n )
#define REFarg(X)   OCLarg( (void*)&X, sizeof(X) )

class OCLarg{
    public:
    int  kind=0;
    int  nbytes=0;
    union{
        float  f;
        int    i;
        void*  ptr;
    };
    inline void setFloat(  float f_ ){ f=f_; kind=OCL_FLOAT; }
    inline void setInt  (  int   i_ ){ i=i_; kind=OCL_INT  ; }
    inline void setBuff (  int   i_ ){ i=i_; kind=OCL_BUFF;  }
    OCLarg(){};
    OCLarg( float f_                ):f(f_), kind(OCL_FLOAT) {}//  printf( "!!!!! OCLarg(f %f kind %i )\n", f, kind ); }
    OCLarg( int   i_,   int kind_   ):i(i_), kind(kind_)     {}//  printf( "!!!!! OCLarg(i %i kind %i )\n", i, kind ); }
    OCLarg( void* ptr_, int nbytes_ ):ptr(ptr_), kind(OCL_PTR), nbytes(nbytes_){} //{  printf( "!!!!! OCLarg(ptr %li kind %i  nbytes %i )\n", (long)ptr, kind, nbytes ); }
};

class OCLtask{
    public:
    OCLsystem  * cl;
    //cl_kernel  * kernel;
    size_t      ikernel   = 0;
    size_t      dim       = 1;
    size_t      global[3] = {0,0,0};
    size_t      local [3] = {0,0,0};

    std::vector<OCLarg> args;

    int useArgs(){
        int err = CL_SUCCESS;
        cl_kernel kernel = cl->kernels[ikernel];
        for(size_t i=0; i<args.size(); i++){
            OCLarg& arg = args[i];
            //printf( "useArgs args[%i].kind: %i\n", i, arg.kind );
            switch(arg.kind){
                case OCL_BUFF:
                    //printf( "buffArg args[%i] ibuff %i p_gpu %i \n", i, arg.i, cl->buffers[arg.i].p_gpu );
                    err |= clSetKernelArg( kernel, i, sizeof(cl_mem), &(cl->buffers[arg.i].p_gpu) );             OCL_checkError(err, "setAsArg"); break;
                //case OCL_BUFF:  err |= cl->buffers[arg.i].setAsArg( kernel, i );                               OCL_checkError(err, "setAsArg"); break;
                case OCL_INT:   err |= clSetKernelArg( kernel, i, sizeof(int)  , &(arg.i) );                     OCL_checkError(err, "setAsArg"); break;
                case OCL_FLOAT: err |= clSetKernelArg( kernel, i, sizeof(float), &(arg.f) );                     OCL_checkError(err, "setAsArg"); break;
                case OCL_LBUFF: err |= clSetKernelArg( kernel, i, arg.i,          NULL    );                     OCL_checkError(err, "setAsArg"); break;
                case OCL_PTR:   err |= clSetKernelArg( kernel, i, arg.nbytes,     arg.ptr );                     OCL_checkError(err, "setAsArg"); break;
            }
        }
        return err;
    }

    inline int enque_raw(  ){
        //printf("enque_raw %i %i (%i,%i,%i) (%i,%i,%i)\n", ikernel, dim, global[0],global[1],global[2], local[0],local[1],local[2] );
        if(local[0]==0){ return clEnqueueNDRangeKernel( cl->commands, cl->kernels[ikernel], dim, NULL, global, NULL,  0, NULL, NULL );   }
        else{            return clEnqueueNDRangeKernel( cl->commands, cl->kernels[ikernel], dim, NULL, global, local, 0, NULL, NULL );   }
    }

    virtual int enque( ){
        int err;
        if( args.size() > 0 ) useArgs();
        err = enque_raw( );  OCL_checkError(err, "enque_raw");
        return err;
    }

    void print_arg_list(){
        printf("kernel[narg=%li]( ", args.size() );
        for(size_t i=0; i<args.size(); i++){
            switch(args[i].kind){
                case OCL_INT:   printf( "[%li]int %i, ",     i, args[i].i ); break;
                case OCL_FLOAT: printf( "[%li]float %g, ",   i, args[i].f ); break;
                case OCL_BUFF:  printf( "[%li]buff[%i]:%s, ",i, args[i].i, cl->buffers[args[i].i].name );   break;
                default:        printf( "[%li]arg(type=%i,val=%i) ",  i, args[i].kind, args[i].i );   break;
            }
        }
        printf(")\n");
    }

    inline void setup( OCLsystem  * cl_, size_t ikernel_, size_t dim_, size_t global_, size_t local_ ){ cl=cl_; ikernel=ikernel_; dim=dim_;  global[0]=global_;global[1]=global_;global[2]=global_; local[0]=local_;local[1]=local_;local[2]=local_; };
    OCLtask          ( OCLsystem  * cl_, size_t ikernel_, size_t dim_, size_t global_, size_t local_ ): cl(cl_),ikernel(ikernel_),dim(dim_){ global[0]=global_;global[1]=global_;global[2]=global_; local[0]=local_;local[1]=local_;local[2]=local_; };
    OCLtask(){};
};

// ========== Helper functions for coverting buffers to OpenCL format




#endif
