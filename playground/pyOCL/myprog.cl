
__kernel void mul(
    const int N,
    __global float* A,
    __global float* B,
    __global float* out
){
    const size_t i = get_global_id(0);
    if(i<N){ 
        out[i] = A[i] * B[i]; 
        //out[i] = sin( i*0.1 ); 
    }
};
