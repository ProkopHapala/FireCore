#!/bin/bash

# from https://www.perplexity.ai/search/learn-how-to-do-proper-profili-h3ulu.hhR1uCc3wakOlp_Q

#python3 test_mmff_cuda_vs_ocl.py

sudo /usr/local/cuda-12.1/bin/ncu --metrics sm__sass_thread_inst_executed_op_fp32_pred_on,sm__cycles_active,dram__bytes_read,dram__bytes_write,l1tex__t_sectors_pipe_lsu_mem_global_op_ld,l1tex__t_sectors_pipe_lsu_mem_global_op_st --target-processes all ./run.sh | tee ncu.OUT

#ncu --metrics sm__sass_thread_inst_executed_op_fp32_pred_on,sm__cycles_active,dram__bytes_read,dram__bytes_write,l1tex__t_sectors_pipe_lsu_mem_global_op_ld,l1tex__t_sectors_pipe_lsu_mem_global_op_st --target-processes all ./run.sh | tee ncu.OUT