#!/usr/bin/env python3
import os, sys, subprocess, argparse, textwrap
import numpy as np


def _get_libasan():
    try:
        return subprocess.check_output(['gcc','-print-file-name=libasan.so'], text=True).strip()
    except Exception:
        return ''


def _default_cpp_build_path():
    here=os.path.dirname(os.path.abspath(__file__))
    cand1=os.path.normpath(os.path.join(here, '../../cpp/Build-opt/libs'))
    cand2=os.path.normpath(os.path.join(here, '../../cpp/Build-dbg/libs'))
    if os.path.isdir(cand1):
        return cand1
    if os.path.isdir(cand2):
        return cand2
    return ''


def worker_cpu(xyz, out_txt):
    from pyBall import eFF as eff
    eff.setVerbosity(5,0)
    Es,_,_ = eff.processXYZ_e(xyz, nstepMax=0, xyz_out=None, fgo_out=None, bOutputs=(1,0,0))
    E=np.array(Es, dtype=float)
    with open(out_txt,'w') as f:
        f.write(f"CPU file={xyz}\n")
        f.write(f"CPU Etot={E[0,0]:.10f} Ek={E[0,1]:.10f} Eee={E[0,2]:.10f} Eae={E[0,3]:.10f} Eaa={E[0,4]:.10f}\n")


def worker_gpu(xyz, out_txt, nloc=32, device_index=0, dbg_pair=False, idbg_sys=0, offload_core=False):
    from pyBall.OCL import eFF_ocl
    ocl = eFF_ocl.EFF_OCL(nloc=nloc, device_index=device_index, bEnergyKernel=True, dbg_pair=dbg_pair, idbg_sys=idbg_sys)
    ocl.load_xyzs(xyz, bPrint=False)
    ocl.realloc_buffers(); ocl.upload_data()
    coreConsts5 = None
    if offload_core:
        coreConsts5 = ocl.compute_core_constants_cpu()
    if dbg_pair:
        ocl.eval_energies(bOffloadCore=offload_core, coreConsts5=coreConsts5)
    Es = np.array(ocl.eval_energies_localmd(bOffloadCore=offload_core, coreConsts5=coreConsts5), dtype=float)
    with open(out_txt,'w') as f:
        f.write(f"GPU file={xyz}\n")
        f.write(f"GPU offload_core={int(offload_core)} dbg_pair={int(dbg_pair)}\n")
        f.write(f"GPU Etot={Es[0,0]:.10f} Ek={Es[0,1]:.10f} Eee={Es[0,2]:.10f} Eae={Es[0,3]:.10f} Eaa={Es[0,4]:.10f}\n")


def run_subprocess(worker, xyz, out_txt, nloc=32, device=0):
    cmd=[sys.executable, '-u', os.path.abspath(__file__), '--worker', worker, '--xyz', xyz, '--out', out_txt, '--nloc', str(nloc), '--device', str(device)]
    env=os.environ.copy()
    if worker=='cpu':
        libasan=_get_libasan()
        if libasan and os.path.isfile(libasan):
            env['LD_PRELOAD']=libasan + ((':'+env['LD_PRELOAD']) if env.get('LD_PRELOAD') else '')
            env.setdefault('ASAN_OPTIONS','detect_leaks=0')
        env.setdefault('CPP_BUILD_PATH', _default_cpp_build_path())
    else:
        env.pop('LD_PRELOAD', None)
    subprocess.check_call(cmd, env=env)


def main():
    ap=argparse.ArgumentParser()
    ap.add_argument('--worker', choices=['cpu','gpu'], default=None)
    ap.add_argument('--xyz', required=True)
    ap.add_argument('--out', required=True)
    ap.add_argument('--nloc', type=int, default=32)
    ap.add_argument('--device', type=int, default=0)
    ap.add_argument('--dbg_pair', action='store_true', default=False)
    ap.add_argument('--idbg_sys', type=int, default=0)
    ap.add_argument('--offload_core', action='store_true', default=False)
    args=ap.parse_args()

    if args.worker=='cpu':
        worker_cpu(args.xyz, args.out); return
    if args.worker=='gpu':
        worker_gpu(args.xyz, args.out, nloc=args.nloc, device_index=args.device, dbg_pair=args.dbg_pair, idbg_sys=args.idbg_sys, offload_core=args.offload_core); return

    # orchestrator mode
    out_base=args.out
    out_cpu=out_base + '__cpu.txt'
    out_gpu=out_base + '__gpu.txt'
    run_subprocess('cpu', args.xyz, out_cpu, nloc=args.nloc, device=args.device)
    cmd=[sys.executable, '-u', os.path.abspath(__file__), '--worker', 'gpu', '--xyz', args.xyz, '--out', out_gpu, '--nloc', str(args.nloc), '--device', str(args.device)]
    if args.dbg_pair: cmd += ['--dbg_pair', '--idbg_sys', str(args.idbg_sys)]
    if args.offload_core: cmd += ['--offload_core']
    env=os.environ.copy(); env.pop('LD_PRELOAD', None)
    subprocess.check_call(cmd, env=env)

    def _parse_E(path):
        with open(path,'r') as f:
            for line in f:
                if line.startswith('CPU Etot=') or line.startswith('GPU Etot='):
                    s=line.split()
                    kv={}
                    for item in s[1:]:
                        if '=' in item:
                            k,v=item.split('=')
                            kv[k]=float(v)
                    return kv
        return {}

    Ec=_parse_E(out_cpu)
    Eg=_parse_E(out_gpu)
    keys=['Etot','Ek','Eee','Eae','Eaa']
    lines=[]
    lines.append(f"XYZ {args.xyz}")
    lines.append("component        CPU               GPU               GPU-CPU")
    for k in keys:
        vc=Ec.get(k,np.nan); vg=Eg.get(k,np.nan)
        lines.append(f"{k:8s} {vc:16.8f} {vg:16.8f} {vg-vc:16.8f}")
    rep='\n'.join(lines)+'\n'
    with open(out_base + '__report.txt','w') as f:
        f.write(rep)
    print(rep)


if __name__ == '__main__':
    main()
