#!/usr/bin/env python3
import os, re, argparse, subprocess, sys

import numpy as np
import matplotlib
import matplotlib.pyplot as plt


def _parse_scan_params_from_comment(line):
    """Parse e.g. '# dist 0.600 Etot -434.750' or '# ang 1.23 dist 0.9' from comment line.
    Returns dict with float values if present.
    """
    if line is None: return {}
    s=line.strip()
    if '|' in s:
        # keep only tail after first '|', because our headers are: 'na,ne,core ... | <aux>'
        s=s.split('|',1)[1].strip()
    if s.startswith('#'): s=s[1:].strip()
    out={}
    for k in ['dist','ang','Etot','E','E_DFT','dft','DFT']:
        m=re.search(rf"\b{k}\b\s*([+-]?[0-9]*\.?[0-9]+)", s)
        if m:
            out[k]=float(m.group(1))
    return out


def extract_xyz_comments(xyz_path):
    """Return list of comment lines (2nd line of each block)."""
    comments=[]
    with open(xyz_path,'r') as f:
        while True:
            line=f.readline()
            if not line: break
            if not line.strip():
                continue
            try:
                n=int(line.strip().split()[0])
            except:
                continue
            c=f.readline()
            comments.append(c.rstrip('\n'))
            for _ in range(n):
                f.readline()
    return comments


def extract_xyz_atom_frames(xyz_path):
    """Parse multi-geometry xyz-with-electrons file; return (na, frames) where frames is list of (elems, pos[na,3]).

    Requires 'na,ne,core <na> <ne> <mode>' in comment line.
    """
    frames=[]
    na=None
    with open(xyz_path,'r') as f:
        while True:
            line=f.readline()
            if not line: break
            if not line.strip():
                continue
            try:
                ntot=int(line.strip().split()[0])
            except:
                continue
            c=f.readline()
            if c is None: break
            m=re.search(r"na,ne,core\s+(\d+)\s+(\d+)\s+(\w)", c)
            if not m:
                # skip block
                for _ in range(ntot):
                    f.readline()
                continue
            na_ = int(m.group(1))
            if na is None:
                na = na_
            elif na != na_:
                raise RuntimeError(f"Inconsistent na in {xyz_path}: {na} vs {na_}")
            elems=[]
            pos=np.zeros((na,3), dtype=float)
            iatom=0
            # read ntot lines, but allow blank lines
            nread=0
            while nread < ntot:
                l=f.readline()
                if not l: break
                s=l.split()
                if len(s)==0 or s[0].startswith('#'):
                    continue
                if iatom < na:
                    elems.append(s[0])
                    pos[iatom,:]=[float(s[1]),float(s[2]),float(s[3])]
                    iatom += 1
                nread += 1
            if iatom != na:
                raise RuntimeError(f"Failed to read {na} atoms in a frame of {xyz_path}, got {iatom}")
            frames.append((elems,pos))
    return na, frames


def compute_scan_coords_from_geom(xyz_path):
    """Compute (dist,ang) arrays from atomic geometry when not present in headers."""
    na, frames = extract_xyz_atom_frames(xyz_path)
    nconf=len(frames)
    dist=np.full(nconf, np.nan, dtype=float)
    ang =np.full(nconf, np.nan, dtype=float)

    if na == 2:
        for i,(es,ps) in enumerate(frames):
            dist[i]=float(np.linalg.norm(ps[1]-ps[0]))
        return dist, ang

    # H2O heuristic: 3 atoms, one O ("O" or "8") and two H ("H" or "1")
    if na == 3:
        for i,(es,ps) in enumerate(frames):
            # map numeric symbols
            es2=[('O' if e=='8' else ('H' if e=='1' else e)) for e in es]
            try:
                iO = es2.index('O')
            except ValueError:
                continue
            iH = [j for j,e in enumerate(es2) if e=='H']
            if len(iH) < 2:
                continue
            pO=ps[iO]; pH1=ps[iH[0]]; pH2=ps[iH[1]]
            v1=pH1-pO; v2=pH2-pO
            c=np.dot(v1,v2)/(np.linalg.norm(v1)*np.linalg.norm(v2))
            c=np.clip(c,-1.0,1.0)
            ang[i]=float(np.arccos(c))
            dist[i]=float(np.linalg.norm(pH1-pO))
        return dist, ang

    return dist, ang


def _get_libasan():
    try:
        return subprocess.check_output(['gcc','-print-file-name=libasan.so'], text=True).strip()
    except Exception:
        return ''


def _run_worker(worker, xyz_path, out_npy, nloc=32, device_index=0):
    cmd=[sys.executable, '-u', os.path.abspath(__file__), xyz_path, '--worker', worker, '--out-npy', out_npy, '--nloc', str(nloc), '--device', str(device_index)]
    env=os.environ.copy()
    if worker=='cpu':
        libasan=_get_libasan()
        if not libasan or not os.path.isfile(libasan):
            raise RuntimeError('CPU worker requested but libasan.so not found via `gcc -print-file-name=libasan.so`')
        env['LD_PRELOAD']=libasan + ((':'+env['LD_PRELOAD']) if 'LD_PRELOAD' in env and env['LD_PRELOAD'] else '')
        env.setdefault('ASAN_OPTIONS','detect_leaks=0')
        if 'CPP_BUILD_PATH' not in env or not env['CPP_BUILD_PATH']:
            # prefer optimized build; fallback to dbg
            cand1=os.path.normpath(os.path.join(os.path.dirname(__file__), '../../cpp/Build-opt/libs'))
            cand2=os.path.normpath(os.path.join(os.path.dirname(__file__), '../../cpp/Build-dbg/libs'))
            if os.path.isdir(cand1):
                env['CPP_BUILD_PATH']=cand1
            elif os.path.isdir(cand2):
                env['CPP_BUILD_PATH']=cand2
    elif worker=='gpu':
        # IMPORTANT: do NOT preload ASan when using PyOpenCL in this process
        env.pop('LD_PRELOAD', None)
    else:
        raise ValueError(worker)
    subprocess.check_call(cmd, env=env)


def _worker_cpu(xyz_path, out_npy):
    from pyBall import eFF as eff
    eff.setVerbosity(0,0)
    Es,_,_ = eff.processXYZ_e(xyz_path, nstepMax=0, xyz_out=None, fgo_out=None, bOutputs=(1,0,0))
    np.save(out_npy, np.array(Es, dtype=float))


def _worker_gpu(xyz_path, out_npy, nloc=32, device_index=0):
    from pyBall.OCL import eFF_ocl
    ocl = eFF_ocl.EFF_OCL(nloc=nloc, device_index=device_index, bEnergyKernel=True)
    ocl.load_xyzs(xyz_path, bPrint=False)
    if len(ocl.systems) == 0:
        raise RuntimeError(f"No systems parsed from {xyz_path}")
    ocl.realloc_buffers()
    ocl.upload_data()
    Es = ocl.eval_energies_localmd()
    np.save(out_npy, np.array(Es, dtype=float))


def _gridify(x, y, z):
    x=np.asarray(x); y=np.asarray(y); z=np.asarray(z)
    m=~(np.isnan(x)|np.isnan(y)|np.isnan(z))
    x=x[m]; y=y[m]; z=z[m]
    xs=np.unique(x)
    ys=np.unique(y)
    grid=np.full((len(xs),len(ys)), np.nan, dtype=float)
    xi={v:i for i,v in enumerate(xs)}
    yi={v:i for i,v in enumerate(ys)}
    for xx,yy,zz in zip(x,y,z):
        grid[xi[xx], yi[yy]] = zz
    return xs, ys, grid


def _stats(arr):
    arr=np.asarray(arr)
    m=~np.isnan(arr)
    if m.sum()==0: return dict(min=np.nan,max=np.nan,mean=np.nan,rms=np.nan)
    a=arr[m]
    return dict(min=float(a.min()), max=float(a.max()), mean=float(a.mean()), rms=float(np.sqrt(np.mean(a*a))))


def plot_1d(x, series, labels, title, xlabel, out_png=None, noshow=False):
    plt.figure(figsize=(8,5))
    def _style(lbl):
        u = lbl.upper()
        if 'CPU' in u:
            return dict(lw=1.5, ls=':', marker='o', ms=3)
        if 'GPU' in u:
            return dict(lw=0.8, ls='-', marker='o', ms=3)
        if 'DFT' in u:
            return dict(lw=1.0, ls='--', marker='o', ms=3)
        return dict(lw=1.0, ls='-', marker='o', ms=3)
    for y,lbl in zip(series,labels):
        if y is None: continue
        plt.plot(x, y, label=lbl, **_style(lbl))
    plt.xlabel(xlabel)
    plt.ylabel('Etot [eV]')
    plt.title(title)
    plt.legend()
    plt.tight_layout()
    if out_png: plt.savefig(out_png, dpi=200)
    if not noshow: plt.show()
    plt.close()


def plot_2d(xs, ys, grids, labels, title, out_png=None, noshow=False):
    # Expect 3 grids: CPU, GPU, Diff; place side-by-side (1 row, 3 cols)
    n=len(grids)
    fig,axs=plt.subplots(1,n, figsize=(5*n,4), squeeze=False)
    axs=axs[0]

    # Determine color scales independently for CPU/GPU; diff centered at 0
    for i,(g,lbl) in enumerate(zip(grids,labels)):
        ax=axs[i]
        if g is None:
            ax.set_axis_off();
            continue
        if lbl.lower().startswith('diff'):
            vmax=np.nanmax(np.abs(g))
            vmin=-vmax
            cmap='coolwarm'
        else:
            vmin=np.nanmin(g)
            vmax=np.nanmax(g)
            cmap='inferno'
        im=ax.imshow(g.T, origin='lower', aspect='auto', extent=[xs.min(),xs.max(),ys.min(),ys.max()], vmin=vmin, vmax=vmax, cmap=cmap)
        ax.set_title(lbl)
        ax.set_xlabel('dist')
        ax.set_ylabel('ang')
        fig.colorbar(im, ax=ax, shrink=0.9)
    fig.suptitle(title)
    fig.tight_layout()
    if out_png: fig.savefig(out_png, dpi=200)
    if not noshow: plt.show()
    plt.close(fig)


def main():
    ap=argparse.ArgumentParser()
    ap.add_argument('xyz', help='electron-augmented multi-geometry xyz scan file')
    ap.add_argument('--worker', choices=['cpu','gpu'], default=None, help='internal: run only CPU/GPU eval and write --out-npy')
    ap.add_argument('--out-npy', default=None)
    ap.add_argument('--mode', choices=['auto','1d','2d'], default='auto')
    ap.add_argument('--outdir', default='export/plots_cpu_gpu')
    ap.add_argument('--noshow', action='store_true')
    ap.add_argument('--png', action='store_true', help='save png into outdir')
    ap.add_argument('--nloc', type=int, default=32)
    ap.add_argument('--device', type=int, default=0)
    args=ap.parse_args()

    if args.worker is not None:
        if not args.out_npy:
            raise RuntimeError('--worker requires --out-npy')
        if args.worker=='cpu':
            _worker_cpu(args.xyz, args.out_npy)
        else:
            _worker_gpu(args.xyz, args.out_npy, nloc=args.nloc, device_index=args.device)
        return

    os.makedirs(args.outdir, exist_ok=True)

    comments = extract_xyz_comments(args.xyz)
    pars = [ _parse_scan_params_from_comment(c) for c in comments ]
    dist = np.array([p.get('dist', np.nan) for p in pars], dtype=float)
    ang  = np.array([p.get('ang' , np.nan) for p in pars], dtype=float)
    dft  = np.array([p.get('E_DFT', p.get('DFT', p.get('dft', np.nan))) for p in pars], dtype=float)

    if (not np.isfinite(dist).any()) and (not np.isfinite(ang).any()):
        dist_g, ang_g = compute_scan_coords_from_geom(args.xyz)
        if np.isfinite(dist_g).any() or np.isfinite(ang_g).any():
            dist = dist_g
            ang  = ang_g

    base=os.path.basename(args.xyz)
    tag=os.path.splitext(base)[0]

    tmp_cpu=os.path.join(args.outdir, f".{tag}__cpu.npy")
    tmp_gpu=os.path.join(args.outdir, f".{tag}__gpu.npy")
    _run_worker('cpu', args.xyz, tmp_cpu, nloc=args.nloc, device_index=args.device)
    _run_worker('gpu', args.xyz, tmp_gpu, nloc=args.nloc, device_index=args.device)
    cpuEs=np.load(tmp_cpu)
    gpuEs=np.load(tmp_gpu)
    cpuE = cpuEs[:,0]
    gpuE = gpuEs[:,0]

    # decide 1D vs 2D
    mode=args.mode
    dist_f = dist[np.isfinite(dist)]
    ang_f  = ang [np.isfinite(ang )]
    n_ud   = len(np.unique(dist_f)) if dist_f.size else 0
    n_ua   = len(np.unique(ang_f )) if ang_f.size  else 0
    if mode=='auto':
        mode = '2d' if (n_ud>1 and n_ua>1) else '1d'

    # stats + table
    rows=[]
    rows.append((tag,'DFT',_stats(dft)))
    rows.append((tag,'CPU',_stats(cpuE)))
    rows.append((tag,'GPU',_stats(gpuE)))
    diff = gpuE - cpuE
    rows.append((tag,'GPU-CPU',_stats(diff)))
    diff_stats = rows[-1][2]

    table_path=os.path.join(args.outdir, f"{tag}__stats.txt")
    with open(table_path,'w') as f:
        f.write("name\tseries\tmin\tmax\tmean\trms\n")
        for name,series,st in rows:
            f.write(f"{name}\t{series}\t{st['min']:.6g}\t{st['max']:.6g}\t{st['mean']:.6g}\t{st['rms']:.6g}\n")
    print('wrote',table_path)

    if mode=='1d':
        if n_ud>1:
            x=dist; xlabel='dist'
        elif n_ua>1:
            x=ang;  xlabel='ang'
        else:
            raise RuntimeError('No dist/ang scan coordinate found in comments')
        out_png=os.path.join(args.outdir, f"{tag}__1d.png") if args.png else None
        plot_1d(x, [dft if np.isfinite(dft).any() else None, cpuE, gpuE], ['DFT','EFF_CPU','EFF_GPU'], tag, xlabel=xlabel, out_png=out_png, noshow=args.noshow)
    else:
        xs, ys, g_cpu = _gridify(dist, ang, cpuE)
        _,  _,  g_gpu = _gridify(dist, ang, gpuE)
        gd = g_gpu - g_cpu
        out_png=os.path.join(args.outdir, f"{tag}__2d.png") if args.png else None
        plot_2d(xs, ys, [g_cpu,g_gpu,gd], ['EFF_CPU','EFF_GPU','GPU-CPU'], tag, out_png=out_png, noshow=args.noshow)

    print(f"{tag}: GPU-CPU diff stats min={diff_stats['min']:.6g} max={diff_stats['max']:.6g} mean={diff_stats['mean']:.6g} rms={diff_stats['rms']:.6g}")

if __name__ == '__main__':
    main()
