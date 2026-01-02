import sys
import numpy as np
import os
import matplotlib.pyplot as plt
import argparse

sys.path.append("../../")

from pyBall.tests import ocl_GridFF_new as gff


# ==== Just plot basis
# t = np.linspace(0.0, 1.0, 10)
# ws = gff.Bspline_basis5(t)
# plt.plot(ws)
# plt.figure()
# plt.imshow( ws )
# plt.show()
# exit()




# =============  Functions

R0 = 3.5
E0 = 0.1
a  = 1.8

Q = 0.4
p0 = [-2.0,-2.0,0.0]


def main():
    parser = argparse.ArgumentParser(description="Run GridFF OCL tests")
    parser.add_argument("--name",       default="NaCl_1x1_L3", help="Base name of xyz file in data/xyz (default: NaCl_1x1_L3)")
    parser.add_argument("--job",        default="MorseFit", choices=["MorseFit","PLQ","PLQ_lin","CG"], help="Job type (default: MorseFit)")
    parser.add_argument("--save-name",  default="double3", help="Save name for outputs (default: double3)")
    parser.add_argument("--use-CG",     type=int,     default=1, help="Use CG solver for MorseFit (default: 1)")
    parser.add_argument("--use-tiled",  type=int,     default=0, help="Use CG solver for MorseFit (default: 1)")
    parser.add_argument("--nmaxiter",   type=int,     default=10000, help="Override nmaxiter for fit3D/fit3D_CG")
    parser.add_argument("--nPerStep",   type=int,     default=10, help="Override nPerStep for fit3D/fit3D_CG")
    parser.add_argument("--damp",       type=float,   default=0.15, help="Override damp for MD solver")
    parser.add_argument("--save-fig",   type=int,     default=1, help="Save convergence figure")
    parser.add_argument("--fig-path",   default=None, help="Path to save convergence figure (default: auto under data/<name>/)")
    args = parser.parse_args()

    name = args.name
    job  = args.job
    save_name = args.save_name

    use_CG = args.use_CG
    # If job=="CG", force CG path
    if job == "CG":
        job = "MorseFit"
        use_CG = True
    # defaults matching previous behavior
    nmaxiter = args.nmaxiter
    nPerStep = args.nPerStep
    damp     = args.damp

    kwargs = {}
    if nmaxiter is not None: kwargs["nmaxiter"] = nmaxiter
    if nPerStep is not None: kwargs["nPerStep"] = nPerStep
    if damp is not None:     kwargs["damp"]     = damp
    if use_CG is not None:   kwargs["use_CG"]   = bool(use_CG)
    kwargs["use_tiled"] = bool(args.use_tiled)
    if args.save_fig:        kwargs["save_fig"] = True
    if args.fig_path is not None: kwargs["fig_path"] = args.fig_path

    gff.test_gridFF_ocl( fname=f"data/xyz/{name}.xyz", save_name=save_name, job=job, **kwargs )


if __name__ == "__main__":
    main()
# plt.figure( figsize=(15,5)); 
# plt.subplot(1,3,1); plt.imshow( VPaul[:,:,1] ); plt.colorbar(); plt.title("VPaul fit GPU")
# plt.subplot(1,3,2); plt.imshow( VLond[:,:,1] ); plt.colorbar(); plt.title("VLond fit GPU")
# plt.subplot(1,3,3); plt.imshow( VCoul[:,:,1] ); plt.colorbar(); plt.title("VCoul fit GPU")
# plt.show()

# ======== Ewald

d=0.6
apos=np.array([
    [-d,.0,0.],
    [+d,.0,0.],
    [0.,-d,0.],
    [0.,+d,0.],
])
qs = [ +1.,+1.,-1.,-1. ]

# d=0.6
# apos=np.array([
#     [0.,.0,-d],
#     [0.,.0,+d],
# ])
# qs = [ +1.,-1. ]

#gff.test_Ewald( apos, qs,  ns=(100,100,100), dg=(0.10,0.10,0.10), order=3, bSlab=True,  nPBC=(100,100,0), bOld=True )
#gff.test_Ewald( apos, qs,  ns=(100,100,100), dg=(0.10,0.10,0.10), order=3, bSlab=True,  nPBC=(100,100,0) )
# gff.test_Ewald( apos, qs,  ns=[100,100,150], dg=[0.10,0.10,0.10], order=3, bSlab=True,  nPBC=[100,100,0] )
# gff.test_Ewald( apos, qs,  ns=[100,100,200], dg=[0.10,0.10,0.10], order=3, bSlab=True,  nPBC=[100,100,0] )
# gff.test_Ewald( apos, qs,  ns=[100,100,300], dg=[0.10,0.10,0.10], order=3, bSlab=True,  nPBC=[100,100,0] )
# gff.test_Ewald( apos, qs,  ns=[100,100,400], dg=[0.10,0.10,0.10], order=3, bSlab=True,  nPBC=[100,100,0] )

plt.show()