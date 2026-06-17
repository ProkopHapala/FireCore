import sys, os
import numpy as np

_THIS_DIR = os.path.abspath(os.path.dirname(__file__))
_SHARED_DIR = os.path.abspath(os.path.join(_THIS_DIR, '..', 'shared'))
for _p in (_THIS_DIR, _SHARED_DIR):
    if _p not in sys.path:
        sys.path.insert(0, _p)

from RRsp3 import RRsp3


def main():
    tmp = os.environ.get('RRSP3_TMPDIR', None)
    if not tmp:
        raise RuntimeError('RRSP3_TMPDIR not set')
    nnode_per_group = int(os.environ.get('RRSP3_NNODE_PER_GROUP', '1'))

    pos = np.load(os.path.join(tmp, 'pos.npy'))
    invm = np.load(os.path.join(tmp, 'invm.npy'))
    quat = np.load(os.path.join(tmp, 'quat.npy'))
    rad = np.load(os.path.join(tmp, 'rad.npy'))
    neighs = np.load(os.path.join(tmp, 'neighs.npy'))
    excl1 = np.load(os.path.join(tmp, 'excl1.npy'))
    excl2 = np.load(os.path.join(tmp, 'excl2.npy'))
    port = np.load(os.path.join(tmp, 'port.npy'))
    K = np.load(os.path.join(tmp, 'K.npy'))
    bkSlots = np.load(os.path.join(tmp, 'bkSlots.npy'))

    build_opts_str = os.environ.get('RRSP3_BUILD_OPTS', '')
    build_opts = [s for s in build_opts_str.split(' ') if s]

    natoms = int(pos.shape[0])
    sim = RRsp3(natoms, group_size=64, prefer_gpu=True, build_options=build_opts)

    sim.upload_state(pos, invm, quat=quat)
    sim.upload_radius(rad)
    sim.upload_neighs_and_exclusions(neighs, excl1, excl2)
    sim.upload_cluster_ports(port, K, nnode_per_group=nnode_per_group)
    sim.upload_bkSlots(bkSlots)

    sim.step_cluster(nnode_per_group=nnode_per_group, dt=0.1, k_coll=50.0, relaxation=0.5, bbox_margin=0.5)


if __name__ == '__main__':
    main()
