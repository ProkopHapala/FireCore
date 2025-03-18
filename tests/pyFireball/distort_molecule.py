import sys
import os
import numpy as np
np.set_printoptions(linewidth=1000)

sys.path.append("../../")
from pyBall import atomicUtils as au
from pyBall.AtomicSystem import AtomicSystem
from pyBall import FireCore as fc


# ======== Body

path = './relaxed_mols/'
#mol = au.AtomicSystem( path+'CH2NH.xyz' )
#mol = au.AtomicSystem( path+'NH2NH2.xyz' )
#mol = au.AtomicSystem( path+'NH2OH.xyz' )
#mol = au.AtomicSystem( path+'CH3NH2.xyz' )



def move_group( mol, gi, gj, l=0.1 ):
    #i,j=ij
    #print( "move_group : groups", groups )
    #gi=groups[i]
    #gj=groups[j]
    #print( "move_group i,j: ", i,j )
    ia = gi[0]
    ja = gj[0]
    #print( "move_group ia,ja: ", ia,ja )
    d  = mol.apos[ja] - mol.apos[ia]
    mol.apos[gj,:] += d*l
    #mol.apos[gi,:] -= d*l

def scan_group_distance( mol, groups, l0=-0.4, dl=0.1 ):
    fname = path+"scan_"+name
    print( "scan_group_distance: ", fname )
    gi,gj = groups[0],groups[1]
    move_group( mol, gi, gj, l0 )
    mol.saveXYZ( fname, mode='w' )
    for i in range(10):
        move_group( mol, gi, gj, dl )
        mol.saveXYZ( fname, mode='a' )

def scan_angle( mol, angs, ia=0, j=1, jbs=[2] ):
    fname = path+"angscan_"+name
    #r1 = mol.apos[j1] - mol.apos[i]
    #r2 = mol.apos[j2] - mol.apos[i]
    ax = au.normalize(mol.apos[j] - mol.apos[ia]) 
    vs = []
    rs = []
    for j in jbs:
        v = mol.apos[j] - mol.apos[ia]
        r = np.linalg.norm(v)
        vs.append(v/r)
        rs.append(r)
    mol.saveXYZ( fname, mode='w' )
    for a in angs:
        ca = np.cos(a)
        sa = np.sin(a)
        for jj,j in enumerate(jbs):
            r = rs[jj]
            mol.apos[j,:] = mol.apos[ia,:] + ax*(r*ca) + vs[jj]*(r*sa)
        mol.saveXYZ( fname, comment=f"ang {a}", mode='a' )

def scan_angle_dist( mol, angs, dists, ia=0, j=1, jbs=[2], bEpairs=False, bFireball=False ):
    fname = path+"angdistscan_"+name
    #r1 = mol.apos[j1] - mol.apos[i]
    #r2 = mol.apos[j2] - mol.apos[i]
    ax = au.normalize(mol.apos[j] - mol.apos[ia]) 
    vs = []
    rs = []
    for j in jbs:
        v = mol.apos[j] - mol.apos[ia]
        r = np.linalg.norm(v)
        vs.append(v/r)
        #rs.append(r)
    if bFireball:
        fc.initialize( atomType=mol.atypes, atomPos=mol.apos, verbosity=3 )
        fc.evalForce( mol.apos, nmax_scf=100 )
        Emap = np.zeros( (len(angs),len(dists)) )

    mol.saveXYZ( fname, mode='w' )
    if bEpairs:
        mol.findBonds()
    for iang,a in enumerate(angs):
        ca = np.cos(a)
        sa = np.sin(a)
        for ir,r in enumerate(dists):
            label = f"# ang {a:.3f} dist {r:.3f}"
            for jj,j in enumerate(jbs):
                mol.apos[j,:] = mol.apos[ia,:] + ax*(r*ca) + vs[jj]*(r*sa)
            if bFireball:
                #print("============= CALL FireBall ", r, a)
                forces, Es = fc.evalForce( mol.apos, nmax_scf=100 )
                #print( "Forces:", forces )
                #print( "Energies:", Es, Es[1:].sum() )
                #print( f"FireBall( a: {a} r: {r} ) Energies:", Es )
                Emap[iang,ir] = Es[0]
                label += f" Etot  {Es[0]:.3f}"
                print( "FireBall " + label )
            if bEpairs:
                mol_ = mol.clonePBC()
                mol_.bonds = mol.bonds
                mol_.neighs()
                mol_.add_electron_pairs()
                mol_.saveXYZ( fname, comment=label, mode='a' )
            else:    
                mol.saveXYZ( fname, comment=label, mode='a' )
    
    if bFireball:
        return Emap

def scan_rot(mol, sel, angs, ax, p0=None, i0=0 ):
    '''
    mol : AtomicSystem
    p0  : center of rotation
    ax  : axis of rotation
    angs : angles to scan
    sel : selected atoms to rotate
    implementation: 
      - we also create rotation matrix initial angle
      - then we run in loop creating rotation matrix for difference between successive angles and multiply the selected atom positions by it (substracting and addind the center of rotation) 
    '''
    if p0 is None: p0 = mol.apos[i0]
    fname = path + "rotscan_" + name
    ax = au.normalize(ax)
    oa = angs[0]
    R0 = au.rotation_matrix(ax, oa)
    for i in sel: mol.apos[i,:] = np.dot(mol.apos[i,:]-p0, R0.T) + p0
    mol.saveXYZ(fname, comment=f"ang {oa}", mode='w')
    for a in angs[1:]:
        dR = au.rotation_matrix(ax, a-oa )
        print( f"{a:.3f}  dR", dR )
        for i in sel: mol.apos[i,:] = np.dot(mol.apos[i,:]-p0, dR.T) + p0
        mol.saveXYZ(fname, comment=f"ang {a}", mode='a')
        oa = a
    return mol

def scan_rot_scale(mol, sel, angs, scales, ax, p0=None, i0=0, bEpairs=False, bFireball=False):
    '''
    Scan rotation and scaling of selected atoms.
    
    mol : AtomicSystem
    p0  : center of rotation
    ax  : axis of rotation
    angs : angles to scan
    scales : scaling factors to scan
    sel : selected atoms to rotate and scale
    bEpairs : whether to add electron pairs
    bFireball : whether to use FireBall for energy calculation
    '''
    if p0 is None: p0 = mol.apos[i0]
    fname = path + "rotscalescan_" + name
    ax = au.normalize(ax)
    oa = angs[0]
    R0 = au.rotation_matrix(ax, oa)
    
    if bFireball:
        fc.initialize(atomType=mol.atypes, atomPos=mol.apos, verbosity=3)
        fc.evalForce(mol.apos, nmax_scf=100)
        Emap = np.zeros((len(angs), len(scales)))
    
    mol.saveXYZ(fname, mode='w')
    if bEpairs:
        mol.findBonds()

    rs = [ np.linalg.norm(mol.apos[i] - p0) for i in sel ]
        
    for iang, a in enumerate(angs):
        dR = au.rotation_matrix(ax, a-oa) if iang > 0 else R0
        for i in sel:
            mol.apos[i,:] = np.dot(mol.apos[i,:]-p0, dR.T) + p0

        for iscale, s in enumerate(scales):
            for ii, i in enumerate(sel):
                vec = mol.apos[i,:] - p0
                vec = au.normalize(vec) * (rs[ii]*s)
                mol.apos[i,:] = p0 + vec
                
            label = f"# ang {a:.3f} scale {s:.3f}"
            
            if bFireball:
                forces, Es = fc.evalForce(mol.apos, nmax_scf=100)
                Emap[iang, iscale] = Es[0]
                label += f" Etot {Es[0]:.3f}"
                print("FireBall " + label)
            if bEpairs:
                mol_ = mol.clonePBC()
                mol_.bonds = mol.bonds
                mol_.neighs()
                mol_.add_electron_pairs()
                mol_.saveXYZ(fname, comment=label, mode='a')
            else:
                mol.saveXYZ(fname, comment=label, mode='a')
        oa = a
    
    if bFireball:
        return Emap
    return mol


# find all molecules in path
#molecules = [f for f in os.listdir(path) if os.path.isfile(os.path.join(path, f))]
#print(molecules)
rad2deg = 180./np.pi

#molecules = ["H2O.xyz"]
molecules = ["CH4.xyz"]

# mol    = AtomicSystem( path+"CH4.xyz" )
# #mol.orient(0, (1,2), (3,4) )
# mol.orient(1, (2,3), (4,5) )
# mol.saveXYZ( path+"CH4_oriented.xyz" )
# exit()

#mol = AtomicSystem( path+molecules[0] )
for name in molecules:
    mol    = AtomicSystem( path+name )
    
    # scan group distance
    # groups = mol.find_groups()
    # groups = list(groups.values())
    # groups.sort( key=lambda g: g[0] )
    # print(  name + ".groups ", groups )
    # if len(groups)<2: continue
    # scan_group_distance( mol, groups, l0=-0.4, dl=0.1 )

    '''
    scales = np.linspace(0.7,2.0,20)
    angs = np.linspace(-np.pi/3, np.pi/3, 20);   print( "angs", angs )
    #scan_rot(mol, [1,2], angs, [0.0,0.0,1.0], i0=0 )
    #scan_rot(mol, [1,2], angs, [1.0,0.0,0.0], i0=0 )
    #scan_rot_scale(mol, [1,2], angs, scales, [1.0,0.0,0.0], p0=None, i0=0, bEpairs=False, bFireball=False)
    Emap = scan_rot_scale(mol, [1,2], angs, scales, [1.0,0.0,0.0], p0=None, i0=0, bEpairs=False, bFireball=True)
    extent = [scales[0]*100,scales[-1]*100,angs[0]*rad2deg,angs[-1]*rad2deg]
    '''
    
    
    # scan angle
    angs = np.linspace(np.pi/4, np.pi*7./8., 20)
    dists = np.linspace(0.7,2.5,20)
    #scan_angle( mol, angs, ia=0, j=1, jbs=[2] )
    #scan_angle_dist( mol, angs, dists, ia=0, j=1, jbs=[2] )
    #scan_angle_dist( mol, angs, dists, ia=0, j=1, jbs=[2], bEpairs=True )
    Emap = scan_angle_dist( mol, angs, dists, ia=0, j=1, jbs=[2], bFireball=True )
    extent = [dists[0]*100,dists[-1]*100,angs[0]*rad2deg,angs[-1]*rad2deg]
    

    import matplotlib.pyplot as plt
    Emin = Emap.min()
    Emax = Emin+5.0
    plt.imshow( Emap, origin='lower', extent=extent, vmin=Emin, vmax=Emax )
    plt.colorbar()
    plt.xlabel( "Distance [pm]" )
    #plt.xlabel( "scale [%]" )
    plt.ylabel( "Angle [deg.]" )
    plt.title( name )
    plt.savefig( name + ".png" )
    plt.show()