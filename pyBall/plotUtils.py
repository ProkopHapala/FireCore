import numpy as np
import matplotlib.pyplot as plt
from   matplotlib import collections  as mc
from . import elements
#from . import utils as ut



def read_gnuplot_2d(fname):
    f = open(fname,'r')
    xs=[]
    ys=[]
    vals=[]
    nx=-1
    il=0
    for iil, l in enumerate(f):
        ws=l.split()
        if len(ws)<3:
            if(nx<0): 
                nx=int(il)
            print( iil, il )
            il=0
        else:
            xs  .append( float(ws[0]) )
            ys  .append( float(ws[1]) )
            vals.append( float(ws[2]) )
            il+=1
    xs   = np.array(xs)  .reshape(-1,nx)
    ys   = np.array(ys)  .reshape(-1,nx)
    vals = np.array(vals).reshape(-1,nx)
    return  vals, xs, ys

def read_dat( fname, ni=0, nf=1, iname=0, toRemove=None ):
    #format:            1  -96.294471702523595       -251.76919147019100       -48.443292828581697       # HHH-hhS1_NNO-hpS1 HHH-hhS1_NNO-hpS1 
    f=open( fname, 'r' )
    ints  =[]
    floats=[]
    names =[] 
    nc = ni+nf+2 # number of columns in the file
    for l in f:
        ws = l.split()
        nw = len(ws)
        if(nw<nc):
            ints_i    = [-1    ]*ni
            floats_i  = [np.nan]*nf
            name      = ws[nw-1]
        else:
            ints_i   = [ int(ws[i])   for i in range(0    ,ni           ) ]
            floats_i = [ float(ws[i]) for i in range(ni   ,ni+nf        ) ]
            name     = ws[ni+nf+1+iname]

        if toRemove is not None:
            if name in toRemove:
                continue

        ints      .append( ints_i   )
        floats    .append( floats_i )
        names     .append( name )
    return ints,floats,names

def plotAtoms( apos=None, es=None, atoms=None, bNumbers=False, labels=None, sizes=100., colors='#808080', marker='o', axes=(0,1), selection=None ):
    ax1,ax2=axes
    if apos is None: apos = np.array([ a[1] for a in atoms ])  #;print(apos)
    if selection is not None: 
        apos                  = apos[selection,:]
        if es is not None: es = es  [selection]
    #print( "apos.shape ", apos.shape )
    #print( "apos ", apos )
    plt.scatter( apos[:,ax1],apos[:,ax2], marker=marker, c=colors, s=sizes, cmap='seismic', zorder=2 ); plt.axis('equal'); #plt.grid()
    bLabels = labels is not None
    if bNumbers or bLabels:
        na = len(apos)
        if not bLabels:
            labels = range(na)
        ax= plt.gca()
        for i in range(na):
            ax.annotate( str(labels[i]), (apos[i,ax1], apos[i,ax2]))

def plotBonds( lps=None, links=None, ps=None, lws=None, axes=(0,1), colors='k', labels=None, ls='solid', fnsz=10, fnclr='k', fontweight=None ):
    ax1,ax2=axes
    ax_inds=[ax1,ax2]
    if lps is None:
        ps_=ps
        if ps_.shape[1]!=2:
            ps_ = ps[:,ax_inds]
        links=np.array(links,dtype=np.int32)
        lps=np.zeros( (len(links),2,2) )
        #print( "lps.shape, ps_.shape ", lps.shape, ps_.shape )
        lps[:,0,:] = ps_[ links[:,0],: ]
        lps[:,1,:] = ps_[ links[:,1],: ]
    #lws = (kek.bondOrder-0.8)*5
    lc = mc.LineCollection(lps, linewidths=lws, colors=colors, linestyle=ls )
    ax= plt.gca()
    ax.add_collection(lc)

    if labels is not None:
        for i, s in enumerate(labels):
            p = (lps[i,0,:]+lps[i,1,:])*0.5
            ax.annotate( str(s), p, size=fnsz, color=fnclr, fontweight=fontweight )
    #ax.autoscale()
    #ax.margins(0.1)

def plotAngles( iangs, angs, ps, axes=(0,1), colors='k', labels=None, bPoly=True, alpha=0.2 ):
    ax1,ax2=axes
    ax_inds=[ax1,ax2]
    if isinstance(colors, str ):
        c = colors
        colors = [ c for i in range(len(iangs)) ]
    ax=plt.gca()
    for i,a in enumerate(iangs):
        iang=iangs[i]
        #print( ps.shape, iang, ax_inds ) 
        pp = ps[iang,:];
        pp = pp[:,ax_inds] ; #print(pp)
        t1 = plt.Polygon( pp, color=colors[i], fill=True, alpha=alpha, lw=0 )
        ax.add_patch(t1)
        p = (ps[iang[0],:] + ps[iang[1],:] + ps[iang[2],:])/3.0
        ax.annotate( "%3.0f˚" %(angs[i]*180.0/np.pi), p[ax_inds], color=colors[i] )


def plotSystem( sys , bBonds=True, colors=None, sizes=None, extent=None, sz=50., RvdwCut=0.5, axes=(0,1), bLabels=True, labels=None, _0=1, HBs=None, bHBlabels=True, bBLabels=False  ):    
    if( bBonds ):
        if sys.bonds is None:
            bs, rbs = sys.findBonds( Rcut=3.0, RvdwCut=RvdwCut )
        rb_labs = None
        if bBLabels:
            rb_labs = [ ("%3.2f" %r) for r in rbs ]
        plotBonds( links=sys.bonds, ps=sys.apos, axes=axes )
    
    enames = [ typ.split('_')[0] for typ in sys.enames ]
    if(colors is None): colors = [ elements.ELEMENT_DICT[e][8]    for e in enames ]
    if(sizes  is None): sizes  = [ elements.ELEMENT_DICT[e][6]*sz for e in enames ]
    if((labels is None) and bLabels): labels=[ "%s%i" %(e,i+_0) for i,e in enumerate(sys.enames) ]
    #print( "labels ", labels)
    plotAtoms( apos=sys.apos, es=sys.enames, sizes=sizes, colors=colors, marker='o', axes=axes, labels=labels )

    # H-Bonds
    if HBs is not None:
        if len(HBs)>0:
            hbs,rhbs = HBs
            rh_labs = None
            if bHBlabels:
                rh_labs = [ ("%3.2f" %r) for r in rhbs ]
            plotBonds ( ps=sys.apos, links=hbs, axes=axes, colors="g", labels=rh_labs )

    if extent is not None:
        plt.xlim(extent[0],extent[1])
        plt.ylim(extent[2],extent[3]) 

def plotTrj( trj, bBonds=True, sz=50., numbers=None, axes=(0,1), extent=None, prefix="mol_", RvdwCut=0.5, figsize=(5,5) ):
    if numbers is None: numbers=range(len(trj))
    for i,sys in enumerate(trj):
        print("plot # ", i)
        fig = plt.figure(figsize=figsize)
        
        if( bBonds ):
            if sys.bonds is None:
                sys.findBonds( Rcut=3.0, RvdwCut=RvdwCut )
                plotBonds( links=sys.bonds, ps=sys.apos, axes=(0,1) )
        
        colors = [ elements.ELEMENT_DICT[e][8]    for e in sys.enames ]
        sizes  = [ elements.ELEMENT_DICT[e][6]*sz for e in sys.enames ]
        plotAtoms( apos=sys.apos, es=sys.enames, sizes=sizes, colors=colors, marker='o', axes=axes )

        if extent is not None:
            plt.xlim(extent[0],extent[1])
            plt.ylim(extent[2],extent[3])
        
        plt.savefig( prefix+("%03i.png" %numbers[i]), bbox_inches='tight' )
        plt.close(fig)

def render_POVray(
    sys, filename, 
    # Atom and bond parameters
    atom_scale=1.0, bond_width=0.2,
    # Camera parameters
    look_at=(0.0, 0.0, 0.0),
    camera_pos  =( 0.0, 0.0, 100.0 ),
    camera_up   =( 0.0, 1.0,   0.0 ),
    camera_right=( 1.0, 0.0,   0.0 ),
    sky=(0.0, 0.0, 1.0),
    zoom=30.0,
    orthographic=True,
    # Image parameters
    width=400, height=400,
    # Lighting parameters
    light_pos       =(10.0, 20.0, 30.0),
    light_color     =(2.5, 2.5, 2.5),
    ambient_light   =(1.0, 1.0, 1.0),
    background_color=(1.0, 1.0, 1.0),
    # Material parameters
    ambient=0.5, diffuse=0.6, specular=0.4, roughness=0.01, metallic=-1.0,  phong=10.0, phong_size=120,
    # Rendering options
    shadows=False,
    bond_clr=(0.5,0.5,0.5),
    z_color_shift=False,
    viewAxis=False,
    ):
    """Export system to POV-Ray file with customizable rendering settings
    
    Args:
        sys: System object with atoms/bonds
        filename: Output .pov file
        
        # Atom and bond parameters
        bond_scale: Bond length relative to sum of atomic radii
        atom_scale: Scale factor for atomic radii
        bond_width: Visual width of bonds
        
        # Camera parameters
        camera_pos: Camera position (x,y,z)
        look_at: Point camera looks at (x,y,z)
        sky: Up vector for camera orientation
        zoom: Camera zoom factor
        orthographic: Use orthographic projection if True, perspective if False
        
        # Image parameters
        width, height: Output image dimensions
        
        # Lighting parameters
        light_pos: Position of main light source
        light_color: RGB color of main light
        ambient_light: RGB color of ambient light
        background_color: RGB color of background
        
        # Material parameters
        ambient: Ambient light reflection
        diffuse: Diffuse light reflection
        specular: Specular highlights intensity
        roughness: Surface roughness
        metallic: Metallic finish if True
        phong, phong_size: Phong shading parameters
        
        # Rendering options
        shadows: Enable shadows if True
        z_color_shift: Enable z-dependent color shifting if True
    """

    def makeFinishString( ambient, diffuse, specular, roughness, phong, phong_size, metallic ):
        return f"""
  finish {{
    ambient    {ambient}
    diffuse    {diffuse}
    specular   {specular}
    roughness  {roughness}
    phong      {phong}
    phong_size {phong_size}
    { f"metallic {metallic}" if metallic>0 else ""}
  }}
  """

    with open(filename, 'w') as pov:
        # Write POV header with customizable parameters
        pov.write(
f'''// ***********************************************
// Camera & other global settings
// ***********************************************

#declare Zoom = {zoom};
#declare Width = {width};
#declare Height = {height};

camera{{
  {"orthographic" if orthographic else ""}
  look_at  <{look_at[0]  }    , {look_at[1]  }    , {look_at[2]} >
  location <{camera_pos[0]}   , {camera_pos[1]}   , {camera_pos[2]}>
  up       <{camera_up[0]*zoom}    , {camera_up[1]*zoom}    , {camera_up[2]*zoom} >
  right    <{camera_right[0]*zoom} , {camera_right[1]*zoom} , {camera_right[2]*zoom}>
  //sky      <{sky[0]}       , {sky[1]}        , {sky[2]} >
  sky      <{camera_up[0]}       , {camera_up[1]}        , {camera_up[2]} >
}}

background      {{ color rgb <{background_color[0]}, {background_color[1]}, {background_color[2]}> }}
light_source    {{ < {light_pos[0]}, {light_pos[1]}, {light_pos[2]}>  rgb <{light_color[0]}, {light_color[1]}, {light_color[2]}> }}
global_settings {{ ambient_light rgb< {ambient_light[0]}, {ambient_light[1]}, {ambient_light[2]}> }}

// ===== macros for common shapes

#macro myFinish()
{makeFinishString( ambient, diffuse, specular, roughness, phong, phong_size, metallic )}
#end

#macro a(X,Y,Z,RADIUS,R,G,B,T)
 sphere{{<X,Y,Z>,RADIUS
  pigment{{rgbt<R,G,B,T>}}
  myFinish()
  {"no_shadow" if not shadows else ""}
 }}
#end

#macro b(X1,Y1,Z1,RADIUS1,X2,Y2,Z2,RADIUS2,R,G,B,T)
 cone{{<X1,Y1,Z1>,RADIUS1,<X2,Y2,Z2>,RADIUS2
  pigment{{rgbt<R,G,B,T>}}
  myFinish()
  {"no_shadow" if not shadows else ""}
 }}
#end

{
f"""
 b(    0.0,    0.0,   0.0,    {bond_width*1.5},   1.0,0.0,0.0,    {bond_width},    1.0,0.0,0.0, 0.0 )  // x-axis
 b(    0.0,    0.0,   0.0,    {bond_width*1.5},   0.0,1.0,0.0,    {bond_width},    0.0,1.0,0.0, 0.0 )  // y-axis
 b(    0.0,    0.0,   0.0,    {bond_width*1.5},   0.0,0.0,1.0,    {bond_width},    0.0,0.0,1.0, 0.0 )  // z-axis
 """ if viewAxis else ""
}

''')

        #return
        # Write atoms
        pov.write('// ------ Atoms\n')
        for i in range(len(sys.apos)):
            e = sys.enames[i]  # Elements are 1-based
            e = e.split('_')[0]
            x, y, z = sys.apos[i]
            rad = elements.ELEMENT_DICT[e][6] * atom_scale
            clr = elements.getColor( e , bFloat=True)            
            pov.write('a( {:10.5f}, {:10.5f}, {:10.5f}, {:10.5f}, {:10.5f}, {:10.5f}, {:10.5f}, 0.0 )\n'.format( x, y, z, rad, clr[0], clr[1], clr[2] ))

        # Write bonds
        if sys.bonds is not None:
            pov.write('\n// ------ Bonds\n')
            for bond in sys.bonds:
                i, j = bond[:2]
                # Get positions
                pos1 = sys.apos[i]
                pos2 = sys.apos[j]                
                pov.write('b( {:10.5f}, {:10.5f}, {:10.5f}, {:10.5f}, {:10.5f}, {:10.5f}, {:10.5f}, {:10.5f}, {:10.5f}, {:10.5f}, {:10.5f}, 0.0 )\n'.format( 
                              pos1[0], pos1[1], pos1[2], bond_width, 
                              pos2[0], pos2[1], pos2[2], bond_width, 
                              bond_clr[0],  bond_clr[1], bond_clr[2]  ))

