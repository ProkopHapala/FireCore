import numpy             as np
import matplotlib.pyplot as plt


'''
def dampFunc_old( vf, vv,ff, damping=0.1 ):
    mask=(vf<0.)
    cf = damping * np.sqrt( vv/ff )
    cv = 1.-damping + 0.*vf
    cv[mask]=0.
    cf[mask]=0.
    return cv,cf

def dampFunc_new( vf, vv,ff, damping=0.1, c0=-0.3 ):
    c  = vf/np.sqrt( vv*ff )
    mask=(c<c0)
    cd = damping + (1-c)/(1-c0)*(1-damping)
    cf = cd * np.sqrt( vv/ff )
    cv = 1.-cd
    cv[mask]=0.
    cf[mask]=0.
    return cv,cf #,c

def dampFunc_c3( vf, vv,ff, damping=0.1, c0=0.0 ):
    c  = vf/np.sqrt( vv*ff )
    mask=(c<c0)
    c_ = ( ( np.abs(c)**0.25)*np.sign(c)  + 1.0)*0.5
    cd = damping + (1.-c_)*(1-damping)
    cf = cd * np.sqrt( vv/ff )
    cv = 1.-cd
    #cv[mask]=0.
    #cf[mask]=0.
    return cv,cf #,c

def plot_with_damp( dampf=None, label0="", ls="-" ):
    cv,cf = dampf( vf,vv,ff, damping=0.1 )
    #cv,cf = dampFunc_new( vf,vv,ff, damping=0.1 )
    vx_ = vx*cv + fx*cf
    vy_ = vy*cv + fy*cf
    plt.plot( phi, cv    ,ls+'k',label=label0+'cv')
    #plt.plot( phi, vx_/vr ,ls+'b',label=label0+'v_||')
    #plt.plot( phi, vy_/vr ,ls+'r',label=label0+'v_T')
    #if(c is not None): plt.plot( phi, c      ,'grey',label='cos(v,f)')

phi = np.linspace(-np.pi,np.pi,100)
vr = 0.66
vx = vr
vy = 0.0

fr = 1.8
fx = np.sin(phi)
fy = np.cos(phi)

vv= vx*vx + vy*vy   + phi*0.
ff= fx*fx + fy*fy
vf= vx*fx + vy*fy


plt.title("FIRE velocity damping")
plt.xlabel("angle between velocity and force [rad]")

plot_with_damp( dampf=dampFunc_old, label0="old ", ls="--"   )
#plot_with_damp( dampf=dampFunc_new, label0="new ", ls="-"  )
plot_with_damp( dampf=dampFunc_c3, label0="c3 ", ls="-"  )

c  = vf/np.sqrt( vv*ff )
plt.plot( phi, c      ,'grey',label='cos(v,f)')
'''

phi = np.linspace(-np.pi,np.pi,1000)
s   = np.sin(phi)
c   = np.cos(phi)
s2 = 1-c*c

mask=(c<0)
cv       = 1-0.5*(c      -1)**4
cv[mask] =   0.5*(c[mask]+1)**4
# cf = d(dv)/dphi
#cf       =   (c      -1)**3*s
#cf[mask] =  -(c[mask]+1)**3*s[mask]
cf       =   s2**8

#f1 = c + (c**3)*-0.35 + (c**5)*-0.1
#f1 = (c + c2*0.3)/(1+0.3)
#f1 = (c + c4*-0.1)/(1-0.1)

cd = 0.1

plt.title("FIRE velocity damping")
plt.plot(phi,c,label='cos(x)')
#plt.plot(phi,s,label='sin(x)')
#plt.plot(phi,c2,label='cos(2x)')
#plt.plot(phi,c4,label='cos(4x)')
plt.plot(phi,cv*(1-cd),label='cv (damp)')
plt.plot(phi,cf*(1-cd) + cd*cv,label='cf (turn)')
#plt.plot(c**5,label='c**5')

plt.legend()
plt.grid()

plt.show()