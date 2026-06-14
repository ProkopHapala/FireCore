### Charge conservation $Q_{ep} + Q_{host} = Q_{resp}$

What was exactly the problem with Host atom of free electron pair ?

Can we **remove**:

```C++
    if(aisep[i]==1){dE_dri+=corr_elec(ir,ir2,Q,d,adirs,i,nep,bs,nj,jpos,j,ps,dE_dQ);}  // if i is in elecron pair subtract the charge transfer
    if(jisep[j]==1){dE_drj-=corr_elec(ir,ir2,Q,d,jdirs,j,nep,bs,0, ps,i,jpos,dE_dQ);}  // if j is in elecron pair subtract the charge transfer    
```

```C++
double corr_elec( double ir, double ir2, double Q, Vec3d d, Vec3d* dirs, int i, int nep, Vec2i* bs, int nj, Vec3d* pos, int j, Vec3d* ps, double &dE_dQ){
    double dE_dr = ir * ir2 * COULOMB_CONST * Q * d.dot(dirs[i]);
    for(int k=0; k<nep; k++){
        if(bs[k].y==nj+i){
            Vec3d dd = pos[j] - ps[bs[k].x-nj];
            dE_dQ -= COULOMB_CONST / sqrt( dd.norm2() );
            break;
        }
    }
    return dE_dr;
}
```

and instead **add**: 

```C++
void acumHostDerivs( int nepair, Vec2i* epairAndHostIndex, int* types, double dEw  ){
    for(int i=0; i<nepair; i++){
        Vec2i ab         = epairAndHostIndex[i];
        int t            = types[ab.i];  // map atom index i to atom type t
        const Quat4i& tt = typToREQ[t];  // get index of degrees of freedom for atom type t
        const Quat4d&  f = fs[ab.j];     // get variation from the host atom
        if(tt.z>=0)fDOFs[tt.z]  -= f.z*dEw;  // we subtract the variational derivative of the host atom charge because the charge is transfered from host to the electron pair
    }
}

// in evalDerivsSamp() ... 
    acumDerivs    ( atoms->natoms,     jtyp, dEw );
    AddedData* data = (AddedData*)atoms->userData;
    if(data && data->nep > 0){
        acumHostDerivs( data->nep, data->epairAndHostIndex, jtyp, dEw );
    }
```

### Variation of H1 depends on Q ?

Why we need to care about REQ.z = Q = Charge when calculating variation wrt H1 correction ?

```C++
            if(aisep[i]==1){
                fsi.z = dE_dQ * jq[j];                              // dEtot/dQi
                fsi.w = dE_dri;                                     // dEtot/dri
            }else{
                fsi.z = dE_dH1 * H1 / (REQi.z+sign(REQi.z)*1e-300); // dEtot/dH1i
            }
```