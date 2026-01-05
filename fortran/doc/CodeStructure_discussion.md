
## Lets define minimal necessary core:

I think the ***services*** around like `READFILES`, `LOOPS`,  `SOLVESH`, `FORM_RHO`,  `INITIALIZERS`, `ALLOCATIONS`, `NEIGHBORS`, `MODULES`  etc. are secondary … most of this can be better replaced by some functions from numpy/scipy/ASE  … the rest can be easily adjusted after we define the ***core***. 

What I mean by the ***core*** is the minimal set of functions you need to **assemble Fireball Hamiltonian** (because this is the part special for Fireball, which cannot be found in some general numerical package). (And then corresponding forces, since I want relax molecular geometry).  And there I’m a bit lost, because I don’t know what each of these assemblers does, and which are needed always and which are supplementary ( also it depends on settings like `idogs` and `iharris` `mcWEDA` `iTheory` etc.  …. we should choose just one for the start (can add other options later) )

So if somebody can say which of the files are ***essential***  and which are ***optional/supplementary***?


- `ASSEMBLERS`  / `assemble_*.f90`  …  These are for Hamiltonian
- `DASSEMBLERS` /`Dassemble_*.f90` … These are for derivative of Hamiltonian => Forces … (Right?). Since each `Dassemble_*.f90` is connected to particular `assemble_*.f90` we don’t have to discuss them separately
- `INTERACTIONS` … These are functions called from `ASSEMBLERS` and `DASSEMBLERS` for each matrix element? Right? Are they all essential ?
- `ROTATIONS`    … This is probably all essential (?)  … actually this is one thing which I would like to reuse for many other things (perhaps rewrite it to C/OpenCL)
- `INTERPOLATERS` … perhaps need them all (?) to read form `Fdata` tables.
# There are shortcuts, what they mean:

(also if you can point out if they are general, or special for some theory-setting)
_1c, _2c, _3c  … 1 centre, 2 center, 3 center
_ca_2c, _ca_3c  …. ?
_vdip, _dip   … something about dipole (?)
_eh …. Hartree energy?
_usr    …. ?
_F      … ?
_Ir     … ?
mcweda … mcWeda
KS … KhonSham
_S … overlap matrix
_olsxc,  _snxc, _xczw    …. different kinds of exchange correlation functionals ?
_zw  …. ?
_off, _on  … offsite, onsite (?)
_VNA … neutral atom potential
_VXC  … exchange correlation potential
_VNL  … non-local potential?

# PLEASE SELECT WHICH CAN BE REMOVED:

There ifs should be removed  `getforces_mcweda.f90`
```
          if (itheory_xc .eq. 1) then
           !write (*,*) ' Dassemble on-site SNXC force contributions. '
           if (itheory .eq. 1) then
             call Dassemble_ca_snxc_on (nprocs, iordern)
             call Dassemble_ca_snxc_2c (nprocs, iordern)
           else
             call Dassemble_snxc_on (nprocs, iordern)
             call Dassemble_snxc_2c (nprocs, iordern)
           endif
          end if
          if (itheory_xc .eq. 2) then
           !write (*,*) ' Dassemble on-site OSLXC force contributions. '
           if (itheory .eq. 1) then
             call Dassemble_ca_olsxc_on (nprocs, iordern)
             call Dassemble_ca_olsxc_2c (nprocs, iordern)
           else
             call Dassemble_olsxc_on (nprocs, iordern)
             call Dassemble_olsxc_2c (nprocs, iordern)
           endif
          end if
```