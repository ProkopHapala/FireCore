

# **H-Bond fitting**


## Ideas

* Polarizability by fitting force-constants of springs attached to the electron pairs to host-atom
   * Try non-linear fit with relaxation of the e-pairs and sigma-holes  
   * Should make relaxation with harmonic approximation (hook law dX=k\*F) this will model polarizability. We can analytically derive derivative of Params with respect to distorted geometry.  
   * We simply evaluate hessian for that atom (e-pairs) using our intra-molecular forcefield (MMFF) 
      * Problem may be that equilibrium geometry of that molecule from reference (DFT) is not exactly that of our forcefield. 
      * But as long as we deflect just the epair (not atom) it should not be a problem because epair positions are not present in DFT reference anyway (it is hidden degree of freedom which we choose).  

* Fast-Potentials
  * Approximate $exp(-b r)$ by $A(1-b r/n)^n$
  * Use 1\\r8 rather then 1\\r12 in Lenard Jones (LJ is too steep)  
  * Damping of $1/r$ by attaching parabola $a r^2 + C$ for $r < r_min$
     * this is in fact apprixamation of Boys function $erf(r)/r$
  * Eventually we can use $A(1-b r/n)^n$ also for Coulomb potential with finite cutoff      

* Possible Paper names:
  * *Hydrogen bonds correction with explicit electron pairs with negligible overhead*

* H-bond corrections in GridFF 
   * Hydrogen bond correction in GridFF should be stored separately (i.e. two channels) for donors (Hb\>0) and for acceptors (Hb\<0)  

* Sample Hbonds-type in the x-points of reference E(x) profile  (WTF? What did I mean?)
* Make Monte-Carlo optimizer in Python using sampleHbondType as solver  



## Implementation on GPU

* We need to build bonding topology in python
* It is better if we store that bonding topology into some .mol file or similar format



## Practical problems

* H-bond fitting is in separate branch `fitREQH` (in `/home/prokop/git/FireCore-fitREQH`), we should perhaps merge it into `prokop` (`/home/prokop/git/FireCore`)


* there is already `FitREQ_ocl.h`, `FitREQ.cl` and also backup `FitREQ_bak.h`
* `Directional_Barrier.py` showing effect of e-pair dummy atoms on the barrier profile of interactomic interaction



Merge Conflict

```
prokop@GTX3090:~/git/FireCore$ git merge origin/fitREQH
Auto-merging cpp/common/IO_utils.h
Auto-merging cpp/common/dataStructures/LimitedGraph.h
CONFLICT (content): Merge conflict in cpp/common/dataStructures/LimitedGraph.h
Auto-merging cpp/common/math/Forces.h
Auto-merging cpp/common/math/quaternion.h
Auto-merging cpp/common/molecular/Atoms.h
CONFLICT (content): Merge conflict in cpp/common/molecular/Atoms.h
Auto-merging cpp/common/molecular/FitREQ.h
CONFLICT (content): Merge conflict in cpp/common/molecular/FitREQ.h
Auto-merging cpp/common/molecular/MMFFBuilder.h
CONFLICT (content): Merge conflict in cpp/common/molecular/MMFFBuilder.h
Auto-merging cpp/common/molecular/MMFFparams.h
CONFLICT (distinct types): cpp/common_resources/AtomTypes.dat had different types on each side; renamed one of them so each can be recorded somewhere.
CONFLICT (distinct types): cpp/common_resources/BondTypes.dat had different types on each side; renamed one of them so each can be recorded somewhere.
CONFLICT (distinct types): cpp/common_resources/ElementTypes.dat had different types on each side; renamed one of them so each can be recorded somewhere.
Auto-merging tests/tFitREQ/run.sh
CONFLICT (content): Merge conflict in tests/tFitREQ/run.sh
Automatic merge failed; fix conflicts and then commit the result.
prokop@GTX3090:~/git/FireCore$ 
```