# NBFF.h

## Includes

- `fastmath.h`
- `Vec3.h`
- `quaternion.h`
- `Atoms.h`
- `Buckets.h`
- `Forces.h`
- `ForceField.h`
- `simd.h`


## Free functions

- `fitAABB` - 
- `makePBCshifts_` - 
- `evalPointCoulPBC` - 
- `sampleCoulombPBC` - 


## Types (classes and structs)

### class `NBFF`

[Short description what is the purpose of the class, what its role in the bigger context]

**Inheritance**

- ForceField

#### properties

- `nBBs`:`int` - 
- `BBs`:`Vec6d*` - 
- `pointBBs`:`Buckets` - 
- `drSR`:`double` - 
- `ampSR`:`double` - 
- `alphaMorse`:`double` - 
- `Rdamp`:`double` - 
- `bPBC`:`bool` - 
- `npbc`:`int` - 


#### methods

- `torq` - 
- `bindShifts` - 
- `makePBCshifts` - 
- `evalPLQs` - 
- `makePLQs` - 
- `evalPLQd` - 
- `makePLQd` - 
- `updatePointBBs` - 
- `selectInBox` - 
- `evalSortRange_BBs` - 
- `evalLJQs` - 
- `evalLJQs_ng4_atom` - 
- `evalLJQs_ng4_omp` - 
- `evalLJQs_ng4` - 
- `addMorseQH_PBC_omp` - 
- `getMorseQH_PBC_omp` - 
- `getLJQs_PBC_omp` - 
- `evalLJQs_PBC_atom_omp` - 
- `evalLJQs_PBC_simd` - 
- `evalLJQs_ng4_PBC_atom_omp` - 
- `evalLJQs_ng4_PBC_simd` - 
- `evalLJQs_atom_omp` - 
- `evalLJQs_simd` - 
- `evalLJQs_ng4_atom_omp` - 
- `evalLJQs_ng4_simd` - 
- `evalLJQs_atom_avx` - 
- `evalCollisionDamp_atom_omp` - 
- `evalCollisionDamp_omp` - 
- `evalLJQs_ng4_PBC_atom` - 
- `evalLJQs_ng4_PBC_omp` - 
- `evalLJQs_ng4_PBC` - 
- `evalLJQ` - 
- `evalMorse` - 
- `evalMorsePLQ` - 
- `evalR` - 
- `makePBCshifts` - 
- `print_nonbonded` - 
- `bindOrRealloc` - 
- `dealloc` - 

