
from Makefile_machines import *

OBJECTS_COM = [ 
    'DFTD3', 'INITMPI', 'ORDERN', 'ALLOCATIONS', 'ASSEMBLERS', 
    'DASSEMBLERS', 'INITIALIZERS', 'INTERACTIONS', 'INTERPOLATERS',
    'LOOPS', 'MD',  'NEIGHBORS', 'PRESSURE', 'READFILES', 
    'ROTATIONS', 'SOLVESH_DIAG', 'FORM_RHO', 'UMBRELLA','UTIL', 
    'UTIL_SPARSE', 'VISUALIZATION', 'XC', 'CG', 'DOS', 'THERMOINT', 
    'NEB', 'TRANS', 'GRID', 'TDSE', 'BIAS', 'NAC',
]

'''
$(DFTD3) $(INITMPI) $(ORDERN) $(ALLOCATIONS) $(ASSEMBLERS) \
        $(DASSEMBLERS) $(INITIALIZERS) $(INTERACTIONS) $(INTERPOLATERS) \
        $(LOOPS) $(MD) $(NEIGHBORS) $(PRESSURE) $(READFILES) \
        $(ROTATIONS) $(SOLVESH_DIAG) $(FORM_RHO) $(UMBRELLA) $(UTIL) \
        $(UTIL_SPARSE) $(VISUALIZATION) $(XC) $(CG) $(DOS) $(THERMOINT) \
        $(NEB) $(TRANS) $(GRID) $(TDSE) $(BIAS) $(NAC)
'''

OBJECTS              = [ 'MODULES' ] + OBJECTS_COM + [ 'MAIN' ]
OBJECTS_QMMM         = [ 'MODULES' ] + OBJECTS_COM + [ 'QMMM' ]
OBJECTS_SERVER       = [ 'MODULES' ] + OBJECTS_COM + [ 'MAIN_SERVER' ]
OBJECTS_SERVER_AMBER = [ 'MODULES' ] + OBJECTS_COM + [ 'MAIN_SERVER_AMBER' ]

_gobals_ = globals()

#import inspect

default_obj_name_exclude = set('o')

def retrieve_name(var):
    ret = [key for key,val in _gobals_.iteritems() if (val is var) ]
    if len(ret)>0: 
        return ret[0]
    else:
        raise Exception( "Module "+__name__+" does not contain global variable with value ", var )
        return  None
    
'''
def writeTarget( fout, name, depend, objs, compiler="$(F90)", fflags="$(FFLAGS)", lflags="$(LFLAGS)" ):
    fout.write( name +" : "+depend+"\n" )
    fout.write( "\t"+compiler+" -o "+name+" " )
    fout.write( " "+fflags )
    for o in objs:
        oname = retrieve_name( o )
        fout.write( " $(%s)" %oname )
    fout.write( " "+lflags )
    fout.write( "\n\n" )
'''

# ==================================================

inline_targets = {

".PHONY" :" clean veryclean extraclean",

"clean" : '''
	rm -f -r core *.o .nfs* rii_files fireball.x.ip*  *.mod ldtmp* *.vo *~ *.il''',

"veryclean": ''' clean
	rm -f fireball.x libfireball.a''',

"extraclean": "veryclean",

"all":'''
	make fireball.x''',

"libfireball": '''$(OBJECTS_QMMM)
	ar rv libfireball.a $(OBJECTS_QMMM)
	ranlib libfireball.a''',
}
