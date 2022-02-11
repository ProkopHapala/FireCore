#!/usr/bin/python

import sys
import os

CC_KINDS={
'F90': ( '', '.f90', '$(FFLAGS)',   '$(INCLUDES)', '$(F90)' ),
'F77': ( '', '.f',   '$(FFLAGS77)', '$(INCLUDES)', '$(F77)' ),
'CC' : ( '', '.c',   '$(FFLAGSC)',  '$(INCLUDES)', '$(CC)'  ),
}

def insertGroup( groups, group_names, name, variant, objs ):
    #groups.append( (group_name,variant,objs) )
    if name in groups:
        groups[name][variant] = objs
    else:
        groups[name] = {variant:objs}
        group_names.append(name)

def fromMakeFile( fname  ):
    fin = open(fname,"r")
    group_names = []
    var_names   = [ "" ]
    var_set     = set()
    group_dict  = {}
    group_name  = None
    objs = []
    variant = ""
    key_if    = "findstring"
    key_if_nc = len(key_if)
    for line in fin:
        ws = line.split()
        nw = len(ws)
        #if   (nw>=1):  print ws[0]
        if   (nw>=1) and (ws[0] == 'ifneq' ):
            i0 = line     .find( key_if )+key_if_nc
            i1 = line[i0:].find( "," )
            variant = line[i0:i0+i1].strip()
            if not variant in var_set: 
                var_names.append(variant)
                var_set  .add(variant)
            #print variant, i0, i1
        elif (nw>=1) and (ws[0] == 'endif' ):
            variant=""
        elif (nw>=2) and (ws[1] == '='):
            if group_name is not None:
                insertGroup( group_dict, group_names, group_name, variant, objs )
            objs = []
            group_name = ws[0]
        objs += [ w[:-2] for w in ws if w[-2:]=='.o' ]
    if group_name is not None:
        insertGroup( group_dict, group_names, group_name, variant, objs )
    return group_dict, group_names, var_names

def diffTwoSets(s1,s2):
    s12 = s2 & s1
    return s12, s1-s12, s2-s12

def diffListSet(lst,st):
    common = [ o for o in lst if o     in st ]
    diff   = [ o for o in lst if o not in st ]
    return diff, common

def diffTwoLists( lst1, lst2 ):
    diff1,common = diffListSet( lst1, set(lst2) )
    diff2,common = diffListSet( lst2, set(lst1) )
    return common, diff1, diff2
    #s12,d1,d2 = diffTwoSets( set(lst1), set(lst2) ) 
    #return list(s12),list(d1),list(d2)

def pruneVariants(groups):
    '''
    remove objects files form specialized variants which are already in default variant
    '''
    for name,vardict in groups:
        if "" in vardict:
            default_objs = set(vardict[""])
            for key in vardict:
                if key != "":
                    vardict[key] , _ = diffListSet( vardict[key], default_objs )

def flatten_groups( group_dict, variants ):
    group_dict_flat = {}
    #print " flatten_groups : group_dict ", group_dict
    for k,vs in group_dict.iteritems():
        group_dict_flat[k] = { o for var in variants for o in vs.get(var,[]) }
    return group_dict_flat

def groups2files( group_dict, variants, pre="", post="" ):
    group_dict_flat = {}
    for k,vs in group_dict.iteritems():
        group_dict_flat[k] = { pre+o+post for var in variants for o in vs.get(var,[]) }
    return group_dict_flat

def groups2files_cc( group_dict, variants, pre="", cc_kinds=CC_KINDS, inv_SPECIAL_CC={}, default='F90' ):
    group_dict_flat = {}
    for k,vs in group_dict.iteritems():
        group_dict_flat[k] = { pre+o+CC_KINDS[inv_SPECIAL_CC.get(o,default)][1] for var in variants for o in vs.get(var,[]) }
    return group_dict_flat

def extractVariants( group_dict ):
    return { v for k,vs in group_dict.iteritems() for v in vs }

def toMakefile_list( fout, name, lst, pre=" ", post=".o", nclmax=120 ):
    s   = name+" = "
    ncl = len ( s )
    fout.write( s )
    no  = len(lst)
    for i in range(no):
        o = lst[i]
        s = pre+o+post
        ncl += len(s)
        if ( i<(no-1) ):
            if ncl > nclmax:
                s += " \\\n\t"
                ncl = 0
        fout.write( s )
    fout.write( "\n" )

def select_obj_groups( group_dict, group_names, var_name_set ):
    return [ name for name in group_names if any( v in var_name_set for v in group_dict[name] ) ]

def toMakefile_obj_groups( fout, group_dict, group_names, var_names, build_path="", nclmax=120 ):
    ncl = 0
    for name in group_names:
        varinants = group_dict[name]
        if len(varinants) == 1:
            #toMakefile_obj_group( fout, name, varinants[ varinants.keys()[0] ],    nclmax=nclmax )
            for var in var_names:
                if var in varinants:
                    toMakefile_list( fout, name, varinants[var],    nclmax=nclmax )
                    fout.write( "\n\n\n" )
                    break
        else:
            for var in var_names:
                if var in varinants:
                    #toMakefile_obj_group( fout, name+"_"+var, varinants[var],    nclmax=nclmax )
                    toMakefile_list( fout, name+"_"+var, varinants[var],    nclmax=nclmax )
            fout.write( name+" = " )
            for var in var_names:
                if var in varinants:
                    fout.write( " $(%s_%s)" %(name,var) )
            fout.write( "\n\n\n" )

def toMakefile_cc_obj( fout, o, pre="../", post=".f90", fflags="$(FFLAGS)", Iflags="$(INCLUDES)", compiler="$(F90)" ):
    src = pre+o+post
    obj = o+".o"
    fout.write( obj+" : "+src+"\n" )
    fout.write( "\t" + compiler + " " + fflags +" "+ Iflags + " -c "+src + "\n" )

def toMakefile_cc_kind( fout, o, pre="", cc_kinds=CC_KINDS, inv_SPECIAL_CC={}, default_cc='F90' ):
    cc_name = inv_SPECIAL_CC.get( o, default_cc)
    cc      = cc_kinds[cc_name] 
    toMakefile_cc_obj( fout, o, pre=pre, post=cc[1], fflags=cc[2], Iflags=cc[3], compiler=cc[4] )

#def toMakefile_cc_objs( fout, group_dict, group_names, var_names, special_cc={}, src_path="", ccomment='*', ncomment=50 ):
def toMakefile_cc_objs( fout, group_dict, group_names, var_names, cc_kinds=CC_KINDS, inv_SPECIAL_CC={}, default_cc="F90", src_path="", ccomment='*', ncomment=50 ):
    for name in group_names:
        varinants = group_dict[name]
        fout.write("\n\n")
        cmline = "#"   +(ccomment*ncomment)+"\n"
        fout.write( cmline )
        fout.write( "#   " + name +"\n" )
        fout.write( cmline )
        for var in var_names:
            if var in varinants:
                fout.write("#====== variant : '" + var +"'\n" )
                for o in varinants[var]:
                    #Iflags=""
                    #if o in special_cc:
                    #    Iflags=" ".join(special_cc[o])
                    #toMakefile_cc_obj( fout, o, pre=src_path+name+"/", Iflags=Iflags )
                    #cc_name = inv_SPECIAL_CC.get( o, default_cc)
                    #cc = cc_kinds[cc_name] 
                    #toMakefile_cc_obj( fout, o, pre=src_path+name+"/" )
                    #if o=='lanc':
                    #    print "!!!!!!! lanc !!!!!!! ", cc_name, cc
                    #toMakefile_cc_obj( fout, o, pre=src_path+name+"/", post=cc[1], fflags=cc[2], Iflags=cc[3], compiler=cc[4] )
                    toMakefile_cc_kind( fout, o, pre=src_path+name+"/", cc_kinds=cc_kinds, inv_SPECIAL_CC=inv_SPECIAL_CC, default_cc=default_cc )
        fout.write("\n\n")

def listToPython(fout, name, lst, pre="'", mid="' : [", post="],\n", nclmax=120 ):
    s = pre+ name +mid
    ncl = len(s)
    fout.write( s )
    no = len(lst)
    for i in range(no):
        o = lst[i]
        s = "'"+o+"'"
        ncl += len(s)
        if ( i<no-1 ):
            s += ","
            if ncl > nclmax:
                s += "\n        "
                ncl = 0
        fout.write( s )
    fout.write( post )
    ncl = 0

def toPython_groups( fout, group_dict, group_names, var_names, nclmax=120 ):
    ncl = 0
    #fout.write( "'%s' : {\n" %name )
    listToPython(fout, "group_names",     group_names, pre="", mid=" = [", post="]\n\n", nclmax=nclmax )
    listToPython(fout, "variant_names",   var_names  , pre="", mid=" = [", post="]\n\n", nclmax=nclmax )
    fout.write( "GROUPS = {\n" )
    for name in group_names:
        varinants = group_dict[name]
        fout.write( "'%s' : {\n" %name )
        for var in var_names:
            if var in varinants:
                #fout.write("# ====== variant : " + var )
                listToPython(fout, var, varinants[var] )
        fout.write("}, #END %s\n\n"  %(name)   )
    fout.write("} #END GROUPS\n\n")

def toMakefile_cc_obj_c(fout, o, pre="../MODULES/", post=".c", fflags="$(CFLAGS)", compiler="$(CC)" ):
    toMakefile_cc_obj( fout, o, pre=pre, post=post, fflags=fflags, compiler=compiler )

def toMakefile_list_vars( fout, name, lst, pre=" $(", post=") " ):
    toMakefile_list( fout, name, lst,  pre=pre, post=post )
    fout.write( "\n" )

def toMakefile_name( fout, name, val ):
    fout.write( name + " = " + val + "\n\n" )

def toMakefile_tar_inline_target( fout, body ):
    fout.write( body+"\n\n" )

def toMakefile_target( fout, name, depend, objs, compiler="$(F90)", fflags="$(FFLAGS)", lflags="$(LFLAGS)" ):
    fout.write( name +" : "+depend+"\n" )
    fout.write( "\t"+compiler+" -o "+name+" " )
    fout.write( " "+fflags )
    #for o in objs:
    #    oname = retrieve_name( o )
    #    fout.write( " $(%s)" %oname )
    fout.write( " "+objs    )
    fout.write( " "+lflags )
    fout.write( "\n\n" )

def checkFiles( path, names,  bRedudant=True, bMissing=True ):
    ls  = os.listdir(path)
    #ls_    = set(ls)
    #names_ = set(names)
    shared, dl, dn = diffTwoLists( ls, names )
    if ( bMissing  and (len(dn)>0) ): print path," : missing  : " , dn #, "                checkFiles "
    if ( bRedudant and (len(dl)>0) ): print path," : redudant : " , dl #, "                checkFiles "
    return shared, dl, dn

def checkFilesInPaths( path_names, name_dct, pre="", bRedudant=True, bMissing=True ):
    #print "checkFilesInPaths: path_names ", path_names
    #print "checkFilesInPaths: name_dct ", name_dct
    for name in path_names:
        path = pre+name
        #print "-------------- ", path    #,"     (in checkFilesInPaths)"
        checkFiles( path, name_dct[name],  bRedudant=bRedudant, bMissing=bMissing )

def fomatLinkLog( fname_in, fname_out ):
    with open(fname_in ,"r") as fin: 
        words = [ w.strip() for l in fin for w in l.split() ]
    with open(fname_out,"w") as fout:
        for w in words:
            fout.write(w+"\n")

if __name__ == "__main__":

    '''
    groups, var_names = fromMakeFile( "Makefile.in" )
    pruneVariants(groups)

    # groups_dct = dict(groups)
    # common, diff1, diff2 = diffTwoLists( groups_dct['INTERACTIONS']['SCALAPACK'], groups_dct['INTERACTIONS']['ORDERN'] )
    # print "common ", common
    # print "diff1 ", diff1
    # print "diff2 ", diff2
    # exit(0)

    with open("Makefile",'w') as fout:
        toMakefile_obj_groups( fout, groups, var_names )
        toMakefile_cc_objs   ( fout, groups, var_names )
    with open("Makefile_objects-New.py",'w') as fout:
        toPython_groups( fout, groups, var_names )

    fomatLinkLog( "link_ref.log", "link_ref.log" )
    fomatLinkLog( "link.log", "link.log" )
    '''

    import Makefile_targets
    import Makefile_machines
    from Makefile_targets  import *
    from Makefile_objects  import *
    from Makefile_machines import *

    #import inspect
    #print [item for item in dir(Makefile_targets) if not item.startswith("__")]
    #print dict(inspect.getmembers( Makefile_targets ))["f_globals"]
    #print dir(Makefile_targets)
    #print Makefile_targets["OBJECTS"]
    #print "====================================="
    #for g in Makefile_targets._gobals_:
    #    print "\n ---------: "+g+" : \n", Makefile_targets._gobals_[g]
    #print  Makefile_targets._gobals_
    #exit()

    variants_new = extractVariants( GROUPS ); print " variants_new ", variants_new 

    # see Understanding nested list comprehension syntax in Python : https://spapas.github.io/2016/04/27/python-nested-list-comprehensions/
    inv_SPECIAL_CC = { v:k for k,vs in SPECIAL_CC.iteritems() for v in vs }
    print " inv_SPECIAL_CC : ",  inv_SPECIAL_CC

    #build_path = "build/"
    #src_path_make  = "../SRC/"
    #src_path_check = "SRC/"
    src_path_make  = "../fortran/"
    src_path_check = "fortran/"
    #src_path = "SRC/"
    #MKL_PATH = "/home/prokop/SW/intel"
    #MKL_PATH = "/home/prokop/intel"
    MKL_PATH = "/home/prokop/SW/intel"
    MPI_PATH = "/usr/lib/x86_64-linux-gnu/openmpi"
    #FFLAGS, LFLAGS_, LPATHS = genFlags( ["OPT"], MKL_PATH=MKL_PATH, MPI_PATH=MPI_PATH )
    LFLAGS_, LPATHS   = genLFLAGS( MKL_PATH=MKL_PATH, MPI_PATH=MPI_PATH )
    print "LFLAGS_ ", LFLAGS_
    print "LPATH   ", LPATHS

    _FFLAGS = Makefile_machines._FFLAGS 
    #mode_opt = 'OPT'
    mode_opt = 'DEBUG'
    #mode_opt = 'DEBUGw'
    if len(sys.argv)>1: mode_opt=sys.argv[1]

    #FFLAGS = _FFLAGS[mode_opt] + _FFLAGS['F90']
    #FFLAGS = _FFLAGS[mode_opt] + _FFLAGS['F90']
    #FFLAGS = _FFLAGS[mode_opt] + _FFLAGS['F90']

    INCLUDES = "-I/usr/include/mpich/"

    #OPTIONAL_MODULES = [ 'NAC','QMMM','DFTD3','BIAS','TRANS','TDSE' ]

    my_variant_names = ['','DOUBLE','PROGRAM'] + all_optional_modules
    my_variant_names_set = set(my_variant_names)
    #variant_names_ = ['','DOUBLE','GAMMA']
    #variant_names_ = ['','DOUBLE','MPI-k','QMMM','ORDERN']
    #variant_names_ = ['','DOUBLE','MPI-k','QMMM','ORDERN']

    my_group_names = select_obj_groups( GROUPS, all_group_names, my_variant_names_set )
    print " === my_group_names ", my_group_names
    #exit(-1)

    #groups_flat = flatten_groups( GROUPS, variant_names )
    #groups_flat = groups2files( GROUPS, variant_names, post=".f90" )
    groups_flat = groups2files_cc( GROUPS, all_variant_names, cc_kinds=CC_KINDS, inv_SPECIAL_CC=inv_SPECIAL_CC )

    print "======= CHECK GROUP FOLDERs "; checkFiles       ( src_path_check,  all_group_names,                  bRedudant=True, bMissing=True )
    print "======= CHECK GROUP FILEs   "; checkFilesInPaths( all_group_names, groups_flat, pre=src_path_check, bRedudant=True, bMissing=True )
    
    #print "======= CHECK GROUP FOLDERs "; checkFiles       ( src_path_check, all_group_names,                  bRedudant=False, bMissing=True )
    #print "======= CHECK GROUP FILEs   "; checkFilesInPaths( all_group_names, groups_flat, pre=src_path_check, bRedudant=False, bMissing=True )

    with open("Makefile",'w') as fout:
        toMakefile_obj_groups( fout, GROUPS, all_group_names, my_variant_names )
        #toMakefile_list_vars( fout, "OBJECTS_SERVER", OBJECTS_SERVER, )
        #toMakefile_list_vars( fout, "OBJECTS",        OBJECTS,        )
        toMakefile_list_vars( fout, "OBJECTS", my_group_names )
        #fout.write( "(OBJ)/%.o" + "\n\n" )
        fout.write( "F90 = gfortran\n" )
        #toMakefile_name( fout, "FFLAGS",   FFLAGS   )
        #toMakefile_name( fout, "FFLAGS77", FFLAGS   )
        #toMakefile_name( fout, "LFLAGS_",  LFLAGS_  )
        toMakefile_name( fout, "FFLAGS",   "  -fPIC  "+ _FFLAGS[mode_opt] + _FFLAGS['F90']  )
        toMakefile_name( fout, "FFLAGS77", "  -fPIC  "+ _FFLAGS[mode_opt] + _FFLAGS['F77']  )
        toMakefile_name( fout, "FFLAGSC",  "  -fPIC  "+ _FFLAGS[mode_opt] + _FFLAGS['CC' ]  )
        toMakefile_name( fout, "LFLAGS_",  LFLAGS_  )
        toMakefile_name( fout, "LPATHS",   LPATHS   )
        toMakefile_name( fout, "INCLUDES", INCLUDES )
        #toMakefile_list_vars( fout, "LFLAGS", ["LPATHS","LFLAGS_"] )
        #LFLAGS = " -L/home/prokop/intel/mkl/lib/intel64 -lmkl_scalapack_lp64 -lmkl_blacs_openmpi_lp64 -lmkl_lapack95_lp64 -lmkl_intel_lp64 -lmkl_intel_thread -lmkl_core -liomp5 -lpthread -L/home/prokop/intel/mkl/intel64 -I/home/prokop/intel/mkl/include/fftw -lfftw3xf_gnu -lm -Bdynamic -L/usr/lib/x86_64-linux-gnu/openmpi/lib -I/usr/lib/x86_64-linux-gnu/openmpi/include -lmpi "
        LFLAGS = " -L"+MKL_PATH+"/mkl/lib/intel64 -lmkl_scalapack_lp64 -lmkl_blacs_openmpi_lp64 -lmkl_lapack95_lp64 -lmkl_intel_lp64 -lmkl_intel_thread -lmkl_core -liomp5 -lpthread -L/home/prokop/intel/mkl/intel64 -I/home/prokop/intel/mkl/include/fftw -lfftw3xf_gnu -lm -Bdynamic -L/usr/lib/x86_64-linux-gnu/openmpi/lib -I/usr/lib/x86_64-linux-gnu/openmpi/include -lmpi "
        
	toMakefile_name( fout, "LFLAGS", LFLAGS )
        #toMakefile_name( fout, "FFLAGS", FFLAGS )
        #writeTarget( fout, "fireball.x"       , "$(OBJECTS)", [OBJECTS       ] )
        #writeTarget( fout, "fireball_server.x", "$(OBJECTS)", [OBJECTS_SERVER] )
        for body in inline_targets :
            #writeInlineTarget( fout, key, body )
            toMakefile_tar_inline_target( fout, body )
        toMakefile_target( fout, "fireball.x"    , "$(OBJECTS)", "$(OBJECTS)"        )
        #toMakefile_target( fout, "libfireball.so", "$(OBJECTS)", "$(OBJECTS)", compiler="$(F90)", fflags=" -shared -fPIC $(FFLAGS)", lflags="$(LFLAGS)" )
        toMakefile_target( fout, "libFireCore.so", "$(OBJECTS)", "$(OBJECTS)", compiler="$(F90)", fflags=" -shared -fPIC $(FFLAGS)", lflags="$(LFLAGS)" )
        #writeTarget( fout, "fireball_server.x", "$(OBJECTS)", "OBJECTS_SERVER" )
        #toMakefile_cc_objs   ( fout, GROUPS, group_names, variant_names, special_cc=SPECIAL_CC, src_path=src_path_make )
        toMakefile_cc_objs   ( fout, GROUPS, all_group_names, my_variant_names, src_path=src_path_make, cc_kinds=CC_KINDS, inv_SPECIAL_CC=inv_SPECIAL_CC )
    
    print " ==== gen_makefile.py DONE ==== "




