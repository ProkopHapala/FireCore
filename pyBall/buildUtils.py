#!/usr/bin/python

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
    groups      = {}
    group_name = None
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
                insertGroup( groups, group_names, group_name, variant, objs )
            objs = []
            group_name = ws[0]
        objs += [ w[:-2] for w in ws if w[-2:]=='.o' ]
    if group_name is not None:
        insertGroup( groups, group_names, group_name, variant, objs )
    #groups_dict = dict(groups)
    print " variants ", var_names
    print " groups ", group_names
    return [ (key,groups[key]) for key in group_names ], var_names

def pruneVariants(groups):
    '''
    remove objects files form specialized variants which are already in default variant
    '''
    for name,vardict in groups:
        if "" in vardict:
            default_objs = set(vardict[""])
            for key in vardict:
                if key != "":
                    objs = []
                    for o in vardict[key]:
                        if not o in default_objs:
                            objs.append(o)
                    vardict[key] = objs


def toMakefile_obj_group( fout, name, objs , nclmax=120 ):
    s   = name+" = "
    ncl = len ( s )
    fout.write( s )
    no  = len(objs)
    for i in range(no):
        o = objs[i]
        s = o+".o "
        ncl += len(s)
        #print i,ncl,"   " , no, nclmax
        if ( i<(no-1) ):
            #print " i<no "
            if ncl > nclmax:
                #print " ncl > nclmax "
                s += "\\\n\t"
                ncl = 0
        fout.write( s )
    fout.write("\n")

def toMakefile_obj_groups( fout, groups,var_names, nclmax=120 ):
    ncl = 0
    for name,varinants in groups:
        if len(varinants) == 1:
            toMakefile_obj_group( fout, name, varinants[ varinants.keys()[0] ],    nclmax=nclmax )
        else:
            for var in var_names:
                if var in varinants:
                    toMakefile_obj_group( fout, name+"_"+var, varinants[var],    nclmax=nclmax )
            fout.write( name+" = " )
            for var in var_names:
                if var in varinants:
                    fout.write( "$(%s_%s) " %(name,var) )
        fout.write( "\n\n\n" )

# dimensions.o : MODULES/dimensions.f90
#	$(F90) $(FFLAGS) -c MODULES/dimensions.f90
def toMakefile_cc_objs( fout, groups,var_names, ccomment='*', ncomment=50 ):
    ncl = 0
    for name,varinants in groups:
        fout.write("\n\n")
        cmline = "#"   +(ccomment*ncomment)+"\n"
        fout.write( cmline )
        fout.write( "#   " + name +"\n" )
        fout.write( cmline )
        for var in var_names:
            if var in varinants:
                fout.write("#====== variant : '" + var +"'\n" )
                for o in varinants[var]:
                    src = name+"/"+o+".f90"
                    fout.write( o+".o : "+src+"\n" )
                    fout.write( "\t$(F90) $(FFLAGS) -c "+src+"\n" )
        fout.write("\n\n")

def toPython( fout, groups, var_names, nclmax=120 ):
    ncl = 0
    for name,varinants in groups:
        fout.write( name+" = {\n" )
        for var in var_names:
            if var in varinants:
                #fout.write("# ====== variant : " + var )
                s = "  '"+ var +"' : ["
                ncl += len(s)
                fout.write( s )
                objs = varinants[var]
                no = len(objs)
                for i in range(no):
                    o = objs[i]
                    s = "'"+o+"'"
                    ncl += len(s)
                    if ( i<no-1 ):
                        s += ","
                        if ncl > nclmax:
                            s += "\n        "
                            ncl = 0
                    fout.write( s )
                fout.write("]\n")
                ncl = 0
        fout.write("}\n\n")

if __name__ == "__main__":
    groups,var_names = fromMakeFile( "Makefile.in" )
    pruneVariants(groups)
    with open("Makefile",'w') as fout:
        toMakefile_obj_groups( fout, groups, var_names )
        toMakefile_cc_objs   ( fout, groups, var_names )
    with open("Makefile.py",'w') as fout:
        toPython( fout, groups, var_names )



