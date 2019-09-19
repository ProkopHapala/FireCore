
import os
import sys
import re
import traceback
import string

scope_begins = {
    'if'        :('else','then'),
    'select'    :('case',),
    'do'        :(),
    'program'   :(),
    'module'    :('contains',),
    'type'      :('contains',),    # cause problems - see   MODULES/configuration.f90
    'subroutine':(),
    'function'  :(),
}

scope_ends = { 'enddo':'do', 'endif':'if' } 
scope_mids = { 'else','case','contains' }
re_split  = re.compile( '''!|'.*?'|".*?"|\w+''' )

def checkScope_re( line, stack, il ):
    ws = re_split.findall(line)
    #print ws
    if len(ws)==0: return
    w = ws[0]
    #print w
    if (w == 'end'):
        if len(stack)==0:
            print "DEBUG Line ", line
            raise Exception("line[%i] found 'end' but no scope is open" %il )
        if len(ws)>1:
            if ws[1] != stack[-1]:
                print "DEBUG Line  ", line
                print "DEBUG stack ", stack
                raise Exception("line[%i] colosing other scope ( '%s' not '%s' )" %(il,stack[-1],ws[1])  )
        #print "DEBUG stack.pop ", stack
        stack.pop()
        return -1
    elif w in scope_ends:
        if len(stack)==0:
            print "DEBUG Line ", line
            raise Exception("line[%i] found 'end' but no scope is open" %il )
        wend = scope_ends[w]
        if wend != stack[-1]:
                print "DEBUG Line  ", line
                print "DEBUG stack ", stack
                raise Exception("line[%i] colosing other scope ( '%s' not '%s' )" %(il,stack[-1],wend)  )
        #print "DEBUG stack.pop ", stack
        stack.pop()
        return -1
    elif w in scope_mids:
        if len(stack)==0:
            print "DEBUG Line ", line
            raise Exception("line[%i] found '%s' but no scope is open" %(il,w) )
        try:
            scb = stack[-1]
            tpw = scope_begins[scb]
            smid = tpw[0]
        except  Exception as e:
            raise Exception("line[%i] found '%s' but '%s'-scope is open" %(il,w,scb) )
            #traceback.print_exc()
        if w != smid:
            raise Exception("line[%i] found '%s' but '%s'-scope is open" %(il,w,scb) )
        return -1
    elif w in scope_begins:
        #print " scope_begin ", w
        tpw = scope_begins[w]
        if len(tpw)>1:   # check if 'if' has 'then'
            #print "tpw : ", len(tpw), tpw
            must_contain = tpw[1]
            if any( must_contain==wi for wi in ws[1:] ):
                stack.append(w)
                #print "DEBUG stack.append ", stack
        else:            # other than if
            stack.append(w)
            #print "DEBUG stack.append ", stack
    return 0

def rewriteFile_joinSplitLines( fin, fout, bNoEmptyLines=False, bNoComments=False, bScope=False, tab='    ' ):
    acum         = []
    acum_comment = []
    stack        = []
    for il,line in enumerate( fin ):

        level = len(stack)
        #print "DEBUG_0 :  ", line

        # ----- comments
        if len(line)==0: continue
        c=line[0]
        if c=='C' or c=='c' or c=='*' or c=='!' :
            if ( bNoComments ): continue 
            comment = line.rstrip()
            line=''
        else:
            ic = line.find('!')
            if ic >=0:
                if ( bNoComments ):
                    comment = ''
                else:
                    comment = line[ic:].rstrip()
                line    = line[:ic]
            else:
                comment = ''

        #print "DEBUG_1 :  ", line

        # ----- rsrip
        line = line.rstrip()
        if len(line)==0: 
            if ( (not bNoComments) and (not ( bNoEmptyLines and  (len(comment)==0) ) ) ):
                #print il,len(line),len(comment),bNoComments, bNoEmptyLines, line+comment
                fout.write(comment+'\n')
            continue

        if bScope: dlevel = checkScope_re( line, stack, il )

        #print "DEBUG_2 :  ", line

        if len(acum)>0:
            line = line.lstrip()
            if len(line)==0:  continue
            #if line[0]  != '&' :
            #    print " DEBUG acum : ", acum
            #    print " DEBUG line #"+il+" : '"+line+"'"
            #    raise Exception( "'&' expected to continue line above" )
            if line[0]  == '&' :
                line = line[1:].lstrip()
            #line = line[1:].lstrip()
            if len(line)==0: continue

        #print "DEBUG_3 :  ", line

        # ----- right &
        if line[-1] == '&':
            line = line[:-1].rstrip()
            if not ( bNoEmptyLines and (len(line)==0) ):
                acum.append( line )
        else:
            bLine = not ( bNoEmptyLines and ((len(line)+len(comment))==0) )
            bAcum = len(acum)>0
            if bScope and (bLine or bAcum ):
                fout.write( tab*(level+dlevel) )
                if bAcum:
                    acum[0]=acum[0].lstrip()
                line=line.lstrip()
            if bAcum:
                for l in acum: 
                    fout.write(l)
                    fout.write(' ')
                if not bLine:
                    fout.write('\n')
                acum = []
            if ( bLine ):
                #print il,len(line),len(comment), line+comment
                fout.write(line+" "+comment+'\n')

    if len(stack)>0:
        raise Exception( "File ended but scope.level is non-zero ( %i )" %len(stack) )

    for l in acum: fout.write(l)
    fout.write('\n')

def rewritePath( path_from, path_to, exts=['.f90'], rewriteFile_func=rewriteFile_joinSplitLines, verbose=True, bNoComments=False, bNoEmptyLines=False,  bScope=False, tab='    ' ):
    tab_dir  = "----"
    tab_file = "    "
    exts_set = set(exts)
    for path, dirs, files in os.walk(path_from):
        path_lst = path.split(os.sep)
        level=len(path_lst)
        path_rel = os.path.relpath( path, path_from )
        new_path = os.path.join(path_to,path_rel)
        if not os.path.isdir(new_path):
            os.makedirs (new_path)
        if(verbose):  print tab_dir*(level-1), os.path.basename(path)
        for fname in files:
            name,ext = os.path.splitext(fname)
            if ext in exts_set:
                if(verbose): print tab_file*level, fname
                fname_in  = os.path.join(path_from,path_rel,fname)
                fname_out = os.path.join(path_to  ,path_rel,fname)
                with open(fname_in,'r') as fin, open(fname_out,'w') as fout:
                    rewriteFile_func( fin, fout, bNoComments=bNoComments, bNoEmptyLines=bNoEmptyLines,  bScope=bScope, tab=tab )

def writeTabled(fout,s,ic,nc_min):
    nspace = nc_min-ic
    fout.write( (" "*nspace ) )
    fout.write(   s           )
    return ic+nspace+len(s)

def writeListNewLine(fout,lst, ntab=4 ):
    fout.write("[")
    for item in lst:
        fout.write("\n")
        writeTabled(fout,str(item),0,ntab)
        fout.write(",")
    fout.write("]")

if __name__ == "__main__":
    #path = "../SRC/LOOPS"
    #path = "../fortran"
    #fin  = DASSEMBLERS/Dassemble_lr_dip.f90

    '''
    #fname_in  = '../fortran/DASSEMBLERS/Dassemble_lr_dip.f90'
    #fname_out = '../fortran_tmp/DASSEMBLERS/Dassemble_lr_dip.f90'
    fname_in  = '../fortran/MODULES/wavefunction.f90'
    fname_out = '../fortran_tmp/MODULES/wavefunction.f90'
    with open( fname_in,'r') as fin, open( fname_out,'w') as fout:
        #rewriteFile_joinSplitLines( fin, fout )
        #rewriteFile_joinSplitLines( fin, fout, bNoEmtyLines=True )
        #rewriteFile_joinSplitLines( fin, fout,    bNoEmptyLines=True, bNoComments=True )
        #rewriteFile_joinSplitLines( fin, fout,    bNoEmptyLines=False, bNoComments=False, bScope=True, tab='\t' )
        #rewriteFile_joinSplitLines( fin, fout,    bNoEmptyLines=False, bNoComments=False, bScope=True, tab='  ' )
        rewriteFile_joinSplitLines( fin, fout,    bNoEmptyLines=True, bNoComments=True, bScope=True, tab='  ' )
    '''

    #rewritePath( "../fortran", "../fortran_tmp", exts=['.f90'], bNoEmptyLines=True, bNoComments=True, bScope=True, tab='  ' )
    rewritePath( "../fortran", "../fortran_tmp", exts=['.f90'], bNoEmptyLines=True, bNoComments=True, bScope=False, tab='  ' )



line_replce = {

'''
! Program Description
! ===========================================================================
'''
:"! ================ Program Description",

'''
! Procedure
! ===========================================================================
''',
:"! ================ Body ",

'''
! Local Parameters and Data Declaration
! ===========================================================================
'''
:"! ================ Local Parameters ",

'''
! Local Variable Declaration and Description
! ===========================================================================
'''
:"! ================ Local Variables  ",

'''
! Argument Declaration and Description
! ===========================================================================
'''
:"! ================ Arguments ",

'''
! Format Statements
! ===========================================================================
'''
:"! ================ Formats ",

'''
! Allocate Arrays
! ===========================================================================
''',
:"! ================ Allocations ",



}