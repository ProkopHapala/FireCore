
import os
import sys
import re
import traceback
import string

def insertRecord( dct, name, fname, rec ):
    if name in dct:
        dct_name = dct[name]
        dct_name[fname] = dct_name.get(fname,[]) + [rec]
    else:
        dct[name] = {fname:[rec]}

''''
def insertRecord( dct, name, fname, il, iw ):
    loc = [il]

def insertRecordLine( dct, name, fname, il, l ):
    rec = [(il,l)]
    if name in dct:
        dct_name = dct[name]
        dct_name[fname] = dct_name.get(fname,[]) + rec
    else:
        dct[name] = {fname:rec}
'''

def searchNamesInFile_re( fpath, fname, keyw_dict, re_split, bToLowerCase=True, bLine=True, found=None, verbose=False ):
    if found is None:
        found = [ {} for k in keyw_dict ]
    with open( fpath, 'r' ) as fin:
        for iline,line in enumerate(fin):
            c=line[0]
            if c=='C' or c=='c' or c=='*' or c=='!' : continue     # old-style comments
            ws = re_split.findall(line)
            for iw,w in enumerate(ws):
                if w == '!': break
                w = w.lower()
                ikey = keyw_dict.get(w,-1)
                if ikey >= 0:
                    insertRecord( found[ikey], w, fname, (iline,line.strip()) )
    return found

def searchKeywordsInFile_re( fpath, fname, keyw_dict, re_split, bToLowerCase=True, found=None, verbose=False ):
    if found is None:
        found = [ {} for k in keyw_dict ]
    with open( fpath, 'r' ) as fin:
        for iline,line in enumerate(fin):
            c=line[0]
            if c=='C' or c=='c' or c=='*' or c=='!' : continue     # old-style comments
            ws = re_split.findall(line)
            for iw,w in enumerate(ws):
                if w == '!': break
                w = w.lower()
                ikey = keyw_dict.get(w,-1)
                if ikey >=0:
                    try:
                        w_prev = ws[iw-1].lower()
                        if ((iw==0)or(w_prev!='end')):
                            w_next = ws[iw+1].lower()
                            i = keyw_dict[w]
                            if(verbose): print "found.append( %s )" %str(rec)
                            insertRecord( found[ikey], w_next, fname, iline )
                    except Exception as e:
                        traceback.print_exc()
                        print "DEBUG w{"+w+"} line{"+line+"}"
                        print "DEBUG ws = ",  ws
                        raise Exception("ERROR in searchKeywordsInFile_split [%i,%i] of %s" %(iline,iw,fname)  )
    return found

#def searchKeywordsInPath( root_path, keywords, exts=['.f90','.f'], wchars="a-zA-Z0-9_", tab_dir='----', tab_file='    ', verbose=True, bToLowerCase=True, bStorePath=True ):
def searchInPath( root_path, keywords, exts=['.f90','.f'], tab_dir='----', tab_file='    ', verbose=True, bToLowerCase=True, bStorePath=True, serchFunc=searchKeywordsInFile_re ):
    exts_set = set(exts)
    keyw_dict = { kw:i for i,kw in enumerate(keywords) }
    re_split  = re.compile( '''!|'.*?'|".*?"|\w+''' )
    found = [ {} for k in keyw_dict ]
    nkw = len(keywords)
    for path, dirs, files in os.walk(root_path):
        path_lst     = path.split(os.sep)
        level=len(path_lst)
        path_rel = os.path.relpath( path, root_path )
        if(verbose): print tab_dir*(level-1), os.path.basename(path)
        for fname in files:
            name,ext = os.path.splitext(fname)
            fpath = path+os.sep+fname
            if ext in exts_set:
                if(verbose):
                    nfound0 = [ len(found[i]) for i in xrange(nkw) ]
                fname_=fname
                if bStorePath:fname_= path_rel +os.sep+ fname
                serchFunc( fpath, fname_, keyw_dict, re_split, found=found, bToLowerCase=bToLowerCase )
                if(verbose):
                    nfound  = [ len(found[i])-nfound0[i] for i in xrange(nkw) ]
                    print tab_file*level, fname, "       ", nfound
    return found




"""
def searchKeywordsInPath( root_path, keywords, exts=['.f90','.f'], wchars="a-zA-Z0-9_", tab_dir='----', tab_file='    ', verbose=True, bToLowerCase=True, bStorePath=True ):
    exts_set = set(exts)
    keyw_dict = { kw:i for i,kw in enumerate(keywords) }
    re_split  = re.compile( '''!|'.*?'|".*?"|\w+''' )
    found = [ {} for k in keyw_dict ]
    nkw = len(keywords)
    for path, dirs, files in os.walk(root_path):
        path_lst     = path.split(os.sep)
        level=len(path_lst)
        path_rel = os.path.relpath( path, root_path )
        if(verbose): print tab_dir*(level-1), os.path.basename(path)
        for fname in files:
            name,ext = os.path.splitext(fname)
            fpath = path+os.sep+fname
            if ext in exts_set:
                if(verbose):
                    nfound0 = [ len(found[i]) for i in xrange(nkw) ]
                fname_=fname
                if bStorePath:fname_= path_rel +os.sep+ fname
                searchKeywordsInFile_re( fpath, fname_, keyw_dict, re_split, found=found, bToLowerCase=bToLowerCase )
                if(verbose):
                    nfound  = [ len(found[i])-nfound0[i] for i in xrange(nkw) ]
                    print tab_file*level, fname, "       ", nfound
    return found
"""

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

def writeTags( fout, dct, ntab1=35, ntab2=80, bSortNames=True, bNewLineLine=False, bAsVar=True ):
    names = dct.keys()
    if bSortNames:
        names=sorted(names)
    for name in names:
        ic=0
        recs = dct[name]
        nrec = len(recs)
        if nrec==0:
            print "WARRNING: record for name{"+name+"} is empty" 
        elif nrec==1:
            fname=next(iter(recs))
            if bAsVar:
                ic= writeTabled(fout, "%s" %name, 0,0 )
                ic= writeTabled(fout, "={",   ic,ntab1-5 )
            else:
                ic= writeTabled(fout, "'%s'" %name, 0,0 )
                ic= writeTabled(fout, ":{",   ic,ntab1-5 )
            ic= writeTabled(fout, "'%s'" %str(fname      ), ic, ntab1 )
            if bNewLineLine:
                fout.write(":")
                writeListNewLine(fout, recs[fname], ntab=ntab2 )
            else:
                ic= writeTabled(fout, ":%s"  %str(recs[fname]), ic, ntab2 )
            ic= writeTabled(fout, "},\n", ic, ntab2 )
        else:
            if bAsVar:
                ic= writeTabled(fout, "%s" %name, 0,0 )
                ic= writeTabled(fout, "={", ic,ntab1-5 )
            else:
                ic= writeTabled(fout, "'%s'" %name, 0,0 )
                ic= writeTabled(fout, ":{", ic,ntab1-5 )
            fnames=recs.keys()
            if bSortNames:
                fnames=sorted(fnames)
            for fname in fnames:
                fout.write("\n")
                ic = writeTabled(fout, "'%s'"  %str(fname), 0, ntab1 )
                if bNewLineLine:
                    fout.write(":")
                    writeListNewLine(fout, recs[fname], ntab=ntab2 )
                    fout.write(",")
                else:
                    ic= writeTabled(fout, ":%s,"  %str(recs[fname]), ic, ntab2 )
            ic= writeTabled(fout, "},\n", ic, ntab2 )

def writeFoundTagFile_py( fname, keywords, found, ntab1=4, ntab2=8, bSortNames=True, bNewLineLine=False, bAsVar=True ):
    with open(fname, 'w' ) as fout:
        for ik,key in enumerate(keywords):
            writeTags( fout, found[ik], ntab1=ntab1, ntab2=ntab2, bSortNames=bSortNames, bNewLineLine=bNewLineLine, bAsVar=bAsVar )

def writeFoundTagFiles_py( keywords, found, ntab1=35, ntab2=80, bSortNames=True, bNewLineLine=False, bAsVar=True ):
    for ik,key in enumerate(keywords):
        with open("tags_%s.py" %key, 'w' ) as fout:
            fout.write( "tag_dict_%s ={\n" %key )
            writeTags(fout, found[ik], ntab1=ntab1, ntab2=ntab2, bSortNames=bSortNames, bNewLineLine=bNewLineLine, bAsVar=bAsVar )
            fout.write( "}\n\n" )

def writeDefCallPairs( fout, defs, calls, tab="    ", bSort=True ):
    defs_set = set( defs.keys()  )
    call_set = set( calls.keys() )
    matched_set = call_set.intersection(defs_set)
    matched_lst = list( matched_set )
    d_calls = list( call_set - matched_set )
    d_defs  = list( defs_set - matched_set )
    if bSort:
        d_calls.sort()
        d_defs.sort()
        matched_lst.sort()
    fout.write( "unmatched_calls = "); fout.write( str( d_calls ) ); fout.write("\n" )
    fout.write( "unmatched_defs  = "); fout.write( str( d_defs  ) ); fout.write("\n" )
    for name in matched_lst:
        defs_name  = list( defs[name] .keys() )
        calls_name = list( calls[name].keys() )
        if bSort:
            defs_name.sort ()
            calls_name.sort()
        #print name," : {" 
        #print "'Defs'  : ", tab, defs_name  ,"," 
        #print "'Calls' : ", tab, calls_name ,","
        #print"}"
        fout.write( name); fout.write(" = {\n" )
        fout.write( "'Defs'   : "); fout.write( str(defs_name ) ); fout.write(",\n" )
        fout.write( "'Calls'  : "); fout.write( str(calls_name) ); fout.write("\n"  )
        fout.write( "}\n" )

if __name__ == "__main__":

    #path = "../SRC/LOOPS"
    path = "../fortran"
    #path = "../fortran_tmp"

    names = [
        'itheory','itheory_xc',
        'iks','iqmmm','iordern','ibias',
        'iforce',
        'idipole','igauss','ivdw','iharmonic',
        'iwrtfpieces'
    ]
    found = searchInPath( path, names, serchFunc=searchNamesInFile_re )
    #writeFoundTagFiles_py( names, found, ntab1=4, ntab2=8, bNewLineLine=True )
    writeFoundTagFile_py( "tags_ioption.py",  names, found, ntab1=4, ntab2=40, bSortNames=True, bNewLineLine=False )
    #exit()

    keywords = [ 'function', 'subroutine', 'call', 'use' ]
    found = searchInPath( path, keywords, serchFunc=searchKeywordsInFile_re )
    writeFoundTagFiles_py( keywords, found, bAsVar=False )
    
    from tags_subroutine import *
    from tags_call       import *
    with open('tags_subroutine2call.py','w') as fout:
        fout.write("\n\n\n\n\n\n")
        writeDefCallPairs( fout, tag_dict_subroutine, tag_dict_call  )

