
import os
import sys
import re
import traceback
import string

def insertRecord( dct, name, fname, il, iw ):
    loc = [il]
    if name in dct:
        dct_name = dct[name]
        dct_name[fname] = dct_name.get(fname,[]) + loc
    else:
        dct[name] = {fname:loc}

def searchKeywordsInFile_re( fpath, fname, keyw_dict, re_split, bToLowerCase=True, bFname=True, found=None, verbose=False ):
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
                            insertRecord( found[ikey], w_next, fname, iline, iw )
                    except Exception as e:
                        traceback.print_exc()
                        print "DEBUG w{"+w+"} line{"+line+"}"
                        print "DEBUG ws = ",  ws
                        raise Exception("ERROR in searchKeywordsInFile_split [%i,%i] of %s" %(iline,iw,fname)  )
    return found

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

def writeTabled(fout,s,ic,nc_min):
    nspace = nc_min-ic
    fout.write( (" "*nspace ) )
    fout.write( s                  )
    return ic+nspace+len(s)

def writeTags( fout, dct, tab="    ",  ntab1=35, ntab2=80, bSortNames=True ):
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
            rec=next(iter(recs))
            ic= writeTabled(fout, "'%s'" %name, 0,0 )
            ic= writeTabled(fout, ":{", ic,ntab1-5 )
            ic= writeTabled(fout, "'%s'"  %str(rec), ic, ntab1 )
            ic= writeTabled(fout, ":%s"  %str(recs[rec]), ic, ntab2 )
            ic= writeTabled(fout, "},\n", ic, ntab2 )
        else:
            ic= writeTabled(fout, "'%s'" %name, 0,0 )
            ic= writeTabled(fout, ":{", ic,ntab1-5 )
            fnames=recs.keys()
            if bSortNames:
                fnames=sorted(fnames)
            for rec in fnames:
                fout.write("\n")
                ic = writeTabled(fout, "'%s'"  %str(rec), 0, ntab1 )
                ic = writeTabled(fout, ":%s,"  %str(recs[rec]), ic, ntab2 )
            ic= writeTabled(fout, "},\n", ic, ntab2 )

def writeFoundTagFiles_py( keywords, found ):
    for ik,key in enumerate(keywords):
        with open("tags_%s.py" %key, 'w' ) as fout:
            fout.write( "tag_dict_%s ={\n" %key )
            writeTags(fout, found[ik] )
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

    
    keywords = [ 'function', 'subroutine', 'call', 'use' ]
    #path = "../SRC/LOOPS"
    path = "../fortran"
    found = searchKeywordsInPath( path, keywords )
    writeFoundTagFiles_py( keywords, found )
    
    from tags_subroutine import *
    from tags_call       import *
    with open('tags_subroutine2call.py','w') as fout:
        fout.write("\n\n\n\n\n\n")
        writeDefCallPairs( fout, tag_dict_subroutine, tag_dict_call )

