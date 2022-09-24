
#ifndef SMILESparser_h
#define SMILESparser_h

#include <vector>
#include <unordered_map>
#include <cstdio>
#include <cstring>

#include "MMFFBuilder.h"


//char charCheck[] = "                                 !"#$%&'()*+,-./0123456789:;<=>?@ABCDEFGHIJKLMNOPQRSTUVWXYZ[\]^_`abcdefghijklmnopqrstuvwxyz{|}~";
char charTypes[]   = "                                 -'B---'()--,B--0000000000---B---AAAAAAAAAAAAAAAAAAAAAAAAAA(-)--'aaaaaaaaaaaaaaaaaaaaaaaaaa(-)-";

class SMILESparser{ public:
    //std::vector<ParserItem> items;
    MM::Builder* builder;

    int ich,nch;
    int itok;
    char * str;

    // params
    char cOPEN  = '(';
    char cCLOSE = ')';
    char cNEXT  = ';';

    constexpr static int ntokmax=16;
    char tok[ntokmax];

    std::unordered_map<int,int> anchors;

    // ============== Functions

    int getBondOrder( char bond ){
        int order=-1;
        if(bond=='-'){ order=1; }
        if(bond=='='){ order=2; }
        if(bond=='#'){ order=3; }
        return order;
    }

    void insertAnchor(int ia, char bond){
        tok[itok]='\0';
        int ianch = atoi(tok);
        auto found = anchors.find(ianch);
        if( found==anchors.end() ){     // NEW
            anchors.insert({ianch,ia});
        }else{                         // OLD
            int ia2 = found->second;
            int order = getBondOrder( bond );
            builder->insertBond( {ia,ia2}, order );
        }
    };

    int insertAtom(int ia, char bond){
        tok[itok]='\0';
        printf("insertAtom ia=%i bond=%c tok=%s\n", ia, bond, tok );
        int ia2 = builder->insertAtom( tok, 0,0,0,0 );
        if(ia>=0){
            int order = getBondOrder( bond );
            builder->insertBond( {ia,ia2}, order );
        }
        return ia2;
    };   

    void parseString( int nch_, char * str_ ){
        str = str_;
        nch = nch_;
        ich = 0;
        //printf( "%s\n", str );
        parse( 0, 0, -1 );
    }

    int parse( int ich0, int level, int parent ){
        ich         = ich0;
        int ia      = -1; 
        char bond   = '-';
        //printf( "%s\n", str );
        printf( "ich0 %i level %i \n", ich0, level );
        char tokType=' ';
        while( ich < nch ){
            char ch = str[ich];
            printf( "parse[%i](%c)\n", ich, ch );
            if (ch==cOPEN) {  // start sub-chain
                ich = parse( ich+1, level+1, ia );
            }else if (ch==cCLOSE) { // end sub-chain
                return ich;
            }else{ // token
                char t=charTypes[ch];
                printf( "tokType(%c):t(%c)\n", tokType, t );
                // proces previous token
                if(t!=tokType){ // close token
                    if      (tokType=='0'){ // anchor point
                        insertAnchor(ia,bond);         
                        itok=0;
                    }else if(tokType=='a'){ // atom
                        ia = insertAtom(ia,bond);
                        bond='-'; 
                        itok=0;
                    }
                }
                tokType=t;
                if     (t=='A'){ tokType='a'; }
                if(t=='B'){ 
                    bond=ch;
                }else{
                    tok[itok]=ch; 
                    itok++;     
                    if(itok>256){printf("ERROR: itok(%i)>tok.size(%i)", itok, ntokmax ); exit(0);}
                }
                if(ch=='\0'){ nch=0; break;} // Terminate on null character
            }
            ich++;
        }
        return ich;
    }

};

#endif


