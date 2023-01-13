
#ifndef LatticeMatch2D_h
#define LatticeMatch2D_h

#include "Vec2.h";



void findLatticeR( double Rmin, double Rmax, Vec2d a, Vec2d b ){
    double R0 = (Rmin-Rmax)*0.5;
    double ra = a.norm();
    double rb = b.norm();
    double dna = R0/ra; int na=(int)(dna+0.5); dna-=na; 
    double dnb = R0/rb; int nb=(int)(dnb+0.5); dnb-=nb;

    double bonus = (a.x*b.x + a.y*b.y)*2;

    n.x

}



class LatticeMatch2D{
    Vec2d lat0[2];
    Vec2d lat1[2];



}

#endif



