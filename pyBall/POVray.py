
import os,sys
import numpy as np

def write_Camera( fout, angle=16.0, up=(500,0,0), right=(500,0,0),  location=(250.0,250.0,0), look_at=(250.0,250.0,1000) ):
    '''
        perspective
        angle 16.260204
        right < 500, 0, 0>
        up < 0, -500, 0 >
        sky < 0, -1, 0 >
        location < 250.0, 250.0, 0>
        look_at < 250.0, 250.0, 1000 >
    '''
    right =

    fout.write( f"""camera{{
    perspective
    angle {angle}
    right <{right[0]},{right[1]},{right[2]}>
    up <{up[0]},{up[1]},{up[2]}>
    sky < 0, -1, 0 >
    location <{location[0]},{location[1]},{location[2]}>
    look_at  <{look_at[0]},{look_at[1]},{look_at[2]}>
    }}""")

def write_material( fout, ambient=0.45, diffuse=0.84, specular=0.22, roughness=.00001, metallic=0.0, phong=0.9, phong_size=120. ):
    '''
    ambient    0.45
    diffuse    0.84
    specular   0.22
    roughness .00001
    metallic   shineFactor
    phong      0.9*shineFactor
    phong_size 120*shineFactor
    '''
    S = f"""
    #macro atomMaterial(T)
    finish {{
    ambient    {ambient}
    diffuse    {diffuse}
    specular   {specular}
    roughness  {roughness}
    metallic   {metallic}
    phong      {phong}
    phong_size {phong_size}
    }}#end
    """
    fout.write( S )

def write_basic_macros():
    S = '''
    #macro a(X,Y,Z,RADIUS,R,G,B,T)
        sphere{<X,Y,Z>,RADIUS
        pigment{rgbt<R,G,B,T>}
        atomMaterial(T)
    }#end

    #macro b(X1,Y1,Z1,R1, X2,Y2,Z2, R,G,B,T )
        cylinder{<X1,Y1,Z1>,<X2,Y2,Z2>,R
        pigment{rgbt<R,G,B,T>}
    }#end

    #macro c(X1,Y1,Z1,R1, X2,Y2,Z2,R2,  R,G,B,T)
        cone{<X1,Y1,Z1>,R1,<X2,Y2,Z2>,R2
        pigment{rgbt<R,G,B,T>}
    }#end
'''