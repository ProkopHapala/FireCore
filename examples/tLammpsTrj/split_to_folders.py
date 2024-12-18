#!/usr/bin/python
import os,sys
import BasePairUtils as bpu

'''
mols_gs=[
("D",["H-p","H-h_1","H-h_2","H-hh",]),
("A",["N-h","N-hh","O-p","O-h",]),
("DD",["HH-p_1","HH-p_2","HH-h_1","HH-h_2","HH-pp","HH-hp","HH-h-p","HH-hh","HH-hh-p",]),
("DA",["HN-h","HN-pp","HN-hp_1","HN-hp_2","HN-h-p","HN-hh","HO-h",]),
("AA",["NN-pp","NN-hp","NN-hh","NO-p","NO-h","OO-h",]),   #"NO-h-p",
("DDD",["HHH-p","HHH-h","HHH-hp","HHH-h-p","HHH-hh",]),
("DDA",["HHN-hp","HHN-hh","HHO-p","HHO-h","HHO-hp","HHO-h-p_1","HHO-h-p_2","HHO-hh",]),
("DAD",["HNH-p","HNH-h","HNH-hp","HNH-h-p","HNH-hh",]),
("DAA",["HNN-hp","HNN-hh","HNO-p","HNO-h","HNO-hp","HNO-h-p","HNO-hh",]),
("ADA",["NHO-hp","NHO-hh","OHO-p","OHO-h_1","OHO-h_2","OHO-h-p"]),
("AAA",["NNN-hhh","NNO-hh_1","NNO-hh_2",]), #"NNO-hp","ONO-p"
]
'''


#dirs=[ 'H','HH','HHH',  'e','ee','eee',   'He', 'HHe', 'Hee', 'HeH', 'eHe' ]


work_dir='/home/prokop/Desktop/CARBSIS/Mithun/B3LYP_Barbora/all'

# got to work_dir
os.chdir(work_dir)
# for d in dirs:
#     try:
#         os.mkdir(d)
#     except:
#         print("WARNING: directory", d, "already exists")

# list all .xyz files in directory
files  = [ f for f in os.listdir(work_dir) if os.path.isfile(f) and f.endswith(".xyz") ]

known_folders = set()
for f in files:
    # split extension
    name = os.path.splitext(f)[0]
    pair = bpu.split_pair_with_S( name )
    folder = bpu.name_to_class( pair[0] )+"_"+ bpu.name_to_class( pair[1] )

    print( folder, "      ", f  )
    if folder not in known_folders:
        known_folders.add( folder )
        print( "make dir: ", folder )
        try:
            os.mkdir( folder )
        except:
            print("WARNING: directory", folder, "already exists")

    # copy file to appropriate directory
    os.system('cp '+f+' '+folder+'/'+f)




