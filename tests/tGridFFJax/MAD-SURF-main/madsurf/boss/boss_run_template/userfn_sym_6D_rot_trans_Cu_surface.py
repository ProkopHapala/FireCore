import os
#import decimal
import numpy as np

def f(X):
    "User function calling an external script"

    # write variables to file
    f = open('energy_calc/variables.in', 'w')
    for i in range(len(X[0])):
        f.write('%15.7f\n'%(X[0][i]))
    f.close()
    
    # Sets of equivalent variables going into BOSS search
   #decimal.getcontext().prec = 7 # Precision set to 7 decimals
    # X1 = decimal.Decimal(np.copy(X))
    X1 = np.copy(X)
    X2 = np.copy(X)
    X3 = np.copy(X)
    X4 = np.copy(X)  
    x1 = X[0, 0] # x-variable
    x2 = X[0, 1] # y-variable
    x6 = X[0, 5] # Rz-variable
    
    # Define symmetries
    #X2 = np.zeros((1, 2))
    x1s = 0 # x low bound
    x2s = 0 # y low bound
    x1h = 1.284 # halfway-point of x-coordinate
    x2h = 2.224 # halfway-point of y-coordinate
    #x1e = 2.568 # x high bound
    x2e = 4.448 # y high bound
    
    # Rounded variables for comparison-purposes:
    x2sym1 = round((x2e/3), 7)
    x2sym2 = round((x2e*(5/6)), 7)
    x2sym3 = round(((x2e/6)+x2h), 7)
    x2sym4 = round(((x2e/6)), 7)
    
    # Translational symmetry (two-fold translational symmetry on Cu111):
    if x1 < x1h:
        if x2 > x2h:
            X2[0, 0] = x1 + x1h
            X2[0, 1] = x2 - x2h
        elif x2 < x2h:
            X2[0, 0] = x1 + x1h
            X2[0, 1] = x2 + x2h
    elif x1 > x1h:
        if x2 > x2h:
            X2[0, 0] = x1 - x1h
            X2[0, 1] = x2 - x2h
        elif x2 < x2h:
            X2[0, 0] = x1 - x1h
            X2[0, 1] = x2 + x2h
    else:
        X2[0, 0] = 0
        X2[0, 1] = 0
    
    # Rotational symmetry in Rz (three-fold at high-symmetry points with Cu111):
    if x1 == x1s and x2 == x2s: # First layer
        if x6 <= -60.0:       
            X3[0, 5] = x6 + 240.0
            X4[0, 5] = x6 + 120.0
        elif x6 >= 60.0:
            X3[0, 5] = x6 - 240.0
            X4[0, 5] = x6 - 120.0
        elif x6 < 60.0 and x6 > -60.0: 
            X3[0, 5] = x6 - 120.0
            X4[0, 5] = x6 + 120.0
    elif x1 == x1h and x2 == x2h:        
        if x6 <= -60.0:       
            X3[0, 5] = x6 + 240.0
            X4[0, 5] = x6 + 120.0
        elif x6 >= 60.0:
            X3[0, 5] = x6 - 240.0
            X4[0, 5] = x6 - 120.0
        elif x6 < 60.0 and x6 > -60.0: 
            X3[0, 5] = x6 - 120.0
            X4[0, 5] = x6 + 120.0
    elif x1 == x1s and x2 == x2sym1: # Second layer
        if x6 <= -60.0:       
            X3[0, 5] = x6 + 240.0
            X4[0, 5] = x6 + 120.0
        elif x6 >= 60.0:
            X3[0, 5] = x6 - 240.0
            X4[0, 5] = x6 - 120.0
        elif x6 < 60.0 and x6 > -60.0: 
            X3[0, 5] = x6 - 120.0
            X4[0, 5] = x6 + 120.0
    elif x1 == x1h and x2 == x2sym2:        
        if x6 <= -60.0:       
            X3[0, 5] = x6 + 240.0
            X4[0, 5] = x6 + 120.0
        elif x6 >= 60.0:
            X3[0, 5] = x6 - 240.0
            X4[0, 5] = x6 - 120.0
        elif x6 < 60.0 and x6 > -60.0: 
            X3[0, 5] = x6 - 120.0
            X4[0, 5] = x6 + 120.0
    elif x1 == x1s and x2 == x2sym3: # Third layer
        if x6 <= -60.0:       
            X3[0, 5] = x6 + 240.0
            X4[0, 5] = x6 + 120.0
        elif x6 >= 60.0:
            X3[0, 5] = x6 - 240.0
            X4[0, 5] = x6 - 120.0
        elif x6 < 60.0 and x6 > -60.0: 
            X3[0, 5] = x6 - 120.0
            X4[0, 5] = x6 + 120.0
    elif x1 == x1h and x2 == x2sym4:        
        if x6 <= -60.0:       
            X3[0, 5] = x6 + 240.0
            X4[0, 5] = x6 + 120.0
        elif x6 >= 60.0:
            X3[0, 5] = x6 - 240.0
            X4[0, 5] = x6 - 120.0
        elif x6 < 60.0 and x6 > -60.0: 
            X3[0, 5] = x6 - 120.0
            X4[0, 5] = x6 + 120.0
    
    # Defining new arrays with rotationally equivalent variables for the translated x,y-variables: 
    X5 = np.copy(X2) 
    X6 = np.copy(X2)
    X5[0, 5] = X3[0, 5]
    X6[0, 5] = X4[0, 5]
           
    # Call bash script for AIMS simulation
    os.system('./run_6D_rot_trans_Cu_surf.sh')
    
    # Read energy from file and return it
    f = open('energy_calc/energy.out')
    e = float(f.readline())
    f.close()
    
    # Adjust the energy if it is above a given threshold, i.e. lower the energy of higher energy regions:
    if e > 1.0:
    	e = np.log10(e) + 1.0
    
    # Check for changes in rotational coordinate Rz. If none, return only translated variables:
    if np.array_equal(X1, X3): 
        # Remove duplicates:
        unique_coords = []
        for coord in np.concatenate((X1, X2), axis=0):
            if coord.tolist() not in unique_coords:
                unique_coords.append(coord.tolist())
        X = np.array(unique_coords)
        # Returning the x-locations as rows in a 2D array.
       # X = np.concatenate((X1, X2), axis=0)
        # The corresponding evaluations must be returned in a 1D array.
        #E = np.array([e, e])
    else:
        # Remove duplicates:
        unique_coords = []
        for coord in np.concatenate((X1, X2, X3, X4, X5, X6), axis=0):
            if coord.tolist() not in unique_coords:
                unique_coords.append(coord.tolist())
        X = np.array(unique_coords)
        #X = np.concatenate((X1, X2, X3, X4, X5, X6), axis=0)
        #E = np.array([e, e, e, e, e, e])  

    # Returning the x-locations as rows in a 2D array.
    E = np.array([e]*len(X))   
    
    return X, E

