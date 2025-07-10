from ctypes import wstring_at
import sys
import numpy as np
import os
import matplotlib.pyplot as plt
import time
# import torch
import tensorflow as tf
import scipy.optimize as opt 

sys.path.append("../../")
from scipy.optimize import dual_annealing
from pyBall import eFF as eff

elementPath_e = "export/scan_data/angdistscan_CH4_e2.xyz"
elementPath = "export/scan_data/angdistscan_CH4.xyz"
fileToSavePath = "results/AI/result_dual_anneal.txt"
fileToSaveProcess = "processXYZ.xyz"
numOfGetVar = 0


# class Angdist_sim:
def count_mask_lines( fgo_file, mask='#iconf' ):
    line = None
    nline = 0
    with open(fgo_file) as f:
        for l in f:
            if mask in l:
                line = l
                nline += 1
    return nline, line
                

def extract_nae( xyz_file ):
    with open(xyz_file) as f:
        ws = f.readline().strip().split()
        na= int(ws[0])
        ne= int(ws[1])
        return na, ne

def extract_blocks(xyz_file):
    """Extract parameters from XYZ file comments (lines starting with #)
    Returns:
        dict: Dictionary of extracted parameters with NaNs for missing values
        (e.g. {'ang': [...], 'dist': [...], 'Etot': [...]})
    """
    all_keys = set()
    records = []
    # First pass: collect all keys and raw records
    with open(xyz_file) as f:
        for line in f:
            if line.startswith('#'):
                parts = line[1:].strip().split()
                record = {}
                for i in range(0, len(parts)-1, 2):
                    key = parts[i]
                    val = float(parts[i+1])
                    record[key] = val
                    all_keys.add(key)
                records.append(record)
    # Initialize params with NaN-filled arrays
    nrec= len(records)
    params = {key: np.full(nrec, np.nan) for key in all_keys}
    # Second pass: fill values
    for i, record in enumerate(records):
        for key, val in record.items():
            params[key][i] = val
    return params, nrec

def plot_energy_landscape( Xs, Ys, Es, Espan=None):
    """Plot energy landscape from XYZ file using imshow (simple and robust)."""
    #params, nrec = extract_blocks(xyz_file)
    xs, idx = np.unique ( Xs, return_inverse=True )
    ys, idy  = np.unique( Ys, return_inverse=True )
    energy_grid = np.full((len(xs), len(ys)), np.nan)
    # Fill the grid
    for i in range(len(Xs)):
        energy_grid[ idx[i], idy[i]] = Es[i]
    # Plot
    plt.figure(figsize=(10,8))
    if Espan is not None: 
        vmin = np.nanmin(energy_grid)
        vmax = vmin + Espan
    else:
        vmin = None
        vmax = None
    plt.imshow(energy_grid, extent=[ys.min(), ys.max(), xs.max(), xs.min()], aspect='auto', cmap='inferno', vmin=vmin, vmax=vmax)
    plt.colorbar(label='Total Energy ')
    plt.xlabel('Distance (Å)')
    plt.ylabel('Angle (rad)')
    #plt.title('Potential Energy Surface')     


def save_simulation(angleArr, distArr, minTheta ,fileToSavePath, allEtot=None, flexVar=None, variance=None):
    if os.path.exists(fileToSavePath):
        os.remove(fileToSavePath)   
    
    with open(fileToSavePath, "w") as f:
        # Convert and flatten arrays (to ensure 1D), then write manually
        angleArrnp = np.array(angleArr).flatten()
        distArrnp = np.array(distArr).flatten()
        
        f.write(" ".join(f"{val:.6f}" for val in minTheta) + "\n")

        if (allEtot is not None) and (flexVar is not None) and (variance is not None):
            # Write angle and dist arrays in one line each
            f.write(" ".join(f"{val:.6f}" for val in angleArrnp) + "\n")
            f.write(" ".join(f"{val:.6f}" for val in distArrnp) + "\n")
            f.write("\n\n\n")  # 3 blank lines

            for i in range(len(flexVar)):
                # Ensure each input is 1D
                flex_row = np.array(flexVar[i]).flatten()
                var_row = np.array(variance[i]).flatten()
                etot_row = np.array(allEtot[i]).flatten()

                # Write each on a single line
                f.write(" ".join(f"{val:.6f}" for val in flex_row) + "\n")
                f.write(" ".join(f"{val:.6f}" for val in var_row) + "\n")
                f.write(" ".join(f"{val:.6f}" for val in etot_row) + "\n")
                f.write("\n")  # blank line between each block

def get_variance_from_KRSrho(theta): #Function we need to minimize
    global numOfGetVar
    outEs = np.zeros((nrec,5))
    eff.setKRSrho(theta)
    eff.processXYZ_e(elementPath_e, outEs=outEs, nstepMax=100, dt=0.005, Fconv=1e-3) # FOR NORMAL PURPOSES USE nstepMax=10000
    outEsdiff = np.zeros((nrec,5))
    outEsdiff[:,0] = params["Etot"] - outEs[:,0]
    outEsdiff[:,0] -= sum(outEsdiff[:,0])/len(outEsdiff[:,0])
    outEsdiff[:,0][params['dist'] == 0.7] = 0 #getting rid of unwanted valuess 
    var = np.var(outEsdiff[:,0])
    numOfGetVar += 1
    print("Number of get variance", numOfGetVar)
    print("Flexible variable: ", theta)
    print("Variance: ", var)
    print("\n")
    return var #, outEsdiff[:,0]

def get_variance_from_KRSrho_tf(theta): #Function we need to minimize
    global numOfGetVar
    outEs = np.zeros((nrec,5))
    eff.setKRSrho(theta)
    eff.processXYZ_e(elementPath_e, outEs=outEs, nstepMax=100, dt=0.005, Fconv=1e-3) # FOR NORMAL PURPOSES USE nstepMax=10000
    outEsdiff = np.zeros((nrec,5))
    outEsdiff[:,0] = params["Etot"] - outEs[:,0]
    outEsdiff[:,0] -= sum(outEsdiff[:,0])/len(outEsdiff[:,0])
    outEsdiff[:,0][params['dist'] == 0.7] = 0 #getting rid of unwanted valuess 
    var = np.var(outEsdiff[:,0])
    numOfGetVar += 1
    print("Number of get variance", numOfGetVar)
    print("Flexible variable: ", theta)
    print("Variance: ", var)
    print("\n")
    return var #, outEsdiff[:,0]

def minimizeBySteps(func, x0):
    allEtot = []
    flexVar = [] # Flexible variable, AKA theta
    variance = [] # vraiance we want to minimize
    iteretions = 3 #FOR NORMAL PURPUSES USE 50

    theta0 = x0.copy()
    steps = np.full_like(x0, 0.1)
    # steps = np.array([0.05, 0.05, 0.05])    
    # steps = np.array([0.01, 0.01, 0.01])

    precision = 0.001
    deltaTheta = np.full_like(x0, 0.0)

    normalConst = np.sqrt(len(theta0))/len(theta0) #such that each vector of movement is on ELIPSE that enclose each point, calculated through the step method

    for j in range(iteretions):
        var0, outEsTot0 = func(theta0)
        for i, step in enumerate(steps): 
            var, outEsTot = var0.copy(), outEsTot0.copy()
            theta1 = theta0.copy()
            theta1[i] = theta1[i] + step
            theta2 = theta0.copy()
            theta2[i] = theta2[i] - step
            
            var1, outEsTot1 = func(theta1)
            var2, outEsTot2 = func(theta2)
            
            
            newStep = step
            if var1 < var2:
                var = var1
                outEsTot = outEsTot1
                theta = theta1.copy()
                deltaTheta[i] = step

                if var1 > var0:
                    newStep = step/2
            else:
                var = var2
                outEsTot = outEsTot2
                theta = theta2.copy()
                deltaTheta[i] = -step
                
                if var2 > var0:
                    newStep = step/2
                
            steps[i] = newStep

            allEtot.append(outEsTot.copy())
            flexVar.append(theta.copy())
            variance.append(var)
            print("----------------------------------")
            print("Theta1: ", theta1)
            print("Theta0: ", theta0)
            print("Theta2: ", theta2)
            print("Var1: ", var1)
            print("Var0: ", var0)
            print("Var2: ", var2)
            print("Theta: ", theta)
            print("Iteretion: ", j)
            print("Value index: ", i)
            print("Step: ", newStep)
            print("Variance: ", var)
            # print("--------------------")
            # input("Press Enter to continue...")

        print(deltaTheta)
        print(theta0)
        theta0 = theta0 + deltaTheta.copy()*normalConst
        print("======================================================================================================")
        print("Theta: ", theta)
        print("Theta0: ", theta0)
        print("Steps: ", steps)
        print("Variance: ", var)
        print("======================================================================================================")

        if ((steps < precision) & (steps > -precision)).all():
            print("PRECISSION MET")
            break

    return theta, steps, var, allEtot, flexVar, variance

@tf.custom_gradient
def get_variance_wrapper(theta):
    # Use tf.py_function to wrap the numpy-based black-box function.
    # This allows it to be part of the TF graph execution. It tells TensorFlow:
    # "Run this Python code; I will provide the gradient for it myself."
    loss = tf.py_function(func=get_variance_from_KRSrho, inp=[theta], Tout=tf.float64)

    # Define the backward pass (the gradient function)
    def grad_fn(dy):
        # dy is the upstream gradient. For a scalar loss, it's usually 1.0.
        
        # We must now compute the gradient of `loss` w.r.t. `theta` numerically
        # using the finite difference method.
        theta_np = theta.numpy()
        grads_np = np.zeros_like(theta_np)
        epsilon = 1e-6  # A small step for numerical differentiation
        loss0 = get_variance_from_KRSrho(theta_np)

        # Calculate partial derivative for each element of theta
        for i in range(len(theta_np)):
            # Perturb theta[i] positively
            theta_plus = theta_np.copy()
            theta_plus[i] += epsilon
            loss_plus = get_variance_from_KRSrho(theta_plus)

            # Central difference formula for the partial derivative
            partial_grad = (loss_plus - loss0) / (2 * epsilon)
            grads_np[i] = partial_grad
            
        # The gradient must be a TensorFlow tensor of the same dtype as theta.
        # We also multiply by the upstream gradient `dy`.
        final_grad = tf.convert_to_tensor(grads_np, dtype=theta.dtype) * dy
        return final_grad

    # The custom_gradient function must return the result of the forward pass
    # and the gradient function.
    return loss, grad_fn


if __name__ == "__main__":
    print("#=========== RUN /home/gabriel/git/FireCore/tests/tEFF/AI_rough.py")
    if os.path.exists(fileToSaveProcess):
        os.remove(fileToSaveProcess) # deleting useless information
    eff.setVerbosity(0,0)
    atomParams = np.array([
    #  Q   sQ   sP   cP
    [ 0.,  1.0, 1.00, 0.0 ], # 0
    [ 1.,  0.0, 0.00, 0.0 ], # 1 H
    [ 0.,  1.0, 1.00, 1.0 ], # 2 He
    [ 1.,  0.0, 0.10, 1.0 ], # 3 Li
    [ 2.,  0.0, 0.10, 1.0 ], # 4 Be
    [ 3.,  0.0, 0.10, 1.0 ], # 5 B
    [ 4.,  0.0, 0.10, 1.0 ], # 6 C
    [ 5.,  0.0, 0.10, 1.0 ], # 7 N
    [ 6.,  0.0, 0.2, 1.0 ], # 8 O
    [ 7.,  0.0, 0.10, 1.0 ], # 9 F
    ], dtype=np.float64)
    eff.setAtomParams( atomParams )
    # print("set atom par")

    params, nrec = extract_blocks(elementPath)
    eff.initOpt( dt=0.005, damping=0.005, f_limit=1000.0)
    # print("init opt")
    bCoreElectrons = False
    eff.setSwitches( coreCoul=1 )
    eff.preAllocateXYZ(elementPath, Rfac=-1.35, bCoreElectrons=bCoreElectrons )
    eff.getBuffs()
    # print("get buffs")

    eff.info()
    eff.esize[:]=0.7
    
    theta0 = np.array([1.125, 0.9, -0.2])
    # theta0 = np.array([1.1, 0.7, -0.4])

    # Minimize by steps method
    # minTheta, steps, var, allEtot, flexVar, variance = minimizeBySteps(get_varianceFromKRSrho, theta0)

    # Dual annealing method
    # bounds = [
    # (0.1, 10.0),
    # (0.1, 5.0),
    # (-5.0, -0.1),]

    # results = dual_annealing(get_variance_from_KRSrho, bounds=bounds, maxiter=1000000, maxfun=1000)
    # minTheta = results.x
    # numOfGetVar = results.nfev
    # minVar = results.fun
    # funEval = results.nfev
    # print("======================================================================================================")
    # print(results.message)

    # ADAM optimizer method
    theta = tf.Variable(theta0, dtype=tf.float64)
    optimizer = tf.keras.optimizers.Adam(learning_rate=0.01)

    for step in range(1000):
        with tf.GradientTape() as tape:
            # Call the wrapper function which has the custom gradient defined
            loss = get_variance_wrapper(theta)
        grads = tape.gradient(loss, [theta])  # Compute gradients
        optimizer.apply_gradients(zip(grads, [theta]))  # Apply gradients

    if step % 10 == 0:
        print(f"Step {step}: x = {theta.numpy()}, loss = {loss.numpy():.6f}")
    minTheta = theta.numpy()
    minVar = loss.numpy()
    uncertainty = tf.abs(grads[0])


    eps = 1e-6
    grad = opt.approx_fprime(minTheta, get_variance_from_KRSrho, eps)
    grad_norm = np.linalg.norm(grad)

    # print("Grad vector: ", uncertainty)
    print("Grad vector 2: ", grad)


    angleArr = params['ang']
    distArr = params['dist']
    save_simulation(angleArr, distArr, minTheta, fileToSavePath)

    print("#===========")
    print("Theta final: ", minTheta)
    # print("Steps final: ", steps)
    print("Variance final: ", minVar)
    print("Number of get var: ", numOfGetVar)
    # print("Number of function evaluation", funEval)
    print("Uncertainty of final variable: ", grad_norm)
    print("#=========== DONE /home/gabriel/git/FireCore/tests/tEFF/AI_angdist_sim.py")