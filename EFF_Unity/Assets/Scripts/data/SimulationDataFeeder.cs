using System;
using System.Runtime.InteropServices;
using UnityEngine;


public class SimulationDataFeeder : IDataFeeder
{
    private Func<(int electronCount, int atomCount)> countFetcher;
    private int electronCount => countFetcher().electronCount;
    private int atomCount => countFetcher().atomCount;
    
    public (int electronCount, int atomCount, int[] spins) Init(string fileName, Func<(int electronCount, int atomCount)> countFetcher)
    {
        this.countFetcher = countFetcher;
        var dataPtr = SimulationManager.unityInit(fileName);
        // First get the electron and atom counts
        int eCount = Marshal.ReadInt32(dataPtr, 0);
        int aCount = Marshal.ReadInt32(dataPtr, sizeof(int));
    
        // Create and populate the espins array
        var spins = new int[eCount];
        for (int i = 0; i < eCount; i++) {
            spins[i] = Marshal.ReadInt32(dataPtr, (2 + i) * sizeof(int));
        }
    
        // Free the unmanaged memory
        SimulationManager.cleanupInitData(dataPtr);

        return (eCount, aCount, spins);
    }
    
    public (Vector3[] positions, float[] sizes) GetNextFrame()
    {
        int size;
        IntPtr positionsPtr = SimulationManager.unityNextFrame(out size);
    
        // Convert to managed array
        float[] rawPositions = new float[size];

        Marshal.Copy(positionsPtr, rawPositions, 0, size);
    
        // Clean up native memory
        SimulationManager.cleanupPositions(positionsPtr);
    
        Vector3[] positions = new Vector3[atomCount + electronCount];
        float[] sizes = new float[electronCount];
    
        int idx = 0;
        // First, extract electron positions
        for (int i = 0; i < electronCount; i++)
        {
            positions[i] = new Vector3(
                rawPositions[idx++],
                rawPositions[idx++],
                rawPositions[idx++]
            );
        }
    
        // Then extract atom positions
        for (int i = 0; i < atomCount; i++)
        {
            positions[electronCount + i] = new Vector3(
                rawPositions[idx++],
                rawPositions[idx++],
                rawPositions[idx++]
            );
        }
    
        // Finally extract electron sizes
        for (int i = 0; i < electronCount; i++) {
            sizes[i] = rawPositions[idx++];
        }
    
    
        return (positions, sizes);
    }
}
