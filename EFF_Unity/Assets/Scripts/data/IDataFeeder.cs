using System;
using UnityEngine;



public interface IDataFeeder
{
    
    (int electronCount, int atomCount, int[] spins) Init(string fileName, Func<(int electronCount, int atomCount)> countFetcher);
    (Vector3[] positions, float[] sizes) GetNextFrame();
    
}
