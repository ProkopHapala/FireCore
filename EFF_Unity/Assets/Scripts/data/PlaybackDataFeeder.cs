using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using System.Text;
using UnityEngine;


public class PlaybackDataFeeder : IDataFeeder
{
    
    private Func<(int electronCount, int atomCount)> countFetcher;
    private int electronCount => countFetcher().electronCount;
    private int atomCount => countFetcher().atomCount;
    
    private StreamReader reader;
    private int frameLength;

    private int frameCounter = 0;
    private int skipFrames = 5;
    private int skipFramesCounter;
    
    private (Vector3[] positions, float[] sizes) lastFrame;
    
    public (int electronCount, int atomCount, int[] spins) Init(string fileName, Func<(int electronCount, int atomCount)> countFetcher)
    {
        this.countFetcher = countFetcher;
        reader = new StreamReader(fileName);
        reader.ReadLine(); // skip #playback
        IEnumerable<string> lines = File.ReadLines(fileName).Skip(1);

        string[] lineSplit = SplitLine(lines.First(), 2);
        int aCount = int.Parse(lineSplit[0]);
        int eCount = int.Parse(lineSplit[1]);
        frameLength = aCount + eCount + 1;
        
        string[] electrons = lines.Skip(1 + aCount).Take(eCount).Select(e => e.Trim()).ToArray();
        List<int> spins = new List<int>();
        Array.ForEach(electrons, e => spins.Add(e.Substring(e.LastIndexOf(' ') + 1) == "-1" ? -1 : 1));

        skipFramesCounter = skipFrames;
        
        return (eCount, aCount, spins.ToArray());
    }

    public (Vector3[] positions, float[] sizes) GetNextFrame()
    {
        // Debug.Log(++frameCounter);
        if (skipFramesCounter < skipFrames) {
            skipFramesCounter++;
            return lastFrame;
        }
        skipFramesCounter = 0;
        
        reader.ReadLine(); // skip counts
        Vector3[] positions = new Vector3[electronCount + atomCount];
        float[] sizes = new float[electronCount];
        
        for (int i = 0; i < atomCount; i++)
        {
            float[] line = SplitLine(reader.ReadLine(), 3).Select(s => float.Parse(s)).ToArray();
            positions[electronCount + i] = new Vector3(line[0], line[1], line[2]);
        }
        
        for (int i = 0; i < electronCount; i++)
        {
            float[] line = SplitLine(reader.ReadLine(), 4).Select(s => float.Parse(s)).ToArray();
            positions[i] = new Vector3(line[0], line[1], line[2]);
            sizes[i] = line[3];
        }
        lastFrame = (positions, sizes);
        return lastFrame;

    }

    private string[] SplitLine(string line, int cap = -1) {
        StringBuilder sb = new();
        List<string> output = new();
        
        foreach (char c in line) {
            if (c == ' ') {
                if (sb.Length == 0) {
                    continue;
                }
                output.Add(sb.ToString());
                if (output.Count == cap && cap != -1) {
                    break;
                }
                sb = new();
            }
            else {
                sb.Append(c);
            }
        }
        return output.ToArray();
    }
}
