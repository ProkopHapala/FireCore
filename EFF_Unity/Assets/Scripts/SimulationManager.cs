using System;
using System.Runtime.InteropServices;

public class SimulationManager
{
    public const string PATH_TO_EFF_APP = "/home/perry/FireCore/cpp/Build/apps/EFF/libEFFapp_console.so";

    [DllImport(PATH_TO_EFF_APP, CallingConvention = CallingConvention.Cdecl, CharSet = CharSet.Ansi)]
    public static extern IntPtr unityInit(string fileName);

    [DllImport(PATH_TO_EFF_APP, CallingConvention = CallingConvention.Cdecl, CharSet = CharSet.Ansi)]
    public static extern IntPtr unityNextFrame(out int size);

    [DllImport(PATH_TO_EFF_APP, CallingConvention = CallingConvention.Cdecl, CharSet = CharSet.Ansi)]
    public static extern void cleanupPositions(IntPtr positions);

    [DllImport(PATH_TO_EFF_APP, CallingConvention = CallingConvention.Cdecl, CharSet = CharSet.Ansi)]
    public static extern void cleanupInitData(IntPtr data);

    [DllImport(PATH_TO_EFF_APP, CallingConvention = CallingConvention.Cdecl, CharSet = CharSet.Ansi)]
    public static extern void unitySetElectronPosition(int electronIndex, float x, float y, float z, float size);

    [DllImport(PATH_TO_EFF_APP, CallingConvention = CallingConvention.Cdecl, CharSet = CharSet.Ansi)]
    public static extern void unitySetAtomPosition(int atomIndex, float x, float y, float z);

    [DllImport(PATH_TO_EFF_APP, CallingConvention = CallingConvention.Cdecl, CharSet = CharSet.Ansi)]
    public static extern void fixAtom(int atomIndex);

    [DllImport(PATH_TO_EFF_APP, CallingConvention = CallingConvention.Cdecl, CharSet = CharSet.Ansi)]
    public static extern void fixElectron(int atomIndex);

    [DllImport(PATH_TO_EFF_APP, CallingConvention = CallingConvention.Cdecl, CharSet = CharSet.Ansi)]
    public static extern void unfixAtom(int atomIndex);

    [DllImport(PATH_TO_EFF_APP, CallingConvention = CallingConvention.Cdecl, CharSet = CharSet.Ansi)]
    public static extern void unfixElectron(int atomIndex);

    [DllImport(PATH_TO_EFF_APP, CallingConvention = CallingConvention.Cdecl, CharSet = CharSet.Ansi)]
    public static extern IntPtr unityGetProtonNumbers();

    [DllImport(PATH_TO_EFF_APP, CallingConvention = CallingConvention.Cdecl, CharSet = CharSet.Ansi)]
    public static extern void cleanupProtonNumbers(IntPtr protonNumbers);

}