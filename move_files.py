#!/usr/bin/env python3
import os, shutil, argparse
from pathlib import Path

def find_files(src_dir, suffixes, ext_check=['.h','.hpp','.c','.cpp']):
    """Return list of paths ending with any suffix in suffixes and with a matching code file in same dir"""
    matches = []
    bCheckExt = ext_check is not None
    for root, dirs, files in os.walk(src_dir):
        for f in files:
            #print(f)
            for suf in suffixes:
                if f.endswith(suf):
                    #print(f)
                    base = f[:-len(suf)]  # remove suffix to get base name
                    if bCheckExt:
                        # check for corresponding code file
                        if not any((Path(root)/(base+ext)).exists() for ext in ext_check):
                            continue
                    matches.append(Path(root)/f)
                    break
    return matches

def move_files(files, src_dir, dest_dir):
    """Move each file to dest_dir, recreating directory structure"""
    for src_path in files:
        rel = src_path.relative_to(src_dir)
        dest_path = Path(dest_dir)/rel
        dest_path.parent.mkdir(parents=True, exist_ok=True)
        shutil.move(str(src_path), str(dest_path))
        print(f"Moved {src_path} -> {dest_path}")

def main():
    p = argparse.ArgumentParser(description='List or move markdown skeleton files matching code files')
    p.add_argument('--src',      default='.',                type=str, help='Root directory to search')
    p.add_argument('--suffixes', default='.skeleton.md,.md', type=str, help='Comma-separated suffixes to match')
    p.add_argument('--dest',     default=None,               type=str, help='If set, move found files to this directory')
    args = p.parse_args()

    suffixes = [s.strip() for s in args.suffixes.split(',') if s.strip()]
    files = find_files(args.src, suffixes, ext_check=None )
    if not files:
        print('No matching files found.')
        return
    print('Found files:')
    for f in files:
        print(f)
    if args.dest:
        move_files(files, args.src, args.dest)

if __name__=='__main__':
    main()
