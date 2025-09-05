#!/usr/bin/env python3
"""
xyz_simplify_symbols.py

Process .xyz trajectory/movie files where atom names carry type suffixes like "C_2", "H_OH", etc.
For each atom line, replace the first column with the plain element symbol (the substring before the first underscore),
leaving all other content and whitespace unchanged. Non-atom lines (e.g., frame size, comment starting with '#') are left unchanged.

- Minimal dependencies: stdlib only
- Does not parse floats to avoid any precision/formatting changes; works with strings and preserves whitespace

Usage examples:
  # Process a single file to an output directory
  python scripts/xyz_simplify_symbols.py -i path/to/input.xyz -o out_dir

  # Process all .xyz files in a directory (non-recursive)
  python scripts/xyz_simplify_symbols.py -i path/to/dir -o out_dir

  # Recursive with custom glob
  python scripts/xyz_simplify_symbols.py -i path/to/dir -o out_dir -r --glob "**/*.xyz"

  # Overwrite in place (be careful!)
  python scripts/xyz_simplify_symbols.py -i path/to/dir --in-place
"""
from __future__ import annotations

import argparse
import os
import sys
from pathlib import Path
import contextlib
from typing import Iterable, Tuple, List, Optional


def simplify_line(line: str) -> str:
    """Replace the first token like "C_2" -> "C" only if the token contains an underscore.

    Preserve leading and intra-token whitespace and the remainder of the line exactly.
    If the line is blank or starts with a non-atom token (e.g., '#', '59'), return unchanged.
    """
    if not line:
        return line

    # Find leading whitespace
    i = 0
    L = len(line)
    while i < L and line[i].isspace():
        i += 1

    # If the line is all whitespace, or empty after whitespace
    if i >= L:
        return line

    # Identify first token (non-whitespace span)
    j = i
    while j < L and not line[j].isspace():
        j += 1

    first = line[i:j]

    # Leave comment or size lines untouched (e.g., starts with '#', or digit)
    if not first:
        return line
    if first[0] == '#':
        return line
    if first[0].isdigit():
        return line

    # Only act if token contains an underscore and starts with an alpha (element-like)
    if '_' in first and first[0].isalpha():
        base = first.split('_', 1)[0]
        # Reconstruct preserving whitespace exactly
        return f"{line[:i]}{base}{line[j:]}"

    return line


def process_stream(instream: Iterable[str], outstream) -> Tuple[int, int]:
    """Process lines from instream, write to outstream.

    Returns (lines_total, lines_modified).
    """
    total = 0
    modified = 0
    for line in instream:
        total += 1
        new_line = simplify_line(line)
        if new_line is not line and new_line != line:
            modified += 1
        outstream.write(new_line)
    return total, modified


def parse_int_token(s: str) -> Optional[int]:
    """Parse an integer from a line that is expected to be the atom count.
    Returns None if parsing fails."""
    try:
        return int(s.strip())
    except Exception:
        return None


def iterate_frames(instream: Iterable[str]) -> Iterable[Tuple[str, str, List[str]]]:
    """Yield XYZ frames as (count_line, comment_line, atom_lines).

    - Preserves lines exactly as read
    - Attempts to follow the conventional order: count line, comment line, N atom lines
    - If structure is inconsistent, it will try to resynchronize by scanning forward
    """
    it = iter(instream)
    for line in it:
        # Seek count line
        n = parse_int_token(line)
        if n is None:
            # Not a valid frame start, skip this line (write responsibility is with caller)
            continue
        count_line = line
        # Expect comment line next; if EOF, yield incomplete frame conservatively
        try:
            comment_line = next(it)
        except StopIteration:
            yield count_line, "", []
            break
        # Collect atom lines
        atoms: List[str] = []
        for _ in range(n):
            try:
                atoms.append(next(it))
            except StopIteration:
                # Incomplete frame; yield what we have
                break
        yield count_line, comment_line, atoms


def detect_direction(comment_line: str) -> Optional[str]:
    """Detect scan direction token ('y' or 'z') from the comment line.
    The token is expected as a standalone word in whitespace-separated tokens.
    Returns 'y', 'z', or None if not found.
    """
    tokens = comment_line.strip().split()
    for t in tokens:
        if t == 'y' or t == 'Y':
            return 'y'
        if t == 'z' or t == 'Z':
            return 'z'
    return None


def iter_xyz_files(root: Path, recursive: bool, glob: str) -> Iterable[Path]:
    if root.is_file():
        yield root
        return
    if recursive:
        yield from root.rglob(glob)
    else:
        yield from root.glob(glob)


def ensure_out_path(in_file: Path, out_dir: Path | None, suffix: str | None, in_place: bool) -> Path:
    if in_place:
        return in_file
    if out_dir is None:
        # Default to sibling directory named "converted" if not provided
        out_dir = in_file.parent / "converted"
    out_dir.mkdir(parents=True, exist_ok=True)
    name = in_file.name
    if suffix:
        stem = in_file.stem
        ext = ''.join(in_file.suffixes)  # keeps compound suffixes if any
        name = f"{stem}{suffix}{ext}"
    return out_dir / name


def out_paths_for_split(in_file: Path, out_dir: Optional[Path], suffix: str) -> Tuple[Path, Path]:
    """Return output paths for -y and -z split files based on input filename.
    Filenames are of the form '<stem>-y{suffix}<ext>' and '<stem>-z{suffix}<ext>'.
    """
    if out_dir is None:
        out_dir = in_file.parent / "converted"
    out_dir.mkdir(parents=True, exist_ok=True)
    stem = in_file.stem
    ext = ''.join(in_file.suffixes)
    yname = f"{stem}-y{suffix}{ext}"
    zname = f"{stem}-z{suffix}{ext}"
    return out_dir / yname, out_dir / zname


def main(argv: Iterable[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description="Utilities for .xyz movies: simplify atom types (e.g., C_2 -> C) and/or split frames by scan direction (y/z). All formatting preserved.")
    ap.add_argument("-i", "--input", required=True, help="Input .xyz file or directory containing .xyz files")
    ap.add_argument("-o", "--output", help="Output directory. If not provided with --in-place, defaults to '<input_dir>/converted'")
    ap.add_argument("-r", "--recursive", action="store_true", help="Recurse into subdirectories when input is a directory")
    ap.add_argument("--glob", default="*.xyz", help="Glob to match files inside the input directory (default: *.xyz). Use with -r for patterns like **/*.xyz")
    ap.add_argument("--suffix", default="", help="Optional suffix to append to output filenames before extension (e.g., '_elm')")
    ap.add_argument("--in-place", action="store_true", help="Overwrite input files in place (use with care)")
    ap.add_argument("--dry-run", action="store_true", help="Do not write files, just report what would change")
    ap.add_argument("-v", "--verbose", action="store_true", help="Verbose logging")
    ap.add_argument("--simplify", dest="simplify", action="store_true", help="Simplify atom types in the first column (e.g., H_OH -> H)")
    ap.add_argument("--no-simplify", dest="simplify", action="store_false", help="Do not simplify atom types")
    # Tri-state: None means not specified; we'll derive default behavior after parsing
    ap.set_defaults(simplify=None)
    ap.add_argument("--split", action="store_true", help="Split frames into separate '-y' and '-z' files based on the comment line token")

    args = ap.parse_args(list(argv) if argv is not None else None)

    in_path = Path(args.input)
    out_dir = Path(args.output).resolve() if args.output else None

    if args.in_place and out_dir is not None:
        print("[WARN] --in-place specified; --output will be ignored", file=sys.stderr)
    if args.split and args.in_place:
        print("[ERROR] --split cannot be combined with --in-place (splitting writes multiple outputs)", file=sys.stderr)
        return 2

    # Decide final behavior: if simplify not explicitly set, default is:
    #  - simplify when not splitting
    #  - do not simplify when splitting (user can add --simplify to enable)
    do_simplify = (args.simplify if args.simplify is not None else (not args.split))

    # Collect files
    if not in_path.exists():
        print(f"[ERROR] Input path does not exist: {in_path}", file=sys.stderr)
        return 2

    files = list(iter_xyz_files(in_path, args.recursive, args.glob))
    if not files:
        print(f"[WARN] No files matched. input={in_path} recursive={args.recursive} glob='{args.glob}'", file=sys.stderr)
        return 1

    total_files = 0
    total_lines = 0
    total_mod = 0

    for f in files:
        if not f.is_file():
            continue
        total_files += 1
        if args.split:
            y_path, z_path = out_paths_for_split(f, out_dir, args.suffix)
            if args.verbose or args.dry_run:
                print(f"Processing (split): {f} -> {y_path} | {z_path}")

            lines = 0
            mods = 0
            counts = {"y": 0, "z": 0, None: 0}
            if args.dry_run:
                with f.open('r', encoding='utf-8', errors='replace') as fin:
                    for cnt, com, atoms in iterate_frames(fin):
                        lines += 2 + len(atoms)
                        direction = detect_direction(com)
                        counts[direction] = counts.get(direction, 0) + 1
                        if do_simplify:
                            # Count potential modifications in atom lines
                            for al in atoms:
                                if simplify_line(al) != al:
                                    mods += 1
                total_lines += lines
                total_mod += mods
                if args.verbose:
                    print(f"  frames_y={counts.get('y',0)} frames_z={counts.get('z',0)} frames_other={counts.get(None,0)}")
                    print(f"  lines={lines} modified_atom_lines={mods}")
                continue

            # Write to two files
            with contextlib.ExitStack() as stack:
                fin = stack.enter_context(f.open('r', encoding='utf-8', errors='replace'))
                fy = stack.enter_context(y_path.open('w', encoding='utf-8', errors='strict'))
                fz = stack.enter_context(z_path.open('w', encoding='utf-8', errors='strict'))
                for cnt, com, atoms in iterate_frames(fin):
                    direction = detect_direction(com)
                    target = fy if direction == 'y' else fz if direction == 'z' else None
                    if target is None:
                        # Skip frames with no detectable direction
                        continue
                    # Write frame preserving text; optionally simplify atom lines
                    target.write(cnt)
                    target.write(com)
                    for al in atoms:
                        new_al = simplify_line(al) if do_simplify else al
                        if new_al is not al and new_al != al:
                            mods += 1
                        target.write(new_al)
                        lines += 1
                    # account for header lines
                    lines += 2
            total_lines += lines
            total_mod += mods
            if args.verbose:
                print(f"  lines={lines} modified_atom_lines={mods}")
        else:
            out_path = ensure_out_path(f, out_dir, args.suffix, args.in_place)
            if args.verbose or args.dry_run:
                target = "(in-place)" if args.in_place else str(out_path)
                print(f"Processing: {f} -> {target}")

            if not args.dry_run:
                out_path.parent.mkdir(parents=True, exist_ok=True)

            with f.open('r', encoding='utf-8', errors='replace') as fin:
                if args.dry_run:
                    # Count modifications without writing output
                    lines = 0
                    mods = 0
                    if do_simplify:
                        for line in fin:
                            lines += 1
                            if simplify_line(line) != line:
                                mods += 1
                    else:
                        for line in fin:
                            lines += 1
                    total_lines += lines
                    total_mod += mods
                    if args.verbose:
                        if do_simplify:
                            print(f"  lines={lines} modified={mods}")
                        else:
                            print(f"  lines={lines} (no simplification)")
                    continue

                # Write transformed or identical content
                with out_path.open('w', encoding='utf-8', errors='strict') as fout:
                    if do_simplify:
                        lines, mods = process_stream(fin, fout)
                    else:
                        # Copy as-is
                        lines = 0
                        mods = 0
                        for line in fin:
                            lines += 1
                            fout.write(line)
                total_lines += lines
                total_mod += mods
                if args.verbose:
                    if do_simplify:
                        print(f"  lines={lines} modified={mods}")
                    else:
                        print(f"  lines={lines} (no simplification)")

    print(f"Done. files={total_files} lines={total_lines} modified_lines={total_mod}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
