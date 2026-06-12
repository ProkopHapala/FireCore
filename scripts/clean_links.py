#!/usr/bin/env python3
"""Clean up markdown links in topical audit files.

Converts IDE-style links to clean repo-relative paths with optional line ranges:
  [text](cci:1://file:///home/prokop/.../path:230:0-284:17)  -> [text](path:230-284)  (default)
  [text](cci:1://file:///home/prokop/.../path:230:0-284:17)  -> [text](path)          (--no-lines)
  [text](cci:7://file:///home/prokop/.../path:0:0-0:0)       -> [text](path)
  (1234 lines)                                              -> (removed)
"""

import argparse
import re
import sys
import pathlib


def clean_line(line: str, keep_lines: bool = True) -> str:
    """Rewrite a single line of markdown."""

    # 1. Strip leading absolute-path prefix from cci: links, preserving line ranges
    def repl_link(m: re.Match) -> str:
        text = m.group(1)
        path = m.group(2)           # repo-relative path
        start_line = m.group(3)     # e.g. '230'
        end_line   = m.group(5)     # e.g. '284'

        # If line numbers are 0:0-0:0, drop them entirely
        if start_line == '0' and end_line == '0':
            return f'[{text}]({path})'
        if keep_lines:
            return f'[{text}]({path}:{start_line}-{end_line})'
        return f'[{text}]({path})'

    # Match: [anything](cci:<num>://file:///home/prokophapala/git/FireCore/<path>:<sl>:<sc>-<el>:<ec>)
    line = re.sub(
        r'\[([^\]]+)\]\(cci:[0-9]+://file:///home/prokophapala/git/FireCore/([^:]+):([0-9]+):([0-9]+)-([0-9]+):([0-9]+)\)',
        repl_link,
        line,
    )

    # 2. Also catch the shorter form without line ranges (cci:9://...:0:0-0:0 used for dirs/images)
    line = re.sub(
        r'\[([^\]]+)\]\(cci:[0-9]+://file:///home/prokophapala/git/FireCore/([^:]+)\)',
        r'[\1](\2)',
        line,
    )

    # 3. Remove (N lines) annotations
    line = re.sub(r' \([0-9]+ lines\)', '', line)

    return line


def main():
    parser = argparse.ArgumentParser(
        description='Clean IDE-style markdown links to repo-relative paths.'
    )
    parser.add_argument('file', type=pathlib.Path, help='Markdown file to clean in-place')
    parser.add_argument('--no-lines', action='store_true',
                        help='Strip line-number ranges from links (default: keep them)')
    parser.add_argument('--dry-run', action='store_true',
                        help='Print cleaned output instead of writing file')
    args = parser.parse_args()

    if not args.file.exists():
        print(f"File not found: {args.file}")
        sys.exit(1)

    text = args.file.read_text()
    cleaned = ''.join(clean_line(l, keep_lines=not args.no_lines)
                      for l in text.splitlines(keepends=True))

    if args.dry_run:
        print(cleaned, end='')
    else:
        args.file.write_text(cleaned)
        print(f"Cleaned {args.file} (keep_lines={not args.no_lines})")


if __name__ == '__main__':
    main()
