import os
import subprocess
import argparse
from datetime import datetime

def get_git_info(filepath):
    """Fetches the last commit date for a specific file."""
    try:
        # Get last commit date in ISO format
        cmd = ["git", "log", "-1", "--format=%ai", "--", filepath]
        result = subprocess.check_output(cmd, stderr=subprocess.DEVNULL).decode('utf-8').strip()
        if not result:
            return "Not Tracked"
        # Return only YYYY-MM-DD HH:MM
        return result[:16]
    except subprocess.CalledProcessError:
        return "Not in Repo"

def get_file_stats(filepath):
    """Fetches line count, character count, and system birth time."""
    # Line and Char count
    line_count = 0
    char_count = 0
    try:
        with open(filepath, 'rb') as f:
            content = f.read()
            char_count = len(content)
            line_count = content.count(b'\n')
    except Exception:
        pass

    # Creation Time (Birth time)
    try:
        stat = os.stat(filepath)
        # st_birthtime is available on macOS/BSD; on Linux it's harder to get via Python
        # so we fallback to st_mtime (modification) or a placeholder
        creation_ts = getattr(stat, 'st_birthtime', stat.st_ctime)
        creation_date = datetime.fromtimestamp(creation_ts).strftime('%Y-%m-%d %H:%M')
    except Exception:
        creation_date = "Unknown"

    return line_count, char_count, creation_date

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="List file modification dates, git history, and sizes.")
    parser.add_argument("path", nargs="?", default=".", help="Directory to scan (default: current)")
    args = parser.parse_args()

    target_dir = args.path
    
    # Header
    print(f"{'Filename':<25} | {'Created (Sys)':<17} | {'Last Git Edit':<17} | {'Lines':<7} | {'Chars':<8}")
    print("-" * 85)

    try:
        files = sorted([f for f in os.listdir(target_dir) if os.path.isfile(os.path.join(target_dir, f))])
    except OSError as e:
        print(f"Error: {e}")
        exit(1)

    for filename in files:
        full_path = os.path.join(target_dir, filename)
        
        git_date = get_git_info(full_path)
        lines, chars, created = get_file_stats(full_path)

        print(f"{filename[:25]:<25} | {created:<17} | {git_date:<17} | {lines:<7} | {chars:<8}")