#!/usr/bin/env python3
"""
Scan a directory of Fortran files for non-ASCII characters.

Usage:
    python find_nonascii_dir.py /path/to/fortran/src
"""

import os
import sys

def find_nonascii_in_file(filepath):
    """Return a list of (line_no, col_no, char) for non-ASCII characters."""
    results = []
    with open(filepath, "r", encoding="utf-8", errors="replace") as f:
        for lineno, line in enumerate(f, start=1):
            for colno, char in enumerate(line, start=1):
                if ord(char) > 127:
                    results.append((lineno, colno, char))
    return results

def scan_directory(root_dir, extension=".f90"):
    """Scan all files with the given extension in a directory."""
    files_with_nonascii = 0
    # changed: gather all files first
    all_files = []
    for dirpath, _, filenames in os.walk(root_dir):
        for fname in filenames:
            if fname.lower().endswith(extension):
                all_files.append(os.path.join(dirpath, fname))  # changed

    # changed: sort files alphabetically
    all_files.sort()  # changed

    for full_path in all_files:
        nonascii = find_nonascii_in_file(full_path)
        if nonascii:
            files_with_nonascii += 1
            print(f"\nFile: {full_path}")
            for lineno, colno, char in nonascii:
                print(f"  Line {lineno}, Col {colno}: {repr(char)}")
    print("\nSummary:")
    print(f"  Total files with non-ASCII characters: {files_with_nonascii}")

if __name__ == "__main__":
    if len(sys.argv) != 2:
        print("Usage: python find_nonascii_dir.py /path/to/fortran/src")
        sys.exit(1)

    src_dir = sys.argv[1]
    scan_directory(src_dir)
