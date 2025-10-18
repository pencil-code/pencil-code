#!/usr/bin/env python3
import os
import glob
from rstcloth import RstCloth

def process_file(file) -> tuple[str, str]:
    """
    Process the given f90 file, and extract the module comment from it.

    The script looks at comments both before and after the 'module' declaration,
    but stops at the first Fortran instruction that is not 'module' (e.g., 'use', 'implicit', etc.).

    The script disregards initial blank lines, CVS $Id$ strings, separator lines such as "------", "!!!!!!", etc.
    BUT once it has detected at least one valid comment line, it will stop at the first blank/separator line.
    In other words, only the first paragraph in the comments is returned.

    If there are lines tagged with MODULE_DOC, it will return those and only those lines.

    Lines starting with ** or with CPARAM are ignored, and all lines following those are ignored as well.

    :file: Path to the f90 file.
    :return: A tuple containing the module name (assumed to be the same as the file, without the '.f90' extension)
        and the module description string extracted from the file. If the file does not contain a 'module'
        declaration, the comment will be 'Not a module'.
    """
    module_comments = [] # Lines of module comments
    module_started = False # We encountered the "module NAME" declaration
    has_module_doc = False
    in_comments = True # We are still adding lines to the module_comments list
    with open(file, "r") as fp:
        for line_raw in fp:
            line = line_raw.strip()
            if line.lower().startswith("module "):
                module_started = True
                if not in_comments:
                    break
                continue
            if not line:
                if module_comments:
                    in_comments = False
                continue
            if line.startswith("!!"):
                continue
            if line.startswith("!"):
                if not in_comments:
                    continue
                comment_line = line[1:].replace('`', "'").strip()
                if comment_line.startswith("$Id") or comment_line.startswith("--"):
                    continue
                if comment_line.startswith("MODULE_DOC:") or comment_line.startswith("MOUDLE_DOC"):
                    has_module_doc = True
                    comment_line = comment_line[11:].strip()
                    module_comments.append(comment_line)
                    continue
                if has_module_doc:
                    in_comments = False
                    continue
                if not comment_line:
                    if not module_comments:
                        continue
                    else:
                        in_comments = False
                        continue
                if comment_line.startswith("**") or comment_line.startswith("CPARAM"):
                    in_comments = False
                    continue
                module_comments.append(comment_line)
                continue
            if in_comments:
                in_comments = False
            continue
    if not module_started:
        module_comments = ["Not a module"]

    return ("``" + os.path.basename(file).replace(".f90", "") + "``", " ".join(module_comments).replace("_", r"\_"))

def process_directory(path) -> list:
    """
    Process all f90 files in the given directory.

    :param path: Path to the directory.
    :return: Tuple containing the directory (e.g., ``src``, ``src/experimental``)
        and the list of processed files (see output of ``process_file``).
    """
    table = []
    for file in sorted(glob.glob(os.path.join(path, "*.f90")), key=str.casefold):
        table.append(process_file(file))
    return table

def create_fortran_modules_rst(path_to_src: str,
                               output_file: str) -> None:
    """
    Create the rst file containing the tables with the Fortran modules.

    :param path_to_src: Relative or absolute path to the Pencil Code src/ directory.
    :param output_file: Output rst file. Will be overwritten if it exists.
    """
    modules = {}
    modules["src"] = process_directory(path_to_src)
    for subdir in sorted([p for p in os.scandir(path_to_src) if p.is_dir()], key=lambda x:x.name):
        table = process_directory(subdir.path)
        if not table:
            continue
        modules[subdir.name] = table
    total_modules = sum(len(val) for val in modules.values())

    with open(output_file, "w") as f:
        d = RstCloth(f)
        d.title("Fortran modules")
        d.content(f"Currently, the Pencil Code contains {total_modules} Fortran files.")
        d.newline()
        for dirname, table in modules.items():
            if dirname == "src":
                dirshow = "src"
                d.h2(f"Main source files ({dirshow})")
            else:
                dirshow = f"src/{dirname}"
                d.h2(f"*{dirname}* source files ({dirshow})")
            nmodules = len(table)
            d.content(f"The *{dirshow}* directory contains {nmodules} Fortran file{'s' if nmodules != 1 else ''}.")
            d.table(header=["File", "Description"], data=table)
