#!/usr/bin/env python3
import os
import re
import glob
import shutil
import subprocess
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

    return (os.path.basename(file).replace(".f90", ""), " ".join(module_comments).replace("_", r"\_"))

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

FILES_THAT_DONT_WORK = [
]

def create_fortran_modules_rst(path_to_src: str) -> list[str]:
    """
    Create the rst files containing the tables with the Fortran modules.

    :param path_to_src: Relative or absolute path to the Pencil Code src/ directory.
    :param output_file: Output rst file. Will be overwritten if it exists.
    :return: List of f90 files, to be used in conf.py as fortran_src.
    """
    retval = []
    modules = {}
    modules["src"] = process_directory(path_to_src)
    for subdir in sorted([p for p in os.scandir(path_to_src) if p.is_dir()], key=lambda x:x.name):
        table = process_directory(subdir.path)
        if not table:
            continue
        modules[subdir.name] = table
    total_modules = sum(len(val) for val in modules.values())

    dirshow = lambda module: "src" if module == "src" else f"src/{module}"
    dirtitle = lambda module: f"Main source files ({dirshow(module)})" if module == "src" else f"{module} source files ({dirshow(module)})"

    # Create the root directory for sourceFortran
    DOCSROOT = os.path.dirname(__file__)
    F90ROOT = os.path.join(DOCSROOT, "code/sourceFortran")
    shutil.rmtree(F90ROOT, ignore_errors=True)
    os.makedirs(F90ROOT, exist_ok=True)
    # Make the automodule files
    parsed_f90 = []
    for dirname, table in modules.items():
        os.makedirs(os.path.join(F90ROOT, dirname), exist_ok=True)
        for module, _ in table:
            # Modules that won't compile - skip them
            # Changed to only do this files for testing.
            if f"{dirname}/{module}.f90" in FILES_THAT_DONT_WORK:
                continue
            with open(os.path.join(F90ROOT, dirname, f"{module}.rst"), "w") as f:
                d = RstCloth(f)
                d.title(module)
                d.newline()
                d.directive("f:autosrcfile", arg=f"{module}.f90")
                actual_path = os.path.join(os.path.dirname(path_to_src), dirshow(dirname), f"{module}.f90")
                if not os.path.isfile(actual_path):
                    import sys
                    print(f"Attention: file {actual_path} does not exist; this is probably a bug in fortran_rst_generator.py!")
                    sys.exit(1)
                retval.append(actual_path)
                parsed_f90.append(f"{dirname}/{module}")
    # Main index
    with open(os.path.join(F90ROOT, "index.rst"), "w") as f:
        d = RstCloth(f)
        d.ref_target('fortran_modules')
        d.title("Fortran modules")
        d.content(f"Currently, the Pencil Code contains {total_modules} Fortran files.")
        d.newline()
        d.directive("raw", "html", content="<div>Filter: <input type='text' id='custommodsearch' /></div><br/>")
        d.newline()
        for dirname, table in modules.items():
            d.ref_target(dirname)
            d.newline()
            d.h2(dirtitle(dirname))
            d.newline()
            nmodules = len(table)
            d.content(f"The *{dirshow(dirname)}* directory contains {nmodules} Fortran file{'s' if nmodules != 1 else ''}.")
            table_content = []
            for it in table:
                if f"{dirname}/{it[0]}" in parsed_f90:
                    table_content.append((f":doc:`{dirname}/{it[0]}`", it[1]))
                else:
                    table_content.append(it)
            d.newline()
            d.table_list(headers=["File", "Description"], data=table_content, widths=[25,75])
            d.directive("toctree", fields=[("hidden", ""), ("maxdepth", "1")], content=[f"{dirname}/{it[0]}" for it in table])
            d.newline()
        d.newline()
        d.directive("toctree", fields=[("hidden", ""), ("titlesonly", ""), ("maxdepth", "2")], content=[f":ref:`{key}`" for key in modules.keys()])

    return retval



def process_diag(marker: str, in_file: str):
    files = [it.replace("../../src/", "") for it in glob.glob("../../src/**/*.f90", recursive=True)]
    files = sorted(files, key=lambda x: (x.count("/"), x))
    vars = {}
    total_vars = 0
    empty_vars = 0
    unique_vars = set()
    used_times = {}
    for file in files:
        if not os.path.isfile(f"../../src/{file}"):
            continue
        try:
            lines = subprocess.check_output(["grep", "-R", f" {marker}", f"../../src/{file}"]).decode().splitlines()
        except subprocess.CalledProcessError:
            continue
        filevars = {}
        current_var = None
        current_comment = ""
        for line in lines:
            if matches := re.findall(r":: idiag_(.*?)\s*=", line):
                if current_var:
                    filevars[current_var] = current_comment
                    unique_vars.add(current_var)
                    total_vars += 1
                    if current_var in used_times:
                        used_times[current_var] += 1
                    else:
                        used_times[current_var] = 1
                    if not current_comment.strip():
                        empty_vars += 1
                    current_var = None
                    current_comment = ""
                current_var = matches[0]
            if matches := re.findall(rf"!\s*{marker}: (.*)", line):
                current_comment += matches[0].strip()
        if current_var:
            filevars[current_var] = current_comment
            unique_vars.add(current_var)
            total_vars += 1
            if not current_comment.strip():
                empty_vars += 1
            if current_var in used_times:
                used_times[current_var] += 1
            else:
                used_times[current_var] = 1
            current_var = None
            current_comment = ""
        if filevars:
            vars[file] = filevars
    frac_undoc = int(empty_vars*100/total_vars)
    if empty_vars > 0 and frac_undoc == 0:
        frac_undoc = 1 # Don't show 0 if there is at least 1
    with open(f"code/tables/{in_file}.rst", "w") as f:
        d = RstCloth(f, line_width=5000)
        d.title(f"List of parameters for ``{in_file}``")
        d.content(f"This page lists {total_vars} variables distributed into {len(vars)} files. Of these:")
        d.newline()
        d.li(f"{total_vars - empty_vars} ({100 - frac_undoc}%) are documented;")
        d.li(f"{empty_vars} ({frac_undoc}%) are undocumented.")
        d.newline()
        d.content(f"Some of the variable names are shared amongst modules, so there are {len(unique_vars)} unique names:")
        d.newline()
        value_counts = {}
        for v in used_times.values():
            value_counts[v] = value_counts.get(v, 0) + 1
        # Sort by the value (the dict key)
        value_counts = dict(sorted(value_counts.items()))
        nvc = len(value_counts)
        for i, (k, v) in enumerate(value_counts.items()):
            d.li(f"{v} variable{'s' if v != 1 else ''} appear{'s' if v == 1 else ''} {k} time{'s' if k != 1 else ''}{';' if i < nvc - 1 else '.'}")
        d.newline()
        d.directive("raw", "html", content="<div>Filter: <input type='text' id='customvarsearch' /></div><br/>")
        d.newline()
        for file, filevars in vars.items():
            d.h3(f"Module *{file}*")
            table = []
            for k, v in filevars.items():
                # MathJax is more picky than LaTeX with spaces and formatting. Fix some common issues...
                vv = re.sub(r"\$(.*?)\$", r":math:`\1`", v)
                vv = vv.replace("\\rm", "\\mathrm")
                vv = vv.replace("\\quad", ":math:`\\quad` ")
                vv = re.sub(r":math:`(.*?)`", r":math:`\1` ", vv)
                vv = re.sub(r":math:` ", r":math:`", vv)
                vv = re.sub(r"(\w):math:", r"\1 :math:", vv)
                vv = vv.replace(r"} `", r"}`")
                table.append((f"*{k}*", vv))
            d.table_list(["Variable", "Meaning"], data=table, widths=[25, 75])

def process_all_idiag():
    diag_list = [
        ("DIAG_DOC", "print.in"),
        ("PHIAVG_DOC", "phiaver.in"),
        ("XYAVG_DOC", "xyaver.in"),
        ("XZAVG_DOC", "xzaver.in"),
        ("YZAVG_DOC", "yzaver.in"),
        ("YAVG_DOC", "yaver.in"),
        ("ZAVG_DOC", "zaver.in")
    ]
    os.makedirs("code/tables", exist_ok=True)
    with open(f"code/tables/index.rst", "w") as f:
        d = RstCloth(f, line_width=5000)
        d.title(f"Startup and run-time parameters")
        for pars in diag_list:
            process_diag(*pars)
        d.directive("toctree", fields=[("maxdepth", "1")], content=[pars[1] for pars in diag_list])
