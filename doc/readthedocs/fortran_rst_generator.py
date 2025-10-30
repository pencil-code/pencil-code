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
    files = sorted(files, key=lambda x: (x.count("/") > 0, x.lower()))
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
                if current_comment:
                    current_comment += " "
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

def process_boundary_conditions():
    search_for = {
        "BCX_DOC": "Boundary condition bcx",
        "BCY_DOC": "Boundary condition bcy",
        "BCZ_DOC": "Boundary condition bcz"
    }
    search_keys = list(search_for.keys())
    estring = [f"-e {it}" for it in search_keys]
    files = subprocess.check_output(
        f"grep -rl {' '.join(estring)} ../../src/*.f90",
        shell=True,
        text=True
    ).splitlines()
    files = [it.replace("../../src/", "") for it in files]
    files = sorted(files, key=lambda x: (x.count("/") > 0, x.lower()))
    vars = {}
    current_var = None
    current_comment = ""
    current_block = None
    prev_block = None
    for file in files:
        with open(f"../../src/{file}", "r") as fh:
            lines = [it.strip() for it in fh.readlines()]
        nlines = len(lines)
        i = 0
        while i < nlines:
            for key in search_keys:
                if lines[i].startswith(f"! {key}:"):
                    # Find the beginning of the SELECT statement
                    j = i-1
                    while j >= 0 and not lines[j].startswith("select case "): j -= 1
                    # Loop until the end of the SELECT statement
                    while j < nlines and (not lines[j].startswith("endselect")) and (not lines[j].startswith("end select")):
                        if lines[j].startswith("case ("):
                            current_var = re.findall(r"case \((.*)\)", lines[j])[0]
                            if "," in current_var:
                                lst = [it.strip("'") for it in current_var.split(",")]
                                lst = ["''" if not it else it for it in lst]
                                current_var = ", ".join(lst)
                            else:
                                current_var = current_var.strip("'")
                        else:
                            found = False
                            for key in search_keys:
                                if lines[j].startswith(f"! {key}:"):
                                    found = True
                                    if current_block is None:
                                        current_block = key.split("_")[0]
                                    prev_block = current_block
                                    if matches := re.findall(rf"!\s*{key}: (.*)", lines[j]):
                                        if current_comment:
                                            current_comment += " "
                                        current_comment += matches[0].strip()
                            if not found:
                                if current_var:
                                    if current_block is None:
                                        current_block = prev_block # uncommented variables
                                    if current_block not in vars:
                                        vars[current_block] = {}
                                    if file not in vars[current_block]:
                                        vars[current_block][file] = {}
                                    vars[current_block][file][current_var] = current_comment
                                    current_var = None
                                    current_comment = ""
                                    current_block = None
                        j += 1
                    if current_var:
                        if current_block is None:
                            current_block = prev_block
                        if current_block not in vars:
                            vars[current_block] = {}
                        if file not in vars[current_block]:
                            vars[current_block][file] = {}
                        vars[current_block][file][current_var] = current_comment
                        current_var = None
                        current_comment = ""
                        current_block = None
                    i = j
                    continue
            i += 1
    with open(f"code/tables/boundary.rst", "w") as f:
        d = RstCloth(f, line_width=5000)
        d.title(f"Boundary conditions")
        d.content("The following tables list all possible boundary condition labels.")
        d.newline()
        d.directive("raw", "html", content="<div>Filter: <input type='text' id='custombcsearch' /></div><br/>")
        d.newline()
        for block, blockitems in vars.items():
            d.h2(f"Boundary conditions *{block.lower()}*")
            for file, fileitems in blockitems.items():
                d.h3(f"Module *{file}*")
                table = []
                for k, v in fileitems.items():
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

def separate_fortran_comment(line: str) -> str:
    """
    Remove Fortran comments while preserving literal '!' inside quotes.
    Supports single and double-quoted strings.
    """
    in_single = False
    in_double = False

    for i, ch in enumerate(line):
        if ch == "'" and not in_double:
            in_single = not in_single
        elif ch == '"' and not in_single:
            in_double = not in_double
        elif ch == "!" and not in_single and not in_double:
            # True comment start → stop reading further
            return line[:i].rstrip(), line[i+1:].lstrip()

    return line.rstrip(), ""

def preprocess_fortran(lines):
    cleaned = []
    original = []
    comments = []
    buffer = ""

    for line in lines:
        # Remove comments
        line, comment = separate_fortran_comment(line)

        # Trim whitespace around continuation
        lstripped = line.lstrip()

        # If this line continues the previous with a leading &
        if lstripped.startswith('&'):
            # Remove leading &, join to existing buffer
            continuation = lstripped[1:].strip()
            buffer += continuation
            original.append(line)
            comments.append(line)
            continue

        # If this line ends with &
        if line.endswith('&'):
            buffer += line[:-1].rstrip()
            original.append(line)
            comments.append(comment)
            continue

        # If previous line was a continuation
        if buffer:
            buffer += " " + lstripped
            original.append(line)
            comments.append(comment)
            cleaned.append((buffer, original, comments))
            buffer = ""
            original = []
            comments = []
        else:
            if lstripped:
                cleaned.append((lstripped, [line], [comment]))

    # Append leftover if exists
    if buffer:
        cleaned.append((buffer, original, comments))

    return cleaned

# identifier: starts with letter, followed by letters/digits/underscore
IDENT_RE = r'[A-Za-z][A-Za-z0-9_]*'   # Fortran identifier (starts with letter)

def extract_var_info(clean_lines, query_vars, filename, known_cparams = [], debug = True):
    """
    Scan cleaned Fortran declaration lines and return info for each variable in query_vars.

    Returns a dict keyed by the original query name (preserving case), e.g.
      { "FI_mixfrac_pdf2d": {"type": "...", "value": "...", "decl": "...", "parsed_name": "fi_mixfrac_pdf2d"} }

    Behavior:
    - Matching is case-insensitive.
    - Tokens must start with a letter and can include letters, digits, underscore.
    - Splits commas not inside parentheses.
    """
    # map lowercase -> original query name(s) (in case user asked same name twice with different case)
    lc_to_query = {}
    for q in query_vars:
        lc = q.lower()
        lc_to_query.setdefault(lc, []).append(q)

    results = {}
    need = set(lc_to_query.keys())  # lowercase names we still need
    get_all = not known_cparams

    # patterns
    decl_re = re.compile(r'^\s*(?P<type>[^:]+?)\s*::\s*(?P<vars>.+)$', flags=re.IGNORECASE)
    tok_re  = re.compile(rf'^\s*(?P<name>{IDENT_RE})\s*(?P<vdims>\([^)]*\))?\s*(?:=\s*(?P<val>.+))?\s*$', flags=re.IGNORECASE)
    split_commas_not_in_parens = re.compile(r',(?![^()]*\))')

    for line, original, comments in clean_lines:
        if not get_all and not need:
            break  # all found

        mdecl = decl_re.match(line)
        if not mdecl:
            continue

        type_part = mdecl.group('type').strip()
        varlist = mdecl.group('vars').strip()
        # safe split
        tokens = [t.strip() for t in split_commas_not_in_parens.split(varlist) if t.strip()]

        for tok in tokens:
            mt = tok_re.match(tok)
            if not mt:
                # token didn't parse; skip
                if debug:
                    print("  FAILED TOKEN PARSE:", repr(tok))
                continue

            parsed_name = mt.group('name')
            lc_name = parsed_name.lower()
            vdims = mt.group('vdims') or ""
            val = mt.group('val').strip() if mt.group('val') is not None else None

            if lc_name in need or get_all:
                comment = ""
                # extract comment, if available
                for i, origline in enumerate(original):
                    if re.findall(rf"\b{lc_name}\b", origline, flags=re.IGNORECASE):
                        comment = comments[i]
                        break
                # for each original query spelling for this name, store result keyed by that original
                valref = val.strip().replace(".", "\\.") if val is not None else ""
                if known_cparams and valref:
                    pattern = r"\b(" + "|".join(map(re.escape, known_cparams)) + r")\b"
                    valref = re.sub(pattern, rf"\ :f:var:`~cparam/\1`", valref)
                if not get_all:
                    for orig_query in lc_to_query[lc_name]:
                        results[orig_query] = {
                            "type": (type_part + vdims).strip(),
                            "value": valref,
                            "comment": comment,
                            #"decl": line,
                            #"parsed_name": parsed_name
                        }
                else:
                    results[lc_name] = {
                        "type": (type_part + vdims).strip(),
                        "value": valref,
                        "comment": comment
                    }
                # mark as satisfied
                need.discard(lc_name)

    # Optionally warn about queries not found
    if debug and need:
        for missing_lc in need:
            #print("NOT FOUND:", missing_lc, "requested as", lc_to_query.get(missing_lc), " in file ", filename)
            for orig_query in lc_to_query[missing_lc]:
                results[orig_query] = {
                    "type": "",
                    "value": "",
                    "comment": "(not found)"
                }
    return results

def get_init_and_run_pars(filename, known_cparams):
    """
    Extract the init and run parameters from the corresponding namelists
    in the given file.
    """
    with open(f"../../src/{filename}", "r") as f:
        lines = [it.strip() for it in f.readlines()]
    initparsname = ""
    runparsname = ""
    initpars = []
    runpars = []
    ininit = False
    inrun = False
    def _clean_line(line):
        if "!" in line:
            # If the line accidentally contains a comment (it happens!), remove it
            line = line[:line.index("!")].rstrip()
        # Remove the continuation character
        return line.strip("&")
    # Loop through all the lines in the file
    for line in lines:
        # If it's the start of a block, start appending lines
        # and mark that we are inside that block
        if matches := re.findall(r"namelist /(.*?_?init_pars)/", line):
            initparsname = matches[0]
            initpars.append(_clean_line(line))
            ininit = True
        elif matches := re.findall(r"namelist /(.*?_?run_pars)/", line):
            runparsname = matches[0]
            runpars.append(_clean_line(line))
            inrun = True
        # If we are inside a block, keep appending lines
        elif ininit:
            initpars.append(_clean_line(line))
            # but as soon as we no longer encounter the continuation character, the block has ended
            if not line.endswith("&"):
                ininit = False
        # same for run_pars
        elif inrun:
            runpars.append(_clean_line(line))
            if not line.endswith("&"):
                inrun = False
        # If we are in neither block, but both lists have elements, it means we have passed both blocks.
        # No need to continue going through the lines.
        elif initpars and runpars:
            break
    # Now put all lines together, and clean them up
    initstr = ''.join(initpars).replace(f"namelist /{initparsname}/", "")
    runstr = ''.join(runpars).replace(f"namelist /{runparsname}/", "")
    # Split the lines by comma, strip all blank spaces, and form the final list of parameters
    init_pars = dict.fromkeys(sorted(list(filter(None, map(str.strip, initstr.split(",")))), key=str.lower))
    run_pars = dict.fromkeys(sorted(list(filter(None, map(str.strip, runstr.split(",")))), key=str.lower))
    # Extract type and default value
    if filename == "param_io.f90":
        with open("../../src/cdata.f90", "r") as f:
            lines.extend([it.strip() for it in f.readlines()])
    clean_lines = preprocess_fortran(lines)
    init_info = extract_var_info(clean_lines, init_pars.keys(), filename, known_cparams)
    for k, v in init_info.items():
        assert k in init_pars
        init_pars[k] = v
    run_info = extract_var_info(clean_lines, run_pars.keys(), filename, known_cparams)
    for k, v in run_info.items():
        assert k in run_pars
        run_pars[k] = v

    return initparsname, init_pars, runparsname, run_pars

def process_init_run_pars():
    """
    Extract the init and run parameters from the corresponding namelists
    from all the files.
    """

    # Get all constants from cparam; we'll link them
    with open("../../src/cparam.f90") as f:
        cparamlines = f.readlines()
    clean_cparam_lines = preprocess_fortran(cparamlines)
    cparam_info = extract_var_info(clean_cparam_lines, [], "cparam.f90")
    known_cparams = list(cparam_info.keys())

    files = subprocess.check_output(
        f'grep -rl -e "namelist /.*init_pars/" -e "namelist /.*run_pars/" ../../src/',
        shell=True,
        text=True
    ).splitlines()
    files = [it for it in files if it.endswith(".f90") and not it.endswith("/cdata.f90") and not it.endswith("/param_io.f90")]
    files = sorted([it.replace("../../src/", "") for it in files], key=lambda x: (x.count("/") > 0, x.lower()))
    files.insert(0, "param_io.f90")
    init_pars = {}
    run_pars = {}
    for f in files:
        print(f)
        iname, ipars, rname, rpars = get_init_and_run_pars(f, known_cparams)
        if iname and ipars:
            init_pars[f] = (iname, ipars)
        if rname and rpars:
            run_pars[f] = (rname, rpars)

    stats_commented = 0
    stats_uncommented = 0
    stats_not_found = 0
    total = 0
    for vars in init_pars.values():
        for v in vars[1].values():
            if not v['comment']:
                stats_uncommented += 1
            elif v['comment'] == '(not found)':
                stats_not_found += 1
            else:
                stats_commented += 1
            total += 1
    frac_commented = round(stats_commented * 100 / total, 1)
    if frac_commented == 0 and stats_commented > 0:
        frac_commented = 0.1
    frac_not_found = round(stats_not_found * 100 / total)
    if frac_not_found == 0 and stats_not_found > 0:
        frac_not_found = 0.1
    frac_uncommented = 100 - frac_commented - frac_not_found

    with open(f"code/tables/init.rst", "w") as f:
        d = RstCloth(f, line_width=5000)
        d.title(f"Startup parameters for ``start.in``")
        d.content(f"This page lists {total} variables distributed into {len(init_pars)} files. Of these:")
        d.newline()
        d.li(f"{stats_commented} ({frac_commented:.1f}%) are documented;")
        d.li(f"{stats_uncommented} ({frac_uncommented:.1f}%) are undocumented;")
        d.li(f"{stats_not_found} ({frac_not_found:.1f}%) could not be found (they appear in the namelist, but the declaration could not be found).")
        d.newline()
        d.directive("raw", "html", content="<div>Filter: <input type='text' id='customvarsearch' /></div><br/>")
        d.newline()
        for f, vars in init_pars.items():
            if f == "param_io.f90":
                d.h3("Module *cdata.f90* / *param_io.f90*")
            else:
                d.h3(f"Module *{f}*")
            d.content(f"The following variables are part of the *{vars[0]}* namelist:")
            d.newline()
            table = []
            for k, v in vars[1].items():
                if v:
                    table.append((f"*{k}*", v['type'], v['value'], v['comment']))
                else:
                    table.append((f"*{k}*", "", "", ""))
            d.table_list(["Variable", "Type", "Default", "Meaning"], data=table, widths=[20, 10, 10, 60])

    stats_commented = 0
    stats_uncommented = 0
    stats_not_found = 0
    total = 0
    for vars in run_pars.values():
        for v in vars[1].values():
            if not v['comment']:
                stats_uncommented += 1
            elif v['comment'] == '(not found)':
                stats_not_found += 1
            else:
                stats_commented += 1
            total += 1
    frac_commented = round(stats_commented * 100 / total, 1)
    if frac_commented == 0 and stats_commented > 0:
        frac_commented = 0.1
    frac_not_found = round(stats_not_found * 100 / total)
    if frac_not_found == 0 and stats_not_found > 0:
        frac_not_found = 0.1
    frac_uncommented = 100 - frac_commented - frac_not_found

    with open(f"code/tables/run.rst", "w") as f:
        d = RstCloth(f, line_width=5000)
        d.title(f"Runtime parameters for ``run.in``")
        d.content(f"This page lists {total} variables distributed into {len(run_pars)} files. Of these:")
        d.newline()
        d.li(f"{stats_commented} ({frac_commented:.1f}%) are documented;")
        d.li(f"{stats_uncommented} ({frac_uncommented:.1f}%) are undocumented;")
        d.li(f"{stats_not_found} ({frac_not_found:.1f}%) could not be found (they appear in the namelist, but the declaration could not be found).")
        d.newline()
        d.directive("raw", "html", content="<div>Filter: <input type='text' id='customvarsearch' /></div><br/>")
        d.newline()
        for f, vars in run_pars.items():
            if f == "param_io.f90":
                d.h3("Module *cdata.f90* / *param_io.f90*")
            else:
                d.h3(f"Module *{f}*")
            d.content(f"The following variables are part of the *{vars[0]}* namelist:")
            d.newline()
            table = []
            for k, v in vars[1].items():
                if v:
                    table.append((f"*{k}*", v['type'], v['value'], v['comment']))
                else:
                    table.append((f"*{k}*", "", "", ""))
            d.table_list(["Variable", "Type", "Default", "Meaning"], data=table, widths=[20, 10, 10, 60])

def process_all_pcparam():
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
        process_init_run_pars()
        for pars in diag_list:
            process_diag(*pars)
        process_boundary_conditions()
        toclist = ["init", "run"]
        toclist.extend([pars[1] for pars in diag_list])
        toclist.append("boundary")
        d.directive("toctree", fields=[("maxdepth", "1")], content=toclist)