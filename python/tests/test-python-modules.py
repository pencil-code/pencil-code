#!/usr/bin/python3
# -*- coding: utf-8 -*-   vim: set fileencoding=utf-8 :

"""Run unit tests for Pencil Code Python modules."""

import argparse
import json
import os
import pathlib
import re
import subprocess
import sys
import warnings

try:
    import coverage
except ImportError:
    coverage_pkg_present = False
else:
    coverage_pkg_present = True

pytest_fast_flags = ["-x", "--failed-first"]
pytest_coverage_flags = [
    "--cov=pencil",
    "--cov-report=",
    "--cov-context=test",
    "--cov-append",
    "--script-test-coverage",
    ]

def call_pytest(fast=False):
    #Keep this import here so that call_tox works without pytest installed.
    import pytest

    if fast:
        fast_flags = pytest_fast_flags
    else:
        fast_flags = []

    sys.exit(pytest.main([
        '-c',
        str(pathlib.Path(__file__).parent/"pytest.ini"),
        "-m",
        "not integration",
        *fast_flags
        ]))

class get_repo_version:
    def __init__(self, repo_loc):
        self.commit = "Unknown"
        self.branch = "Unknown"

        kwargs = {
            'capture_output': True,
            'text': True,
            'check': True,
            'cwd': repo_loc,
            }
        try:
            p = subprocess.run(["git", "rev-parse", "HEAD"], **kwargs)
            self.commit = p.stdout.strip()

            p = subprocess.run(["git", "rev-parse", "--abbrev-ref", "HEAD"], **kwargs)
            self.branch = p.stdout.strip()
        except Exception as e:
            warnings.warn(f"getting git repository details failed: {e}")

def call_tox(output_dir, report_coverage=True, fast=False):
    """
    output_dir: pathlib.Path instance
    """
    output_dir = output_dir.absolute()
    if not output_dir.exists():
        output_dir.mkdir()

    json_filename = output_dir/"report.json"
    html_filename = output_dir/"index.html"
    py_tests_dir = pathlib.Path(__file__).parent

    if report_coverage:
        subprocess.run(["coverage", "erase"], check=True, cwd=py_tests_dir)
        coverage_flags = " ".join(pytest_coverage_flags)
    else:
        coverage_flags = ""

    if fast:
        fast_flags = " ".join(pytest_fast_flags)
    else:
        fast_flags = ""

    p = subprocess.Popen("bash", stdin=subprocess.PIPE, text=True)
    _, _ = p.communicate(
        f"""
        source {py_tests_dir/"../../sourceme.sh"}
        tox run --conf "{py_tests_dir}/tox.ini" \
            --result-json "{json_filename}" \
            --colored no \
            --override "testenv.setenv=PYTEST_ADDOPTS='--color=no'" \
            -- \
            {coverage_flags} {fast_flags}
        """
        )

    if report_coverage:
        htmlcov_dir = output_dir/"htmlcov"
        if not htmlcov_dir.exists():
            htmlcov_dir.mkdir()

        subprocess.run(["coverage", "html", f"--directory={htmlcov_dir}"], check=True, cwd=py_tests_dir)

    json_to_html(
        json_filename,
        html_filename,
        repo_info=get_repo_version(py_tests_dir),
        )
    print(f"Saved HTML test report in {html_filename}")
    sys.exit(p.returncode)

_ansi_escape = re.compile(r'''
    \x1B  # ESC
    (?:
        # OSC: any sequence of characters, terminated by \x1B\\
        \]
        .*
        \x1B\\ #string terminator
    |   #pytest seems to use the following for hyperlinks
        \]
        8;;[^\x1B\x07]*\x07
    |   # 7-bit C1 Fe (except CSI)
        [@-Z\\-_]
    |   # or [ for CSI, followed by a control sequence
        \[
        [0-?]*  # Parameter bytes
        [ -/]*  # Intermediate bytes
        [@-~]   # Final byte
    )
''', re.VERBOSE)

def strip_ansi(text):
    """
    Strip ANSI control sequences from text. This is needed because pytest uses
    OSC sequences in its output even when --color=no is passed.

    Copied from https://stackoverflow.com/questions/14693701/how-can-i-remove-the-ansi-escape-sequences-from-a-string-in-python/14693789#14693789
    """
    return _ansi_escape.sub('', text)

def json_to_html(
    json_filename,
    html_filename,
    repo_info = None,
    add_coverage_link = True,
    ):
    """
    repo_info: get_repo_version instance
    """

    with open(json_filename) as f:
        report = json.load(f)

    if report['reportversion'] != "1":
        raise NotImplementedError(f"for tox report version {report['reportversion']}")

    css = """
        h1 {
            font-size: 1.25rem;
            }
        pre {
            font-family: monospace;
            overflow-x: auto;
            }
        footer {
            text-align: center;
            background: #f8f8f8;
            border-top: 1px solid #ccc;
            }
        
        th, td {
            border-spacing: 1em 0em;
            text-align: left;
            vertical-align: top;
            }

        .success_message {
            color: green;
            }
        .skip_message {
            color: gold;
            }
        .failure_message {
            color: red;
            }

        #test_summary {
            width: fit-content;
            padding: 2rem 1rem;
            margin: 0 auto;
            background: #f8f8f8;
            border: 1px solid #ccc;
            border-radius: 1em;

            h1 {
                text-align: center;
                }
            }
        """
    html_head = f"""
        <!DOCTYPE html>
        <html>
        <head>
        <style>{css}</style>
        <title>Tox test report</title>
        </head>
        <body>
        """
    html_foot = f"""
        <footer>
            Generated by tox {report['toxversion']} on {report['host']}.
        </footer>
        </body>
        </html>"""

    #These will be built in the loop that follows
    html_body = ""
    html_summary = ""

    if repo_info is not None:
        html_body += f"""
    <p>Tested git commit: {repo_info.commit} (branch {repo_info.branch})</p>
    """

    if add_coverage_link:
        html_body += f"""
        <p><a href="htmlcov/index.html">click to view code coverage report</a></p>
        """

    all_succeeded = True
    for name, testenv in report['testenvs'].items():
        if 'test' not in testenv:
            status = "<span class='skip_message'>skipped</span> (interpreter not available?)"
            skipped = True
        elif testenv['result']['success']:
            status = "<span class='success_message'>succeeded</span>"
            skipped = False
        else:
            status = "<span class='failure_message'>failed</span>"
            all_succeeded = False
            skipped = False

        html_summary += f"""<tr>
            <td><a href='#section_results_{name}'>{name}</a></td>
            <td>{status}</td>
            </tr>"""

        html_body += f"""
            <h1>Test details for <a id='section_results_{name}'>{name}</a></h1>
            <p>status: {status}</p>
            """

        if not skipped:
            html_body += f"""
                <p>Python version: {testenv['python']['version']}</p>
                """

            for heading, k, k_o, show_empty in [
                ("Setup", 'setup', 'output', True),
                ("Tests", 'test', 'output', True),
                ("Errors", 'test', 'err', False),
                ]:
                details = ""
                for d in testenv[k]:
                    if len(details) > 0:
                        details += "\n"

                    output = strip_ansi(d[k_o])
                    if len(output) > 0:
                        details += f"$ {" ".join(d['command'])}\n{output}"
                    elif show_empty:
                        details += f"$ {" ".join(d['command'])}"

                if len(details) > 0:
                    html_body += f"""
                        <details>
                            <summary>{heading}</summary>
                            <pre>{details}</pre>
                        </details>
                        """

    if all_succeeded:
        overall_status = "<span class='success_message'>All succeeded</span>"
    else:
        overall_status = "<span class='failure_message'>Some failures</span>"

    html_summary = f"""
        <div id="test_summary">
        <h1>{overall_status}</h1>
        <table>
            <tr>
                <th>Environment</th><th>Status</th>
            </tr>
            {html_summary}
        </table>
        </div>
        """

    with open(html_filename, 'w') as f:
        f.write("\n".join([html_head, html_summary, html_body, html_foot]))

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--outputdir",
        help = "location to store HTML output in (only has an effect with `--tox`)",
        default = "."
        )
    parser.add_argument(
        "--fast",
        help = "exit immediately on the first failure, rather than running all the tests regardless. Also, on repeated runs, run the previously failed tests first.",
        default = False,
        action = 'store_true',
        )
    parser.add_argument(
        "--tox",
        help = "Run the full set of tests by installing the required Python modules in an isolated environment. Also generates a HTML report of the test status",
        default = False,
        action = 'store_true',
        )
    parser.add_argument(
        "--coverage",
        help = "Generate the code coverage report (only has an effect with `--tox`).",
        default = False,
        action = 'store_true',
        )

    args = parser.parse_args()

    if args.coverage and not coverage_pkg_present:
        raise RuntimeError("`coverage` (Python package) must be installed to use the --coverage option.")

    if args.tox:
        call_tox(
            output_dir = pathlib.Path(args.outputdir),
            fast = args.fast,
            report_coverage = args.coverage,
            )
    else:
        call_pytest(fast=args.fast)
