# Copyright 2012-2013 (C) Daniel Watkins <daniel@daniel-watkins.co.uk>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.


from datetime import datetime
import os
import re
import sys
import subprocess
import six
from docutils import nodes
from docutils.parsers.rst import Directive, directives
from git import Repo
import importlib


# pylint: disable=too-few-public-methods
class GitChangelog(Directive):
    default_sha_length = 7

    def run(self):
        commits = self._commits_to_display()
        markup = self._build_markup(commits)
        return markup

    def _commits_to_display(self):
        env = self.state.document.settings.env
        self.file = f"{env.temp_data['docname']}.rst"
        self.repo = Repo(env.srcdir, search_parent_directories=True)
        commits = self._filter_commits(self.repo)
        return commits

    def _filter_commits(self, repo):
        if self.file == None:
            return []
        hashes = (
            subprocess.check_output(
                ["git", "log", "--follow", "--pretty=format:%H", self.file]
            )
            .decode()
            .splitlines()
        )

        commits = list(repo.iter_commits())
        return self._filter_commits_on_filenames(commits, hashes)

    def _filter_commits_on_filenames(self, commits, hashes):
        filtered_commits = []
        for commit in commits:
            if commit.hexsha in hashes:
                filtered_commits.append(commit)
        return filtered_commits

    def _build_markup(self, commits):
        parent = nodes.container()
        list_node = nodes.bullet_list()
        for commit in commits:
            # Construct a backlink to BitBucket
            url = f"https://github.com/pencil-code/pencil-code/commit/{commit.hexsha}"
            ref = nodes.emphasis(
                "",
                "",
                nodes.reference(
                    "", commit.hexsha[: self.default_sha_length], refuri=url
                ),
            )

            # Other commit details
            date_str = datetime.fromtimestamp(commit.authored_date)

            item = nodes.list_item()
            par = nodes.paragraph()
            par += [nodes.inline(text=commit.message)]

            author = six.text_type(commit.author)
            par += [
                nodes.inline(text=" by ", classes=["gray"]),
                nodes.emphasis(text=author, classes=["gray"]),
            ]
            par += [
                nodes.inline(text=" at ", classes=["gray"]),
                nodes.emphasis(text=str(date_str), classes=["gray"]),
            ]
            par += nodes.inline(text=", ")
            par += ref
            item.append(par)

            list_node.append(item)
        parent.update_basic_atts({"classes": ["collapsable"]})
        parent.append(list_node)
        return [parent]


def add_changelog(app, docname, source):
    """
    Add revision history to all non-generated rst files.
    """
    if "code/sourcePython" in docname or docname.startswith("code/py"):
        return
    source[0] += "\n\n----------\n\n.. rubric:: Revision history\n.. docs_changelog::\n"


def setup(app):
    app.add_directive("docs_changelog", GitChangelog)
    # Add revision history to all rst files
    app.connect("source-read", add_changelog)
