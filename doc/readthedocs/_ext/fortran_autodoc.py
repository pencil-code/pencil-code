# -*- coding: utf-8 -*-
"""Sphinx extension for autodocumenting fortran codes.


"""
# Copyright or © or Copr. Actimar/IFREMER (2010-2019)
#
# This software is a computer program whose purpose is to provide
# utilities for handling oceanographic and atmospheric data,
# with the ultimate goal of validating the MARS model from IFREMER.
#
# This software is governed by the CeCILL license under French law and
# abiding by the rules of distribution of free software.  You can  use,
# modify and/ or redistribute the software under the terms of the CeCILL
# license as circulated by CEA, CNRS and INRIA at the following URL
# "http://www.cecill.info".
#
# As a counterpart to the access to the source code and  rights to copy,
# modify and redistribute granted by the license, users are provided only
# with a limited warranty  and the software's author,  the holder of the
# economic rights,  and the successive licensors  have only  limited
# liability.
#
# In this respect, the user's attention is drawn to the risks associated
# with loading,  using,  modifying and/or developing or reproducing the
# software by the user in light of its specific status of free software,
# that may mean  that it is complicated to manipulate,  and  that  also
# therefore means  that it is reserved for developers  and  experienced
# professionals having in-depth computer knowledge. Users are therefore
# encouraged to load and test the software's suitability as regards their
# requirements in conditions enabling the security of their systems and/or
# data to be ensured and,  more generally, to use and operate it in the
# same conditions as regards security.
#
# The fact that you are presently reading this means that you have had
# knowledge of the CeCILL license and that you accept its terms.
#
from __future__ import print_function
from builtins import str
from builtins import range
from builtins import object
import six
import re
import os
from collections import OrderedDict
from operator import itemgetter

from docutils.parsers.rst import Directive
from docutils.parsers.rst.directives import unchanged
from docutils.statemachine import string2lines
from sphinx.util.console import bold
from glob import glob
from numpy.f2py.crackfortran import crackfortran, fortrantypes
from sphinx.util import logging

from sphinxfortran.fortran_domain import FortranDomain


logger = logging.getLogger(__name__)


# Fortran parser and formatter
# ----------------------------

class F90toRstException(Exception):
    pass


class F90toRst(object):
    '''Fortran 90 parser and restructeredtext formatter

    :Parameters:

         - **ffiles**: Fortran files (glob expression allowed) or dir (or list of)

    :Options:

        - **ic**: Indentation char.
        - **ulc**: Underline char for titles.
        - **sst**: Subsection type.
        - **vl**: Verbose level (0=quiet).
    '''
    _re_unended_match = re.compile(r'(.*)&\s*', re.I).match
    _re_unstarted_match = re.compile(r'\s*&(.*)', re.I).match
    _re_comment_match = re.compile(r'\s*!(.*)', re.I).match
    _re_space_prefix_match = re.compile(r'^(\s*).*$', re.I).match
#    _re_sub_numvardesc_match = re.compile(r'\$(\d+)',re.I).match
#    _re_vardesc_match = re.compile(r'\s*(?P<varname>\w+)[\s,&]*!\s*(?P<vardesc>.+)',re.I).match
#    _re_type_var_findall = re.compile(r'(\w+)',re.I).match
    _fmt_vardesc = ':%(role)s %(vtype)s %(vname)s%(vdim)s%(vattr)s: %(vdesc)s'
    #_fmt_vardesc = '**%(vname)s** %(vdim)s :: %(vtype)s %(vattr)s: %(vdesc)s'
    _fmt_vattr = " [%(vattr)s]"
    _fmt_fvardesc = '%(vtype)s%(vdim)s %(vattr)s%(vdesc)s'
    #_fmt_fvardim = ', **dim%(vdim)s**'
#    _re_call_subr_search = re.compile(r'\bcall\s+(\w+)\(',re.I).search
#    _re_call_func_findall = re.compile(r'\b(\w+)\(',re.I).search

    def __init__(
            self,
            ffiles,
            ic='\t',
            ulc='-',
            vl=0,
            encoding='utf8',
            sst='rubric'):
        # Be sure to have a list
        if not isinstance(ffiles, list):
            ffiles = [ffiles]

        # Read and store them
        self.src = OrderedDict()
        self.ffiles = ffiles
        for ff in ffiles:
            f = open(ff)
            self.src[ff] = []
            for l in f.readlines():
                try:
                    self.src[ff].append(l[:-1].decode(encoding))
                except AttributeError:
                    self.src[ff].append(l[:-1])
                except BaseException:
                    raise F90toRstException(
                        'Encoding error\n  file = %s\n  line = %s' %
                        (ff, l))
            #self.src[ff] = [l.decode(encoding) for l in f.readlines()]
            #self.src[ff] = [l.decode(encoding) for l in f.read().split('\n')]

            f.close()

        # Crack files
        global verbose, quiet
        import numpy.f2py.crackfortran
        self._verbose = numpy.f2py.crackfortran.verbose = verbose = vl
        numpy.f2py.crackfortran.quiet = quiet = 1 - verbose
        self.crack = []
        for ff in ffiles:
            self.crack.extend(crackfortran(ff))

        # Build index
        self.build_index()

        # Scan all files to get description, etc
        self.scan()

        # Add function 'call from' to index
        self.build_callfrom_index()

        # Other inits
        self.rst = OrderedDict()
        self._ic = ic
        self._ulc = ulc
        self._sst = sst

    # Indexing ---

    def build_index(self):
        """Register modules, functions, types and module variables for quick access

        Index constituents are:

            .. attribute:: modules

                Dictionary where each key is a module name, and each value is the cracked block.

            .. attribute:: routines

                Module specific functions and subroutines

            .. attribute:: types

                Module specific types

            .. attribute:: variables

                Module specific variables
        """
        # Containers
        self.modules = OrderedDict()
        self.types = OrderedDict()
        self.routines = self.functions = self.subroutines = OrderedDict()
        self.variables = OrderedDict()
        self.programs = OrderedDict()

        # Loop on all blocks and subblocks
        for block in self.crack:

            # Modules
            if block['block'] == 'module':

                # Index module
                module = block['name']
                self.modules[module] = block

                # Loop inside module
                for subblock in block['body']:

                    # Types and routines as subblock
                    if subblock['block'] in ['function', 'type', 'subroutine']:

                        # Index
                        container = getattr(self, subblock['block'] + 's')
                        container[subblock['name']] = subblock
                        subblock['module'] = module

                        # Variables
                        varnames = subblock['sortvars']
                        if subblock['block'] == 'type':
                            varnames.sort()
                        for varname in varnames:
                            subblock['vars'][varname]['name'] = varname

                # Function aliases from "use only" (rescan)
                for bfunc in list(self.routines.values()):
                    bfunc['aliases'] = []
                for subblock in block['body']:
                    if not subblock['block'] == 'use':
                        continue
                    for monly in list(block['use'].values()):
                        if not monly:
                            continue
                        for fname, falias in list(monly['map'].items()):
                            self.routines[falias] = self.routines[fname]
                            if falias not in self.routines[fname]['aliases']:
                                self.routines[fname]['aliases'].append(falias)

                # Module variables
                for varname in sorted(block['sortvars']):
                    bvar = block['vars'][varname]
                    if varname not in self.routines:
                        self.variables[varname] = bvar
                        # self.variables.pop(varname)
                        bvar['name'] = varname
                        bvar['module'] = module

            # Local functions, subroutines and programs
            elif block['block'] in ['function', 'subroutine', 'program']:

                # Index
                container = getattr(self, block['block'] + 's')
                container[block['name']] = block

                # Variables
                for varname in block['sortvars']:
                    block['vars'][varname]['name'] = varname

        # Regular expression for fast search
        # - function calls
        subs = [block['name'].lower() for block in list(
            self.routines.values()) if block['block'] == 'subroutine']
        self._re_callsub_findall = subs and re.compile(
            r'call\s+\b(%s)\b' %
            ('|'.join(subs)),
            re.I).findall or (
            lambda line: [])
        funcs = [block['name'].lower() for block in list(
            self.routines.values()) if block['block'] == 'function']
        self._re_callfunc_findall = funcs and re.compile(
            r'\b(%s)\b\s*\(' %
            ('|'.join(funcs)),
            re.I).findall or (
            lambda line: [])
        # - function variables
        for block in list(self.routines.values())+list(self.types.values()):
            vars = (r'|').join(block['sortvars']) + r'|\$(?P<varnum>\d+)'
            sreg = r'[\s\*\-:]*(?:@param\w*)?\b(?P<varname>%s)\b\W+(?P<vardesc>.*)' % vars
            block['vardescmatch'] = re.compile(sreg).match
        # - variables with description
#        for block in self.types.values()+self.modules.values():
        for block in list(self.types.values()) + \
                list(self.modules.values()) + list(self.routines.values()):
            #sreg = r'\b(?P<varname>%s)\b[\W\d]*!\s*(?P<vardesc>.*)'%'|'.join(block['sortvars'])
            #sreg = r'[\W\(\),\b\*=\-\&]*?:?:[ \t\&]*(?P<varname>%s)\b[\w\s\(\)\*,_=]*!\s*(?P<vardesc>.*)'%'|'.join(block['sortvars'])
            #sreg = r'.*[\W\(\),\b\*=\-\&]*?:?:[ \t\&]*(?P<varname>%s)\b[\w\s\(\)\*,_=\.]*!\s*(?P<vardesc>.*)'%'|'.join(block['sortvars'])
            # reversed+sorted is a hack to avoid conflicts when variables share
            # the same prefix
            if block['sortvars']:
                sreg = r'.*\b(?P<varname>%s)\b\s*(?P<dims>\([\*:,\w]+\))?[^!\)]*!\s*(?P<vardesc>.*)\s*' % '|'.join(
                    reversed(sorted(block['sortvars'])))
                block['vardescsearch'] = re.compile(sreg, re.I).search
            else:
                block['vardescsearch'] = lambda x: None

    def build_callfrom_index(self):
        """For each function, index which function call it"""
        for bfunc in list(self.routines.values()):
            bfunc['callfrom'] = []
            for bfuncall in list(self.routines.values()) + \
                    list(self.programs.values()):
                if bfunc['name'] in bfuncall['callto']:
                    bfunc['callfrom'].append(bfuncall['name'])

    def filter_by_srcfile(self, sfile, mode=None, objtype=None, **kwargs):
        """Search for subblocks according to origin file

        :Params:
            - **sfile**: Source file name.
            - **mode**, optional: Mode for searching for sfile.
              If ``"strict"``, exact match is needed, else only basename.
            - **objtype**, optional: Restrict search to one or a list of object types
              (i.e. ``"function"``, ``"program"``, etc).
        """
        if mode is None:
            mode = 'basename'
        if objtype and not isinstance(objtype, list):
            objtype = [objtype]
        bb = []
        if mode != 'strict':
            sfile = os.path.basename(sfile)
        for b in self.crack:
            if objtype and objtype != 'all' and b['block'] not in objtype:
                continue
            bfile = b['from'].split(':')[0]  # remove module name
            if mode != 'strict':
                bfile = os.path.basename(bfile)
            if sfile == bfile:
                bb.append(b)

        return bb

    def scan(self):
        """Scan """
        # Loop on all blocks
        for block in self.crack:

            # Modules
            if block['block'] == 'module':

                # Get source lines
                modsrc = self.get_blocksrc(block)

                # Get module description
                block['desc'] = self.get_comment(modsrc, aslist=True)

                # Scan types and routines
                for subblock in block['body']:

                    self.scan_container(subblock, insrc=modsrc)

                # Scan module variables
                self.strip_blocksrc(
                    block, [
                        'type', 'function', 'subroutine'], src=modsrc)
                if modsrc:
                    for line in modsrc:
                        if line.strip().startswith('!'):
                            continue
                        m = block['vardescsearch'](line)
                        if m:
                            block['vars'][m.group('varname').lower()]['desc'] = m.group(
                                'vardesc')
                for bvar in list(block['vars'].values()):
                    bvar.setdefault('desc', '')

            # Routines
            elif block['block'] in ['function', 'subroutine', 'program']:

                self.scan_container(block)

    def scan_container(self, block, insrc=None):
        """Scan a block of program, routines or type"""

        # Check block type
        if block['block'] not in ['type', 'function', 'subroutine', 'program']:
            return

        # Source lines
        subsrc = self.get_blocksrc(block, insrc)

        # Comment
        block['desc'] = self.get_comment(subsrc, aslist=True)

        # Scan comments to find descriptions
        if block['desc'] and block['block'] in [
                'function', 'subroutine', 'type']:
            varname = None
            for iline, line in enumerate(block['desc']):
                if 'vardescmatch' in block:
                    m = block['vardescmatch'](line)
                    if m:  # There is a variable and its description

                        # Name of variable
                        varname = m.group('varname')

                        # Numeric name
                        if m.group('varnum'):  # $[1-9]+
                            ivar = int(m.group('varnum'))
                            ivar = ivar - 1
                            if ivar < 0 or ivar >= len(block['args']):
                                continue
                            block['desc'][iline] = line.replace(
                                varname, block['args'][ivar])
                            varname = block['args'][ivar]

                        # Store description
                        ifirst = len(line) - len(line.strip())
                        if varname != '':
                            block['vars'][varname]['desc'] = m.group(
                                'vardesc')

                elif line.strip() and varname is not None and \
                        (len(line) - len(line.strip())) > ifirst:
                    # Description continuation?

                    block['vars'][varname]['desc'].append(
                        ' ' + line.strip())

                else:
                    varname = None

        # Callable objects
        if block['block'] in ['function', 'subroutine', 'program']:

            # Index calls
            block['callto'] = []
            if subsrc is not None:
                self.join_src(subsrc)
                for line in subsrc[1:-1]:
                    if line.strip().startswith('!'):
                        continue
                    line = line.lower()
                    for fn in self._re_callsub_findall(
                            line) + self._re_callfunc_findall(line):
                        if fn not in block['callto']:
                            block['callto'].append(fn)
                            pass

        # Get description of variables from inline comment that overwrites
        # description inferred from header comment
        if block['block'] in [
            'function',
            'subroutine',
                'type'] and subsrc is not None:
            for line in subsrc:
                if line.strip().startswith('!'):
                    continue
                if 'vardescsearch' in block:
                    m = block['vardescsearch'](line)
                    if m:
                        block['vars'][m.group('varname').lower()]['desc'] = m.group(
                            'vardesc')

        # Fill empty descriptions
        for bvar in list(block['vars'].values()):
            bvar.setdefault('desc', '')

        del subsrc

    # Getting info ---

    def get_module(self, block):
        """Get the name of the current module"""
        while block['block'] != 'module':
            if block['parentblock'] == 'unknown':
                break
            block = block['parentblock']
        return block['name']

    def get_src(self, block):
        """Get the source lines of the file including this block"""
        srcfile = block['from'].split(':')[0]
        return self.src[srcfile]

    def join_src(self, src):
        """Join unended lines that does not finish with a comment"""
        for iline, line in enumerate(src):
            m = self._re_unended_match(line)
            if m:
                thisline = m.group(1)
                m = self._re_unstarted_match(src[iline + 1])
                nextline = m.group(1) if m else src[iline + 1]
                src[iline] = thisline + nextline
                del src[iline + 1]
        return src

    def get_blocksrc(
            self,
            block,
            src=None,
            istart=0,
            getidx=False,
            stopmatch=None,
            exclude=None):
        """Extract an identified block of source code

        :Parameters:

            - *block*: Cracked block

        :Options:

            - *src*: List of source line including the block
            - *istart*: Start searching from this line number

        :Return:

            ``None`` or a list of lines
        """
        # Set up
        if src is None:
            src = self.get_src(block)
        blocktype = block['block'].lower()
        blockname = block['name'].lower()
        ftypes = r'(?:(?:%s).*\s+)?' % fortrantypes if blocktype == 'function' else ''
        rstart = re.compile(
            r"^\s*%s%s\s+%s\b.*$" %
            (ftypes, blocktype, blockname), re.I).match
        rend = re.compile(r"^\s*end\s+%s\b.*$" % blocktype, re.I).match
        if isinstance(stopmatch, str):
            stopmatch = re.compile(stopmatch).match

        # Beginning
        for ifirst in range(istart, len(src)):
            # Simple stop on match
            if stopmatch and stopmatch(src[ifirst]):
                return
            # Ok, now check
            if rstart(src[ifirst]):
                break
        else:
            return

        # End
        for ilast in range(ifirst, len(src)):
            if stopmatch and stopmatch(src[ilast]):
                break
            if rend(src[ilast].lower()):
                break

        # Extraction
        mysrc = list(src[ifirst:ilast + 1])

        # Block exclusions
        self.strip_blocksrc(block, exclude, src=mysrc)

        if getidx:
            return mysrc, (ifirst, ilast)
        return mysrc

    def strip_blocksrc(self, block, exc, src=None):
        """Strip blocks from source lines

        :Parameters:

            - *block*:
            - *exc* list of block type to remove

        :Options:

            - *src*: list of source lines

        :Example:

        >>> obj.strip_blocksrc(lines, 'type')
        >>> obj.strip_blocksrc(lines, ['function', 'type']
        """
        if src is None:
            self.get_blocksrc(block)
        if exc is None:
            return
        if not isinstance(exc, list):
            exc = [exc]
        for subblock in block['body']:
            if subblock['block'] in exc:
                subsrc = self.get_blocksrc(subblock, src=src, getidx=True)
                if subsrc is None:
                    continue  # Already stripped
                del src[subsrc[1][0]:subsrc[1][1]]
                del subsrc

    def get_comment(
            self,
            src,
            iline=1,
            aslist=False,
            stripped=False,
            getilast=False,
            rightafter=True):
        """Search for and return the comment starting after ``iline`` in ``src``

        :Params:

            - **src**: A list of lines.
            - **iline**, optional: Index of first line to read.
            - **aslist**, optional: Return the comment as a list.
            - **stripped**, optional: Strip each line of comment.
            - **getilast**, optional: Get also index of last line of comment.
            - **rightafter**, optional: Suppose the comment right after
              the signature line. If True, it prevents from reading a comment
              that is not a description of the routine.

        :Return:

            - ``scomment``: string or list
            - OR ``scomment,ilast``: if ``getilast is True``
        """
        scomment = []
        if src:
            in_a_breaked_line = src[0].strip().endswith('&')
            for iline in range(iline, len(src)):
                line = src[iline].strip()

                # Breaked line
                if line.startswith('&'):
                    continue

                # Manage no comment line
                m = self._re_comment_match(line)
                if m is None:

                    if not scomment:

                        # Not the end of signature
                        if line.endswith('&'):
                            in_a_breaked_line = True
                            continue
                        if in_a_breaked_line:
                            in_a_breaked_line = False
                            continue

                        # Empty line but we continue searching
                        if not rightafter and not line:
                            continue

                    # Stop searching
                    break

                # Get comment part
                comment = m.group(1)

                # Strip?
                if stripped:
                    comment = comment.strip()

                # Load and remove space prefix
                if not scomment:
                    prefix = self._re_space_prefix_match(comment).group(1)
                if comment.startswith(prefix):
                    comment = comment[len(prefix):]

                # Save comment
                scomment.append(comment)

        if not aslist:
            scomment = self.format_lines(scomment, nlc=' ')
        if getilast:
            return scomment, iline
        return scomment

    def get_synopsis(self, block, nmax=2):
        """Get the first ``nmax`` non empty lines of the function, type or module comment as 1 line.

        If the header has more than ``nmax`` lines, the first one is taken and appended of '...'.
        If description if empty, it returns an empty string.
        """
        sd = []
        for line in block['desc']:
            line = line.strip()
            if not line:
                if not sd:
                    continue
                break
            sd.append(line)
            if len(sd) > nmax:
                if sd[-1].endswith('.'):
                    sd[-1] += '...'
                break
        if not sd:
            return ''
        sd = ' '.join(sd)
        return sd

    def get_blocklist(self, choice, module, sort=True):
        """Get the list of types, variables or function of a module"""
        choice = choice.lower()
        if not choice.endswith('s'):
            choice += 's'
        assert choice in ['types', 'variables', 'functions',
                          'subroutines'], "Wrong type of declaration"
        module = module.lower()
        assert module in self.modules, "Wrong module name"
        baselist = list(getattr(self, choice).values())
        sellist = [v for v in baselist if 'module' in v and v['module']
                   == module.lower()]
        if sort:
            sellist.sort(key=itemgetter('name'))
        return sellist

    # Formating ---

    def set_ulc(self, ulc):
        """Set the underline character for title inside module description"""
        self._ulc = ulc

    def get_ulc(self):
        """Get the underline character for title inside module description"""
        return self._ulc
    ulc = property(
        get_ulc,
        set_ulc,
        doc='Underline character for title inside module description')

    def set_ic(self, ic):
        """Set the indentation character"""
        self._ic = ic

    def get_ic(self):
        """Get the indentation character"""
        return self._ic
    ic = property(get_ic, set_ic, doc='Indentation character')

    def set_sst(self, sst):
        """Set the subsection type"""
        self._sst = sst

    def get_sst(self):
        """Get the subsection type"""
        return self._sst
    sst = property(
        get_sst,
        set_sst,
        doc='Subsection type ("title" or "rubric")')

    def indent(self, n):
        """Get a proper indentation"""
        return n * self.ic

    def format_lines(
            self,
            lines,
            indent=0,
            bullet=None,
            nlc='\n',
            strip=False):
        """Convert a list of lines to text"""
        if not lines:
            return ''

        # Bullet
        if bullet is True:
            bullet = '-'
        bullet = (str(bullet) + ' ') if bullet else ''

        # Get current indentation for reduction
        if isinstance(lines, six.string_types):
            lines = [lines]

        # Split lines
        tmp = []
        for line in lines:
            if not line:
                tmp.append(line)
            else:
                tmp.extend(line.splitlines())
        lines = tmp
        del tmp

        # Expand tabs
        lines = [line.expandtabs(4) for line in lines]

        # Strip
        if strip:
            tmp = []
            for line in lines:
                if not tmp and not line:
                    continue
                tmp.append(line)
            lines = tmp
            del tmp
        if not lines:
            return ''

        # Column of first non space car
        goodlines = [(len(line) - len(line.lstrip()))
                     for line in lines if line.expandtabs().strip()]
        firstchar = goodlines and min(goodlines) or 0
        del goodlines

        # Redent
        mylines = [self.indent(indent) + bullet + line[firstchar:]
                   for line in lines]

        # Create text block
        text = nlc.join(mylines) + nlc
        del mylines
        return text

    def format_title(self, text, ulc=None, indent=0):
        """Create a simple rst titlec with indentation

        :Parameters:

            - *text*: text of the title

        :Options:

            - *ulc*: underline character (default to attribute :attr:`ucl`)

        :Example:

            >>> print(o.format_title('My title', '-'))
            My title
            --------
        """
        if ulc is None:
            ulc = self.ulc
        return self.format_lines([text, ulc * len(text)], indent=indent) + '\n'

    def format_rubric(self, text, indent=0):
        """Create a simple rst rubric with indentation

        :Parameters:

            - *text*: text of the rubric

        :Example:

            >>> print(o.format_rubric('My title', '-'))
            .. rubric:: My rubric
        """
        return self.format_lines('.. rubric:: ' + text, indent=indent) + '\n'

    def format_subsection(self, text, indent=indent, **kwargs):
        """Format a subsection for describing list of subroutines, types, etc"""
        if self.sst == 'title':
            return self.format_title(text, indent=indent, **kwargs)
        return self.format_rubric(text, indent=indent, **kwargs)

    def format_declaration(
            self,
            dectype,
            name,
            description=None,
            indent=0,
            bullet=None,
            options=None):
        """Create an simple rst declaration

        :Example:

        >>> print(format_declaration('var', 'myvar', 'my description', indent=1, bullet='-'))
            - .. f:var:: myvar

                my description
        """
        declaration = self.format_lines(
            '.. f:%(dectype)s:: %(name)s' %
            locals(), bullet=bullet, indent=indent)
        if options:
            declaration += self.format_options(options, indent=indent + 1)
        declaration += '\n'
        if description:
            declaration += self.format_lines(description, indent=indent + 1)
        return declaration + '\n'

    def format_options(self, options, indent=0):
        """Format directive options"""
        options = [
            ':%s: %s' %
            option for option in list(
                options.items()) if option[1] is not None]
        return self.format_lines(options, indent=indent)

    def format_funcref(
            self,
            fname,
            current_module=None,
            aliasof=None,
            module=None):
        """Format the reference link to a module function

        Formatting may vary depending on if function is local
        and is an alias.

        :Example:

        >>> print(obj.format_type('myfunc'))
        :f:func:`~mymodule.myfunc`
        """
        # Alias?
        fname = fname.lower()
        if aliasof is not None:
            falias = fname
            fname = aliasof
        if fname in self.routines and fname in self.routines[fname].get('aliases', [
        ]):
            falias = fname
            fname = self.routines[fname]['name']
        else:
            falias = None

        # Local reference ?
        if module is None and fname in self.routines:
            module = self.routines[fname].get('module')
        if module is None or (
                current_module is not None and module == current_module):
            if falias:
                return ':f:func:`%(falias)s<%(fname)s>`' % locals()
            return ':f:func:`%(fname)s`' % locals()

        # Remote reference
        from sphinxfortran.fortran_domain import f_sep
        if falias:
            ':f:func:`%(falias)s<~%(module)s%(f_sep)s%(fname)s>`' % locals()
        return ':f:func:`~%(module)s%(f_sep)s%(fname)s`' % locals()

    def format_use(self, block, indent=0, short=False):
        """Format use statement

        :Parameters:

            - *block*: a module block
        """
        # TODO: format aliases bug
        use = ''
        if 'use' in block:
            if short:
                use = ':use: '
            else:
                use = self.format_subsection('Needed modules', indent=indent)
            lines = []
            for mname, monly in list(block['use'].items()):

                # Reference to the module
                line = (self.indent(indent) if not short else '') + \
                    ':f:mod:`%s`' % mname

                # Reference to the routines
                if monly:
                    funcs = []
                    for fname, falias in list(monly['map'].items()):
                        func = self.format_funcref(fname, module=mname)
                        if fname != falias:
                            falias = self.format_funcref(
                                falias, module=mname, aliasof=fname)
                            func = '%s => %s' % (falias, func)
                        funcs.append(func)
                    line += ' (%s)' % ', '.join(funcs)

                # Short description
                if mname in self.modules and not short:
                    sdesc = self.get_synopsis(self.modules[mname])
                    if sdesc:
                        line += ': ' + sdesc

                # Append
                lines.append(line)
            if short:
                use += ', '.join(lines)
                use = self.format_lines(use, indent)
            else:
                use += self.format_lines(lines, indent, bullet='-') + '\n'
            del lines
        return use

    # def format_arithm(self, expr):
        #"""Format an arithmetic expression"""
        #ops = re.findall(r'(\W+)',expr)
        #nums = re.split(r'\W+', expr)
        # if expr.startswith(ops[0]):
        #nums = ['']+nums
        #if len(nums)!=len(ops): ops.append('')
        #newexpr = ''
        # for num, op in zip(nums, ops):
        # if num.isalpha():
        #num = ' :f:var:`%s` '%num
        # if '*' in op:
        #op = op.replace('*', '\*')
        #newexpr += num+op
        # return newexpr

    def format_argdim(self, block):
        """Format the dimension of a variable

         :Parameters:

            - *block*: a variable block
        """
        if 'dimension' in block:
            return'(%s)' % (','.join([s.strip('()') for s in block['dimension']]))
            # return'(%s)'%(','.join([s.replace(':','\:').strip('()') for s in
            # block['dimension']]))
        return ''

    def format_argattr(self, block):
        """Filter and format the attributes (optional, in/out/inout, etc) of a variable

         :Parameters:

            - *block*: a variable block
        """
        vattr = []
        if 'intent' in block and block['intent']:
            vattr.append('/'.join(block['intent']))
        if 'attrspec' in block and block['attrspec']:
            newattrs = []
            default_value = (block['='] if '=' in block and
                             not block['='].startswith('shape(') else None)
            for attr in block['attrspec']:
                #                if '=' in block:
                #                    if attr=='optional':
                #                        continue
                #                    elif attr=='parameter':
                #                        attr += '='+block['=']
                #                        #attr += '='+self.format_arithm(block['='])
                #                if attr in []:
                #                    attr = attr.upper()
                if attr == 'optional' and default_value is None:
                    continue
                newattrs.append(attr)
            if default_value is not None:
                newattrs.append('default=' + default_value)
            if 'private' in newattrs and 'public' in newattrs:
                newattrs.remove('private')
            block['attrspec'] = newattrs
            vattr.append('/'.join(block['attrspec']))
        if not vattr:
            return ''
        vattr = ','.join(vattr)
        return self._fmt_vattr % locals() if vattr else ''

    def format_argtype(self, block):
        if 'typespec' not in block:
            return ''
        vtype = block['typespec']
        if vtype == 'type':
            vtype = block['typename']
        return vtype

    def format_argfield(self, blockvar, role=None, block=None):
        """Format the description of a variable

         :Parameters:

            - *block*: a variable block
        """
        vname = blockvar['name']
        vtype = self.format_argtype(blockvar)  # ['typespec']
        #if vtype=='type': vtype = block['typename']

        vdim = self.format_argdim(blockvar)
        if ':' in vdim:
            vdim = vdim.replace(':', '*')
        vattr = self.format_argattr(blockvar)
        vdesc = blockvar['desc'] if 'desc' in blockvar else ''
        optional = 'attrspec' in blockvar and 'optional' in blockvar['attrspec']
        if not role:
            if block and vname in [block['name'], block.get('result')]:
                role = 'r'
            else:
                role = 'o' if optional else 'p'
        return self._fmt_vardesc % locals()

    def format_type(self, block, indent=0, bullet=True):
        """Format the description of a module type

         :Parameters:

            - *block*: block of the type
        """
        # Declaration and description
        declaration = self.format_declaration(
            'type',
            block['name'],
            block['desc'],
            indent=indent,
            bullet=bullet) + '\n'

        # Variables
        vlines = []
        for varname in sorted(block['sortvars']):
            bvar = block['vars'][varname]
            vlines.append(self.format_argfield(bvar, role='f'))
        variables = self.format_lines(vlines, indent=indent + 1) + '\n'
        del vlines

        return declaration + variables

    def get_varopts(self, block):
        """Get options for variable declaration as a dict"""
        options = OrderedDict()
        vdim = self.format_argdim(block)
        #if vdim!='': vdim = self._fmt_fvardim%locals()
        if vdim:
            options['shape'] = vdim
        options['type'] = self.format_argtype(block)
        vattr = self.format_argattr(block).strip(' []')
        if vattr:
            options['attrs'] = vattr
        return options

    # def format_vardesc(self, block):
        #"""Format the specification of a variable at top of its declaration content"""
        #specs = self.format_varspecs(block)
        #vdesc = block.get('desc', '')
        # return specs+'\n'+vdesc+'\n'

    def format_var(self, block, indent=0, bullet=True):
        """Format the description of a module type

         :Parameters:

            - *block*: block of the variable
        """
        # Description of the variable
        options = self.get_varopts(block)
        description = block.get('desc', None)
        if 'name' in block:
            declaration = self.format_declaration(
                'variable',
                block['name'],
                description=description,
                options=options,
                indent=indent,
                bullet=bullet)
        else:
            declaration = ''

        # Description of the sub-variables
        # if block['typespec']=='type':
            # xxxx
            #variables = []
            #btype = self.types[block['typename']]
            # for bvar in btype['vars'].values():
            #vname = '%s%%%s'%(block['name'], bvar['name'])
            #desc = self.format_vardesc(bvar, indent=0, vname=vname)
            #variables.append(self.format_declaration('var', vname, desc, indent=indent+1, bullet=True)[:-1])
            #variables = '\n'.join(variables)
        # else:
            #variables = ''

        return declaration  # +variables

    def format_signature(self, block):
        signature = ''
        nopt = 0
        for i, var in enumerate(block['args']):
            optional = 'optional' in block['vars'][var]['attrspec'] and '=' not in block['vars'][var] \
                if 'attrspec' in block['vars'][var] else False
            signature += '[' if optional else ''
            signature += ', ' if i else ''
            if optional:
                nopt += 1
            signature += var
        signature += nopt * ']'
        return signature

    def format_routine(self, block, indent=0):
        """Format the description of a function, a subroutine or a program"""
        # Declaration of a subroutine or function
        if isinstance(block, six.string_types):
            if block not in list(self.programs.keys()) + \
                    list(self.routines.keys()):
                raise F90toRstException(
                    'Unknown function, subroutine or program: %s' %
                    block)
            if block in self.programs:
                block = self.programs[block]
            else:
                block = self.routines[block]
        elif block['name'] not in list(self.modules.keys()) + list(self.routines.keys()) + list(self.programs.keys()):
            raise F90toRstException(
                'Unknown %s: %s' %
                (block['block'], block['name']))

        name = block['name']
        blocktype = block['block']
        signature = '(%s)' % self.format_signature(
            block) if blocktype != 'program' else ''
        declaration = self.format_declaration(
            blocktype, '%(name)s%(signature)s' %
            locals(), indent=indent)
        #declaration = self.indent(indent)+'.. f:%(blocktype)s:: %(name)s%(signature)s\n\n'%locals()

        # Treat variables in comment (subroutines and functions only)
        comments = list(block['desc']) + ['']
        if blocktype != 'program':
            found = []
            for iline in range(len(comments)):
                if 'vardescmatch' in block:
                    m = block['vardescmatch'](comments[iline])
                    if m:
                        varname = m.group('varname')
                        found.append(varname)
                        if varname != '':
                            comments[iline] = self.format_argfield(
                                block['vars'][varname], block=block)
            for varname in block['args'] + block['sortvars']:
                if varname not in found:
                    comments.append(
                        self.format_argfield(
                            block['vars'][varname],
                            block=block))
                    found.append(varname)

        # Description
        description = self.format_lines(comments, indent + 1)

        # Add use of modules
        use = self.format_use(block, indent=indent + 1, short=True)

        # Add calls
        calls = []
        module = block.get('module')
        # - call froms
        if blocktype in ['function', 'subroutine']:
            if 'callfrom' in block and block['callfrom']:
                callfrom = []

                for fromname in block['callfrom']:
                    if fromname in self.routines:
                        cf = self.format_funcref(fromname, module)
                    else:
                        cf = ':f:prog:`%s`' % fromname
                    callfrom.append(cf)

                # callfrom += ', '.join([self.format_funcref(getattr(self,
                # routines[fn]['name'], module) for fn in block['callfrom']])
                callfrom = ':from: ' + ', '.join(callfrom)

                calls.append(callfrom)
        # - call tos
        if block['callto']:
            callto = ', '.join([self.format_funcref(fn, module)
                                for fn in block['callto']])
            #callto = ', '.join([self.format_funcref(self.routines[fn]['name'], module) for fn in block['callto']])
            if callto == '':
                callto = 'None'
            callto = ':to: ' + callto
            calls.append(callto)
        calls = '\n' + self.format_lines(calls, indent=indent + 1)
        return declaration + description + use + calls + '\n\n'

    format_function = format_routine
    format_subroutine = format_routine

    def format_quickaccess(self, module, indent=indent):
        """Format an abstract of all types, variables and routines of a module"""
        if not isinstance(module, six.string_types):
            module = module['name']

        # Title
        title = self.format_subsection('Quick access', indent=indent) + '\n'

        # Types
        decs = []
        tlist = sorted(self.get_blocklist('types', module))
        if tlist:
            decs.append(':Types: ' +
                        ', '.join([':f:type:`%s`' %
                                   tt['name'] for tt in tlist]))

        # Variables
        vlist = self.get_blocklist('variables', module)
        if vlist:
            decs.append(':Variables: ' +
                        ', '.join([':f:var:`%s`' %
                                   vv['name'] for vv in vlist]))

        # Functions and subroutines
        flist = self.get_blocklist('functions', module)
        if flist:
            decs.append(':Routines: ' +
                        ', '.join([':f:func:`~%s/%s`' %
                                   (module, ff['name']) for ff in flist]))

        if decs:
            return self.format_lines(title + '\n'.join(decs)) + '\n\n'
        return ''

    def format_types(self, block, indent=0):
        """Format the description of all fortran types"""
        types = []
        for subblock in block['body']:
            if subblock['block'] == 'type':
                types.append(self.format_type(subblock, indent=indent))
        if types:
            types = self.format_subsection(
                'Types', indent=indent) + '\n'.join(types)
        else:
            types = ''
        return types

    def format_variables(self, block, indent=0):
        """Format the description of all variables (global or module)"""
        variables = ''
        if block['vars']:
            varnames = block['sortvars']
            if block['block'] == 'module':
                varnames.sort()
            for varname in varnames:
                bvar = block['vars'][varname]
                variables += self.format_var(bvar, indent=indent)
            variables = self.format_subsection(
                'Variables', indent=indent) + variables + '\n\n'
        return variables

    def format_description(self, block, indent=0):
        """Format the description of an object"""
        description = ''
        if block['desc']:
            description = self.format_subsection('Description', indent=indent)
            description += self.format_lines(
                block['desc'], indent=indent, strip=True) + '\n'
        return description

    def format_routines(self, block, indent=0):
        """Format the list of all subroutines and functions"""
        routines = ''
        blocks = block if isinstance(block, list) else block['body']
        fdecs = []
        for subblock in blocks:  # block['body']:
            if subblock['block'] in ['function', 'subroutine']:
                fdecs.append(self.format_routine(subblock, indent))
        if fdecs:
            fdecs = '\n'.join(fdecs)
            routines = self.format_subsection(
                'Subroutines and functions', indent=indent) + fdecs
        return routines

    def format_module(self, block, indent=0):
        """Recursively format a module and its declarations"""

        # Declaration of the module
        if isinstance(block, six.string_types):
            if block not in self.modules:
                raise F90toRstException('Unknown module: %s' % block)
            block = self.modules[block]
        elif block['name'] not in self.modules:
            raise F90toRstException('Unknown module: %' % block['name'])

        modname = block['name']
        declaration = self.format_declaration(
            'module', modname, indent=indent, options=dict(
                synopsis=self.get_synopsis(block).strip() or None))

        # Description
        description = self.format_description(block, indent=indent)

        # Quick access
        quickaccess = self.format_quickaccess(modname, indent=indent)

        # Use of other modules
        use = self.format_use(block, indent=indent)

        # Types
        types = self.format_types(block, indent=indent)

        # Variables
        variables = self.format_variables(block, indent=indent)

        # Subroutines and functions
        routines = self.format_routines(block, indent=indent)

        return declaration + description + quickaccess + \
            use + types + variables + routines

    def format_srcfile(
            self,
            srcfile,
            indent=0,
            objtype=None,
            search_mode='basename',
            **kwargs):
        """Format all declaration of a file, except modules"""
        rst = ''
        if objtype is not None and not isinstance(objtype, (list, tuple)):
            objtype = [objtype]

        # Programs
        if objtype is None or 'program' in objtype:
            bprog = self.filter_by_srcfile(
                srcfile, objtype='program', mode=search_mode)
            if bprog:
                rst += self.format_subsection('Program', indent=indent) + '\n'
                rst += self.format_routine(bprog[0], indent=indent) + '\n'

        # Modules
        if objtype is None or 'module' in objtype:
            bmod = self.filter_by_srcfile(
                srcfile, objtype='module', mode=search_mode)
            if bmod:
                rst += self.format_subsection('Module', indent=indent) + '\n'
                rst += self.format_module(bmod[0], indent=indent) + '\n'

        # Functions and subroutines
        oal = ['function', 'subroutine']
        oo = [o for o in oal if o in objtype] if objtype is not None else oal
        if oo:
            brouts = self.filter_by_srcfile(
                srcfile, objtype=oo, mode=search_mode)
            rst += self.format_routines(brouts, indent=indent) + '\n'

        return rst

    def __getitem__(self, module):
        return self.format_module(self.modules[module])


# Sphinx directive ---

def list_files(fortran_src, exts=['f', 'f90', 'f95'], absolute=True):
    """Get the list of fortran files"""

    # Extensions (suffixes)
    if not isinstance(exts, list):
        exts = list(exts)
    for e in exts:
        if e.lower() not in exts:
            exts.append(e.lower())
        if e.upper() not in exts:
            exts.append(e.upper())
    exts = list(set(exts))

    # List the files using globs
    ffiles = []
    for fg in fortran_src:
        if not isinstance(fg, six.string_types):
            continue
        if os.path.isdir(fg):
            for ext in exts:
                ffiles.extend(glob(os.path.join(fg, '*.' + ext)))
        else:
            ffiles.extend(glob(fg))
    if absolute:
        ffiles = [os.path.abspath(ffile) for ffile in ffiles]
    ffiles.sort()
    return ffiles


def fortran_parse(app):
    env = app.builder.env
    if isinstance(app.config.fortran_src, (str, list)):
        logger.info(bold('parsing fortran sources...'), nonl=True)

        # Sources a list
        if not isinstance(app.config.fortran_src, list):
            app.config.fortran_src = [app.config.fortran_src]

        # All files
        ffiles = list_files(app.config.fortran_src, app.config.fortran_ext)

        # Parse files
        if not ffiles:
            logger.info(" no fortran files found")
            app.config._f90torst = None
        else:
            app.config.fortran_indent = fmt_indent(app.config.fortran_indent)
            app.config._f90torst = F90toRst(
                ffiles,
                ic=app.config.fortran_indent,
                ulc=app.config.fortran_title_underline,
                encoding=app.config.fortran_encoding)
            logger.info(' done')
        app._status.flush()

    else:
        logger.warning(
            "wrong list of fortran 90 source specifications: " + str(app.config.fortran_src))
        app.config._f90torst = None
#    app.config._f90files = []


def fmt_indent(string):
    if string is None:
        return
    if isinstance(string, int):
        string = ' ' * string
    if string == 'tab':
        string = '\t'
    elif string == 'space':
        string = ' '
    return string


# class fortran_module(nodes.General, nodes.Element):
    # pass


class FortranAutoModuleDirective(Directive):
    has_content = True
    option_spec = dict(title_underline=unchanged, indent=fmt_indent,
                       subsection_type=unchanged)
    required_arguments = 1
    optional_arguments = 0

    def run(self):

        # Get environment
        f90torst = self.state.document.settings.env.config._f90torst
        if f90torst is None:
            return []

        # Check module name
        module = self.arguments[0]
        if module not in f90torst.modules:
            self.state_machine.reporter.warning(
                'Wrong fortran module name: ' + module, line=self.lineno)
#            self.warn('Wrong fortran module name: '+module)

        # Options
        ic = f90torst.ic
        ulc = f90torst.ulc
        sst = f90torst.sst
        if self.options.get('indent'):
            f90torst.ic = self.options['indent']
        if self.options.get('title_underline'):
            f90torst.ulc = self.options['title_underline']
        if self.options.get('subsection_type'):
            f90torst.sst = self.options['subsection_type']

        # Get rst
        raw_text = f90torst.format_module(module)

        # Insert it
        source = self.state_machine.input_lines.source(
            self.lineno - self.state_machine.input_offset - 1)
        include_lines = string2lines(raw_text, convert_whitespace=1)
        self.state_machine.insert_input(include_lines, source)

        # Restore defaults
        if 'indent' in self.options:
            f90torst.ic = ic
        if 'title_underline' in self.options:
            f90torst.ulc = ulc
        if 'subsection_type' in self.options:
            f90torst.sst = sst

        return []


class FortranAutoObjectDirective(Directive):
    """Generic directive for fortran object auto-documentation

    Redefine :attr:`_warning` and :attr:`_objtype` attribute when subcassling.

    .. attribute:: _warning

        Warning message when object is not found, like:

        >>> _warning = 'Wrong function or subroutine name: %s'

    .. attribute:: _objtype

        Type of fortran object.
        If "toto" is set as object type, then :class:`F90toRst` must have
        attribute :attr:`totos` containg index of all related fortran objects,
        and method :meth:`format_totos` for formatting the object.

    """
    has_content = False
    option_spec = OrderedDict()
    required_arguments = 1
    optional_arguments = 0
    _warning = 'Wrong routine name: %s'
    _objtype = 'routine'

    def run(self):

        # Get environment
        f90torst = self.state.document.settings.env.config._f90torst
        if f90torst is None:
            return []

        # Check object name
        objname = self.arguments[0].lower()
        from sphinxfortran.fortran_domain import f_sep
        if f_sep in objname:
            objname = objname.split(f_sep)[-1]  # remove module name
        objects = getattr(f90torst, self._objtype + 's')
        if objname not in objects:
#            print(self._warning % objname)
            self.state_machine.reporter.warning(
                self._warning %
                objname, line=self.lineno)
#            self.warn(self._warning%objname)

        # Get rst
        raw_text = getattr(f90torst, 'format_' + self._objtype)(objname)

        # Check if inside module
        b = objects[objname]
        if 'parent_block' in b:
            curmod_text = '.. f:currentmodule:: %s\n\n' % b['parent_block']['name']
            raw_text = curmod_text + raw_text

        # Insert it
        source = self.state_machine.input_lines.source(
            self.lineno - self.state_machine.input_offset - 1)
        include_lines = string2lines(raw_text, convert_whitespace=1)
        self.state_machine.insert_input(include_lines, source)

        return []


class FortranAutoFunctionDirective(FortranAutoObjectDirective):
    _warning = 'Wrong function name: %s'
    _objtype = 'function'


class FortranAutoSubroutineDirective(FortranAutoObjectDirective):
    _warning = 'Wrong subroutine name: %s'
    _objtype = 'subroutine'


class FortranAutoTypeDirective(FortranAutoObjectDirective):
    _warning = 'Wrong type name: %s'
    _objtype = 'type'


class FortranAutoVariableDirective(FortranAutoObjectDirective):
    _warning = 'Wrong variable name: %s'
    _objtype = 'variable'


class FortranAutoProgramDirective(Directive):
    has_content = False
    option_spec = OrderedDict()
    required_arguments = 1
    optional_arguments = 0

    def run(self):

        # Get environment
        f90torst = self.state.document.settings.env.config._f90torst
        if f90torst is None:
            return []

        # Check routine name
        program = self.arguments[0].lower()
        if program not in f90torst.programs:
#            print('Wrong program name: ' + program)
            self.state_machine.reporter.warning(
                'Wrong program name: ' + program, line=self.lineno)
#            self.warning('Wrong program name: '+program)

        # Get rst
        raw_text = f90torst.format_routine(program)

        # Insert it
        source = self.state_machine.input_lines.source(
            self.lineno - self.state_machine.input_offset - 1)
        include_lines = string2lines(raw_text, convert_whitespace=1)
        self.state_machine.insert_input(include_lines, source)

        return []


class FortranAutoSrcfileDirective(Directive):
    has_content = False
    option_spec = dict(search_mode=unchanged, objtype=unchanged)
    required_arguments = 1
    optional_arguments = 0

    def run(self):

        # Get environment
        f90torst = self.state.document.settings.env.config._f90torst
        if f90torst is None:
            return []

        # Options
        search_mode = self.options.get('search_mode')
        objtype = self.options.get('objtype')
        if objtype:
            objtype = objtype.split(' ,')

        # Get rst
        srcfile = self.arguments[0].lower()
        raw_text = f90torst.format_srcfile(
            srcfile, search_mode=search_mode, objtype=objtype)
        if not raw_text:
            msg = 'No valid content found for file: ' + srcfile
#            print(msg)
            self.state_machine.reporter.warning(msg, line=self.lineno)
#            self.warning('No valid content found for file: '+srcfile)

        # Insert it
        source = self.state_machine.input_lines.source(
            self.lineno - self.state_machine.input_offset - 1)
        include_lines = string2lines(raw_text, convert_whitespace=1)
        self.state_machine.insert_input(include_lines, source)

        return []


def setup(app):

    app.add_object_type(
        'ftype',
        'ftype',
        indextemplate='pair: %s; Fortran type',
    )
    app.add_object_type(
        'fvar',
        'fvar',
        indextemplate='pair: %s; Fortran variable',
    )

    app.add_config_value('fortran_title_underline', '-', False)
    app.add_config_value('fortran_indent', 4, False)
    app.add_config_value('fortran_subsection_type', 'rubric', False)
    app.add_config_value('fortran_src', ['.'], False)
    app.add_config_value('fortran_ext', ['f90', 'f95'], False)
    app.add_config_value('fortran_encoding', 'utf8', False)

    #app.add_directive('fortran_module', IncludeFortranDirective)
    FortranDomain.directives.update(
        automodule=FortranAutoModuleDirective,
        autoroutine=FortranAutoObjectDirective,
        autofunction=FortranAutoFunctionDirective,
        autosubroutine=FortranAutoSubroutineDirective,
        autoprogram=FortranAutoProgramDirective,
        autosrcfile=FortranAutoSrcfileDirective,
    )
    app.connect('builder-inited', fortran_parse)
