# -*- coding: utf-8 -*-
"""
A fortran domain for sphinx

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
from builtins import zip
from builtins import str
from builtins import object
import re
from collections import OrderedDict

from docutils import nodes
from docutils.parsers.rst import Directive, directives

from sphinx import addnodes, version_info
from sphinx.roles import XRefRole
from sphinx.locale import _
from sphinx.domains import Domain, ObjType, Index
from sphinx.directives import ObjectDescription
from sphinx.util.nodes import make_refnode
from sphinx.util.docfields import Field, GroupedField, TypedField, DocFieldTransformer, _is_single_paragraph

import six


# FIXME: surlignage en jaune de la recherche inactive si "/" dans target

# Utilities

def convert_arithm(node, expr, modname=None, nodefmt=nodes.Text):
    """Format an arithmetic expression for a node"""
    ops = re.findall(r'(\W+)', expr)
    nums = re.split(r'\W+', expr)
    if len(nums) != len(ops):
        ops.append('')
    for num, op in zip(nums, ops):
        if num:
            if num[0].isalpha():
                refnode = addnodes.pending_xref(
                    '', refdomain='f', reftype='var', reftarget=num,
                    modname=modname)
                refnode += nodefmt(num, num)
                node += refnode
            else:
                node += nodefmt(num, num)
        if op:
            op = op.replace(':', '*')
            node += nodefmt(op, op)


def parse_shape(shape):
    if not shape:
        return
    if not shape.startswith('('):
        shape = '(' + shape
    if not shape.endswith(')'):
        shape += ')'
    return shape


def add_shape(node, shape, modname=None, nodefmt=nodes.Text):
    """Format a shape expression for a node"""
    dims = re.split(r'\s*,\s*', shape.strip('( )'))
    node += nodefmt(' (', ' (')
    convert_arithm(node, shape.strip('( )'), modname=modname, nodefmt=nodefmt)
    node += nodefmt(')', ')')

#class fortranfield(nodes.Admonition, nodes.TextElement): pass

# Doc fields


re_name_shape = re.compile(r'(\w+)(\(.+\))?')

re_fieldname_match = re.compile(
    r'(?P<type>\b\w+\b(?P<kind>\s*\(.*\))?)?\s*(?P<name>\b\w+\b)\s*(?P<shape>\(.*\))?\s*(?P<sattrs>\[.+\])?').match


class FortranField(Field):
    def make_xref(self, rolename, domain, target, innernode=nodes.emphasis,
                  modname=None, typename=None):
        if not rolename:
            return innernode(target, target)
        refnode = addnodes.pending_xref(
            '',
            refdomain=domain,
            refexplicit=False,
            reftype=rolename,
            reftarget=target,
            modname=modname,
            typename=typename)
        refnode += innernode(target, target)
        return refnode


class FortranCallField(FortranField):
    is_grouped = True

    def __init__(self, name, names=(), label=None, rolename=None):
        Field.__init__(self, name, names, label, True, rolename)

    def make_field(self, types, domain, items, **kwargs):

        fieldname = nodes.field_name('', self.label)
        #par = Field.make_field(self, types, domain, items[0])
        par = nodes.paragraph()
        for i, item in enumerate(items):
            if i:
                par += nodes.Text(' ')
            par += item[1]  # Field.make_field(self, types, domain, item)
        fieldbody = nodes.field_body('', par)
        return nodes.field('', fieldname, fieldbody)


class FortranCompleteField(FortranField, GroupedField):
    """
    A doc field that is grouped and has type information for the arguments.  It
    always has an argument.  The argument can be linked using the given
    *rolename*, the type using the given *typerolename*.

    Two uses are possible: either parameter and type description are given
    separately, using a field from *names* and one from *typenames*,
    respectively, or both are given using a field from *names*, see the example.

    Example::

       :param foo: description of parameter foo
       :type foo:  SomeClass

       -- or --

       :param SomeClass foo: description of parameter foo
    """
    is_typed = 2

    def __init__(self, name, names=(), typenames=(), label=None,
                 rolename=None, typerolename=None,
                 shapenames=None, attrnames=None,
                 prefix=None,
                 strong=True,
                 can_collapse=False):
        GroupedField.__init__(self, name, names, label, rolename, can_collapse)
        self.typenames = typenames
        self.typerolename = typerolename
        self.shapenames = shapenames
        self.attrnames = attrnames
        self.prefix = prefix
        if strong:
            self.namefmt = nodes.strong
        else:
            self.namefmt = addnodes.desc_name

    def make_field(self, types, domain, items, shapes=None, attrs=None,
                   modname=None, typename=None):
        def handle_item(fieldarg, content):
            par = nodes.paragraph()
            if self.prefix:
                par += self.namefmt(self.prefix, self.prefix)

            par += self.make_xref(self.rolename,
                                  domain,
                                  fieldarg,
                                  self.namefmt,
                                  modname=modname,
                                  typename=typename)
            #par += self.namefmt(fieldarg, fieldarg)

            fieldtype = types.pop(fieldarg, None)
            fieldshape = shapes and shapes.pop(fieldarg, None)
            fieldattrs = attrs and attrs.pop(fieldarg, None)
            if fieldshape:
                shape = parse_shape(fieldshape[0].astext())
                #par += nodes.Text(' %s'%shape)
                add_shape(par, shape, modname=modname)
            if fieldtype or fieldattrs:
                par += nodes.emphasis(' [', ' [')
            if fieldtype:
                if len(fieldtype) == 1 and isinstance(
                        fieldtype[0], nodes.Text):
                    thistypename = fieldtype[0].astext()
                    #typename = u''.join(n.astext() for n in fieldtype)
                    par += self.make_xref(self.typerolename,
                                          domain,
                                          thistypename,
                                          modname=modname,
                                          typename=typename)
                else:
                    par += fieldtype
            if fieldattrs:
                if fieldtype:
                    par += nodes.emphasis(',', ',')
                par += fieldattrs
            if fieldtype or fieldattrs:
                par += nodes.emphasis(']', ']')
            if content:
                par += nodes.Text(' :: ')
                par += content
            return par

        if len(items) == 1 and self.can_collapse:
            fieldarg, content = items[0]
            bodynode = handle_item(fieldarg, content)
        else:
            bodynode = self.list_type()
            for fieldarg, content in items:
                bodynode += nodes.list_item('', handle_item(fieldarg, content))
        label = self.label or ''
        fieldname = nodes.field_name('', label)
        fieldbody = nodes.field_body('', bodynode)
        return nodes.field('', fieldname, fieldbody)


class FortranDocFieldTransformer(DocFieldTransformer):
    """
    Transforms field lists in "doc field" syntax into better-looking
    equivalents, using the field type definitions given on a domain.
    """

    def __init__(self, directive, modname=None, typename=None):
        self.domain = directive.domain
        if '_doc_field_type_map' not in directive.__class__.__dict__:
            directive.__class__._doc_field_type_map = \
                self.preprocess_fieldtypes(directive.__class__.doc_field_types)
        self.typemap = directive._doc_field_type_map
        self.modname = modname
        self.typename = typename

    def preprocess_fieldtypes(self, types):
        typemap = OrderedDict()
        for fieldtype in types:
            for name in fieldtype.names:
                typemap[name] = fieldtype, False
            if fieldtype.is_typed:
                for name in fieldtype.typenames:
                    typemap[name] = fieldtype, 'types'
                for name in fieldtype.shapenames:
                    typemap[name] = fieldtype, 'shapes'
                for name in fieldtype.attrnames:
                    typemap[name] = fieldtype, 'attrs'
        return typemap

    def scan_fieldarg(self, fieldname):
        """Extract type, name, shape and attributes from a field name.

        :Some possible syntaxes:

            - ``p name``
            - ``p type name(shape) [attr1,attr2]``
            - ``p type name``
            - ``p name [attr1, attr2]``

        :Returns: ``name, shape, type, list of attributes``.
            if no shape is specified, it is set to ``None``,
        """
        m = re_fieldname_match(fieldname.strip())
        if not m:
            raise ValueError(
                'Wrong field (%s). It must have at least one parameter name and one argument' %
                fieldname)
        ftype, kind, name, shape, attrs = m.groups()
        attrs = attrs and attrs[1:-1]
        # if attrs:
        #attrs = [a.strip() for a in attrs[1:-1].split(',')]
        # else:
        #attrs = []
        return name, shape, ftype, attrs

    def transform(self, node):
        """Transform a single field list *node*."""
        typemap = self.typemap
        fmodname = self.modname
        ftypename = self.typename

        entries = []
        groupindices = OrderedDict()
        types = OrderedDict()
        shapes = OrderedDict()
        attrs = OrderedDict()

        # step 1: traverse all fields and collect field types and content
        for field in node:

            fieldname, fieldbody = field
            try:
                # split into field type and argument
                fieldtype, fieldarg = fieldname.astext().split(None, 1)
            except ValueError:
                # maybe an argument-less field type?
                fieldtype, fieldarg = fieldname.astext(), ''
            typedesc, is_typefield = typemap.get(fieldtype, (None, None))

            # sort out unknown fields
            if typedesc is None:  # or typedesc.has_arg != bool(fieldarg):
                # either the field name is unknown, or the argument doesn't
                # match the spec; capitalize field name and be done with it
                new_fieldname = fieldtype.capitalize() + ' ' + fieldarg
                fieldname[0] = nodes.Text(new_fieldname)
                entries.append(field)
                continue

            typename = typedesc.name

            # collect the content, trying not to keep unnecessary paragraphs
            if _is_single_paragraph(fieldbody):
                content = fieldbody.children[0].children
            else:
                content = fieldbody.children

            # if the field specifies a type, put it in the types collection
            if is_typefield:
                # filter out only inline nodes; others will result in invalid
                # markup being written out
                content = [n for n in content if isinstance(n, nodes.Inline) or
                           isinstance(n, nodes.Text)]
                if content:
                    eval(is_typefield).setdefault(
                        typename, OrderedDict())[fieldarg] = content
                continue

            # also support syntax like ``:param type name [attrs]:``
            if typedesc.is_typed == 2:
                argname, argshape, argtype, argattrs = self.scan_fieldarg(
                    fieldarg)
                if argtype:
                    types.setdefault(typename, OrderedDict())[argname] = \
                        [nodes.Text(argtype)]
                if argshape:
                    shapes.setdefault(typename, OrderedDict())[argname] = \
                        [nodes.Text(argshape)]
                if argattrs:
                    attrs.setdefault(typename, OrderedDict())[argname] = \
                        [nodes.emphasis(argattrs, argattrs)]
                fieldarg = argname
            elif typedesc.is_typed:
                try:
                    argtype, argname = fieldarg.split(None, 1)
                except ValueError:
                    pass
                else:
                    types.setdefault(typename, OrderedDict())[argname] = \
                        [nodes.Text(argtype)]
                    fieldarg = argname

            # grouped entries need to be collected in one entry, while others
            # get one entry per field
            if typedesc.is_grouped:
                if typename in groupindices:
                    group = entries[groupindices[typename]]
                else:
                    groupindices[typename] = len(entries)
                    group = [typedesc, []]
                    entries.append(group)
                group[1].append(typedesc.make_entry(fieldarg, content))
            else:
                entries.append([typedesc,
                                typedesc.make_entry(fieldarg, content)])

        # step 2: all entries are collected, construct the new field list
        new_list = nodes.field_list()
        for entry in entries:
            if isinstance(entry, nodes.field):
                # pass-through old field
                new_list += entry
            else:
                fieldtype, content = entry
                fieldtypes = types.get(fieldtype.name, OrderedDict())
                fieldshapes = shapes.get(fieldtype.name, OrderedDict())
                fieldattrs = attrs.get(fieldtype.name, OrderedDict())
                new_list += fieldtype.make_field(fieldtypes,
                                                 self.domain,
                                                 content,
                                                 shapes=fieldshapes,
                                                 attrs=fieldattrs,
                                                 modname=fmodname,
                                                 typename=ftypename)

        node.replace_self(new_list)


# REs for Fortran signatures
f_sep = '/'
f_sig_re = re.compile(
    r'''^ (\w+(?:[^%%%(f_sep)s]%(f_sep)s\w+))? \s*          # type
          (\b(?:subroutine|function))?  \s*             # objtype
          (\b\w+%(f_sep)s)?              # module name
          (\b\w+%%)?              # type name
          (\b\w+)  \s*             # thing name
          (?: \((.*)\))?           # optional: arguments
           $                   # and nothing more
          ''' % dict(f_sep=f_sep), re.VERBOSE + re.I)


# Directives

# RE to split at word boundaries
wsplit_re = re.compile(r'(\W+)')
f_type_re = re.compile(r'^([\w]+).*$')

f_paramlist_re = re.compile(r'([\[\],])')  # split at '[', ']' and ','


# def fortran_rsplit(fullname):
#    items = [item for item in f_separator.findall(fullname)]
#    return ''.join(items[:-2]), items[-1]

def _pseudo_parse_arglist(signode, arglist):
    """"Parse" a list of arguments separated by commas.

    Arguments can have "optional" annotations given by enclosing them in
    brackets.  Currently, this will split at any comma, even if it's inside a
    string literal (e.g. default argument value).
    """
    paramlist = addnodes.desc_parameterlist()
    stack = [paramlist]
    try:
        for argument in arglist.split(','):
            argument = argument.strip()
            ends_open = ends_close = 0
            while argument.startswith('['):
                stack.append(addnodes.desc_optional())
                stack[-2] += stack[-1]
                argument = argument[1:].strip()
            while argument.startswith(']'):
                stack.pop()
                argument = argument[1:].strip()
            while argument.endswith(']'):
                ends_close += 1
                argument = argument[:-1].strip()
            while argument.endswith('['):
                ends_open += 1
                argument = argument[:-1].strip()
            if argument:
                stack[-1] += addnodes.desc_parameter(argument, argument)
            while ends_open:
                stack.append(addnodes.desc_optional())
                stack[-2] += stack[-1]
                ends_open -= 1
            while ends_close:
                stack.pop()
                ends_close -= 1
        if len(stack) != 1:
            raise IndexError
    except IndexError:
        # if there are too few or too many elements on the stack, just give up
        # and treat the whole argument list as one argument, discarding the
        # already partially populated paramlist node
        signode += addnodes.desc_parameterlist()
        signode[-1] += addnodes.desc_parameter(arglist, arglist)
    else:
        signode += paramlist


class FortranObject(ObjectDescription):
    """
    Description of a general Fortran object.
    """
    option_spec = {
        'noindex': directives.flag,
        'module': directives.unchanged,
        'type': directives.unchanged,
        'shape': parse_shape,
        'attrs': directives.unchanged,
    }

    doc_field_types = [
        FortranCompleteField('parameter', label=_('Parameters'),
                             names=(
            'p', 'param', 'parameter', 'a', 'arg', 'argument'),
            # rolename='var',
            typerolename='type',
            typenames=('paramtype', 'type', 'ptype'),
            shapenames=('shape', 'pshape'),
            attrnames=('attrs', 'pattrs', 'attr'),
            can_collapse=True),
        FortranCompleteField('optional', label=_('Options'),
                             names=(
            'o', 'optional', 'opt', 'keyword', 'option'),
            # rolename='var',
            typerolename='type',
            typenames=('optparamtype', 'otype'),
            shapenames=('oshape',),
            attrnames=('oattrs', 'oattr'),
            can_collapse=True),
        FortranCompleteField('typefield', label=_('Type fields'),
                             names=('f', 'field', 'typef', 'typefield'),
                             # rolename='typef',
                             typerolename='type',
                             typenames=('fieldtype', 'ftype'),
                             shapenames=('fshape',),
                             attrnames=('fattrs', 'fattr'),
                             prefix='% ',
                             strong=False,
                             can_collapse=False),
        FortranCompleteField('return', label=_('Return'),
                             names=('r', 'return', 'returns'),
                             typerolename='type',
                             typenames=('returntype', 'rtype'),
                             shapenames=('rshape',),
                             attrnames=('rattrs', 'rattr'),
                             can_collapse=True),
        FortranCallField('calledfrom', label=_('Called from'),
                         names=('calledfrom', 'from')),
        FortranCallField('callto', label=_('Call to'),
                         names=('callto', 'to')),
    ]

    # These Fortran types aren't described anywhere, so don't try to create
    # a cross-reference to them
    stopwords = set(('float', 'integer', 'character', 'double', 'long'))

    _parens = ''

#    def _parse_type(self, node, ftype):
#        m = f_type_re.match(ftype)
#        tnode = nodes.Text(ftype, ftype)
#        modname = self.options.get(
#            'module', self.env.temp_data.get('f:module'))
#        if m :
#            ftype = m.groups(0)
#            if ftype not in self.stopwords:
#                pnode = addnodes.pending_xref(
#                    '', refdomain='f', reftype='type', reftarget=ftype,
#                    modname=modname)
#                pnode += tnode
#                node += pnode
#            else:
#                node += tnode
#        else:
#            node += tnode

    def get_signature_prefix(self, sig):
        """
        May return a prefix to put before the object name in the signature.
        """
        return ''

    def needs_arglist(self):
        """
        May return true if an empty argument list is to be generated even if
        the document contains none.
        """
        return False

    def handle_signature(self, sig, signode):
        """
        Transform a Fortran signature into RST nodes.
        Returns (fully qualified name of the thing, classname if any).

        If inside a class, the current class name is handled intelligently:
        * it is stripped from the displayed name if present
        * it is added to the full name (return value) if not present
        """
        m = f_sig_re.match(sig)
        if m is None:
            raise ValueError
        ftype, objtype, modname, typename, name, arglist = m.groups()
        if not typename:
            typename = ""

        # determine module, type, shape and attributes
        modname = (modname and modname[:-1]) or self.options.get(
            'module', self.env.temp_data.get('f:module'))
        if typename:
            name = typename[:-1]
        attrs = self.options.get('attrs')
        shape = parse_shape(self.options.get('shape'))
        ftype = ftype or self.options.get('type')
        if self.objtype == 'typefield' and not typename:
            raise ValueError

        #if typename: name = typename+'%'+name

        #fullname = name
        # if modname:
        if self.objtype == 'program':
            fullname = name
        else:
            fullname = (modname or '_') + f_sep + name

        signode['module'] = modname
        signode['type'] = typename
        signode['fullname'] = fullname

        # Add "function" or "subroutine" tag
        sig_prefix = self.get_signature_prefix(sig)
        if objtype or sig_prefix:
            objtype = objtype or sig_prefix
            signode += addnodes.desc_annotation(objtype + ' ', objtype + ' ')

        # Add module
        if self.env.config.add_module_names and modname and self.objtype != 'typefield':
            nodetext = modname + f_sep
            signode += addnodes.desc_addname(nodetext, nodetext)

        # Add name
        signode += addnodes.desc_name(name, name)

        # In the parenthesis
        if self.needs_arglist():  # call for functions and subroutines
            if arglist:  # Calling arguments
                _pseudo_parse_arglist(signode, arglist)
            elif self.needs_arglist():  # for callables, add an empty parameter list
                signode += addnodes.desc_parameterlist()
        # Declare shape instead of arguments (variables)
        elif arglist and not shape:
            shape = arglist

        # Add remaining
        self.add_shape_and_attrs(signode, modname, ftype, shape, attrs)

        return fullname, ftype

    def add_shape_and_attrs(self, signode, modname, ftype, shape, attrs):
        # add shape
        if shape:
            add_shape(signode, shape, modname=modname)
            #signode += nodes.Text(' '+shape)
        # add type ('float', 'interger', etc)
        if ftype or attrs:
            signode += nodes.emphasis(' [', ' [')
        if ftype:
            refnode = addnodes.pending_xref(
                '', refdomain='f', reftype='type', reftarget=ftype,
                modname=modname,)
            refnode += nodes.emphasis(ftype, ftype)
            signode += refnode
            #tnode = addnodes.desc_type(ftype, ftype)
            # tnode +=
            #signode += addnodes.desc_type(ftype, ftype)
            # signode +=
    #        signode += addnodes.desc_type('', '')
    #        self._parse_type(signode[-1], ftype)
        if attrs:
            if ftype:
                signode += nodes.emphasis(',', ',')
            for iatt, att in enumerate(re.split(r'\s*,\s*', attrs)):
                if iatt:
                    signode += nodes.emphasis(',', ',')
                if att.startswith('parameter'):
                    value = att.split('=')[1]
                    signode += nodes.emphasis('parameter=', 'parameter=')
                    convert_arithm(signode, value, modname=modname)
                else:
                    signode += nodes.emphasis(att, att)
            #signode += nodes.emphasis(attrs, attrs)

        if ftype or attrs:
            signode += nodes.emphasis(']', ']')

    def add_target_and_index(self, name, sig, signode):
        # modname = self.options.get(
            # 'module', self.env.temp_data.get('f:module'))
        modname = signode.get(
            'module', self.env.temp_data.get('f:module'))
#        fullname = (modname and modname + '/' or '') + name[0]
        fullname = 'f' + f_sep + name[0]

        # note target
        if fullname not in self.state.document.ids:
            signode['names'].append(fullname)
            signode['ids'].append(fullname)
            signode['first'] = (not self.names)
            self.state.document.note_explicit_target(signode)
            objects = self.env.domaindata['f']['objects']
            if fullname in objects:
                self.env.warn(
                    self.env.docname,
                    'duplicate object description of %s, ' % fullname +
                    'other instance in ' +
                    self.env.doc2path(objects[fullname][0]),
                    self.lineno)
            objects[fullname] = (self.env.docname, self.objtype)
        indextext = self.get_index_text(modname, fullname)
        if indextext:
            # self.indexnode['entries'].append(('single', indextext,
                                              # fullname, fullname,None))
            self.indexnode['entries'].append(
                FortranCreateIndexEntry(
                    indextext, fullname, fullname))

    def before_content(self):
        # needed for automatic qualification of fields (reset in subclasses)
        self.typename_set = False

    def after_content(self):
        if self.typename_set:
            self.env.temp_data['f:type'] = None

    def get_index_text(self, modname, name):
        add_modules = self.env.config.add_module_names
        if name.startswith('f' + f_sep):
            name = name[2:]
        mn = modname or '_'
        sobj = ''
        if name.startswith(mn + f_sep):
            name = name[len(mn) + 1:]
        if self.objtype == 'type':
            sobj = _('fortran type')
        if self.objtype == 'typefield':
            sobj = _('fortran type field')
        elif self.objtype == 'variable':
            sobj = _('fortran variable')
        elif self.objtype == 'subroutine':
            sobj = _('fortran subroutine')
        elif self.objtype == 'function':
            sobj = _('fortran function')
        elif self.objtype == 'module':
            sobj = _('fortran module')
            modname = ''
        elif self.objtype == 'program':
            sobj = _('fortran program')
            modname = ''
        sinmodule = (
            _(' in module %s') %
            modname) if modname and add_modules else ''
        return '%s%s (%s%s)' % (name, self._parens, sobj, sinmodule)


class FortranSpecial(object):
    def get_signature_prefix(self, sig):
        """
        May return a prefix to put before the object name in the signature.
        """
        return self.objtype + ' '


class WithFortranDocFieldTransformer(object):
    def run(self):
        """Same as :meth:`sphinx.directives.ObjectDescription`
        but using :class:`FortranDocFieldTransformer`"""
        if ':' in self.name:
            self.domain, self.objtype = self.name.split(':', 1)
        else:
            self.domain, self.objtype = '', self.name
        if not hasattr(self, 'env'):
            self.env = self.state.document.settings.env
        self.indexnode = addnodes.index(entries=[])

        node = addnodes.desc()
        node.document = self.state.document
        node['domain'] = self.domain
        # 'desctype' is a backwards compatible attribute
        node['objtype'] = node['desctype'] = self.objtype
        node['noindex'] = noindex = ('noindex' in self.options)

        self.names = []
        signatures = self.get_signatures()
        for i, sig in enumerate(signatures):
            # add a signature node for each signature in the current unit
            # and add a reference target for it
            signode = addnodes.desc_signature(sig, '')
            signode['first'] = False
            node.append(signode)
            try:
                # name can also be a tuple, e.g. (classname, objname);
                # this is strictly domain-specific (i.e. no assumptions may
                # be made in this base class)
                name = self.handle_signature(sig, signode)
            except ValueError:
                # signature parsing failed
                signode.clear()
                signode += addnodes.desc_name(sig, sig)
                continue  # we don't want an index entry here
            if not isinstance(name[0], six.string_types):
                name = (str(name), name[1])
            if not noindex and name not in self.names:
                # only add target and index entry if this is the first
                # description of the object with this name in this desc block
                self.names.append(name)
                self.add_target_and_index(name, sig, signode)

        modname = signode.get('module')
        typename = signode.get('type')
        contentnode = addnodes.desc_content()
        node.append(contentnode)
        if self.names:
            # needed for association of version{added,changed} directives
            self.env.temp_data['object'] = self.names[0]
        self.before_content()
        self.state.nested_parse(self.content, self.content_offset, contentnode)
        FortranDocFieldTransformer(
            self,
            modname=modname,
            typename=typename).transform_all(contentnode)
        self.env.temp_data['object'] = None
        self.after_content()
        return [self.indexnode, node]


class FortranType(
        FortranSpecial,
        WithFortranDocFieldTransformer,
        FortranObject):
    def before_content(self):
        FortranObject.before_content(self)
        if self.names:
            self.env.temp_data['f:type'] = self.names[0][0].split(f_sep)[-1]
            self.typename_set = True


class FortranTypeField(FortranObject):
    # def handle_signature(self, sig, signode):
        # """
        # Transform a Fortran signature into RST nodes.
        # Returns (fully qualified name of the thing, classname if any).

        # If inside a class, the current class name is handled intelligently:
        # * it is stripped from the displayed name if present
        # * it is added to the full name (return value) if not present
        # """
        #m = f_sig_re.match(sig)
        # if m is None:
            #raise ValueError
        #ftype, objtype, modname, typename, name, arglist = m.groups()
        # print 'handle_signature', ftype, objtype, modname, typename, name, arglist
        #if not typename: typename = ""

        # determine module and type
        # modname = (modname and modname[:-1]) or self.options.get(
            # 'module', self.env.temp_data.get('f:module'))
        # typename = (typename and typename[:-1]) or self.options.get(
            # 'type', self.env.temp_data.get('f:type'))
        # print ' mod type', modname, typename
        # print self.objtype
        # if self.objtype=='typefield' and not typename:
            #raise ValueError

        #if typename: name = typename+'%'+name

        #fullname = name
        # if modname:
            #fullname = modname + f_sep + name

        #signode['module'] = modname
        #signode['type'] = typename
        #signode['fullname'] = fullname

        # Fill node
        #signode += addnodes.desc_name(name, name)
        #shape = self.options.get('shape')
        #if shape: signode += nodes.Text(shape, shape)
        #ftype = self.options.get('type', ftype)
        #attr= self.options.get('attr')
        # if ftype or attr:
            #signode += nodes.Text(' :: ', ' :: ')
            #if ftype: signode += nodes.emphasis('', ftype)
            #if attr: signode += nodes.literal('', '['+attr+']')
        # if self.content:
            #signode += nodes.Text(': ', ': ')
            #argnodes, msgs = self.state.inline_text(' '.join(self.content), self.lineno)
            #signode += argnodes
            #signode += msgs

        # return fullname, ftype

    def before_content(self):
        FortranObject.before_content(self)
        lastname = self.names and self.names[-1][1]
        if lastname and not self.env.temp_data.get('f:type'):
            self.env.temp_data['f:type'] = lastname.split(f_sep)[-1]
            self.typename_set = True


class FortranProgram(
        FortranSpecial,
        WithFortranDocFieldTransformer,
        FortranObject):
    pass


class FortranWithSig(
        FortranSpecial,
        WithFortranDocFieldTransformer,
        FortranObject):
    """
    Description of a function of subroutine
    """
    _parens = '()'

    def needs_arglist(self):
        return True

    def get_signature_prefix(self, sig):
        """
        May return a prefix to put before the object name in the signature.
        """
        return self.objtype + ' '


class FortranField(Directive):
    """
    Directive to describe a change/addition/deprecation in a specific version.
    """

    has_content = True
    required_arguments = 1
    optional_arguments = 0
    final_argument_whitespace = True
    option_spec = {
        'type': directives.unchanged,
        'shape': parse_shape,
        'attrs': directives.unchanged,
    }

    def run(self):
        from docutils import nodes
        node = nodes.paragraph()
        node += addnodes.desc_name(self.arguments[0], self.arguments[0])
        shape = self.options.get('shape')
        if shape:
            #node += nodes.Text(shape, shape)
            add_shape(node, shape)
        type = self.options.get('type')
        attrs = self.options.get('attrs')
        if type or attrs:
            node += nodes.Text(' :: ', ' :: ')
            if type:
                node += nodes.emphasis('', type)
            if attr:
                node += nodes.literal('', '[' + attr + ']')
        if self.content:
            node += nodes.Text(': ', ': ')
            argnodes, msgs = self.state.inline_text(
                ' '.join(self.content), self.lineno)
            node += argnodes
            node += msgs
        ret = [node]
#        env = self.state.document.settings.env
#        env.note_versionchange(node['type'], node['version'], node, self.lineno)
        return ret


class FortranModule(Directive):
    """
    Directive to mark description of a new module.
    """

    has_content = False
    required_arguments = 1
    optional_arguments = 0
    final_argument_whitespace = False
    option_spec = {
        'platform': lambda x: x,
        'synopsis': lambda x: x,
        'noindex': directives.flag,
        'deprecated': directives.flag,
    }

    def run(self):
        env = self.state.document.settings.env
        modname = self.arguments[0].strip()
        noindex = 'noindex' in self.options
        env.temp_data['f:module'] = modname
        env.domaindata['f']['modules'][modname] = \
            (env.docname, self.options.get('synopsis', ''),
             self.options.get('platform', ''), 'deprecated' in self.options)
        env.domaindata['f']['objects']['f' + \
            f_sep + modname] = (env.docname, 'module')
        #targetnode = nodes.target('', '', ids=['module-' + modname], ismod=True)
        targetnode = nodes.target(
            '', '', ids=[
                'f' + f_sep + modname], ismod=True)
        self.state.document.note_explicit_target(targetnode)
        ret = [targetnode]
        # XXX this behavior of the module directive is a mess...
        if 'platform' in self.options:
            platform = self.options['platform']
            node = nodes.paragraph()
            node += nodes.emphasis('', _('Platforms: '))
            node += nodes.Text(platform, platform)
            ret.append(node)
        # the synopsis isn't printed; in fact, it is only used in the
        # modindex currently
        if not noindex:
            indextext = _('%s (module)') % modname
            # inode = addnodes.index(entries=[('single', indextext,
            # 'module-' + modname, modname)])
            # inode = addnodes.index(entries=[('single', indextext,
            # 'f' + f_sep + modname, modname)])
            inode = addnodes.index(
                entries=[
                    FortranCreateIndexEntry(
                        indextext,
                        'f' + f_sep + modname,
                        modname)])
            ret.append(inode)
        return ret


def FortranCreateIndexEntry(indextext, fullname, modname):
    # See https://github.com/sphinx-doc/sphinx/issues/2673
    if version_info < (1, 4):
        return ('single', indextext, fullname, modname)
    else:
        return ('single', indextext, fullname, modname, None)


class FortranCurrentModule(Directive):
    """
    This directive is just to tell Sphinx that we're documenting
    stuff in module foo, but links to module foo won't lead here.
    """

    has_content = False
    required_arguments = 0
    optional_arguments = 1
    final_argument_whitespace = False
    option_spec = OrderedDict()

    def run(self):
        env = self.state.document.settings.env
        modname = self.arguments and (
            self.arguments[0] or self.arguments[0].strip()) or None
        if modname:
            env.temp_data['f:module'] = None
        else:
            env.temp_data['f:module'] = modname
        return []


class FortranXRefRole(XRefRole):
    def process_link(self, env, refnode, has_explicit_title, title, target):
        refnode['f:module'] = env.temp_data.get('f:module')
        refnode['f:type'] = env.temp_data.get('f:type')
        if not has_explicit_title:
            title = title.lstrip('.')   # only has a meaning for the target
            target = target.lstrip('~')  # only has a meaning for the title
            # if the first character is a tilde, don't display the module/class
            # parts of the contents
            if title[0:1] == '~':
                title = title[1:].split(f_sep)[-1]
            # elif '%' not in title and refnode['f:type']:
                #title = '%s%%%s'%(refnode['f:type'], title)
        # if the first character is a sep, search more specific namespaces first
        # else search builtins first
        if target.startswith(f_sep):
            target = target[1:]
            refnode['refspecific'] = True
        return title, target


class FortranModuleIndex(Index):
    """
    Index subclass to provide the Fortran module index.
    """

    name = 'modindex'
    localname = _('Fortran Module Index')
    shortname = _('fortran modules')

    def generate(self, docnames=None):
        content = OrderedDict()
        # list of prefixes to ignore
        ignores = self.domain.env.config['modindex_common_prefix']
        ignores = sorted(ignores, key=len, reverse=True)
        # list of all modules, sorted by module name
        modules = sorted(six.iteritems(self.domain.data['modules']),
                         key=lambda x: x[0].lower())
        # sort out collapsable modules
        prev_modname = ''
        num_toplevels = 0
        for modname, (docname, synopsis, platforms, deprecated) in modules:
            if docnames and docname not in docnames:
                continue

            for ignore in ignores:
                if modname.startswith(ignore):
                    modname = modname[len(ignore):]
                    stripped = ignore
                    break
            else:
                stripped = ''

            # we stripped the whole module name?
            if not modname:
                modname, stripped = stripped, ''

            entries = content.setdefault(modname[0].lower(), [])

            package = modname.split(f_sep)[0]
            if package != modname:
                # it's a submodule
                if prev_modname == package:
                    # first submodule - make parent a group head
                    entries[-1][1] = 1
                elif not prev_modname.startswith(package):
                    # submodule without parent in list, add dummy entry
                    entries.append([stripped + package, 1, '', '', '', '', ''])
                subtype = 2
            else:
                num_toplevels += 1
                subtype = 0

            qualifier = deprecated and _('Deprecated') or ''
            # entries.append([stripped + modname, subtype, docname,
            #'module-' + stripped + modname, platforms,
            # qualifier, synopsis])
            entries.append([stripped + modname, subtype, docname,
                            'f' + f_sep + stripped + modname, platforms,
                            qualifier, synopsis or ''])
            prev_modname = modname

        # apply heuristics when to collapse modindex at page load:
        # only collapse if number of toplevel modules is larger than
        # number of submodules
        collapse = len(modules) - num_toplevels < num_toplevels

        # sort by first letter
        content = sorted(six.iteritems(content))

        return content, collapse


class FortranDomain(Domain):
    """Fortran language domain."""
    name = 'f'
    label = 'Fortran'
    object_types = {
        'program': ObjType(_('program'), 'prog'),
        'type': ObjType(_('type'), 'type'),
        'variable': ObjType(_('variable'), 'var'),
        'function': ObjType(_('function'), 'func'),
        'subroutine': ObjType(_('subroutine'), 'func', 'subr'),
        'module': ObjType(_('module'), 'mod'),
    }

    directives = {
        'program': FortranProgram,
        'type': FortranType,
        'variable': FortranObject,
        'function': FortranWithSig,
        'subroutine': FortranWithSig,
        'module': FortranModule,
        'currentmodule': FortranCurrentModule,
    }

    roles = {
        'prog': FortranXRefRole(),
        'type': FortranXRefRole(),
        'var': FortranXRefRole(),
        'func': FortranXRefRole(fix_parens=True),
        'subr': FortranXRefRole(fix_parens=True),
        'mod': FortranXRefRole(),
    }
    initial_data = {
        'objects': OrderedDict(),  # fullname -> docname, objtype
        'modules': OrderedDict(),  # modname -> docname, synopsis, platform, deprecated
    }
    indices = [
        FortranModuleIndex,
    ]

    def clear_doc(self, docname):
        for fullname, (fn, _) in list(self.data['objects'].items()):
            if fn == docname:
                del self.data['objects'][fullname]
        for modname, (fn, _, _, _) in list(self.data['modules'].items()):
            if fn == docname:
                del self.data['modules'][modname]

    def find_obj(self, env, modname, name, role, searchorder=0):
        """
        Find a Fortran object for "name", perhaps using the given module and/or
        typename.

        :Params:

            - **searchorder**, optional: Start using relative search
        """
        # skip parens
        if name.endswith('()'):
            name = name[:-2]

        if not name:
            return None, None
        if f_sep in name:
            modname, name = name.split(f_sep)
        #modname = modname or '_'
        if '%' in name:
            name, tmp = name.split('%')

        objects = self.data['objects']
        newname = None
        matches = []
        objtypes = self.objtypes_for_role(role)
        if searchorder == 1:  # :role:`/toto`
            if role in ['mod', 'prog']:
                if 'f' + f_sep + name not in objects:  # exact match
                    return []
                newname = 'f' + f_sep + name
            elif modname and 'f' + f_sep + modname + f_sep + name in objects and \
                    objects['f' + f_sep + modname + f_sep + name][1] in objtypes:
                newname = 'f' + f_sep + modname + f_sep + name
            elif 'f' + f_sep + '_' + f_sep + name in objects and \
                    objects['f' + f_sep + '_' + f_sep + name][1] in objtypes:
                newname = 'f' + f_sep + '_' + f_sep + name
            elif 'f' + f_sep + name in objects and \
                    objects['f' + f_sep + name][1] in objtypes:
                newname = 'f' + f_sep + name
            elif name in objects and \
                    objects[name][1] in objtypes:
                newname = name

        else:  # :role:`toto`
            # NOTE: searching for exact match, object type is not considered
            if 'f' + f_sep + name in objects:
                newname = 'f' + f_sep + name
            elif role in ['mod', 'prog']:
                # only exact matches allowed for modules
                return []
            elif 'f' + f_sep + '_' + f_sep + name in objects:
                newname = 'f' + f_sep + '_' + f_sep + name
            elif modname and 'f' + f_sep + modname + f_sep + name in objects:
                newname = 'f' + f_sep + modname + f_sep + name

        # Last chance: fuzzy search
        if newname is None:
            matches = [(oname, objects[oname]) for oname in objects
                       if oname.endswith(f_sep + name)
                       and objects[oname][1] in objtypes]
        else:
            matches.append((newname, objects[newname]))
        return matches

    def resolve_xref(self, env, fromdocname, builder,
                     type, target, node, contnode):
        modname = node.get('f:module', node.get('modname'))
        typename = node.get('f:type', node.get('typename'))
        searchorder = node.hasattr('refspecific') and 1 or 0
        matches = self.find_obj(env, modname, target, type, searchorder)
        if not matches:
            return None
        elif len(matches) > 1:
            env.warn(fromdocname,
                     'more than one target found for cross-reference '
                     '%r: %s' % (target,
                                 ', '.join(match[0] for match in matches)),
                     node.line)
        name, obj = matches[0]

        if obj[1] == 'module':
            # get additional info for modules
            docname, synopsis, platform, deprecated = self.data['modules'][name[1 + len(
                f_sep):]]
            assert docname == obj[0]
            title = name
            if synopsis:
                title += ': ' + synopsis
            if deprecated:
                title += _(' (deprecated)')
            # return make_refnode(builder, fromdocname, docname,
                # 'module-' + name, contnode, title)
            return make_refnode(builder, fromdocname, docname,
                                name, contnode, title)
        else:
            return make_refnode(builder, fromdocname, obj[0], name,
                                contnode, name)

    def get_objects(self):
        for modname, info in six.iteritems(self.data['modules']):
            yield (modname, modname, 'module', info[0], 'module-' + modname, 0)
        for refname, (docname, type) in six.iteritems(self.data['objects']):
            yield (refname, refname, type, docname, refname, 1)


def setup(app):
    app.add_domain(FortranDomain)
