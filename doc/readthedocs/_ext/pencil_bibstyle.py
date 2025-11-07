# my_bibstyle.py
"""
Custom pybtex formatting style for sphinxcontrib-bibtex.
"""
from pybtex.style.formatting.unsrt import Style as UnsrtStyle
from pybtex.style.formatting.unsrt import pages, date
from pybtex.style.labels.alpha import LabelStyle as AlphaLabelStyle
from pybtex.style.sorting.author_year_title import SortingStyle as AuthorYearTitleSortingStyle
from pybtex.database import Person
from pybtex.style.formatting import BaseStyle, toplevel
from pybtex.style.template import (
    join, words, field, optional, first_of,
    names, sentence, tag, optional_field, href
)
from pybtex.richtext import Text, Symbol
import latexcodec
import re
from dataclasses import dataclass, field as dfield
import sphinxcontrib.bibtex.plugin
from sphinxcontrib.bibtex.style.referencing import BracketStyle, PersonStyle
from sphinxcontrib.bibtex.style.referencing.author_year import AuthorYearReferenceStyle
from sphinxcontrib.bibtex import bibfile
from sphinxcontrib.bibtex.bibfile import BibData, BibFile, get_mtime
from pybtex.database.input.bibtex import Parser
old_parse_bibdata = bibfile.parse_bibdata
from sphinx.util.logging import getLogger
logger = getLogger(__name__)
from pybtex.database import BibliographyDataError

# --- Helper: format a Person as "Last"; we won't use initials in this project --------------
def person_to_last(person: Person) -> str:
    """
    Convert a pybtex.database.Person to "Last", e.g. "John Quincy Adams" -> "Adams"
    """
    # join last names (handles compound last names)
    last = " ".join(person.last_names) if person.last_names else ""
    last = last.encode().decode("latex").replace("{", "").replace("}", "")
    return last

# --- Helper: format author list as "A., B., & C." with Oxford styling -------
def format_authors(persons, together=False):
    """
    Given list of pybtex Person objects, produce:
      - Single author: "Last"
      - Two authors: "Last1 & Last2"
      - Three+ authors: "Last1, Last2 & Last3"
    If together=True, do not add spaces (for use in sorting, year duplicate detection, etc.)
    """
    names = [person_to_last(p) for p in persons]
    if not names:
        return ""
    if len(names) == 1:
        return names[0]
    if len(names) == 2:
        return f"{names[0]} & {names[1]}" if not together else f"{names[0]}{names[1]}"
    # 3 or more: comma separated with " & " before last
    return ", ".join(names[:-1]) + " & " + names[-1] if not together else f"{names[0]}-etal"

def handle_years(entries):
    """
    Go through all the entries, and detected entries where the author-year combination is the same.
    Then, for those entries, modify the year field to add 'a', 'b', etc. as needed.
    We need to also extract the last name of the authors, and build a key such as
    Author1LastName,Author2LastName,...-Year to detect duplicates.

    A typical entry contains keys such as:
       persons={'author': [Person('{Freytag}, B.'), Person('{Steffen}, M.'), Person('{Ludwig}, H.-G.'), Person('{Wedemeyer-B{\\"o}hm}, S.'), Person('{Schaffenberger}, W.'), Person('{Steiner}, O.')]})), ('2013oss..prop...80H', Entry('misc',
       fields=[
         ('title', '{High Temperature Mineral Formation by Short Circuits in Protoplanetary Disks}'),
         ('howpublished', 'NASA Proposal \\#13-OSS13-80'),
         ('year', '2013'),
         ('adsurl', 'http://adsabs.harvard.edu/abs/2013oss..prop...80H'),
         ('adsnote', 'Provided by the SAO/NASA Astrophysics Data System'),
         ('pcsection', 'Planet formation')],
    """
    author_year_dict = {}
    for key, entry in entries.items():
        if 'author' not in entry.persons or 'year' not in entry.fields:
            continue
        authors_key = format_authors(entry.persons['author'], together=True)
        year = entry.fields['year']
        author_year_key = f"{authors_key}-{year}"
        if author_year_key not in author_year_dict:
            author_year_dict[author_year_key] = []
        author_year_dict[author_year_key].append(entry)

    # Now go through the dict and assign suffixes as needed
    for entries_list in author_year_dict.values():
        if len(entries_list) > 1:
            # Sort entries by key to ensure consistent ordering
            entries_list.sort(key=lambda e: e.key)
            for i, entry in enumerate(entries_list):
                suffix = chr(ord('a') + i)
                entry.fields['year'] += suffix

    return entries

def parse_bibdata(bibfilenames, encoding):
    """Parse *bibfilenames* with given *encoding*, and return parsed data."""
    logger.info("CUSTOM PARSER")
    parser = Parser(encoding)
    bibfiles = {}
    keys = {}
    for filename in bibfilenames:
        logger.info("parsing bibtex file {0}... ".format(filename), nonl=True)
        if not filename.is_file():
            logger.warning(
                "could not open bibtex file {0}.".format(filename),
                type="bibtex",
                subtype="bibfile_error",
            )
            new_keys = {}
        else:
            try:
                parser.parse_file(filename)
            except BibliographyDataError as exc:
                logger.warning(
                    "bibliography data error in {0}: {1}".format(filename, exc),
                    type="bibtex",
                    subtype="bibfile_data_error",
                )
            keys, old_keys = dict.fromkeys(parser.data.entries.keys()), keys
            parser.data.entries = handle_years(parser.data.entries)
            assert all(key in keys for key in old_keys)
            new_keys = dict.fromkeys(key for key in keys if key not in old_keys)
            logger.info("parsed {0} entries".format(len(new_keys)))
        bibfiles[filename] = BibFile(mtime=get_mtime(filename), keys=new_keys)
    logger.info("ENDING CUSTOM PARSER")
    return BibData(encoding=encoding, bibfiles=bibfiles, data=parser.data)
bibfile.parse_bibdata = parse_bibdata


def bracket_style() -> BracketStyle:
    return BracketStyle(
        left='(',
        right=')',
    )

def person_style() -> PersonStyle:
    return PersonStyle(
        sep2=' & ',
        last_sep=' & '
    )

def format_journal(j):
    if j == "\\nat":
        return "Nature"
    if j == "\\nar":
        return "New Astron. Rev."
    if j == "\\apss":
        return "Astrophys. Space Sci."
    if j == "\\araa":
        return "Ann. Rev. Astron. Astrophys."
    if j == "\\prd":
        return "Phys. Rev. D"
    if j == "\\pre":
        return "Phys. Rev. E"
    if j == "\\prl":
        return "Phys. Rev. Lett."
    if j == "\\aj":
        return "Astron. J."
    if j == "\\apj":
        return "Astrophys. J."
    if j == "\\apjl":
        return "Astrophys. J. Lett."
    if j == "\\apjs":
        return "Astrophys. J. Supp."
    if j == "\\mnras":
        return "Month. Not. Roy. Astron. Soc."
    if j == "\\physrep":
        return "Phys. Rep."
    if j == "\\aap":
        return "Astron. Astrophys."
    if j == "\\jgr":
        return "J. Geophys. Res."
    if j == "\\grl":
        return "Geophys. Res. Lett."
    if j == "\\solphys":
        return "Sol. Phys."
    if j == "\\ssr":
        return "Space Sci. Ref."
    if j == "\\memsai":
        return "Mem. Soc. Astr. Ital."
    if j == "\\physscr":
        return "Phys. Scr."
    if j == "\\pasj":
        return "Pub. Astron. Soc. Japan"
    if j == "\\jcap":
        return "J. Cosmol. Astropart. Phys."
    if j == "\\psj":
        return "Planet. Sci. J."

    return j.replace("\\&", "&")

@dataclass
class MyReferenceStyle(AuthorYearReferenceStyle):
    bracket_parenthetical: BracketStyle = dfield(default_factory=bracket_style)
    bracket_textual: BracketStyle = dfield(default_factory=bracket_style)
    bracket_author: BracketStyle = dfield(default_factory=bracket_style)
    bracket_label: BracketStyle = dfield(default_factory=bracket_style)
    bracket_year: BracketStyle = dfield(default_factory=bracket_style)
    person: PersonStyle = dfield(default_factory=person_style)

class MyPencilAlpha(AlphaLabelStyle):
    def format_label(self, entry):
        return entry.key
    
class MyPencilSorting(AuthorYearTitleSortingStyle):
    def persons_key(self, persons):
        return format_authors(persons, together=True)

# --- Main formatting style ------------------------------------------------
class MyPencilStyle(UnsrtStyle):
    name = 'pencilstyle'
    default_label_style = 'pencilalpha'
    default_sorting_style = 'pencilsorting'

    def format_address_organization_publisher(
        self, e, include_organization=True):
        """Format address, organization, publisher.
        Everything is optional.
        """
        # small difference from unsrt.bst here: unsrt.bst
        # starts a new sentence only if the address is missing;
        # for simplicity here we always start a new sentence
        if include_organization:
            organization = optional_field('organization')
        else:
            organization = None
        return first_of[
            # this will be rendered if there is an address
            optional [
                join(sep=' ') [
                    sentence[
                        field('address')
                    ],
                    sentence[
                        organization,
                        optional_field('publisher'),
                    ],
                ],
            ],
            # if there is no address then we have this
            sentence[
                organization,
                optional_field('publisher')
            ],
        ]

    def get_inbook_template(self, e):
        template = toplevel [
            format_authors(e.persons['author']),
            f"({e.fields['year']}).",
            sentence [
                self.format_btitle(e, 'title', as_sentence=False),
                self.format_chapter_and_pages(e),
            ],
            self.format_volume_and_series(e),
            sentence [
                field('publisher'),
                optional_field('address'),
                optional [
                    words [field('edition'), 'edition']
                ],
                optional_field('note'),
            ],
            self.format_web_refs(e),
        ]
        return template

    def get_phdthesis_template(self, e):
        template = toplevel [
            format_authors(e.persons['author']),
            f"({e.fields['year']}).",
            optional[ tag('em') [field('title')] ],
            sentence[
                first_of [
                    optional_field('type'),
                    'PhD thesis',
                ],
                field('school'),
                optional_field('address'),
            ],
            sentence [ optional_field('note') ],
            self.format_web_refs(e),
        ]
        return template

    def get_misc_template(self, e):
        template = toplevel [
            format_authors(e.persons['author']),
            f"({e.fields['year']}).",
            optional[ tag('em') [field('title')] ],
            sentence[
                optional[ field('howpublished') ],
            ],
            sentence [ optional_field('note') ],
            self.format_web_refs(e),
        ]
        return template

    def get_inproceedings_template(self, e):
        template = toplevel [
            format_authors(e.persons['author']),
            f"({e.fields['year']}).",
            sentence [
                tag('em') [field('title')]
            ],
            words [
                'In',
                sentence [
                    optional[ self.format_editor(e, as_sentence=False) ],
                    self.format_btitle(e, 'booktitle', as_sentence=False),
                    self.format_volume_and_series(e, as_sentence=False),
                    optional[ pages ],
                ],
                self.format_address_organization_publisher(e),
            ],
            sentence [ optional_field('note') ],
            self.format_web_refs(e),
        ]
        return template

    def get_article_template(self, e):
        volume_and_pages = first_of [
            # volume and pages, with optional issue number
            optional [
                join [
                    tag('b') [field('volume')],
                    optional['(', field('number'),')'],
                    ', ', pages
                ],
            ],
            # pages only
            e.fields['pages'] if "pages" in e.fields and not e.fields['pages'].startswith("arXiv") else "",
        ]
        template = toplevel [
            format_authors(e.persons['author']),
            f"({e.fields['year']}).",
            sentence [
                tag('em') [field('title')]
            ],
            sentence [
                format_journal(e.fields['journal']),
                optional[ volume_and_pages ],
                date],
            sentence [ optional_field('note') ],
            self.format_web_refs(e),
        ]
        return template

    def format_web_refs(self, e):
        # based on urlbst output.web.refs
        return sentence(capfirst=False) [
            optional [ self.format_doi(e) ],
            optional [ self.format_eprint(e) ],
            optional [ self.format_pubmed(e) ],
            optional [ self.format_ads(e) ],
            optional [ self.format_url(e) ]
            ]
    
    def format_ads(self, e):
        if "adsurl" not in e.fields:
            return ""
        # based on urlbst format.doi
        return href [
            field('adsurl'),
            join [
                'ads:',
                e.fields["adsurl"].split('/')[-1]
                ]
            ]

# Register MyPencilStyle as a pybtex formatting plugin
from pybtex.plugin import register_plugin
register_plugin('pybtex.style.labels', 'pencilalpha', MyPencilAlpha)
register_plugin('pybtex.style.formatting', 'pencilstyle', MyPencilStyle)
register_plugin('pybtex.style.sorting', 'pencilsorting', MyPencilSorting)
sphinxcontrib.bibtex.plugin.register_plugin('sphinxcontrib.bibtex.style.referencing', 'author_year_round', MyReferenceStyle)