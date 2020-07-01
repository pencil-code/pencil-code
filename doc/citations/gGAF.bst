 %-----------------------------------------
 % For Geophysical and Astrophysical Fluid Dynamics `gGAF.bst',
 % Taylor & Francis Limited.
 % generated with the docstrip utility.
 % Developed by SR Nova, Bangalore, India
 % ----------------------------------------
 % *** Geophysical and Astrophysical Fluid Dynamics (gGAF) ***
 % 
 % This bibliography style file is intended for texts in ENGLISH
 % This is an author-year citation style bibliography. As such, it is
 % non-standard LaTeX, and requires a special package file to function properly.
 % Such a package is    natbib.sty   by Patrick W. Daly
 % The form of the \bibitem entries are
 %   \bibitem[\protect\citeauthoryear{Jones}{1990}]{key}...
 %   \bibitem[\protect\citeauthoryear{Jones and Baker}{1990}]{key}...
 %   \bibitem[\protect\citeauthoryear{Jones {\itshape{et~al.}}}{1990}]{key}...
 % The essential feature is that the label (the part in brackets) consists
 % of the author names, as they should appear in the citation, with the year
 % in parentheses following. There must be no space before the opening
 % parenthesis!
 % In natbib.sty, it is possible to define the type of enclosures that is
 % really wanted (brackets or parentheses), but in either case, there must
 % be parentheses in the label.
 % The \cite command functions as follows:
 %   \citet{key} ==>>                Jones et al. (1990)
 %   \citep{key} ==>>                (Jones et al. 1990)
 %   \citep[chap. 2]{key} ==>>       (Jones et al. 1990, chap. 2)
 %   \citep[e.g.][]{key} ==>>        (e.g. Jones et al. 1990)
 %   \citep[e.g.][p. 32]{key} ==>>   (e.g. Jones et al. 1990, p. 32)
 %   \citep{key,key} ==>>            (Smith and Jones 1990, Patel et al. 2003)
 %   \citeauthor{key} ==>>           Jones et al.
 %   \citeyear{key} ==>>             1990
 %---------------------------------------------------------------------

ENTRY
  { address
    author
    booktitle
    chapter
    chaptertitle
    englishtitle
    edition
    editor
    eid
    howpublished
    institution
    journal
    key
    month
    note
    number
    organization
    pages
    publisher
    location
    school
    series
    seriesnumber
    title
    type
    volume
    year
  }
  {}
  { label extra.label sort.label short.list} 
INTEGERS { output.state before.all mid.sentence after.sentence after.block }
FUNCTION {init.state.consts}
{ #0 'before.all :=
  #1 'mid.sentence :=
  #2 'after.sentence :=
  #3 'after.block :=
}
STRINGS { s t}
FUNCTION {output.nonnull}
{ 's :=
  output.state mid.sentence =
    { ", " * write$ }
    { output.state after.block =
        { add.period$ write$
          newline$
          "\newblock " write$
        }
        { output.state before.all =
            'write$
            { add.period$ " " * write$ }
          if$
        }
      if$
      mid.sentence 'output.state :=
    }
  if$
  s
}
FUNCTION {output}
{ duplicate$ empty$
    'pop$
    'output.nonnull
  if$
}
FUNCTION {output.check}
{ 't :=
  duplicate$ empty$
    { pop$ "empty " t * " in " * cite$ * warning$ }
    'output.nonnull
  if$
}
FUNCTION {fin.entry}
{ add.period$
  write$
  newline$
}

FUNCTION {new.block}
{ output.state before.all =
    'skip$
    { after.block 'output.state := }
  if$
}
FUNCTION {new.sentence}
{ output.state after.block =
    'skip$
    { output.state before.all =
        'skip$
        { after.sentence 'output.state := }
      if$
    }
  if$
}
FUNCTION {add.blank}
{  " " * before.all 'output.state :=
}

FUNCTION {add.semicolon}
{  "; " * before.all 'output.state :=
}

FUNCTION {date.block}
{
  skip$
}

FUNCTION {not}
{   { #0 }
    { #1 }
  if$
}
FUNCTION {and}
{   'skip$
    { pop$ #0 }
  if$
}
FUNCTION {or}
{   { pop$ #1 }
    'skip$
  if$
}
FUNCTION {non.stop}
{ duplicate$
   "}" * add.period$
   #-1 #1 substring$ "." =
}

FUNCTION {new.block.checkb}
{ empty$
  swap$ empty$
  and
    'skip$
    'new.block
  if$
}
FUNCTION {field.or.null}
{ duplicate$ empty$
    { pop$ "" }
    'skip$
  if$
}
FUNCTION {emphasize}
{ duplicate$ empty$
    { pop$ "" }
    { "{\itshape " swap$ * "}" * }
  if$
}

FUNCTION {bold}
{ duplicate$ empty$
    { pop$ "" }
    { "\textbf{" swap$ * "}" * }
  if$
}

FUNCTION {etalemphasize}
{ duplicate$ empty$
    { pop$ "" }
    { "{\itshape{" swap$ * "}}" * }
  if$
}

FUNCTION {tie.or.space.prefix}
{ duplicate$ text.length$ #3 <
    { "~" }
    { " " }
  if$
  swap$
}

FUNCTION {capitalize}
{ "u" change.case$ "t" change.case$ }

FUNCTION {space.word}
{ " " swap$ * " " * }
 % Here are the language-specific definitions for explicit words.
 % Each function has a name bbl.xxx where xxx is the English word.
 % The language selected here is ENGLISH
FUNCTION {bbl.and}
{ "and"}

FUNCTION {bbl.etal}
{ "et~al." }

FUNCTION {bbl.editors}
{ "Eds" }

FUNCTION {bbl.editor}
{ "ed." }

FUNCTION {bbl.edby}
{ "edited by " }

FUNCTION {bbl.edition}
{ "edn" }

FUNCTION {bbl.volume}
{ "Vol. " }

FUNCTION {bbl.of}
{ "of" }

FUNCTION {bbl.number}
{ "no." }

FUNCTION {bbl.nr}
{ "no." }

FUNCTION {bbl.in}
{ "In" }

FUNCTION {bbl.plainin}
{ "in" }

FUNCTION {bbl.pages}
{ "pp." }

FUNCTION {bbl.page}
{ "p." }

FUNCTION {bbl.chapter}
{ "chap." }

FUNCTION {bbl.chaptertitle}
{ "Chapter title." }

FUNCTION {bbl.techrep}
{ "Technical report" }

FUNCTION {bbl.mthesis}
{ "Master's thesis" }

FUNCTION {bbl.phdthesis}
{ "Ph.D. Thesis" }

MACRO {jan} {"Jan."}

MACRO {feb} {"Feb."}

MACRO {mar} {"Mar."}

MACRO {apr} {"Apr."}

MACRO {may} {"May"}

MACRO {jun} {"Jun."}

MACRO {jul} {"Jul."}

MACRO {aug} {"Aug."}

MACRO {sep} {"Sep."}

MACRO {oct} {"Oct."}

MACRO {nov} {"Nov."}

MACRO {dec} {"Dec."}

MACRO {acmcs} {"ACM Computing Surveys"}

MACRO {acta} {"Acta Informatica"}

MACRO {cacm} {"Communications of the ACM"}

MACRO {ibmjrd} {"IBM Journal of Research and Development"}

MACRO {ibmsj} {"IBM Systems Journal"}

MACRO {ieeese} {"IEEE Transactions on Software Engineering"}

MACRO {ieeetc} {"IEEE Transactions on Computers"}

MACRO {ieeetcad}
 {"IEEE Transactions on Computer-Aided Design of Integrated Circuits"}

MACRO {ipl} {"Information Processing Letters"}

MACRO {jacm} {"Journal of the ACM"}

MACRO {jcss} {"Journal of Computer and System Sciences"}

MACRO {scp} {"Science of Computer Programming"}

MACRO {sicomp} {"SIAM Journal on Computing"}

MACRO {tocs} {"ACM Transactions on Computer Systems"}

MACRO {tods} {"ACM Transactions on Database Systems"}

MACRO {tog} {"ACM Transactions on Graphics"}

MACRO {toms} {"ACM Transactions on Mathematical Software"}

MACRO {toois} {"ACM Transactions on Office Information Systems"}

MACRO {toplas} {"ACM Transactions on Programming Languages and Systems"}

MACRO {tcs} {"Theoretical Computer Science"}

FUNCTION {bibinfo.check}
{ swap$
  duplicate$ missing$
    {
      pop$ pop$
      ""
    }
    { duplicate$ empty$
        {
          swap$ pop$
        }
        { swap$
          pop$
        }
      if$
    }
  if$
}
FUNCTION {bibinfo.warn}
{ swap$
  duplicate$ missing$
    {
      swap$ "missing " swap$ * " in " * cite$ * warning$ pop$
      ""
    }
    { duplicate$ empty$
        {
          swap$ "empty " swap$ * " in " * cite$ * warning$
        }
        { swap$
          pop$
        }
      if$
    }
  if$
}
INTEGERS { nameptr namesleft numnames }


STRINGS  { bibinfo}

FUNCTION {format.names}
{ 'bibinfo :=
  duplicate$ empty$ 'skip$ {
  's :=
  "" 't :=
  #1 'nameptr :=
  s num.names$ 'numnames :=
  numnames 'namesleft :=
    { namesleft #0 > }
    { s nameptr
      "{vv~}{ll}{, jj}{, f{.}.}"
      format.name$
      bibinfo bibinfo.check
      't :=
      nameptr #1 >
        {
          namesleft #1 >
            { ", " * t * }
            {
              s nameptr "{ll}" format.name$ duplicate$ "others" =
                { 't := }
                { pop$ }
              if$
              t "others" =
                {
                  " " * bbl.etal emphasize *
                }
                {
                  bbl.and
                  space.word * t *
                }
              if$
            }
          if$
        }
        't
      if$
      nameptr #1 + 'nameptr :=
      namesleft #1 - 'namesleft :=
    }
  while$
  } if$
}
FUNCTION {format.names.ed}
{
  'bibinfo :=
  duplicate$ empty$ 'skip$ {
  's :=
  "" 't :=
  #1 'nameptr :=
  s num.names$ 'numnames :=
  numnames 'namesleft :=
    { namesleft #0 > }
    { s nameptr
      "{f{.}.~}{vv~}{ll}{ jj}"
      format.name$
      bibinfo bibinfo.check
      't :=
      nameptr #1 >
        {
          namesleft #1 >
            { ", " * t * }
            {
              s nameptr "{ll}" format.name$ duplicate$ "others" =
                { 't := }
                { pop$ }
              if$
              t "others" =
                {

                  " " * bbl.etal emphasize *
                }
                {
                  bbl.and
                  space.word * t *
                }
              if$
            }
          if$
        }
        't
      if$
      nameptr #1 + 'nameptr :=
      namesleft #1 - 'namesleft :=
    }
  while$
  } if$
}
FUNCTION {format.key}
{ empty$
    { key field.or.null }
    { "" }
  if$
}

FUNCTION {format.authors}
{ author "author" format.names
}
FUNCTION {get.bbl.editor}
{ editor num.names$ #1 > 'bbl.editors 'bbl.editor if$ }

FUNCTION {format.editors}
{ editor "editor" format.names duplicate$ empty$ 'skip$
    {
      "" *
      " " *
      get.bbl.editor
      capitalize
   "(" swap$ * ")" *
      *
    }
  if$
}

FUNCTION {format.plaineditors}
{ editor "editor" format.names empty$ 'skip$
    {
      bbl.edby swap$ *
%%      " " *
%%      get.bbl.editor
%%    capitalize
%%   "(" swap$ * ")" *
%%      *
    }
  if$
}

FUNCTION {format.note}
{
 note empty$
    { "" }
    { add.blank
      note #1 #1 substring$
      duplicate$ "{" =
        'skip$
        { output.state mid.sentence =
          { "l" }
          { "u" }
        if$
        change.case$
        }
      if$
      note #2 global.max$ substring$ * "note" bibinfo.check
    }
  if$
}

FUNCTION {format.booknote}
{
 note empty$
    { ")" }
    { add.blank
      note #1 #1 substring$
      duplicate$ "{" =
        'skip$
        { output.state mid.sentence =
          { "l" }
          { "u" }
        if$
        change.case$
        }
      if$
      note #2 global.max$ substring$ * "note" bibinfo.check
    }
  if$
}

FUNCTION {format.title}
{ title
  "title" bibinfo.check
  duplicate$ empty$ 'skip$
    {
      "" swap$ *
      ", " *
    }
  if$
}
FUNCTION {format.trntitle}
{ title
  "title" bibinfo.check
  duplicate$ empty$ 'skip$
    { emphasize
      ", " swap$ *
      ", " *
    }
  if$
}


FUNCTION {format.englishtitle}
{ englishtitle
  "englishtitle" bibinfo.check
  duplicate$ empty$ 'skip$
    { emphasize
      "[" swap$ *
      "], " *
    }
  if$
}
FUNCTION {format.plaintitle}
{ title
  "title" bibinfo.check
  duplicate$ empty$ 'skip$
    {
      "" swap$ *
      ". " *
    }
  if$
}
FUNCTION {format.title.dot}
{ title
  "title" bibinfo.check
  duplicate$ empty$ 'skip$
    {
      "``" swap$ *
      "''. " *
    }
  if$
}
FUNCTION {end.quote.title}
{ title empty$
    'skip$
    { before.all 'output.state := }
  if$
}
FUNCTION {format.full.names}
{'s :=
 "" 't :=
  #1 'nameptr :=
  s num.names$ 'numnames :=
  numnames 'namesleft :=
    { namesleft #0 > }
    { s nameptr
      "{vv~}{ll}" format.name$
      't :=
      nameptr #1 >
        {
          namesleft #1 >
            { ", " * t * }
            {
              s nameptr "{ll}" format.name$ duplicate$ "others" =
                { 't := }
                { pop$ }
              if$
              t "others" =
                {
                  " " * bbl.etal emphasize *
                }
                {
                  bbl.and
                  space.word * t *
                }
              if$
            }
          if$
        }
        't
      if$
      nameptr #1 + 'nameptr :=
      namesleft #1 - 'namesleft :=
    }
  while$
}

FUNCTION {author.editor.key.full}
{ author empty$
    { editor empty$
        { key empty$
            { cite$ #1 #3 substring$ }
            'key
          if$
        }
        { editor format.full.names }
      if$
    }
    { author format.full.names }
  if$
}

FUNCTION {author.key.full}
{ author empty$
    { key empty$
         { cite$ #1 #3 substring$ }
          'key
      if$
    }
    { author format.full.names }
  if$
}

FUNCTION {editor.key.full}
{ editor empty$
    { key empty$
         { cite$ #1 #3 substring$ }
          'key
      if$
    }
    { editor format.full.names }
  if$
}

FUNCTION {make.full.names}
{ type$ "book" =
  type$ "inbook" =
  or
    'author.editor.key.full
    { type$ "proceedings" =
        'editor.key.full
        'author.key.full
      if$
    }
  if$
}

FUNCTION {output.bibitem}
{ newline$
  "\bibitem[\protect\citeauthoryear{" write$
  label write$
  "}]{" write$
  cite$ write$
  "}" write$
  newline$
  ""
  before.all 'output.state :=
}

FUNCTION {n.dashify}
{
  't :=
  ""
    { t empty$ not }
    { t #1 #1 substring$ "-" =
        { t #1 #2 substring$ "--" = not
            { "--" *
              t #2 global.max$ substring$ 't :=
            }
            {   { t #1 #1 substring$ "-" = }
                { "-" *
                  t #2 global.max$ substring$ 't :=
                }
              while$
            }
          if$
        }
        { t #1 #1 substring$ *
          t #2 global.max$ substring$ 't :=
        }
      if$
    }
  while$
}

FUNCTION {word.in}
{ bbl.in
  " " * }

FUNCTION {word.plainin}
{ bbl.plainin
  " " * }

FUNCTION {format.date}
{ year "year" bibinfo.check duplicate$ empty$
    {
      "empty year in " cite$ * "; set to ????" * warning$
       pop$ "????"
    }
    'skip$
  if$
  extra.label *
  before.all 'output.state :=
  " " swap$ * "" *
}

FUNCTION {format.unpubdate}
{ year "year" bibinfo.check duplicate$ empty$
    {
      "empty year in " cite$ * "; set to ????" * warning$
       pop$ "????"
    }
    'skip$
  if$
  extra.label *
  before.all 'output.state :=
  " (" swap$ * ")" *
}

FUNCTION {format.plaindate}
{ year "year" bibinfo.check duplicate$ empty$
    {
      "empty year in " cite$ * "; set to ????" * warning$
       pop$ "????"
    }
    'skip$
  if$
  extra.label *
  before.all 'output.state :=
  ", " swap$ * "" *
}

FUNCTION {format.iptitle}
{ title "title" bibinfo.check duplicate$ empty$
    {
      "empty title in " cite$ * "; set to ????" * warning$
       pop$ "????"
    }
    'skip$
  if$
  extra.label *
  before.all 'output.state :=
  ", " swap$ * "; " *
}

FUNCTION {format.ibpublisher}
{ publisher "publisher" bibinfo.check duplicate$ empty$
    {
      "empty publisher in " cite$ * "; set to ????" * warning$
       pop$ "????"
    }
    'skip$
  if$
  extra.label *
  before.all 'output.state :=
  " (" swap$ * "" *
}

FUNCTION {format.ibaddress}
{ address "address" bibinfo.check duplicate$ empty$
    {
      "empty address in " cite$ * "; set to ????" * warning$
       pop$ "????"
    }
    'skip$
  if$
  extra.label *
  before.all 'output.state :=
  ": " swap$ * ")" *
}


FUNCTION {format.tbpublisher}
{ publisher "publisher" bibinfo.check duplicate$ empty$
    {
      "empty publisher in " cite$ * "; set to ????" * warning$
       pop$ "????"
    }
    'skip$
  if$
  extra.label *
  before.all 'output.state :=
  " (" swap$ * "" *
}

FUNCTION {format.tbaddress}
{ address "address" bibinfo.check duplicate$ empty$
    {
      "empty address in " cite$ * "; set to ????" * warning$
       pop$ "????"
    }
    'skip$
  if$
  extra.label *
  before.all 'output.state :=
  ", " swap$ * ")" *
}


FUNCTION {format.btitle}
{ title empty$
    { "" }
    {bbl.chaptertitle output add.blank}
  if$
  title "title" bibinfo.check
  duplicate$ empty$ 'skip$
    {
      emphasize
    }
  if$
}
FUNCTION {format.booktitle}
{ title empty$
    { "" }
    {title "title" bibinfo.check
    }
  if$
  duplicate$ empty$ 'skip$
    {
      emphasize
    }
  if$
}
FUNCTION {format.chaptertitle}
{ title empty$
    { "" }
    {chaptertitle "chaptertitle" bibinfo.check}
  if$
  duplicate$ empty$ 'skip$
    {
      "" swap$ *
      ". " *
    }
  if$
}

FUNCTION {format.inbtitle}
{ title "title" bibinfo.check
  emphasize
  word.in swap$ *
}

FUNCTION {either.or.check}
{ empty$
    'pop$
    { "can't use both " swap$ * " fields in " * cite$ * warning$ }
  if$
}
FUNCTION {format.bvolume}
{ volume empty$
    { "" }
    { bbl.volume volume tie.or.space.prefix
      "volume" bibinfo.check * *
      series "series" bibinfo.check
      duplicate$ empty$ 'pop$
        { swap$ bbl.of space.word * swap$
          emphasize * }
      if$
      "volume and number" number either.or.check
    }
  if$
}

FUNCTION {format.bookvolume}
{ volume empty$
    { ", " }
    { ", " bbl.volume tie.or.space.prefix 
      volume "volume" bibinfo.check * *  
      duplicate$ empty$ 'pop$
        {* ", " 
          * }
      if$
      "volume and number" number either.or.check
    }
  if$
}
FUNCTION {format.number.series}
{ seriesnumber empty$
        { series field.or.null }
        { series empty$
            { seriesnumber "seriesnumber" bibinfo.check }
            { output.state mid.sentence =
                { }
                { }
              if$
              series "series" bibinfo.check *
              seriesnumber tie.or.space.prefix "seriesnumber" bibinfo.check ", " * *
            }
          if$
        }
      if$
    }


FUNCTION {format.number.plainseries}
{ volume empty$
    { number empty$
        { series field.or.null }
        { series empty$
            { number "number" bibinfo.check }
            { output.state mid.sentence =
                { bbl.number }
                { bbl.number capitalize }
              if$
              number tie.or.space.prefix "number" bibinfo.check * *
              bbl.plainin space.word *
              series "series" bibinfo.check *
            }
          if$
        }
      if$
    }
    { "" }
  if$
}

FUNCTION {format.edition}
{ edition duplicate$ empty$ 'skip$
    {
      output.state mid.sentence =
        { "l" }
        { "t" }
      if$ change.case$
      "edition" bibinfo.check
      add.blank
   " (" swap$ " " * bbl.edition * ")" *
      *
    }
  if$
}

FUNCTION {format.bookedition}
{ edition duplicate$ empty$ 'skip$
    { "edition" bibinfo.check
      " edn, " * 
    }
  if$
}

FUNCTION {format.trnbookedition}
{ edition duplicate$ empty$ 'skip$
    { "edition" bibinfo.check
      "translated by " swap$ * 
    }
  if$
}
INTEGERS { multiresult }
FUNCTION {multi.page.check}
{ 't :=
  #0 'multiresult :=
    { multiresult not
      t empty$ not
      and
    }
    { t #1 #1 substring$
      duplicate$ "-" =
      swap$ duplicate$ "," =
      swap$ "+" =
      or or
        { #1 'multiresult := }
        { t #2 global.max$ substring$ 't := }
      if$
    }
  while$
  multiresult
}
FUNCTION {format.pages}
{ pages duplicate$ empty$ 'skip$
    { duplicate$ multi.page.check
        {
          bbl.pages swap$
          n.dashify
        }
        {
          bbl.page swap$
        }
      if$
      tie.or.space.prefix
      "pages" bibinfo.check
      * *
    }
  if$
}


FUNCTION {format.journal.pages}
{ pages empty$
    'skip$
    { duplicate$ empty$
        { pop$ format.pages }
        {
          ", " *
          pages n.dashify *
        }
      if$
    }
  if$
}


FUNCTION {format.journal.eid}
{ eid "eid" bibinfo.check
  duplicate$ empty$ 'pop$
    { swap$ duplicate$ empty$ 'skip$
      {
          ", " *
      }
      if$
      swap$ *
    }
  if$
}


FUNCTION {format.vol.pages}
{ volume field.or.null
  duplicate$ empty$ 'skip$
    {
      "volume" bibinfo.check
      bold
    }
  if$ 
  eid empty$
    { format.journal.pages }
    { format.journal.eid }
  if$
}

FUNCTION {format.vol.num.pages}
{ volume field.or.null
  duplicate$ empty$ 'skip$
    {
      bold
    }
  if$
  number "number" bibinfo.check duplicate$ empty$ 'skip$
    {
      swap$ duplicate$ empty$
        { "there's a number but no volume in " cite$ * warning$ }
        'skip$
      if$
      swap$
      "(" swap$ * ")" *
    }
  if$ * 
  eid empty$
    { format.journal.pages }
    { format.journal.eid }
  if$
}

FUNCTION {format.chapter.pages}
{ chapter empty$
    'format.pages
    { type empty$
        { bbl.chapter }
        { type "l" change.case$
          "type" bibinfo.check
        }
      if$
      chapter tie.or.space.prefix
      "chapter" bibinfo.check
      * *
      pages empty$
        'skip$
        { ", " * format.pages * }
      if$
    }
  if$
}


FUNCTION {format.chapter}
{ chapter empty$
    'format.pages
    { type empty$
        { bbl.chapter }
        { type "l" change.case$
          "type" bibinfo.check
        }
      if$
      chapter tie.or.space.prefix
      "chapter" bibinfo.check
      * *
    }
  if$
}


FUNCTION {format.inchapter}
{ chapter
  "chapter" bibinfo.check
  duplicate$ empty$ 'skip$
    {
      "" swap$ *
      ". " *
    }
  if$
}

FUNCTION {format.edbooktitle}
{
  booktitle "booktitle" bibinfo.check
  emphasize
}
FUNCTION {format.in.ed.booktitle}
{ format.edbooktitle duplicate$ empty$ 'skip$
    {
      editor "editor" format.names.ed duplicate$ empty$ 'pop$
        {
          get.bbl.editor
          " " * swap$ * 
          "" swap$ * "" *
          swap$
          ", " * swap$
          * 
          }
      if$
      word.in swap$ *
    }
  if$
}


FUNCTION {format.in.booktitle}
{
  booktitle "booktitle" bibinfo.check
  emphasize
      word.in swap$ *
}

FUNCTION {format.inproc.booktitle}
{
  booktitle "booktitle" bibinfo.check
       "" word.plainin swap$ emphasize * swap$ emphasize *
}

FUNCTION {format.plainin.booktitle}
{
  booktitle "booktitle" bibinfo.check
}

FUNCTION {format.thesis.type}
{ type duplicate$ empty$
    'pop$
    { swap$ pop$
      "t" change.case$ "type" bibinfo.check
    }
  if$
}
FUNCTION {format.tr.number}
{ number "number" bibinfo.check
  type duplicate$ empty$
    { pop$ bbl.techrep }
    'skip$
  if$
  "type" bibinfo.check
  swap$ duplicate$ empty$
    { pop$ "t" change.case$ }
    { tie.or.space.prefix * * }
  if$
}
FUNCTION {format.article.crossref}
{
  word.in
  " \cite{" * crossref * "}" *
}
FUNCTION {format.book.crossref}
{ volume duplicate$ empty$
    { "empty volume in " cite$ * "'s crossref of " * crossref * warning$
      pop$ word.in
    }
    { bbl.volume
      swap$ tie.or.space.prefix "volume" bibinfo.check * * bbl.of space.word *
    }
  if$
  " \cite{" * crossref * "}" *
}
FUNCTION {format.incoll.inproc.crossref}
{
  word.in
  " \cite{" * crossref * "}" *
}


FUNCTION {format.org.or.pubwithparenths}
{ 't :=
  address empty$ t empty$ and
    'skip$
    { " ("
      address "address" bibinfo.check *
      t empty$
        'skip$
        { address empty$
            'skip$
            {": " * " " bibinfo.check  "" * }
          if$
          t * 
        }
      if$
    }
  if$  
}

FUNCTION {format.org.or.pub}
{ 't :=
  ""
  address empty$ t empty$ and
    'skip$
    {
      t empty$
        { address "address" bibinfo.check *
        }
        { t *
          address empty$
            'skip$
            { ", " * address "address" bibinfo.check * }
          if$
        }
      if$
    }
  if$
}



FUNCTION {format.address.publisher}
{ publisher "publisher" bibinfo.warn format.org.or.pubwithparenths
}

 
FUNCTION {noformat.org.or.pub}
{ 't :=
  ""
  address empty$ t empty$ and
    'skip$
    {
      address "address" bibinfo.check *
      t empty$
        'skip$
        { address empty$
            'skip$
            { ": " * }
          if$
          t * 
        }
      if$
    }
  if$  
}


FUNCTION {format.organization.address}
{ organization "organization" bibinfo.check noformat.org.or.pub
}

FUNCTION {format.location}
{ location empty$
    { "" }
    { location "location" bibinfo.check}
  if$  
}

FUNCTION {format.month}
{ month empty$
    { "" }
    { month "month" bibinfo.check}
  if$  
}


FUNCTION {article}
{ output.bibitem
  format.authors "author" output.check
  author format.key output
  format.plaintitle "title" output.check
  end.quote.title
  crossref missing$
    {
      journal
      "journal" bibinfo.check
      emphasize
      "journal" output.check
      format.plaindate "year" output.check
      format.vol.pages output
    }
    { format.article.crossref output.nonnull
      format.plaindate "year" output.check
      format.pages output
    }
  if$
  format.note output
  fin.entry
}
FUNCTION {book}
{ output.bibitem
  author empty$
    { format.editors "author and editor" output.check
      editor format.key output
      add.blank
    }
    { format.authors output.nonnull
      crossref missing$
        { "author and editor" editor either.or.check }
        'skip$
      if$
      }
  if$
  format.booktitle "title" output.check
  end.quote.title
  crossref missing$
    { format.bookvolume output
      end.quote.title
      format.number.series output
      end.quote.title
      format.bookedition output
      end.quote.title 
      format.date "year" output.check
      end.quote.title
      format.address.publisher output
      end.quote.title
    }
    {
      format.book.crossref output.nonnull
    }
  if$ 
  format.booknote output
  fin.entry
}


FUNCTION {translatedbook}
{ output.bibitem
  author empty$
    { format.editors "author and editor" output.check
      editor format.key output
      add.blank
    }
    { format.authors output.nonnull
      crossref missing$
        { "author and editor" editor either.or.check }
        'skip$
      if$
      }
  if$
  end.quote.title
  format.trntitle "title" output.check
  end.quote.title
  format.trnbookedition output
  crossref missing$
    { format.bvolume output
      format.number.series output 
      format.tbpublisher "publisher" output.check
      format.tbaddress "address" output.check
    }
    {
      format.book.crossref output.nonnull
    }
  if$
  end.quote.title
  format.plaindate "year" output.check
  format.booknote output
  fin.entry
}
FUNCTION {booklet}
{ output.bibitem
  format.authors output
  author format.key output
  format.title "title" output.check
  end.quote.title
  howpublished "howpublished" bibinfo.check output
  address "address" bibinfo.check output
  format.note output
  format.date "year" output.check
  fin.entry
}


FUNCTION {inbook}
{ output.bibitem
  author empty$
    { format.editors "author and editor" output.check
      editor format.key output
    }
    { format.authors output.nonnull
      crossref missing$
        { "author and editor" editor either.or.check }
        'skip$
      if$
    }
  if$
  crossref missing$
    {
      format.bvolume output
      format.number.series output
      format.inchapter "chapter" output.check
      end.quote.title
    }
    {
      format.inchapter "chapter" output.check
      new.block
      format.book.crossref output.nonnull
    }
  if$
  end.quote.title
  format.inbtitle "title" output.check
  format.pages output
  format.edition output
  format.plaindate "year" output.check
  end.quote.title
  format.ibpublisher "publisher" output.check
  format.ibaddress "address" output.check
  format.note output
  fin.entry
}


FUNCTION {incollection}
{ output.bibitem
  format.authors "author" output.check
  author format.key output
  format.plaintitle "title" output.check
  end.quote.title
  crossref missing$
    { format.in.booktitle "booktitle" output.check
      format.number.series output
      format.edition output
      editor "editor" format.names.ed
      format.plaineditors "author and editor" output.check
      format.bvolume output
    }
    { format.incoll.inproc.crossref output.nonnull
    }
  if$
  format.note output
  format.chapter.pages output
  format.plaindate "year" output.check
  end.quote.title
  format.ibpublisher "publisher" output.check
  format.ibaddress "address" output.check
  format.note output  
  fin.entry
}


FUNCTION {inproceedings}
{ output.bibitem
  format.authors "author" output.check
  author format.key output
  format.iptitle "title" output.check
  end.quote.title
  crossref missing$
    { format.inproc.booktitle "booktitle" output.check
      format.edition output
      editor empty$
        { }
        { editor "editor" format.names.ed 
          format.plaineditors "author and editor" output.check
        }
      if$
      format.bvolume output
      format.number.plainseries output
    }
    { format.incoll.inproc.crossref output.nonnull
    }
  if$
  format.location output
  format.month output
  format.note output
  format.plaindate "year" output.check
  format.chapter.pages output
  fin.entry
}

FUNCTION {conference} { inproceedings }

FUNCTION {manual}
{ output.bibitem
  format.authors output
  author format.key output
  format.btitle "title" output.check
  organization "organization" bibinfo.check output
  address "address" bibinfo.check output
  format.edition output
  format.note output
  format.date "year" output.check
  fin.entry
}

FUNCTION {mastersthesis}
{ output.bibitem
  format.authors "author" output.check
  author format.key output
  format.title
  "title" output.check
  end.quote.title
  bbl.mthesis format.thesis.type output.nonnull
  school "school" bibinfo.warn output
  address "address" bibinfo.check output
  format.note output
  format.date "year" output.check
  fin.entry
}

FUNCTION {misc}
{ output.bibitem
  format.authors output
  author format.key output
  format.title output
  end.quote.title
  howpublished "howpublished" bibinfo.check output
  format.note output
  format.date "year" output.check
  fin.entry
}


FUNCTION {phdthesis}
{ output.bibitem
  format.authors "author" output.check
  author format.key output
  format.plaintitle
  "title" output.check
  end.quote.title
  bbl.phdthesis output.nonnull 
  school "school" bibinfo.warn output
  address "address" bibinfo.check output
  format.note output
  format.plaindate "year" output.check
  fin.entry
}

FUNCTION {proceedings}
{ output.bibitem
  format.editors output
  editor format.key output
  format.btitle "title" output.check
  format.bvolume output
  format.number.plainseries output
  publisher empty$
    { format.organization.address output }
    { organization "organization" bibinfo.check output
      format.address.publisher output
    }
  if$
  format.note output
  format.date "year" output.check
  fin.entry
}

FUNCTION {techreport}
{ output.bibitem
  format.authors "author" output.check
  author format.key output
  format.plaintitle
  "title" output.check
  end.quote.title
  format.date "year" output.check
  format.tr.number output.nonnull
  institution "institution" bibinfo.warn output
  address "address" bibinfo.check output
  format.note output
  fin.entry
}

FUNCTION {unpublished}
{ output.bibitem
  format.authors "author" output.check
  author format.key output
  format.plaintitle "title" output.check
  end.quote.title
  format.note "note" output.check
  format.unpubdate output
  fin.entry
}

FUNCTION {default.type} { misc }
READ
FUNCTION {sortify}
{ purify$
  "l" change.case$
}
INTEGERS { len }
FUNCTION {chop.word}
{ 's :=
  'len :=
  s #1 len substring$ =
    { s len #1 + global.max$ substring$ }
    's
  if$
}
FUNCTION {format.lab.names}
{ 's :=
  "" 't :=
  s #1 "{vv~}{ll}" format.name$
  s num.names$ duplicate$
  #2 >
    { pop$
      " " * bbl.etal etalemphasize * 
    }
    { #2 <
        'skip$
        { s #2 "{ff }{vv }{ll}{ jj}" format.name$ "others" =
            {
              " " * bbl.etal etalemphasize *
            }
            { bbl.and space.word * s #2 "{vv~}{ll}" format.name$
              * }
          if$
        }
      if$
    }
  if$
}

FUNCTION {author.key.label}
{ author empty$
    { key empty$
        { cite$ #1 #3 substring$ }
        'key
      if$
    }
    { author format.lab.names }
  if$
}

FUNCTION {author.editor.key.label}
{ author empty$
    { editor empty$
        { key empty$
            { cite$ #1 #3 substring$ }
            'key
          if$
        }
        { editor format.lab.names }
      if$
    }
    { author format.lab.names }
  if$
}

FUNCTION {editor.key.label}
{ editor empty$
    { key empty$
        { cite$ #1 #3 substring$ }
        'key
      if$
    }
    { editor format.lab.names }
  if$
}

FUNCTION {calc.short.authors}
{ type$ "book" =
  type$ "inbook" =
  or
    'author.editor.key.label
    { type$ "proceedings" =
        'editor.key.label
        'author.key.label
      if$
    }
  if$
   'short.list :=
}

FUNCTION {calc.label}
{ calc.short.authors
  short.list
  "}{"
  *
  year duplicate$ empty$
     { pop$ "????" }
     { purify$ #-1 #4 substring$ }
  if$
  *
  'label :=
}

FUNCTION {sort.format.names}
{ 's :=
  #1 'nameptr :=
  ""
  s num.names$ 'numnames :=
  numnames 'namesleft :=
    { namesleft #0 > }
    { s nameptr
      "{vv{ } }{ll{ }}{  ff{ }}{  jj{ }}"
      format.name$ 't :=
      nameptr #1 >
        {
          "   "  *
          namesleft #1 = t "others" = and
            { "zzzzz" * }
            { t sortify * }
          if$
        }
        { t sortify * }
      if$
      nameptr #1 + 'nameptr :=
      namesleft #1 - 'namesleft :=
    }
  while$
}

FUNCTION {sort.format.title}
{ 't :=
  "A " #2
    "An " #3
      "The " #4 t chop.word
    chop.word
  chop.word
  sortify
  #1 global.max$ substring$
}
FUNCTION {author.sort}
{ author empty$
    { key empty$
        { "to sort, need author or key in " cite$ * warning$
          ""
        }
        { key sortify }
      if$
    }
    { author sort.format.names }
  if$
}
FUNCTION {author.editor.sort}
{ author empty$
    { editor empty$
        { key empty$
            { "to sort, need author, editor, or key in " cite$ * warning$
              ""
            }
            { key sortify }
          if$
        }
        { editor sort.format.names }
      if$
    }
    { author sort.format.names }
  if$
}
FUNCTION {editor.sort}
{ editor empty$
    { key empty$
        { "to sort, need editor or key in " cite$ * warning$
          ""
        }
        { key sortify }
      if$
    }
    { editor sort.format.names }
  if$
}
FUNCTION {presort}
{ calc.label
  label sortify
  "    "
  *
  type$ "book" =
  type$ "inbook" =
  or
    'author.editor.sort
    { type$ "proceedings" =
        'editor.sort
        'author.sort
      if$
    }
  if$
  #1 entry.max$ substring$
  'sort.label :=
  sort.label
  *
  "    "
  *
  title field.or.null
  sort.format.title
  *
  #1 entry.max$ substring$
  'sort.key$ :=
}

ITERATE {presort}
SORT
STRINGS { last.label next.extra }
INTEGERS { last.extra.num number.label }
FUNCTION {initialize.extra.label.stuff}
{ #0 int.to.chr$ 'last.label :=
  "" 'next.extra :=
  #0 'last.extra.num :=
  #0 'number.label :=
}
FUNCTION {forward.pass}
{ last.label label =
    { last.extra.num #1 + 'last.extra.num :=
      last.extra.num int.to.chr$ 'extra.label :=
    }
    { "a" chr.to.int$ 'last.extra.num :=
      "" 'extra.label :=
      label 'last.label :=
    }
  if$
  number.label #1 + 'number.label :=
}
FUNCTION {reverse.pass}
{ next.extra "b" =
    { "a" 'extra.label := }
    'skip$
  if$
  extra.label 'next.extra :=
  extra.label
  duplicate$ empty$
    'skip$
    { "{\natexlab{" swap$ * "}}" * }
  if$
  'extra.label :=
  label extra.label * 'label :=
}
EXECUTE {initialize.extra.label.stuff}
ITERATE {forward.pass}
REVERSE {reverse.pass}
FUNCTION {bib.sort.order}
{ sort.label
  "    "
  *
  year field.or.null sortify
  *
  "    "
  *
  title field.or.null
  sort.format.title
  *
  #1 entry.max$ substring$
  'sort.key$ :=
}
ITERATE {bib.sort.order}
SORT
FUNCTION {begin.bib}
{ preamble$ empty$
    'skip$
    { preamble$ write$ newline$ }
  if$
  "\begin{thebibliography}{" number.label int.to.str$ * "}" *
  write$ newline$
  "\providecommand{\natexlab}[1]{#1}"
  write$ newline$
}
EXECUTE {begin.bib}
EXECUTE {init.state.consts}
ITERATE {call.type$}
FUNCTION {end.bib}
{ newline$
  "\end{thebibliography}" write$ newline$
}
EXECUTE {end.bib}
%%
%% End of file `gGAF.bst'.
