The Pencil Code
---------------
The Pencil Code is a high-order finite-difference code for compressible hydrodynamic flows with magnetic fields and particles. It is highly modular and can easily be adapted to different types of problems. The code runs efficiently under MPI on massively parallel shared- or distributed-memory computers.

The Pencil Code is available from
[pencil-code.org](https://pencil-code.org/) and is mirrored to
[GitHub](https://github.com/pencil-code/pencil-code).
It was previously hosted at
[Google Code](https://code.google.com/p/pencil-code/).

In order to checkout the code with read-write premissions
[Subversion](https://subversion.apache.org), use the command
```sh
svn checkout https://pencil-code.org/svn/trunk pencil-code --username=<your-username>
```
where `<your-username>` is your GitHub username that you should use identically to register on [account.pencil-code.org](https://account.pencil-code.org/) for write access to the code repository.

For read-only access via SVN, a username is not required:
```sh
svn checkout https://pencil-code.org/svn/trunk pencil-code
```

To get started, run one of the samples:
```sh
unix>  cd pencil-code
unix>  source sourceme.sh
unix>  cd samples/conv-slab
unix>  mkdir data
```
To set up the symbolic links and compile the code:
```sh
unix>  pc_setupsrc
unix>  pc_build  [ -f /path/to/config/file.conf ]
```
To create the initial condition and run the code:
```sh
unix>  pc_start  [ -f /path/to/config/file.conf ]
unix>  pc_run    [ -f /path/to/config/file.conf ]
```

See `pencil-code/config/hosts/*/*.conf` for sample config files. For more
details, see the manual in the `doc/` directory (also available
[here](http://pencil-code.nordita.org/)).

-----------------------------------------------------------------------------

If you are using bash and you do not want to "source sourceme.sh" on each
session, you can insert the following into your .bashrc and/or .bash_profile:
```sh
export PENCIL_HOME=$HOME/pencil-code  [or wherever you have the code]
_sourceme_quiet=1; . $PENCIL_HOME/sourceme.sh; unset _sourceme_quiet
```

## Documentation

* All documentation is linked from our [documentation overview](https://pencil-code.org/doc.php).
* The [manual][manual] is the main source of information around the code.
* There is also a [quick start][quick_start] to help getting started.
* An auto-generated code documentation is available at [ReadTheDocs](https://pencil-code.readthedocs.io/en/latest/index.html).
* Information about [Python with the Pencil Code][PythonForPencil] and the
  [Python Coding Style][PythonCodingStyle] can be found on the [wiki][wiki].
* Updates to the community are provided through the [newsletter][newsletter].
* The [Pencil Code Office Hours](https://pencil-code.org/contact.php) is a regular online meeting.
* The [Pencil Code User Meeting](http://pencil-code.nordita.org/meetings.php) will be held every year.
* See the [Scientific Usage of the Pencil Code][citations] for papers using or discussing the code.

## List of Contributors

* Around 100 people have contributed to various extent during the
  nearly 20 years of Pencil Code history.
* The current [list of contributors][contributors] shows the temporal
  check-in activity of the those who stayed connected with the code
  over the various host changes (Nordita 2001-2007, Google Code 2007-2015,
  and Github since 2015).
  Some additional contributors are also listed in the [manual][manual].

## How to contribute to the Pencil Code

* For all changes to the code, make sure the auto-test still runs

* If you have write access: check in your changes and make sure you can
  fix possible problems emerging on [travis-ci.com][travis] as well as the
  minutely, hourly, and daily [auto-tests][auto-tests].

* If you have only read access: fork this repository and use pull requests to contribute.

## Code of Conduct

* The Pencil Code community adheres to the [Contributor Covenant Code of Conduct][conduct].
  Please familiarize yourself with its details.

## License
test

* The Pencil Code is under the [GNU public license agreement][license].

[travis]: https://www.travis-ci.com/github/pencil-code/pencil-code
[auto-tests]: http://pencil-code.nordita.org/tests.php
[conduct]: https://github.com/pencil-code/pencil-code/blob/master/license/CODE_OF_CONDUCT.md
[manual]: https://github.com/pencil-code/website/raw/master/doc/manual.pdf
[quick_start]: https://github.com/pencil-code/website/raw/master/doc/quick_start.pdf
[license]: https://github.com/pencil-code/pencil-code/blob/master/license/GNU_public_license.txt
[contributors]: https://github.com/pencil-code/pencil-code/graphs/contributors
[wiki]: https://github.com/pencil-code/pencil-code/wiki
[PythonCodingStyle]: https://github.com/pencil-code/pencil-code/wiki/PythonCodingStyle
[PythonForPencil]: https://github.com/pencil-code/pencil-code/wiki/PythonForPencil
[newsletter]: https://github.com/pencil-code/website/blob/master/NewsLetters/
[citations]: https://github.com/pencil-code/website/raw/master/doc/citations.pdf
[PCSC]: http://norlx65.nordita.org/~brandenb/pencil-code/PCSC/
[meetings]: http://pencil-code.nordita.org/meetings.php

