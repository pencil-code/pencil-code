# Makefile for ING Sphinx documentation
export SHELL:=/bin/bash
export SHELLOPTS:=$(if $(SHELLOPTS),$(SHELLOPTS):)pipefail:errexit

# You can set these variables from the command line, and also
# from the environment for the first two.
SPHINXOPTS    ?=
SPHINXBUILD   ?= sphinx-build
SOURCEDIR     = .
BUILDDIR      = _build

# Put it first so that "make" without argument is like "make help".
help:
	@$(SPHINXBUILD) -M help "$(SOURCEDIR)" "$(BUILDDIR)" $(SPHINXOPTS) $(O)

.PHONY: help Makefile

.ONESHELL:

# make "clean": remove build directory and also the stubs
clean: Makefile
	@rm -rf code/sourcePython/*
	@rm -rf code/sourceFortran/*
	@$(SPHINXBUILD) -M $@ "$(SOURCEDIR)" "$(BUILDDIR)" $(SPHINXOPTS) $(O)

# Make "fast" html: do not include autogenerated code docs
fast: Makefile
	@function fast_setup {
		@echo "Commenting out autogenerated docs"
		@sed -i 's/autosummary:/.. .. autosummary:/' toc/modpython.rst
	}
	@function fast_cleanup {
		@echo "Cleaning up"
		@sed -i 's/\.\. \.\. //' toc/modpython.rst
	}
	@trap fast_cleanup EXIT
	@fast_setup
	@export FASTBUILD=TRUE
	@make clean
	$(SPHINXBUILD) -M html "$(SOURCEDIR)" "$(BUILDDIR)" $(SPHINXOPTS) $(O)

html: Makefile
	@$(SPHINXBUILD) -M $@ "$(SOURCEDIR)" "$(BUILDDIR)" $(SPHINXOPTS) $(O)

# Catch-all target: route all unknown targets to Sphinx using the new
# "make mode" option.  $(O) is meant as a shortcut for $(SPHINXOPTS).
%: Makefile
	@$(SPHINXBUILD) -M $@ "$(SOURCEDIR)" "$(BUILDDIR)" $(SPHINXOPTS) $(O)
