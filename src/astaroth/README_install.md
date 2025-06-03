both codes together with a fresh pull:

  - git clone --recurse-submodules https://<username>@github.com/pencil-code/pencil-code.git
  - git clone --recurse-submodules https://<username>@pencil-code.org/git/ pencil-code
    source sourceme.sh
    cd $PENCIL_HOME
    cd src/astaroth/submodule
    git checkout develop

or add Astaroth to an existing PC installation:

  cd $PENCIL_HOME
  git submodule update --init --recursive
  cd src/astaroth/submodule
  git checkout develop
