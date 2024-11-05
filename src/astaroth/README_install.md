both codes together with a fresh pull:

  - git clone -b gputestv6 --recurse-submodules https://<username>@github.com/pencil-code/pencil-code.git
  - git clone -b gputestv6 --recurse-submodules https://<username>@pencil-code.org/git/ pencil-code
    source sourceme.sh
    cd $PENCIL_HOME
    cd src/astaroth/submodule
    git checkout PCinterface_2019-8-12

or add Astaroth to an existing PC installation:

  cd $PENCIL_HOME
  git checkout gputestv6
  git submodule update --init --recursive
  cd src/astaroth/submodule
  git checkout PCinterface_2019-8-12
