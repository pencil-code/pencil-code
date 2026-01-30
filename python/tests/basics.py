#!/usr/bin/python3
# -*- coding: utf-8 -*-   vim: set fileencoding=utf-8 :

"""Test functionality documented under python/tutorials/basics/.

"""

import os
import sys
from typing import Any

from test_utils import assert_equal, _pretty_print, get_rundir, standalone_test

from pencil.sim.simulation import Simulation


def test_neat_short_tricks() -> None:
    """Code from python/tutorials/basics/neat_short-tricks.py .

    Adapted to use existing data directories

    """
    import pencil as pc

    # Obtain sim object
    try:
        import dill  # noqa
    except ModuleNotFoundError as e:
        print(
            "\nTo use get_sim(), you need to have 'dill' installed"
            "\ne.g."
            "\n  sudo apt install python3-dill"
            "\nor"
            "\n  pip3 install dill",
            file=sys.stderr,
        )
        raise Exception(e)
    sim = pc.get_sim(get_rundir("samples/2d-tests/2d_methane_flame/turbulent_field"))

    # Access a value from start.in / run.in
    _assert_sim_parameter(sim, "inituu", "gaussian-noise")
    assert_equal(sim.get_value("ldensity_nolog"), True)
    assert_equal(sim.get_value("iforce"), "helical")
    assert_equal(sim.get_value("nu"), 1.0)

    assert_equal(sim.get_varlist(), [])
    assert_equal(sim.get_varlist(particle=True), [])
    assert_equal(sim.get_varlist(pos="last10"), [])

def test_neat_short_tricks_standalone() -> None:
    """
    Same as test_neat_short_tricks, but in a separate Python process
    """
    standalone_test([
        f"""sim = pc.get_sim('{get_rundir("samples/2d-tests/2d_methane_flame/turbulent_field")}')""",
        'assert sim.get_value("inituu") == "gaussian-noise"',
        'assert sim.get_value("ldensity_nolog") == True',
        'assert sim.get_value("iforce") == "helical"',
        'assert sim.get_value("nu") == 1.0',
        'assert sim.get_varlist() == []',
        'assert sim.get_varlist(particle=True) == []',
        'assert sim.get_varlist(pos="last10") == []',
        ])

def _assert_sim_parameter(
    sim: Simulation, parameter: str, expected: Any
) -> None:
    """Compare a value from start.in/run.in with a reference value."""
    value = sim.get_value(parameter)
    assert_equal(
        value,
        expected,
        "sim.get_value({}) = {} ≠ {}".format(
            parameter, _pretty_print(value), _pretty_print(expected)
        ),
    )
