# Type declarations for proboscis functions.

from typing import Any, Callable, TypeVar

import proboscis_dummy

F = TypeVar('F', bound=Callable[..., Any])

def bare_decorator(func: F) -> F:
    ...

def decorator_args(url: str) -> Callable[[F], F]:
    ...

from typing import Callable


def test(function: F, **kwargs: Any) -> F: ...


# Really silly: the only way to prevent
#   error: Incompatible import of "TestProgram"
#     (imported name has type "Type[proboscis_dummy.TestProgram]",
#     local name has type "Type[proboscis.TestProgram]"
# (caused by our falling back on proboscis_dummy if proboscis is not
# found) is to pretend that proboscis.TestProgram is a
# proboscis_dummy.TestProgram .
TestProgram = proboscis_dummy.TestProgram
