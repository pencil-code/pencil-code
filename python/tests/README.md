$PENCIL_HOME/python/tests
=========================

Test Pencil Code Python modules, so we can feel better when refactoring the Python code.

```sh
$PENCIL_HOME/python/tests/test-python-modules.py
```

These tests are best run with the [_Proboscis_](https://pythonhosted.org/proboscis/) test runner[^1]:
```sh
pip3 install proboscis
```
but will fall back on a minimal mockup implementation of _Proboscis_ if necessary.


[^1]: The main reason for using _Probioscis_ is test discovery by decorator, as opposed to the standard practice of disovering by name.
