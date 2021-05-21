from setuptools import setup

setup(
        name='pencil',
        version='2.0',
        description='Python for pencil postprocessing scripts',
        packages=['pencil', 'pencil.backpack','pencil.calc','pencil.diag','pencil.export','pencil.pio', 'pencil.ism_dyn', 'pencil.math', 'pencil.read', 'pencil.sim', 'pencil.tool_kit', 'pencil.util', 'pencil.visu'],
        package_dir={'': 'python'}
)
