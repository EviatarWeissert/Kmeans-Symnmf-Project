from setuptools import setup, Extension

module = Extension(
    'symnmfmodule',  # Module name
    sources=['symnmfmodule.c', 'symnmf.c'],  # Source file
   
)

setup(
    name='symnmfmodule',
    version='1.0',
    description='A Python C extension for symnmf algorithm',
    ext_modules=[module]
)
