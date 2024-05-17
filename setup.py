from setuptools import setup, Extension
import os

"""
location = os.sep.join([os.environ["HOME"].rstrip(os.sep), ".local", "diversutils"])
if not (os.path.exists(location) and os.path.isdir(location)):
    os.mkdir(location)
"""

setup(
    name="diversutils",
    author="Louis Estève",
    author_email="louis.esteve@universite-paris-saclay.fr",
    description="DiversUtils - Functions to measure diversity",
    version="0.1.1",
    #location=location,
    packages=["diversutils"],
    package_dir={"diversutils": "diversutils"},
    ext_modules = [
        Extension(
            name="_diversutils",
            sources=["diversutils/src/_diversutilsmodule.c"],
            include_dirs=["diversutils/src/include"],
            define_macros=[("ENABLE_AVX256", "0")],
            library_dirs=[],
            libraries=["m", "rt"],
            extra_compile_args=["-g3", "-Wall", "-Wextra", "-pedantic", "-std=c99", "-pthread"]
        )
    ]
)

