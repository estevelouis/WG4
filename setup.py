from setuptools import setup, Extension
import os

ENABLE_UDPIPE = False
ENABLE_LIB = False

if ENABLE_UDPIPE:
	os.system("cd diversutils ; make -B shared_libraries VERBOSE=1 ; make -B update_bashrc")
else:
	os.system("cd diversutils ; make -B ~/.local/lib/diversutils PYTHON_BUILD=1 VERBOSE=1 ; make -B update_bashrc")

if not ENABLE_LIB:
    os.system("cd diversutils ; make -B build/_diversutilsmodule.c PYTHON_BUILD=1 VERBOSE=1")

setup(
    name="diversutils",
    author="Louis Estève",
    author_email="louis.esteve@universite-paris-saclay.fr",
    description="DiversUtils - Functions to measure diversity",
    version="0.2.0",
    packages=["diversutils"],
    package_dir={"diversutils": "diversutils/diversutils"},
    ext_modules = [
        Extension(
            name="_diversutils",
            sources=[f"diversutils/{'src' if ENABLE_LIB else 'build'}/_diversutilsmodule.c"],
            include_dirs=["diversutils/src/include"],
            define_macros=[("ENABLE_AVX256", "0"), ("ENABLE_AVX512", "0"), ("TOKENIZATION_METHOD", "0")],
            library_dirs=[f"{os.environ['HOME']}/.local/lib/diversutils"] if ENABLE_LIB or ENABLE_UDPIPE else [],
            libraries=["m", "rt"] + (["diversutils"] if ENABLE_LIB else []) + (["udpipe"] if ENABLE_UDPIPE else []),
            extra_compile_args=["-g3", "-Wall", "-Wextra", "-std=c99", "-pthread", "-pedantic", "-Werror", "-fstack-protector-all"]
        )
    ]
)

