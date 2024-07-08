from setuptools import setup, Extension
import os

os.system(f"cd diversutils ; make shared_libraries PYTHON_BUILD=1 ; make update_bashrc")

setup(
    name="diversutils",
    author="Louis Estève",
    author_email="louis.esteve@universite-paris-saclay.fr",
    description="DiversUtils - Functions to measure diversity",
    version="0.1.4",
    packages=["diversutils"],
    package_dir={"diversutils": "diversutils/diversutils"},
    ext_modules = [
        Extension(
            name="_diversutils",
            sources=["diversutils/src/_diversutilsmodule.c"],
            include_dirs=["diversutils/src/include"],
            define_macros=[("ENABLE_AVX256", "0"), ("ENABLE_AVX512", "0")],
            library_dirs=[f"{os.environ['HOME']}/.local/lib/diversutils"],
            libraries=["m", "rt", "diversutils", "udpipe"],
            extra_compile_args=["-g3", "-Wall", "-Wextra", "-pedantic", "-Werror", "-std=c99", "-pthread", "-fstack-protector-all"]
        )
    ]
)

