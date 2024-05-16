from setuptools import setup, Extension

setup(
    name="diversutils",
    author="Louis Estève",
    author_email="louis.esteve@universite-paris-saclay.fr",
    description="DiversUtils - Functions to measure diversity",
    version="0.1.0",
    ext_modules = [
        Extension(
            "diversutils",
            sources=["diversutils/src/diversutilsmodule.c"],
            extra_compile_args=["-g3", "-Wall", "-Wextra", "-pedantic", "-Idiversutils/src/include", "-march=native", "-DENABLE_AVX256=0", "-std=c99", "-lm", "-lrt", "-pthread"]
        )
    ]
)

