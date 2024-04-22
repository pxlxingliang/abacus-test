import setuptools
import os

if os.path.isfile("abacustest/version"):
    with open("abacustest/version") as f:
        version = f.read().strip()
else:
    version = "0.0.0"


setuptools.setup(
    name = "abacustest",
    version= version,
    description = "abacus test system",
    packages=["abacustest"],
    package_dir={"abacustest":"abacustest"},
    package_data = {"abacustest":["*.py","*/*.py","*/*/*.py","version"]},
    entry_points = {'console_scripts': ['abacustest = abacustest.main:main'],},
    python_requires='>=3.6',
    install_requires=["pydflow>=1.8.45","numpy","pandas","pyecharts","seekpath"]
)
