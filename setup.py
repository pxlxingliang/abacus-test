import setuptools
import os

def read_version():
    version_file = os.path.join(os.path.dirname(__file__), "abacustest/version")
    if os.path.isfile(version_file):
        with open(version_file) as f:
            return f.read().strip()
    return "0.0.0"


setuptools.setup(
    name = "abacustest",
    version= read_version(),
    description = "abacus test system",
    packages=["abacustest"],
    package_dir={"abacustest":"abacustest"},
    package_data = {"abacustest":["*.py","*/*.py","*/*/*.py","version"]},
    entry_points = {'console_scripts': ['abacustest = abacustest.main:main'],},
    python_requires='>=3.6',
    install_requires=["pydflow>=1.8.45","numpy","pandas","pyecharts","seekpath","phonopy", "dpdata","ase"]
)
