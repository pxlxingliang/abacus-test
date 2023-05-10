import setuptools

with open("version") as f1: lines = f1.readlines()
version = lines[0].split()[1]
git_commit = lines[1].split()[1]

setuptools.setup(
    name = "abacustest",
    version = version+"." + git_commit,
    description = "abacus test system",
    packages=["abacustest"],
    package_dir={"abacustest":"src"},
    package_data = {"abacustest":["*.py","*/*.py","*/*/*.py"]},
    entry_points = {'console_scripts': ['abacustest = abacustest.main:main'],},
    python_requires='>=3.6',
    install_requires=["pydflow>=1.6.84","numpy","pandas"]
)
