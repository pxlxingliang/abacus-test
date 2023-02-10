import setuptools
setuptools.setup(
    name = "abacustest",
    description = "abacus test system",
    packages=["abacustest"],
    package_dir={"abacustest":"src"},
    package_data = {"abacustest":["abacustest.py","main.py","collectdata.py","lib_collectdata/*","myflow/*","lib_collectdata/*/*","outresult.py"]},
    entry_points = {'console_scripts': ['abacustest = abacustest.main:main']}
)
