import setuptools

setuptools.setup(
    name = "abacustest",
    description = "abacus test system",
    packages=["abacustest"],
    package_dir={"abacustest":"abacustest"},
    package_data = {"abacustest":["*.py","*/*.py","*/*/*.py"]},
    entry_points = {'console_scripts': ['abacustest = abacustest.main:main'],},
    python_requires='>=3.6',
    install_requires=["pydflow>=1.7.85","numpy","pandas","pyecharts"]
)
