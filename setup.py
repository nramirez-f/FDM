from setuptools import setup, find_packages

setup(
    name='fdm',
    version='0.0',
    packages=find_packages(where="src"),
    package_dir={"": "src"},
    install_requires=[
    "ncfiles @ git+https://github.com/nramirez-f/NcFiles.git#egg=ncfiles",
    'plotly',
    'scipy',
    ],
    author="Nramirez",
    description="Finite Differences repository",
    long_description=open("README.md").read(),
    long_description_content_type="text/markdown",
    url="https://github.com/nramirez-f/Finite-Differences",
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
)