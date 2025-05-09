from setuptools import setup

setup(
    name="fdm",
    version="0.0",
    py_modules=["__init__", "advection", "burgers"],
    install_requires=["numpy", "scipy", "ncfiles @ git+https://github.com/nramirez-f/NcFiles.git@main#egg=ncfiles"],
    author="Nramirez",
    description="Finite Difference Method module",
    long_description=open("README.md").read(),
    long_description_content_type="text/markdown",
    url="https://github.com/nramirez-f/FDM",
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
)