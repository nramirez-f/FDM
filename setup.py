from setuptools import setup, find_packages

setup(
    name='FDM',
    version='0.0',
    packages=find_packages('src'),
    package_dir={'': 'src'},
    install_requires=[
    'numpy',
    'netCDF4',
    'plotly'
    ],
)