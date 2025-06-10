"""
Setup for pypi releases of wormcat3
"""
from setuptools import setup, find_packages
from pathlib import Path

__version__ = "0.1.6"

this_directory = Path(__file__).parent
long_description = (this_directory / "README.md").read_text()

setup(name='wormcat3',
      version=__version__,
      description='WormCat is a tool for annotating and visualizing gene set enrichment data from C. elegans microarray, RNA seq or RNAi screen data.',
      long_description_content_type="text/markdown",
      long_description=long_description,

      url='https://github.com/DanHUMassMed/wormcat3.git',
      author='Dan Higgins',
      author_email='daniel.higgins@gatech.edu',
      license='MIT',

      packages=find_packages(),
      install_requires=[
        'ipykernel',
        'scipy==1.15.2',
        'pandas==2.2.3',
        'gseapy==1.1.8',
        'plotnine==0.14.5',
        'statsmodels==0.14.4',
        'XlsxWriter==3.2.2',
        'openpyxl==3.1.5'
      ],
      include_package_data=True,
      zip_safe=False)
