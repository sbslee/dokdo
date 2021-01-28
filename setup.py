from setuptools import setup

exec(open("dokdo/version.py").read())

setup(
    name="dokdo",
    version=__version__,
    author='Seung-been "Steven" Lee',
    author_email="sbstevenlee@gmail.com",
    description=("A Python package for microbiome "
                 "sequencing analysis with QIIME 2"),
    url="https://github.com/sbslee/dokdo",
    entry_points={"console_scripts": ["dokdo=dokdo.__main__:main"]}
)
