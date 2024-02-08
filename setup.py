from setuptools import setup, find_packages

VERSION = "0.6.9"
DESCRIPTION = "A python package to quickly get FASTA and pipmaker annotation files from Ensembl or GenBank for use in mVISTA alignment analysis."

with open ("README.md", "r") as f:
    LONG_DESCRIPTION = f.read()

setup(
    name="getvista",
    version=VERSION,
    description=DESCRIPTION,
    long_description=LONG_DESCRIPTION,
    long_description_content_type="text/markdown",
    author="Jake Leyhr",
    author_email="jakeleyhr535@gmail.com",  
    url="",
    packages=find_packages(where="src"),
    package_dir={"": "src"},
    package_data={'getvista': ['config.ini']},
    install_requires=[
        "requests",
        "packaging",
        "biopython",
        "configparser",
    ],
    extras_require={
        "dev": ["pytest", "black", "setuptools", "twine"],
    },
    python_requires=">=3.11",
    classifiers=[
        "Programming Language :: Python :: 3",
        "Operating System :: OS Independent",
    ],
    py_modules=['envistacoords', 'envistagene', 'gbvistacoords', 'gbvistagene', 'gbgenerecord', 'emailaddress'],  # Specify your modules here)
    entry_points={
        'console_scripts': [
            'encoords = getvista.envistacoords:main',
            'engene = getvista.envistagene:main',
            'gbcoords = getvista.gbvistacoords:main',
            'gbgene = getvista.gbvistagene:main',
            'gbrecord = getvista.gbgenerecord:main',
            'gbemail = getvista.emailaddress:main',
        ]
    },
)