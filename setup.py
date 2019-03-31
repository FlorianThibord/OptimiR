import setuptools

with open("README.md", "r") as fh:
    long_description = fh.read()


    
setuptools.setup(
    name="OptimiR",
    version="1.0",
    author="Florian Thibord",
    author_email="florian.thibord@gmail.com",
    description="integrates genetic information to assess the impact of variants on miRNA expression",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url='https://github.com/FlorianThibord/OptimiR',
    packages=setuptools.find_packages(),
    include_package_data=True,
    classifiers=[
        "Programming Language :: Python",
        "License :: OSI Approved :: GNU General Public License v3 (GPLv3)",
        "Operating System :: OS Independent",
    ],
    install_requires=[
        'pysam',
        'biopython'
    ],
    entry_points={
        'console_scripts': ['optimir = optimir.command_line:main'],
    }
)


