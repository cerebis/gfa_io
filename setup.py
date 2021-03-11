from setuptools import setup

setup(
    name='gfa_io',
    version='0.1',
    packages=['gfa_io'],
    url='https://github.com/cerebis/gfa_io',
    license='GNU  GPLv3',
    author='Matthew Z DeMaere',
    author_email='matt.demaere@gmail.com',
    description='A simple package for reading GFA files',

    classifiers=[
        'Programming Language :: Python :: 3.7',
        'License :: OSI Approved :: GNU Affero General Public License v3',
        'Operating System :: POSIX :: Linux',
        'Intended Audience :: Science/Research',
        'Topic :: Scientific/Engineering :: Bio-Informatics',
        'Development Status :: 4 - Beta'
    ],

    install_requires=['biopython',
                      'networkx'],

)
