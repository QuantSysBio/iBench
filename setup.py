
from setuptools import setup

setup(
    name='ibench',
    version=0.1,
    description='Benchmarking for mass spectrometry identifications.',
    author='John Cormican',
    author_email='john.cormican@mpinat.mpg.de',
    packages=[
        'ibench',
    ],
    long_description=open('README.md').read(),
    py_modules=[
        'ibench',
        'ibench.input',
    ],
    entry_points={
        'console_scripts': [
            'ibench=ibench.run:main'
        ]
    },
    install_requires=[
        'biopython==1.79',
        'certifi==2021.10.8',
        'numpy==1.22.3',
        'pandas==1.4.2',
        'plotly==5.8.0',
        'pyopenms==2.7.0',
        'pyteomics==4.5.3',
        'python-dateutil==2.8.2',
        'pytz==2022.1',
        'PyYAML==6.0',
        'six==1.16.0',
        'tenacity==8.0.1',
    ],
)
