from setuptools import setup, find_packages

setup(
    name='rankchem',
    version='0.0.1',
    packages=find_packages(where='src'),
    package_dir={'': 'src'},
    install_requires=[
        'rdkit',
        'morfeus-ml',
        'pyvistaqt',
        'numpy',
        'py3Dmol',
        'streamlit',
        'stmol'
    ],
    python_requires='>=3.8',
    author='Emma Kappeler, Ludovica Fracassi',
    author_email='emma.kappeler@epfl.ch, ludovica.fracassi@epfl.ch',
    description='RankChem project',
    long_description_content_type='text/markdown',
    url='https://github.com/kappelemma/RankChem',
    classifiers=[
        'Development Status :: 3 - Alpha',
        'Intended Audience :: Developers',
        'License :: OSI Approved :: MIT License',
        'Operating System :: OS Independent',
        'Programming Language :: Python',
        'Programming Language :: Python :: 3',
        'Programming Language :: Python :: 3.8',
        'Programming Language :: Python :: 3.9',
        'Programming Language :: Python :: 3.10',
        'Programming Language :: Python :: 3.11',
    ],
    license='MIT',
)
