#!/usr/bin/env python3
"""
ViroClust setup script for package installation.

Usage:
    python setup.py install
    python setup.py develop
    pip install -e .
"""

from setuptools import setup, find_packages

with open('README.md', 'r', encoding='utf-8') as f:
    long_description = f.read()

with open('requirements.txt', 'r', encoding='utf-8') as f:
    requirements = [line.strip() for line in f if line.strip() and not line.startswith('#')]

setup(
    name='ViroClust',
    version='1.0.0',
    description='Auto-align sequences and generate consensus based on cluster files',
    long_description=long_description,
    long_description_content_type='text/markdown',
    author='ViroClust Team',
    author_email='viroclust@example.com',
    url='https://github.com/user/ViroClust',
    license='MIT',
    
    packages=find_packages(exclude=['test*', 'output*']),
    
    package_data={
        'src': ['*.py'],
        'bin': ['*.py', '*.json'],
    },
    
    entry_points={
        'console_scripts': [
            'viroclust=bin.viroclust:main',
        ],
    },
    
    python_requires='>=3.8',
    
    install_requires=requirements,
    
    extras_require={
        'dev': [
            'pytest>=6.0',
            'pytest-cov>=2.0',
            'black>=21.0',
            'flake8>=3.8',
        ],
    },
    
    classifiers=[
        'Development Status :: 4 - Beta',
        'Intended Audience :: Science/Research',
        'Topic :: Scientific/Engineering :: Bio-Informatics',
        'License :: OSI Approved :: MIT License',
        'Programming Language :: Python :: 3',
        'Programming Language :: Python :: 3.8',
        'Programming Language :: Python :: 3.9',
        'Programming Language :: Python :: 3.10',
        'Programming Language :: Python :: 3.11',
    ],
    
    keywords='bioinformatics sequence alignment consensus virology',
)
</content>
<parameter=filePath>
