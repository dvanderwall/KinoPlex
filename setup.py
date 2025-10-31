#!/usr/bin/env python3
"""
Setup script for Structural Context Pipeline
"""
from setuptools import setup, find_packages
import os

# Read the README file for long description
def read_file(filename):
    here = os.path.abspath(os.path.dirname(__file__))
    filepath = os.path.join(here, filename)
    if os.path.exists(filepath):
        with open(filepath, 'r', encoding='utf-8') as f:
            return f.read()
    return ''

setup(
    name='structural-context-pipeline',
    version='1.0.0',
    description='Pipeline for extracting structural context features from AlphaFold protein structures',
    long_description=read_file('README.md'),
    long_description_content_type='text/markdown',
    author='David',
    author_email='your.email@example.com',  # Update with your email
    url='https://github.com/yourusername/structural-context-pipeline',  # Update with your repo URL
    
    packages=find_packages(where='main'),
    package_dir={'': 'main'},
    
    # Include non-Python files
    include_package_data=True,
    
    # Python version requirement
    python_requires='>=3.8',
    
    # Dependencies
    install_requires=[
        'numpy>=1.20.0',
        'pandas>=1.3.0',
        'biopython>=1.79',
        'scipy>=1.7.0',
    ],
    
    # Optional dependencies for development
    extras_require={
        'dev': [
            'pytest>=6.0',
            'pytest-cov>=2.0',
            'black>=21.0',
            'flake8>=3.9',
        ],
    },
    
    # Console scripts - makes the main pipeline directly callable
    entry_points={
        'console_scripts': [
            'structural-context-pipeline=run_structural_context_pipeline:main',
        ],
    },
    
    # Project classification
    classifiers=[
        'Development Status :: 4 - Beta',
        'Intended Audience :: Science/Research',
        'Topic :: Scientific/Engineering :: Bio-Informatics',
        'License :: OSI Approved :: MIT License',  # Update if different
        'Programming Language :: Python :: 3',
        'Programming Language :: Python :: 3.8',
        'Programming Language :: Python :: 3.9',
        'Programming Language :: Python :: 3.10',
        'Programming Language :: Python :: 3.11',
    ],
    
    # Keywords for discoverability
    keywords='structural-biology bioinformatics alphafold protein-structure phosphorylation',
    
    # Project URLs
    project_urls={
        'Bug Reports': 'https://github.com/yourusername/structural-context-pipeline/issues',
        'Source': 'https://github.com/yourusername/structural-context-pipeline',
    },
)