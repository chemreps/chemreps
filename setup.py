'''
Implementation of various machine learning representations for molecules 
'''
import setuptools
from distutils.core import setup


def setup_chemreps():
    setup(
        name='chemreps',
        maintainer='Dakota Folmsbee, Amanda Dumi, Shiv Upadhyay',
        maintainer_email='dfolmsbee@gmail.com, amandaedumi@gmail.com, shivnupadhyay@gmail.com',
        version='0.0.2',
        packages=['chemreps', 'chemreps.utils'],
        license='MIT License',
        url='https://github.com/dlf57/chemreps',
        description='Molecular machine learning representations',
        long_description=open('README.md').read(),
        install_requires=[
            "cclib>=1.5",
            "numpy>=1.12",
        ],
        classifiers=[
            'Development Status :: 3 - Alpha',
            'Intended Audience :: Science/Research',
            'Programming Language :: Python',
            'License :: OSI Approved :: MIT License',
            'Topic :: Scientific/Engineering :: Chemistry',
        ],
    )


if __name__ == '__main__':
    setup_chemreps()
