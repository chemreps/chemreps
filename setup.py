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
        version=open("VERSION").read().strip(),
        packages=['chemreps', 'chemreps.utils'],
        package_data={'chemreps': ['data/*.pkl']},
        license='MIT License',
        url='https://github.com/chemreps/chemreps',
        description='Molecular machine learning representations',
        long_description=open('README.md').read(),
        long_description_content_type='text/markdown',
        install_requires=[
            "cclib>=1.5",
            "numpy>=1.12",
        ],
        classifiers=[
            'Development Status :: 4 - Beta',
            'Intended Audience :: Science/Research',
            'Programming Language :: Python',
            'License :: OSI Approved :: MIT License',
            'Topic :: Scientific/Engineering :: Chemistry',
        ],
    )


if __name__ == '__main__':
    setup_chemreps()
