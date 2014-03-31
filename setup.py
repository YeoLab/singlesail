__author__ = 'olga'

from setuptools import setup
from setuptools import find_packages


def readme():
    with open('README.md') as f:
        return f.read()


setup(name='singlesail',
      version='0.1.0',
      description='Single-cell hub for object-oriented exploratory research',
      long_description='Functions for common single-cell analyses such as '
                       'PCA, clustering, outlier detection, splicing '
                       'modality clustering, bimodal gene expression '
                       'detection.',
      url='http://github.com/YeoLab/singlesail',
      author='Olga Botvinnik',
      author_email='obotvinn@ucsd.edu',
      license='MIT',
      packages=find_packages(),
      install_requires=map(lambda x: x.rstrip(),
                           open('requirements.txt').readlines()),
      dependency_links=
      ['https://github.com/scikit-fuzzy/scikit-fuzzy/tarball/master#egg'
                        '=scikit-fuzzy'],
      zip_safe=False,
      test_suite='nose.collector',
      tests_require=['nose'],
)
