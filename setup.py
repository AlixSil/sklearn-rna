from setuptools import find_packages, setup

setup(
    name='sklearn-rna',
    packages=find_packages(include=['sklearnRNA']),
    version='0.1.0',
    description='Library offering tool for the treatment of RNA-seq data in the context of machine learning',
    author='Alix Silvert',
    license='MIT',
    install_requires=["pandas", "numpy", "scikit-learn"],
    tests_require=['pytest'],
    test_suite='tests',
)