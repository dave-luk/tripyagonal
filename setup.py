from setuptools import setup

setup(
	name='tripyagonal',
	version='0.1a',
	description='Package with algorithms to solve problems about tridiagonal matrices.',
	author='Dave Luk',
	packages=['tripyagonal'],
	install_requires=[
		'sympy',
		'numpy',
	],
	zip_safe=False)