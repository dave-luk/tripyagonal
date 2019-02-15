from tripyagonal import *

# The parameters to the matrix are the following:
# n : the size of of the matrix,
# diag : entries on the diagonal of the matrix, provided in order,
# sub : entries on the sub-diagonal of the matrix in order,
# super: entries on the super-diagonal of the matrix, in order,
# perturbed : a boolean value indicating whether to apply perturbation on the matrix.

# How to construct a floating point matrix:

m = TridiagonalMatrix(n=5, diag=1, sub=1, sup=2)

# Note that the entries are provided as floats

# How to construct an exact (stochastic) matrix with multiple entries per diagonal:

m2 = TridiagonalMatrix(n=4, diag=['1/6'], sub=['1/2', '1/3'], sup=['1/3', '1/2'], perturbed=True)

# The entries here are provided as list of strings, and will be parsed by the sympy parser.
# The entries will be recycled if insufficient entries are supplied, which can be taken advantage of to create
# alternating patterns.

# To compute eigenvalues and eigenvectors, select one of the available algorithms:
# Kouachi
# Cayley-Hamilton
# Matrices that satisfy the conditions imposed by these algorithms can be processed as follows:

k_m = Kouachi(m)
c_m = CayleyHamilton(m2)

# After applying these algorithms, we can request outputs with respective available functions, refer to each individual
# algorithms for more info:

k_m.eigenvalues()
c_m.computeMatrixPower(10)

# There are few ways to print the matrices to the console:

# Pretty print
pprint(m2)

# Regular printing (with repr())
print(m2)

# Outputted as latex string
print(latex(m2))

# Refer to utility.py for more built-in tools

# EXAMPLE:
mat = TridiagonalMatrix(n=6, sub=['3/4', '3/8'], diag=['0'], sup=['1/4', '1/2'], perturbed=True)
c_m = CayleyHamilton(mat)
print(c_m.eigenvalues())
print(c_m.computeMatrixPower(13))
