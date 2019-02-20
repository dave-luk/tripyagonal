from tripyagonal import *
from sympy import *
from symengine import var

def print_Cs(i, j, size):
	mat = TridiagonalMatrix(n=size, sub=['a'], diag=['b'], sup=['c'])
	print(size, i, j)
	c_m = CayleyHamilton(mat)
	k_m = Kouachi(mat)
	ev = k_m.eigenvalues()
	for m in range(size):
		a, b, c = symbols('a b c')
		cv = c_m.C_list()[m][i,j]
		with open('./C' + str(m), 'w') as output:
			output.write(str(cv))

print_Cs(1, 2, 4)

