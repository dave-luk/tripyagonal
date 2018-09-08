from tripyagonal import *

# How to construct a floating point matrix
#m = TridiagonalMatrix(n=5,diag=1, sub=1, super=2)
# m.pprint()

# How to construct an exact (stochastic) matrix with multiple entries per diagonal
#m2 = TridiagonalMatrix(n=4,diag=['1/6'], sub=['1/2','1/3'], super=['1/3','1/2'], perturbed=True)
#pprint(m2)
#m2.pprint()
# With just one entry, not perturbed
m3 = TridiagonalMatrix(n=4 ,diag='1/6', sub=['1/2','7/12'], super=['1/3','1/4'])
c_m3 = CayleyHamilton(m3)
#pprint(c_m3.computeMatrixPower(3))
print(c_m3.computeMatrixPower(10))
#print(c_m3.stateProb(0,1,2))
#print(c_m3.stateProb(0,1,2).evalf())
#print(utility.latex(k_m2.eigenvalues()))
#print(utility.latex(k_m2.eigenvectors()))

#c_m2 = CayleyHamilton(m2)
# print(c_m2.stateProb(1,2,2))
#pprint(c_m2.computeMatrixPower(2))
#pprint(c_m2.computeMatrixPower(5))

#cm = CirculantMatrix(n=5, diag='1/6', sub='1/3', super='1/2')
#cm.pprint()

# k_cm = Kouachi(cm)

# c_cm = CayleyHamilton(cm)
# pprint(c_cm.computePower(6))
