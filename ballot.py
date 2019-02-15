from tripyagonal import *
from sympy import *

def strips(i, j, n):
	high = i + (j - i + n) / 2
	low = i - (i - j + n) / 2
	if low > 0:
		low = 0
	return {'high': high, 'low': low}

def calculate_solution(i, j, steps):
	size = solution_strip(i, j, steps)
	print_Cs(i, j, size)
	size = possible_strip(i, j, steps)
	print_Cs(i, j, size)

def print_Cs(i, j, size):
	mat = TridiagonalMatrix(n=size, sub=['a'], diag=['b'], sup=['c'])
	print(size, i, j)
	#pprint(mat)
	c_m = CayleyHamilton(mat)
	k_m = Kouachi(mat)
	ev = k_m.eigenvalues()
	for m in range(size):
		a, b, c = symbols('a b c')
		cv = simplify(trigsimp(c_m.C_list()[m][i,j]))

		for k in range(size):
			cv.subs(simplify(ev[k]), symbols('l'+str(k)))
			cv.subs((simplify(ev[k]).subs(a, 1-b-c)), symbols('l'+str(k)))
			cv.subs((simplify(ev[k]).subs(b, 1-a-c)), symbols('l'+str(k)))
			cv.subs((simplify(ev[k]).subs(c, 1-a-b)), symbols('l'+str(k)))

		cv = simplify(trigsimp(simplify(cv.subs(a, 1-b-c))))

		for k in range(size):
			cv.subs(simplify(ev[k]), symbols('l'+str(k)))
			cv.subs((simplify(ev[k]).subs(a, 1-b-c)), symbols('l'+str(k)))
			cv.subs((simplify(ev[k]).subs(b, 1-a-c)), symbols('l'+str(k)))
			cv.subs((simplify(ev[k]).subs(c, 1-a-b)), symbols('l'+str(k)))

		cv = simplify(trigsimp(simplify(cv.subs(b, 1-a-c))))

		for k in range(size):
			cv.subs(simplify(ev[k]), symbols('l'+str(k)))
			cv.subs((simplify(ev[k]).subs(a, 1-b-c)), symbols('l'+str(k)))
			cv.subs((simplify(ev[k]).subs(b, 1-a-c)), symbols('l'+str(k)))
			cv.subs((simplify(ev[k]).subs(c, 1-a-b)), symbols('l'+str(k)))

		cv = simplify(trigsimp(simplify(cv.subs(c, 1-a-b))))

		for k in range(size):
			cv.subs(simplify(ev[k]), symbols('l'+str(k)))
			cv.subs((simplify(ev[k]).subs(a, 1-b-c)), symbols('l'+str(k)))
			cv.subs((simplify(ev[k]).subs(b, 1-a-c)), symbols('l'+str(k)))
			cv.subs((simplify(ev[k]).subs(c, 1-a-b)), symbols('l'+str(k)))
		
		cv = simplify(trigsimp(simplify(cv.subs(1-a-b, c))))

		print("\nc" + str(m) + " = " , cv)

for n in range(2, 5):
	i = 1
	j = 2
	steps = (2*n+1) * (j-i)
	result = strips(i, j, steps)
	solution_size = result['high'] 
	possible_size = result['high'] + abs(result['low'])
	offset = abs(result['low']) 
	
	print("\n\nFrom " + str(i) + " to " + str(j) +  " in " + str(steps) + " steps:")
	#print("\nSolution Strip")
	#print_Cs(int(i), int(j), int(solution_size + 1))
	print("\nPossible Strip")
	print_Cs(int(i + offset), int(j + offset), int(possible_size + 1))
