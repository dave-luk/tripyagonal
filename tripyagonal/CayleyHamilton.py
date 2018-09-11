import itertools
# multi threading
from multiprocessing.dummy import Pool as ThreadPool

import numpy as np
from sympy import *
from functools import reduce
from tripyagonal import CirculantMatrix
from tripyagonal import Kouachi
from tripyagonal import TridiagonalMatrix


class CayleyHamilton:

    def __init__(self, mat):

        self.__n = mat.size()
        self.__C_computed = np.zeros((self.__n, self.__n))
        self.__C_list = []
        for i in range(0, self.__n):
            self.__C_list.append(zeros(self.__n, self.__n))

        if isinstance(mat, TridiagonalMatrix):
            k = Kouachi(mat)
            self.__lambda_ks = k.eigenvalues()
        elif isinstance(mat, CirculantMatrix):
            self.__lambda_ks = mat.eigenvalues()
        else:
            raise TypeError("Incompatible Object")

        self.__init_P_list(mat)
        self.__init_V_inv()

    def eigenvalues(self):
        return self.__lambda_ks

    def __init_P_list(self, mat):
        self.__P_list = [mat]
        for i in range(2, self.__n + 1):
            self.__P_list.append(mat ** i)

    def __init_C_list(self):
        # somewhat empirical...
        if np.sum(self.__C_computed) <= (self.__n ** 2) / 4:
            for a in range(0, self.__n):
                C = zeros(self.__n, self.__n)
                for b in range(0, self.__n):
                    C = C + self.__P_list[b] * self.__vandermonde_inv[a, b]
                self.__C_list[a] = C
        else:
            pool = ThreadPool(4)
            for a in range(0, self.__n):
                for b in range(0, self.__n):
                    if self.__C_computed[a, b] == 0:
                        pool.starmap(self.__compute_C_ij, zip(np.where(self.__C_computed == 0)))

    def __init_V_inv(self):
        self.__vandermonde_inv = Matrix.zeros(self.__n, self.__n)

        for i in range(0, self.__n):
            for j in range(0, self.__n):
                # denominator starts with x_i
                denom = self.__lambda_ks[i]

                # producted with difference of (x_i - x_m)
                for m in range(0, self.__n):
                    if m != i:
                        denom *= (self.__lambda_ks[i] - self.__lambda_ks[m])
                # if we are on the last row, just put 1/denominator
                if j == self.__n - 1:
                    self.__vandermonde_inv[i, j] = 1 / denom
                # otherwise, need to compute numerator
                else:
                    # initialize sum to zero.
                    numer = 0
                    # copy the eigenvalue vector
                    copy = self.__lambda_ks[:]
                    # remove the i-th eigenvalue
                    copy.pop(i)
                    # run a combination on vector with n-j-1 length, aka [n-1 C n-j-1]
                    for subset in itertools.combinations(copy, self.__n - j - 1):
                        # product them and add to the numerator
                        numer += prod(subset, 1)
                    # set the entry to be the numerator divided by denominator, with signs dictated by col
                    self.__vandermonde_inv[i, j] = (-1) ** (self.__n - (j + 1)) * (numer / denom)

        self.__vandermonde_inv = simplify(self.__vandermonde_inv)

    def __compute_C_ij(self, i, j):
        for a in range(0, self.__n):
            C = 0
            for b in range(0, self.__n):
                C = C + self.__P_list[b][i, j] * self.__vandermonde_inv[a, b]
            self.__C_list[a][i, j] = C
        self.__C_computed[i, j] = 1

    def stateProb(self, i, j, k):
        if self.__C_computed[i, j] == 0:
            self.__compute_C_ij(i, j)
        return self.__computePower(i, j, k)

    def C_list(self):
        if not np.sum(self.__C_computed) == self.__n ** 2:
            self.__init_C_list()
        return self.__C_list

    def __computePower(self, i, j, p):
        p_comp = 0
        for counter in range(0, self.__n):
            p_comp += self.__C_list[counter][i, j] * (self.__lambda_ks[i] ** p)
        return simplify(p_comp)

    def computeMatrixPower(self, p):
        if not np.sum(self.__C_computed) == self.__n ** 2:
            self.__init_C_list()
        P_comp = zeros(self.__n, self.__n)
        for i in range(0, self.__n):
            P_comp += self.__C_list[i] * (self.__lambda_ks[i] ** p)
        return simplify(P_comp)
