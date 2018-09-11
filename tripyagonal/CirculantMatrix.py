from .AbstractMatrix import AbstractMatrix
from sympy import *


class CirculantMatrix(AbstractMatrix):

    def __init__(self, **kwargs):
        super(CirculantMatrix, self).__init__(perturbed=True, **kwargs)

    def _perturb(self):
        self._data[0, self._size - 1] -= self.alpha()
        self._data[self._size - 1, 0] -= self.beta()

    def eigenvalues(self):
        w_s = []
        for j in range(self._size):
            w_s.append(exp(2 * pi * I * j / self._size))

        l_s = []

        for i in range(self._size):
            lambd = 0

            for j in range(self._size):
                lambd += self._data[0, j] * w_s[i] ** j

            l_s.append(collect(lambd.expand(complex=True), I))

        return l_s
