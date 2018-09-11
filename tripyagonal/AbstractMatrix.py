from sympy import *
from itertools import cycle
from sympy.parsing.sympy_parser import parse_expr


class AbstractMatrix:
    _arg_dict = {'diag': 0, 'super': 1, 'sub': -1}

    def __init__(self, **kwargs):
        self._size = kwargs['n']
        self._data = zeros(self._size, self._size)
        self._kwargs = {}
        self._perturbed = False

        for key in kwargs:
            if key is 'perturbed':
                self._perturbed = kwargs[key]
                continue
            if key is 'n':
                continue
            if isinstance(kwargs[key], list):
                values = kwargs[key]
            else:
                values = [kwargs[key]]

            self._kwargs.update({key: values})
            ele_cycle = cycle(values)

            for i in range(self._size):
                element = next(ele_cycle)
                if isinstance(element, str):
                    element = parse_expr(element)
                if 0 <= i + self._arg_dict[key] < self._size:
                    self._data[i, i + self._arg_dict[key]] = element
                elif key == 'sub':
                    self._alpha = -element
                elif key == 'super':
                    self._beta = -element

        if not self._perturbed:
            self._alpha = self._beta = 0
        self._perturb()

    def _perturb(self):
        pass

    def size(self):
        return self._size

    def alpha(self):
        return self._alpha

    def beta(self):
        return self._beta

    def kwarg(self, key):
        return self._kwargs[key]

    def pretty(self):
        return pretty(self._data)

    def __pow__(self, n):
        return self._data ** n

    def __mul__(self, b):
        return self._data * b

    def __getitem__(self, key):
        return self._data[key]

    def __str__(self):
        return str(self._data)

    def __repr__(self):
        return repr(self._data)

    def latex(self):
        return latex(self._data)

    def evalf(self):
        return self._data.evalf()
