from itertools import *

from sympy import *
from sympy.parsing.sympy_parser import parse_expr

from tripyagonal import TridiagonalMatrix


class Kouachi:

    def __init__(self, mat):

        if not isinstance(mat, TridiagonalMatrix):
            raise TypeError('Object provided is not a TridiagonalMatrix.')

        self.__mat = mat

        self.__n = mat.size()
        self.__m = int(self.__n / 2)

        self.__theta_ks = []
        self.__lambda_ks = []
        self.__alpha = mat.alpha()
        self.__beta = mat.beta()

        self.__ev = Matrix.zeros(self.__n, self.__n)

        __sub = list(islice(cycle(mat.kwarg('sub')), self.__n))
        # drop alpha entry
        del __sub[0]
        __super = list(islice(cycle(mat.kwarg('super')), self.__n - 1))
        __diag = list(islice(cycle(mat.kwarg('diag')), self.__n))

        self.__sub = [parse_expr(x) for x in __sub] if isinstance(__sub[0], str) else __sub
        self.__sup = [parse_expr(x) for x in __super] if isinstance(__super[0], str) else __super
        self.__diag = [parse_expr(x) for x in __diag] if isinstance(__diag[0], str) else __diag

        # check if b's are in b1-b2 or b's
        if len(mat.kwarg('diag')) == 2:
            self.__b_1 = self.__diag[0]
            self.__b_2 = self.__diag[1]
        elif len(mat.kwarg('diag')) == 1:
            self.__b_1 = self.__b_2 = self.__diag[0]
        else:
            raise ValueError('Matrix does not satisfy the assumptions of Kouachi <b>')

        # a_ic_i is constant (2 or less)

        # temp a
        vec_a = [self.__alpha] + self.__sub
        vec_c = self.__sup + [self.__beta]

        __dict = []
        for i in range(self.__n - 1):
            p = vec_c[i] * vec_a[i + 1]
            __dict.append(p)
        __dict_order = __dict
        __dict = list(set(__dict))
        self.__d_size = len(__dict)

        if self.__d_size == 1:
            self.__d_sq = __dict[0]
        elif self.__d_size == 2:
            if __dict[0] != __dict_order[0]:
                __dict[0], __dict[1] = __dict[1], __dict[0]
            self.__d_sq = list(__dict)
            self.__d = list(map(lambda element: sqrt(element), __dict))
        else:
            print(__dict)
            raise ValueError('Matrix does not satisfy the assumptions of Kouachi <d_sq>, (more than 2 unique d_sq)')

    # alpha = beta = 0
    def __case1(self):
        if self.__n % 2 == 0:
            for k in range(1, self.__n + 1):
                if k <= self.__m:
                    theta_k = '{}*pi/{}'.format(k, 2 * self.__m + 1)
                else:
                    theta_k = '{}*pi/{}'.format((k - self.__m), 2 * self.__m + 1)
                self.__theta_ks.append(parse_expr(theta_k))
        else:
            for k in range(1, self.__n + 1):
                if k <= self.__m:
                    theta_k = '{}*pi/{}'.format(k, 2 * self.__m + 2)
                else:
                    theta_k = '{}*pi/{}'.format((k - self.__m), 2 * self.__m + 2)
                self.__theta_ks.append(parse_expr(theta_k))

        for k in range(1, self.__n + 1):
            if k <= self.__m:
                lambda_k = ((self.__b_1 + self.__b_2) - sqrt(
                    (self.__b_1 - self.__b_2) ** 2 + 16 * self.__d_sq * (cos(self.__theta_ks[k - 1]) ** 2))) / 2
            elif k <= 2 * self.__m:
                lambda_k = ((self.__b_1 + self.__b_2) + sqrt(
                    (self.__b_1 - self.__b_2) ** 2 + 16 * self.__d_sq * (cos(self.__theta_ks[k - 1]) ** 2))) / 2
            else:
                lambda_k = self.__b_1
            self.__lambda_ks.append(factor(lambda_k))

    def __case2(self):
        for k in range(1, self.__n + 1):
            if k <= self.__m - 1:
                theta_k = '{}*pi/{}'.format(k, 2 * self.__m)
            elif k <= self.__n - 2:
                theta_k = '{}*pi/{}'.format((k - self.__m + 1), 2 * self.__m)
            else:
                theta_k = '0'
            self.__theta_ks.append(parse_expr(theta_k))

        for k in range(1, self.__n + 1):
            if k <= self.__m - 1:
                lambda_k = ((self.__b_1 + self.__b_2) - sqrt(
                    (self.__b_1 - self.__b_2) ** 2 + 16 * self.__d_sq * cos(self.__theta_ks[k - 1]) ** 2)) / 2
            elif k <= self.__n - 2:
                lambda_k = ((self.__b_1 + self.__b_2) + sqrt(
                    (self.__b_1 - self.__b_2) ** 2 + 16 * self.__d_sq * cos(self.__theta_ks[k - 1]) ** 2)) / 2
            elif k == self.__n - 1:
                lambda_k = ((self.__b_1 + self.__b_2 - self.__alpha - self.__beta) - sqrt(
                    (self.__b_1 - self.__b_2) ** 2 + (self.__alpha + self.__beta) ** 2 - 2 * (
                            self.__b_1 - self.__b_2) * (self.__alpha - self.__beta))) / 2
            else:
                lambda_k = ((self.__b_1 + self.__b_2 - self.__alpha - self.__beta) + sqrt(
                    (self.__b_1 - self.__b_2) ** 2 + (self.__alpha + self.__beta) ** 2 - 2 * (
                            self.__b_1 - self.__b_2) * (self.__alpha - self.__beta))) / 2
            self.__lambda_ks.append(factor(expand(lambda_k)))

        return self.__lambda_ks

    def __case3(self):
        if self.__n % 2 == 0:
            for k in range(1, self.__n + 1):
                if k <= self.__m:
                    theta_k = '{}*pi/{}'.format(k, 2 * self.__m + 1)
                else:
                    theta_k = '{}*pi/{}'.format((k - self.__m), 2 * self.__m + 1)
                self.__theta_ks.append(parse_expr(theta_k))
        else:
            for k in range(1, self.__n + 1):
                if k <= self.__m:
                    theta_k = '{}*pi/{}'.format(k, 2 * self.__m + 2)
                else:
                    theta_k = '{}*pi/{}'.format((k - self.__m), 2 * self.__m + 2)
                self.__theta_ks.append(parse_expr(theta_k))

        for k in range(1, self.__n + 1):
            if k <= self.__m:
                lambda_k = self.__b_1 - 2 * sqrt(self.__sub[0] * self.__sup[0]) * cos(self.__theta_ks[k - 1])
            elif k <= 2 * self.__m:
                lambda_k = self.__b_1 + 2 * sqrt(self.__sub[0] * self.__sup[0]) * cos(self.__theta_ks[k - 1])
            else:
                lambda_k = self.__b_1
            self.__lambda_ks.append(simplify(factor(expand(lambda_k))))

        return self.__lambda_ks

    def __case4(self):
        # if n is odd
        if self.__n % 2 == 1:
            # Thm 3.1
            if self.__alpha == self.__beta == 0:
                for k in range(1, self.__n + 1):
                    if k <= self.__m:
                        theta_k = '{}*pi/{}'.format(2 * k, self.__n + 1)
                    elif k <= self.__n - 2:
                        theta_k = '{}*pi/{}'.format(2 * (k - self.__m), self.__n + 1)
                    else:
                        theta_k = '0'
                    self.__theta_ks.append(parse_expr(theta_k))

                for k in range(1, self.__n + 1):
                    if k <= self.__m:
                        lambda_k = self.__b_1 + sqrt(
                            self.__d_sq[0] + self.__d_sq[1] + 2 * self.__d[0] * self.__d[1] * cos(
                                self.__theta_ks[k - 1]))
                    elif k <= 2 * self.__m:
                        lambda_k = self.__b_1 - sqrt(
                            self.__d_sq[0] + self.__d_sq[1] + 2 * self.__d[0] * self.__d[1] * cos(
                                self.__theta_ks[k - 1]))
                    else:
                        lambda_k = self.__b_1
                    self.__lambda_ks.append(factor(lambda_k))
            # Thm 3.2
            elif self.__alpha == self.__d[1] and self.__beta == self.__d[0] \
                    or self.__alpha == -self.__d[1] and self.__beta == -self.__d[0]:
                for k in range(1, self.__n + 1):
                    if k <= self.__m:
                        theta_k = '{}*pi/{}'.format(2 * k, self.__n)
                    elif k <= self.__n - 2:
                        theta_k = '{}*pi/{}'.format(2 * (k - self.__m), self.__n)
                    else:
                        theta_k = '0'
                    self.__theta_ks.append(parse_expr(theta_k))

                for k in range(1, self.__n + 1):
                    if k <= self.__m:
                        lambda_k = self.__b_1 + sqrt(
                            self.__d_sq[0] + self.__d_sq[1] + 2 * self.__d[0] * self.__d[1] * cos(
                                self.__theta_ks[k - 1]))
                    elif k <= 2 * self.__m:
                        lambda_k = self.__b_1 - sqrt(
                            self.__d_sq[0] + self.__d_sq[1] + 2 * self.__d[0] * self.__d[1] * cos(
                                self.__theta_ks[k - 1]))
                    else:
                        lambda_k = self.__b_1 - (self.__alpha + self.__beta)
                    self.__lambda_ks.append(factor(lambda_k))
            # Thm 3.3
            elif self.__alpha == -self.__d[1] and self.__beta == self.__d[0] \
                    or self.__alpha == self.__d[1] and self.__beta == -self.__d[0]:
                for k in range(1, self.__n + 1):
                    if k <= self.__m:
                        theta_k = '{}*pi/{}'.format(2 * k - 1, self.__n)
                    elif k <= self.__n - 2:
                        theta_k = '{}*pi/{}'.format(2 * (k - self.__m) - 1, self.__n)
                    else:
                        theta_k = '0'
                    self.__theta_ks.append(parse_expr(theta_k))

                for k in range(1, self.__n + 1):
                    if k <= self.__m:
                        lambda_k = self.__b_1 + sqrt(
                            self.__d_sq[0] + self.__d_sq[1] + 2 * self.__d[0] * self.__d[1] * cos(
                                self.__theta_ks[k - 1]))
                    elif k <= 2 * self.__m:
                        lambda_k = self.__b_1 - sqrt(
                            self.__d_sq[0] + self.__d_sq[1] + 2 * self.__d[0] * self.__d[1] * cos(
                                self.__theta_ks[k - 1]))
                    else:
                        lambda_k = self.__b_1 - (self.__alpha + self.__beta)
                    self.__lambda_ks.append(factor(lambda_k))
            else:
                raise ValueError('Tridiagonal Matrix cannot be solved using Kouachi\'s methods')
        # if n is even
        else:
            # if alpha*beta equals d_2^2
            if self.__alpha * self.__beta == self.__d_sq[1]:
                for k in range(1, self.__n + 1):
                    if k <= self.__m - 1:
                        theta_k = '{}*pi/{}'.format(2 * k, self.__n)
                    elif k <= self.__n - 2:
                        theta_k = '{}*pi/{}'.format(2 * (k - self.__m + 1), self.__n)
                    else:
                        theta_k = '0'
                    self.__theta_ks.append(parse_expr(theta_k))

                for k in range(1, self.__n + 1):
                    if k <= self.__m - 1:
                        lambda_k = self.__b_1 + sqrt(
                            self.__d_sq[0] + self.__d_sq[1] + 2 * self.__d[0] * self.__d[1] * cos(
                                self.__theta_ks[k - 1]))
                    elif k <= self.__n - 2:
                        lambda_k = self.__b_1 - sqrt(
                            self.__d_sq[0] + self.__d_sq[1] + 2 * self.__d[0] * self.__d[1] * cos(
                                self.__theta_ks[k - 1]))
                    elif k == self.__n - 1:
                        lambda_k = self.__b_1 + (-(self.__alpha + self.__beta) + sqrt(
                            (self.__alpha - self.__beta) ** 2 + 4 * self.__d_sq[0])) / 2
                    else:
                        lambda_k = self.__b_1 + (-(self.__alpha + self.__beta) - sqrt(
                            (self.__alpha - self.__beta) ** 2 + 4 * self.__d_sq[0])) / 2
                    self.__lambda_ks.append(factor(expand(lambda_k)))
            # if alpha equals -beta and is +/- d_2
            elif self.__alpha == -self.__beta and abs(self.__alpha) == abs(self.__d[1]):
                for k in range(1, self.__n + 1):
                    if k <= self.__m:
                        theta_k = '{}*pi/{}'.format(2 * k - 1, self.__n)
                    else:
                        theta_k = '{}*pi/{}'.format((2 * (k - self.__m) - 1), self.__n)
                    self.__theta_ks.append(parse_expr(theta_k))

                for k in range(1, self.__n + 1):
                    if k <= self.__m:
                        lambda_k = self.__b_1 + sqrt(
                            self.__d_sq[0] + self.__d_sq[1] + 2 * self.__d[0] * self.__d[1] * cos(
                                self.__theta_ks[k - 1]))
                    else:
                        lambda_k = self.__b_1 - sqrt(
                            self.__d_sq[0] + self.__d_sq[1] + 2 * self.__d[0] * self.__d[1] * cos(
                                self.__theta_ks[k - 1]))
                    self.__lambda_ks.append(factor(expand(lambda_k)))
            else:
                raise ValueError('Tridiagonal Matrix cannot be solved using Kouachi\'s methods')
        return self.__lambda_ks

    def eigenvalues(self):
        if not self.__lambda_ks:
            if self.__d_size == 2 and self.__b_1 == self.__b_2:
                print('Case IV')
                self.__case4()
            elif self.__d_size == 1 and abs(self.__alpha) + abs(self.__beta) == 0:
                if self.__sub[0] + self.__diag[0] + self.__sup[0] == 1 \
                and(0,0,0) <= (self.__sub[0],self.__diag[0] <= self.__sup[0]) <= (1, 1, 1):
                    print('Case III')
                    self.__case3()
                else:
                    print('Case I')
                    self.__case1()
            elif self.__d_size == 1 and self.__n % 2 == 0 \
                 and self.__alpha * self.__beta == self.__d_sq if (type(self.__alpha) is float) else \
                    (self.__alpha * self.__beta).equals(self.__d_sq):
                print('Case II')
                self.__case2()
            else:
                raise ValueError('Tridiagonal Matrix cannot be solved using Kouachi\'s methods')
        return self.__lambda_ks

    def eigenvectors(self):
        if self.__lambda_ks is None:
            self.eigenvalues()

        for r in range(1, self.__n + 1):
            rho_r = ((-sqrt(self.__d_sq)) ** (self.__n - r)) * (prod(self.__sub[:(max(1, r - 1))]))
            for c in range(1, self.__n + 1):
                zeta_r_c = [self.__b_1 - self.__lambda_ks[c - 1], self.__b_2 - self.__lambda_ks[c - 1]]
                if r == 1:
                    if c > self.__n - 2:
                        self.__ev[r - 1, c - 1] = 1
                    else:
                        fT = ((zeta_r_c[1] - self.__beta) * sin(self.__n * self.__theta_ks[c - 1]))
                        sT = (self.__beta * sin((self.__n - 2) * self.__theta_ks[c - 1]))
                        self.__ev[r - 1, c - 1] = fT - sT
                        self.__ev[r - 1, c - 1] *= (-sqrt(self.__d_sq)) ** (self.__n - 1)
                elif r % 2 == 1:  # if row is odd
                    if self.__n % 2 == 1:  # if odd matrix
                        self.__ev[r - 1, c - 1] = (self.__d_sq * sin((self.__n - r + 2) * self.__theta_ks[c - 1])) - (
                                (self.__beta * zeta_r_c[1] - self.__d_sq) * sin(
                            (self.__n - r) * self.__theta_ks[c - 1]))
                        self.__ev[r - 1, c - 1] *= rho_r
                    else:  # even matrix
                        self.__ev[r - 1, c - 1] = ((zeta_r_c[1] - self.__beta) * sin(
                            (self.__n - r + 1) * self.__theta_ks[c - 1])) - (self.__beta * sin(
                            (self.__n - r - 1) * self.__theta_ks[c - 1]))
                        self.__ev[r - 1, c - 1] *= rho_r
                    if c > self.__n - 2:
                        if r == 2:
                            self.__ev[r - 1, c - 1] = (-zeta_r_c[0] + self.__alpha) / self.__sup[r - 2]
                        else:
                            self.__ev[r - 1, c - 1] = ((-zeta_r_c[r % 2] * self.__ev[r - 2, c - 1]) - (
                                    self.__sub[r - 3] * self.__ev[r - 3, c - 1])) / self.__sup[r - 2]
                else:  # row is even
                    if self.__n % 2 == 1:  # if odd matrix
                        self.__ev[r - 1, c - 1] = ((zeta_r_c[0] - beta) * sqrt(self.__d_sq) * sin(
                            (self.__n - r + 1) * self.__theta_ks[c - 1])) - (beta * sqrt(self.__d_sq) * sin(
                            (self.__n - r - 1) * self.__theta_ks[c - 1]))
                        self.__ev[r - 1, c - 1] *= rho_r
                    else:  # if even matrix
                        self.__ev[r - 1, c - 1] = sqrt(self.__d_sq) * sin(
                            (self.__n - r + 2) * self.__theta_ks[c - 1]) \
                                                  - (self.__beta * zeta_r_c[0] / sqrt(self.__d_sq) -
                                                     sqrt(self.__d_sq)) * sin((self.__n - r) * self.__theta_ks[c - 1])
                        self.__ev[r - 1, c - 1] *= rho_r
                    if c > self.__n - 2:
                        if r == 2:
                            self.__ev[r - 1, c - 1] = (-zeta_r_c[0] + self.__alpha) / self.__sup[r - 2]
                        else:
                            self.__ev[r - 1, c - 1] = ((-zeta_r_c[r % 2] * self.__ev[r - 2, c - 1]) - (
                                    self.__sub[r - 3] * self.__ev[r - 3, c - 1])) / self.__sup[r - 2]

        for r in range(0, self.__n):
            for c in range(0, self.__n):
                self.__ev[r, c] /= self.__ev[self.__n - 1, c]

        return self.__ev
