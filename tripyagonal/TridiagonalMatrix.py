from .AbstractMatrix import AbstractMatrix


class TridiagonalMatrix(AbstractMatrix):

    def __init__(self, **kwargs):
        super(TridiagonalMatrix, self).__init__(**kwargs)

    def _perturb(self):
        self._data[0, 0] -= self.alpha()
        self._data[self._size - 1, self._size - 1] -= self.beta()
