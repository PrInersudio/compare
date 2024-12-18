"""Реализует кольцо вычетов."""

from __future__ import annotations

import ctypes

from ring import Ring


_mul = ctypes.CDLL('./libmul.so').mul
_mul.restype = ctypes.c_uint
_mul.argtypes = [ctypes.c_uint, ctypes.c_uint, ctypes.c_uint, ctypes.c_uint]

class IntegerModRing(Ring):
    """Кольцо вычетов."""

    def __init__(self, m: int):
        super().__init__()
        self.m = m
        if self.m & 1:
            self.__b1 = (1 << (self.m.bit_length() << 1)) % self.m
            self._mul = lambda a, b: self(_mul(a, b, self.m, self.__b1))
        else: self._mul = lambda a, b: self(a * b)
        self.has_multiplicative_identity = True
        self.is_commutative = True

    def __eq__(self, other: 'IntegerModRing') -> bool:
        if not isinstance(other, IntegerModRing):
            return False
        return self.m == other.m

    def __repr__(self) -> str:
        return f'IntegerModRing({self.m})'

    def __str__(self) -> str:
        return f'Кольцо целых чисел по модулю {self.m}'

    def __call__(self, value: 'IntegerModRing.Element' | int = 0):
        return self.Element(self, int(value))

    def __iter__(self):
        x = self(0)
        for _ in range(self.m):
            yield x
            x += 1


    class Element(Ring.Element):
        """Элемент кольца вычетов."""

        def __init__(self, ring: 'IntegerModRing', value: 'IntegerModRing.Element' | int):
            super().__init__(ring, value)
            self.value %= self.ring.m

        def __pow__(self, other: int) -> 'IntegerModRing.Element':
            return self.ring(pow(self.value, other, self.ring.m))

        def __mul__(self, other: 'IntegerModRing.Element' | int) -> 'IntegerModRing.Element':
            other, check_result = self._check(other)
            if not check_result:
                return NotImplemented
            return self.ring._mul(self.value, other.value)

        def __truediv__(self, other: 'IntegerModRing.Element' | int) -> 'IntegerModRing.Element':
            other, check_result = self._check(other)
            if not check_result:
                return NotImplemented
            return self * (other ** (-1))

        def __repr__(self) -> str:
            return f'IntegerModRing.Element({repr(self.ring)}, {self.value})'

        def __str__(self) -> str:
            return f'[{self.value}]_{self.ring.m}'
