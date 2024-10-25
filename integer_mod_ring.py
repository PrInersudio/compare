from __future__ import annotations
from ring import Ring
import ctypes
import unittest
from types import MethodType

mul = ctypes.CDLL('./libmul.so').mul
mul.restype = ctypes.c_uint
mul.argtypes = [ctypes.c_uint, ctypes.c_uint, ctypes.c_uint]

class IntegerModRing(Ring):

    def __init__(self, m: int):
        super().__init__()
        self.m = m
        if self.m & 1: self._mul = self.__mul_odd_module
        else: self._mul = self.__mul_even_module
        self.has_multiplicative_identity = True

    def __eq__(self, other: 'IntegerModRing') -> bool:
        if not isinstance(other, IntegerModRing): return False
        return self.m == other.m
    
    def __repr__(self) -> str:
        return f'IntegerModRing({self.m})'

    def __str__(self) -> str:
        return f'Кольцо целых чисел по модулю {self.m}'
    
    def __mul_odd_module(self, a: 'IntegerModRing.Element', b: 'IntegerModRing.Element') -> 'IntegerModRing.Element':
            return self(mul(a.value, b.value, self.m))
                    
    def __mul_even_module(self, a: 'IntegerModRing.Element', b: 'IntegerModRing.Element') -> 'IntegerModRing.Element':
        return self(a.value * b.value)
    
    def __call__(self, value: 'IntegerModRing.Element' | int = 0):
        return self.Element(self, int(value))
    
    def __iter__(self):
        x = self(0)
        for _ in range(self.m):
            yield x
            x += 1
    
    class Element(Ring.Element):

        def __init__(self, ring: 'IntegerModRing', value: 'IntegerModRing.Element' | int):
            super().__init__(ring, value)
            self.value %= self.ring.m
        
        def __pow__(self, other: int) -> 'IntegerModRing.Element':
            return self.ring(pow(self.value, other, self.ring.m))
        
        def __mul__(self, other: 'IntegerModRing.Element' | int) -> 'IntegerModRing.Element':
            other, check_result = self._check(other)
            if not check_result: return NotImplemented
            return self.ring._mul(self, other)

        def __truediv__(self, other: 'IntegerModRing.Element' | int) -> 'IntegerModRing.Element':
            other, check_result = self._check(other)
            if not check_result: return NotImplemented
            return self * (other ** (-1))
        
        def __repr__(self) -> str:
            return f'IntegerModRing.Element({repr(self.ring)}, {self.value})'

        def __str__(self) -> str:
            return f'[{self.value}]_{self.ring.m}'

class TestIntegerModRing(unittest.TestCase):

    def setUp(self):
        return super().setUp()
    
    def test_init(self):
        self.assertEqual(IntegerModRing(7)(9), IntegerModRing(7)(2), "Ошибка создания элемента кольца.")

    def test_int(self):
        self.assertEqual(int(IntegerModRing(7)(5)), 5, "Ошибка преобразования элемента кольца в int.")

    def test_bool(self):
        self.assertTrue(IntegerModRing(7)(5), "Ошибка преобразования элемента кольца в bool True.")
        self.assertFalse(IntegerModRing(7)(), "Ошибка преобразования элемента кольца в bool False.")

    def test_check(self):
        with self.assertRaises(ValueError, msg="Ошибка _op_chec. В случае другого кольца не подымает ValueError."):
            IntegerModRing(5)(3)._check(IntegerModRing(7)(3))
        self.assertEqual(IntegerModRing(5)(3)._check(IntegerModRing(5)(7)), IntegerModRing(5)(2),
            "Ошибка _op_chec. В случае такого же кольца не возвращает переданный элемент."
        )
        self.assertEqual(IntegerModRing(5)(3)._check(7), IntegerModRing(5)(2),
            "Ошибка _op_chec. В случае int не возвращает cоответствующий элемент."
        )

    def test_check(self):
        self.assertFalse(
            IntegerModRing(5)(3)._check(IntegerModRing(7)(3))[1],
            "Ошибка _check. В случае другого кольца не возвращеает False."
        )
        self.assertEqual(
            IntegerModRing(5)(3)._check(IntegerModRing(7)())[0],
            0,
            "Ошибка _check. В случае 0."
        )
        self.assertEqual(
            IntegerModRing(5)(3)._check(IntegerModRing(5)(7))[0],
            IntegerModRing(5)(2),
            "Ошибка _check. В случае такого же кольца не возвращает переданный элемент."
        )
        self.assertEqual(
            IntegerModRing(5)(3)._check(7)[0],
            IntegerModRing(5)(2),
            "Ошибка _check. В случае int не возвращает cоответствующий элемент."
        )

    def test_add(self):
        self.assertEqual(IntegerModRing(7)(5) + 4, IntegerModRing(7)(2), "Ошибка сложения элементов кольца.")

    def test_sub(self):
        self.assertEqual(IntegerModRing(7)(2) - 4, IntegerModRing(7)(5), "Ошибка вычитания элементов кольца.")

    def test_pow(self):
        self.assertEqual(IntegerModRing(7)(5) ** 2, IntegerModRing(7)(4), "Ошибка возведения в положительную степень элементов кольца.")
        self.assertEqual(IntegerModRing(7)(5) ** (-2), IntegerModRing(7)(2), "Ошибка возведения в отрицательную степень элементов кольца.")
        with self.assertRaises(ValueError, msg="""Ошибка возведения в отрицательную степень элементов кольца
            (в случае отсуствия обратного должен поднимать ValueError)."""
        ): IntegerModRing(4)(2) ** (-2)

    def test_mul(self):
        self.assertEqual(mul(15, 6, 17), 5, "Ошибка умножения  в случае нечётного модуля (вызов mul напрямую).")
        self.assertEqual(IntegerModRing(17)(15) * 6, IntegerModRing(17)(5), "Ошибка умножения  в случае нечётного модуля.")
        self.assertEqual(IntegerModRing(4)(2) * 2, 0, "Ошибка умножения на int в случае чётного модуля.")

    def test_truediv(self):
        with self.assertRaises(ValueError, msg="Деление на необратимый int должно давать ValueErrror."): IntegerModRing(5)(3) / 0
        self.assertEqual(IntegerModRing(5)(3) / 2, IntegerModRing(5)(4), "Ошибка деления в общем случае.")

    def test_gt(self):
        self.assertTrue(IntegerModRing(5)(3) > IntegerModRing(5)(2), "Ошибка gt в случае True.")
        self.assertFalse(IntegerModRing(5)(2) > IntegerModRing(5)(3), "Ошибка gt в случае False.")

    def test_lt(self):
        self.assertTrue(IntegerModRing(5)(2) < IntegerModRing(5)(3), "Ошибка lt в случае True.")
        self.assertFalse(IntegerModRing(5)(3) < IntegerModRing(5)(2), "Ошибка lt в случае False.")

    def test_ge(self):
        self.assertFalse(IntegerModRing(5)(2) >= IntegerModRing(5)(3), "Ошибка ge в случае False.")
        self.assertTrue(IntegerModRing(5)(3) >= IntegerModRing(5)(2), "Ошибка ge в случае не равно True.")
        self.assertTrue(IntegerModRing(5)(3) >= IntegerModRing(5)(3), "Ошибка ge в случае равно True.")

    def test_le(self):
        self.assertFalse(IntegerModRing(5)(3) <= IntegerModRing(5)(2), "Ошибка le в случае False.")
        self.assertTrue(IntegerModRing(5)(2) <= IntegerModRing(5)(3), "Ошибка le в случае не равно True.")
        self.assertTrue(IntegerModRing(5)(3) <= IntegerModRing(5)(3), "Ошибка le в случае равно True.")

    def test_iter(self):
        ring = IntegerModRing(5)
        self.assertEqual(list(ring), [ring(0), ring(1), ring(2), ring(3), ring(4)], "Ошибка iter.")

if __name__ == '__main__':
    unittest.main()