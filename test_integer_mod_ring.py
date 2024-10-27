"""Тестирование integer_mod_ring."""

import unittest

from integer_mod_ring import IntegerModRing


class TestIntegerModRing(unittest.TestCase):
    """Тестирование integer_mod_ring."""

    def test_init(self):
        """Тестирование инициализации IntegerModRing и IntegerModRing.Element."""
        self.assertEqual(IntegerModRing(7)(9),
                         IntegerModRing(7)(2),
                         "Ошибка создания элемента кольца.")

    def test_int(self):
        """Тестирование преобразования IntegerModRing.Element в int."""
        self.assertEqual(int(IntegerModRing(7)(5)), 5,
                         "Ошибка преобразования элемента кольца в int.")

    def test_bool(self):
        """Тестирование преобразования IntegerModRing.Element в bool."""
        self.assertTrue(IntegerModRing(7)(5), "Ошибка преобразования элемента кольца в bool True.")
        self.assertFalse(IntegerModRing(7)(), "Ошибка преобразования элемента кольца в bool False.")

    # pylint: disable=protected-access
    def test_check(self):
        """Тестирование метода _check,
        который при выполнении сравнений и операций
        проверяет, что элементы принадлежат одному кольцу.
        """
        self.assertFalse(
            IntegerModRing(5)(3)._check(IntegerModRing(7)(3))[1],
            "Ошибка _check. В случае другого кольца не возвращает False."
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
            "Ошибка _check. В случае int не возвращает соответствующий элемент."
        )
    # pylint: enable=protected-access

    def test_add(self):
        """Тестирование сложения элементов кольца."""
        self.assertEqual(IntegerModRing(7)(5) + 4, IntegerModRing(7)(2),
                         "Ошибка сложения элементов кольца.")

    def test_sub(self):
        """Тестирование вычитания элементов кольца."""
        self.assertEqual(IntegerModRing(7)(2) - 4, IntegerModRing(7)(5),
                         "Ошибка вычитания элементов кольца.")

    def test_pow(self):
        """Тестирование возведения в степень элементов кольца."""
        self.assertEqual(IntegerModRing(7)(5) ** 2, IntegerModRing(7)(4),
                         "Ошибка возведения в положительную степень элементов кольца.")
        self.assertEqual(IntegerModRing(7)(5) ** (-2), IntegerModRing(7)(2),
                         "Ошибка возведения в отрицательную степень элементов кольца.")
        with self.assertRaises(ValueError,
                               msg="""Ошибка возведения в отрицательную степень элементов кольца
                               (в случае отсутствия обратного должен поднимать ValueError)."""):
            IntegerModRing(4)(2) ** (-2) # pylint: disable=expression-not-assigned

    def test_mul(self):
        """Тестирование умножения элементов кольца."""
        self.assertEqual(IntegerModRing(17)(15) * 6, IntegerModRing(17)(5),
                         "Ошибка умножения  в случае нечётного модуля.")
        self.assertEqual(IntegerModRing(4)(2) * 2, 0,
                         "Ошибка умножения на int в случае чётного модуля.")

    def test_truediv(self):
        """Тестирование деления элементов кольца."""
        with self.assertRaises(ValueError,
                               msg="Деление на необратимый int должно давать ValueError."):
            IntegerModRing(5)(3) / 0 # pylint: disable=expression-not-assigned
        self.assertEqual(IntegerModRing(5)(3) / 2, IntegerModRing(5)(4),
                         "Ошибка деления в общем случае.")

    def test_gt(self):
        """Тестирование сравнения элементов кольца (больше)."""
        self.assertTrue(IntegerModRing(5)(3) > IntegerModRing(5)(2), "Ошибка gt в случае True.")
        self.assertFalse(IntegerModRing(5)(2) > IntegerModRing(5)(3), "Ошибка gt в случае False.")

    def test_lt(self):
        """Тестирование сравнения элементов кольца (меньше)."""
        self.assertTrue(IntegerModRing(5)(2) < IntegerModRing(5)(3), "Ошибка lt в случае True.")
        self.assertFalse(IntegerModRing(5)(3) < IntegerModRing(5)(2), "Ошибка lt в случае False.")

    def test_ge(self):
        """Тестирование сравнения элементов кольца (больше-равно)."""
        self.assertFalse(IntegerModRing(5)(2) >= IntegerModRing(5)(3), "Ошибка ge в случае False.")
        self.assertTrue(IntegerModRing(5)(3) >= IntegerModRing(5)(2),
                        "Ошибка ge в случае не равно True.")
        self.assertTrue(IntegerModRing(5)(3) >= IntegerModRing(5)(3),
                        "Ошибка ge в случае равно True.")

    def test_le(self):
        """Тестирование сравнения элементов кольца (меньше-равно)."""
        self.assertFalse(IntegerModRing(5)(3) <= IntegerModRing(5)(2), "Ошибка le в случае False.")
        self.assertTrue(IntegerModRing(5)(2) <= IntegerModRing(5)(3),
                        "Ошибка le в случае не равно True.")
        self.assertTrue(IntegerModRing(5)(3) <= IntegerModRing(5)(3),
                        "Ошибка le в случае равно True.")

    def test_iter(self):
        """Тестирование итерирования по элементам кольца."""
        ring = IntegerModRing(5)
        self.assertEqual(list(ring), [ring(0), ring(1), ring(2), ring(3), ring(4)], "Ошибка iter.")


if __name__ == '__main__':
    unittest.main()
