"""Тестирование linear_comparisons."""

import unittest

from linear_comparisons import solve_linear_comparison, crt
from integer_mod_ring import IntegerModRing


class TestLinearComparisons(unittest.TestCase):
    """Тестирование linear_comparisons."""

    def test_solve_linear_comparison(self):
        """Тестирование solve_linear_comparison."""
        self.assertCountEqual(solve_linear_comparison(2, 3, 5),
                              {IntegerModRing(5)(4)},
                              "Ошибка solve_linear_comparison в случае одного решения.")
        self.assertCountEqual(solve_linear_comparison(2, 3, 4), {},
                              "Ошибка solve_linear_comparison в случае отсутствия решений.")
        self.assertCountEqual(solve_linear_comparison(2, 4, 4),
                              {IntegerModRing(4)(0), IntegerModRing(4)(2)},
                              "Ошибка solve_linear_comparison в случае нескольких решений.")
        self.assertCountEqual(solve_linear_comparison(5, 5, 5), set(IntegerModRing(5)),
                              "Ошибка solve_linear_comparison в случае 0x=0.")
        self.assertCountEqual(solve_linear_comparison(2, 1, 1), {IntegerModRing(1)(0)},
                              "Ошибка solve_linear_comparison в случае m=1.")

    def test_crt(self):
        """Тестирование crt."""
        self.assertEqual(crt([IntegerModRing(5)(4)]), IntegerModRing(5)(4),
                         "Ошибка crt в случае одного сравнения.")
        self.assertEqual(
            crt([
                IntegerModRing(12)(2),
                IntegerModRing(17)(16),
                IntegerModRing(25)(0),
                ]),
            IntegerModRing(5100)(50), "Ошибка crt в случае нескольких сравнений."
        )


if __name__ == '__main__':
    unittest.main()
