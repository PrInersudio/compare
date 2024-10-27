"""Тестирование comparisons."""

import unittest

from field import GF
from polynomial import Polynomial
from integer_mod_ring import IntegerModRing
from comparisons import solve_comparisons, parse_equation


class TestComparisons(unittest.TestCase):
    """Тестирование comparisons."""

    def test_solve_comparisons(self):
        """Тест solve_comparisons."""
        system = [
            Polynomial(IntegerModRing(13))("x^2 - 5"),
            Polynomial(IntegerModRing(13))("x - 5"),
            ]
        self.assertEqual(solve_comparisons(system), set(),
                         "Ошибка solve_comparisons в случае отсутствия решения первого уравнения.")
        system = [
            Polynomial(IntegerModRing(13))("x - 5"),
            Polynomial(IntegerModRing(13))("x^2 - 5"),
            ]
        self.assertEqual(solve_comparisons(system), set(),
                         "Ошибка solve_comparisons в случае отсутствия решения второго уравнения.")
        system = [
            Polynomial(IntegerModRing(189))("162x^10+135x^9+162x^8+56x^4+162x^3+162x^2+113x+188"),
            Polynomial(IntegerModRing(16))("x^3+11x^2+2x+8"),
            Polynomial(IntegerModRing(16))("x - 12")
            ]
        ring = IntegerModRing(48384)
        result = {
            ring(23692), ring(37516), ring(48316), ring(34492), ring(14620), ring(29740),
            ring(40540), ring(32764), ring(26716), ring(43564), ring(44860), ring(35788),
            ring(46588), ring(47884), ring(41836), ring(38812), ring(20668), ring(12028),
            ring(2956), ring(4252), ring(15052), ring(9004), ring(5980), ring(21100), ring(1228),
            ring(13324), ring(7276), ring(18076), ring(24124), ring(10300), ring(36220), ring(5548),
            ring(2524), ring(16348), ring(27148), ring(28444), ring(8572), ring(33196),
            ring(11596), ring(19372), ring(45292), ring(39244), ring(25420), ring(30172),
            ring(31468), ring(17644), ring(42268), ring(22396),
        }
        self.assertEqual(solve_comparisons(system), result,
                         "Ошибка solve_comparisons в общем случае.")

    def test_parse_equation(self):
        """Тест parse_equation."""
        self.assertEqual(parse_equation("x^2 + 1 = 1 (mod 13)"),
                         Polynomial(GF(13))("x^2"),
                         "Ошибка parse_equation в случае константы справа.")
        self.assertEqual(parse_equation("1 = x^2 + 1 (mod 13)"),
                         Polynomial(GF(13))("12x^2"),
                         "Ошибка parse_equation в случае константы слева.")
        self.assertEqual(parse_equation("1 = 1 (mod 13)"),
                         Polynomial(GF(13))(),
                         "Ошибка parse_equation в случае констант c обоих сторон.")
        self.assertEqual(parse_equation("x^2 + 1 = x^2 + 1 (mod 13)"),
                         Polynomial(GF(13))(),
                         "Ошибка parse_equation в общем случае.")
        self.assertEqual(parse_equation("x^2+1=x^2+1(mod13)"),
                         Polynomial(GF(13))(),
                         "Ошибка parse_equation в случае слитного написания.")


if __name__ == '__main__':
    unittest.main()
