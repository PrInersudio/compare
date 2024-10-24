import logging
from typing import List, Tuple, Set
from polynomial import Polynomial
from integer_mod_ring import IntegerModRing
from field import GF
import re
import unittest

logger = logging.getLogger(__name__)
logger.setLevel(logging.DEBUG)
handler = logging.FileHandler(f'log/{__name__}.log')
handler.setLevel(logging.DEBUG)
handler.setFormatter(logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s'))
logger.addHandler(handler) 

def solve_comparisons (equations: List[Polynomial.Element]) -> Set[IntegerModRing.Element]:
    logger.debug(f'solve_comparisons. equations={equations}')
    result = equations[0].roots()
    for poly in equations[1:]:
        logger.debug(f'solve_comparisons. equations={equations}, result={result}')
        new_result = set()
        for a in result:
            new_poly_ring = Polynomial(poly.ring.basic_ring)
            new_poly = new_poly_ring()
            for exp, coeff in poly.monomials():
                new_poly += coeff * new_poly_ring.binomial_theorem(a.ring.m,a.value,exp)
                logger.debug(f'solve_comparisons. result={result}, new_poly={new_poly}')
            for a1 in new_poly.roots():
                new_result.add(IntegerModRing(a.ring.m*a1.ring.m)(a.ring.m*a1.value+a.value))
        result = new_result
    logger.debug(f'solve_comparisons. equations={equations}, result={result}')
    return result

def parse_equation(equation: str) -> Polynomial.Element:
    logger.debug(f'parse_equation. equation={equation}')
    pattern = r'(.+?)\s*=\s*(.+?)\s*\(\s*mod\s*(\d+)\s*\)'
    match = re.match(pattern, equation)
    if not match:
        logger.critical(f'parse_equation. equation={equation}. Неверный формат уравнения.')
        raise ValueError("Неверный формат строки. Ожидается формат 'P(x) = Q(x) (mod m)'")
    ring = Polynomial(IntegerModRing(int(match.group(3).strip())))
    poly = ring(match.group(1).strip()) - ring(match.group(2).strip())
    logger.debug(f'parse_equation. equation={equation}, poly={poly}')
    return poly

class TestComparisons(unittest.TestCase):
    
    def setUp(self):
        return super().setUp()
    
    def test_solve_comparisons(self):
        system = [Polynomial(IntegerModRing(13))("x^2 - 5"), Polynomial(IntegerModRing(13))("x - 5")]
        self.assertEqual(solve_comparisons(system), set(), "Ошибка solve_comparisons в случае отсуствия решения первого уравнения.")
        system = [Polynomial(IntegerModRing(13))("x - 5"), Polynomial(IntegerModRing(13))("x^2 - 5")]
        self.assertEqual(solve_comparisons(system), set(), "Ошибка solve_comparisons в случае отсуствия решения второго уравнения.")
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
            ring(2524), ring(16348), ring(27148), ring(28444), ring(8572), ring(33196), ring(11596),
            ring(19372), ring(45292), ring(39244), ring(25420), ring(30172), ring(31468), ring(17644),
            ring(42268), ring(22396)
        }
        self.assertEqual(solve_comparisons(system), result, "Ошибка solve_comparisons в общем случае.")

    def test_parse_equation(self):
        self.assertEqual(parse_equation("x^2 + 1 = 1 (mod 13)"), Polynomial(GF(13))("x^2"), "Ошибка parse_equation в случае константы справа.")
        self.assertEqual(parse_equation("1 = x^2 + 1 (mod 13)"), Polynomial(GF(13))("12x^2"), "Ошибка parse_equation в случае константы слева.")
        self.assertEqual(parse_equation("1 = 1 (mod 13)"), Polynomial(GF(13))(), "Ошибка parse_equation в случае констант c обоих сторон.")
        self.assertEqual(parse_equation("x^2 + 1 = x^2 + 1 (mod 13)"), Polynomial(GF(13))(), "Ошибка parse_equation в общем случае.")
        self.assertEqual(parse_equation("x^2+1=x^2+1(mod13)"), Polynomial(GF(13))(), "Ошибка parse_equation в случае слитного написания.")

if __name__ == '__main__':
    unittest.main()