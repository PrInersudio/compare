from typing import List, Tuple, Set
from polynomial import *
import re
import unittest

logger = logging.getLogger(__name__)
logger.setLevel(logging.DEBUG)
handler = logging.FileHandler(f'log/{__name__}.log')
handler.setLevel(logging.DEBUG)
handler.setFormatter(logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s'))
logger.addHandler(handler) 

def solve_comparisons (equations: List[Tuple[Polynomial, int]]) -> Set[Tuple[int,int]]:
    logger.debug(f'solve_comparisons. equations={equations}')
    result = equations[0][0].roots(equations[0][1])
    for poly, module in equations[1:]:
        logger.debug(f'solve_comparisons. equations={equations}, result={result}')
        new_result = set()
        for a,m in result:
            new_poly = Polynomial()
            for exp, coeff in poly.monomials():
                new_poly += coeff * binomial_theorem(m,a,exp)
                logger.debug(f'solve_comparisons. result={result}, new_poly={new_poly}')
            for a1,m1 in new_poly.roots(module):
                new_result.add((m*a1+a,m*m1))
        result = new_result
    logger.debug(f'solve_comparisons. equations={equations}, result={result}')
    return result

def parse_equation(equation: str) -> Tuple[Polynomial, int]:
    logger.debug(f'parse_equation. equation={equation}')
    pattern = r'(.+?)\s*=\s*(.+?)\s*\(\s*mod\s*(\d+)\s*\)'
    match = re.match(pattern, equation)
    if not match:
        logger.critical(f'parse_equation. equation={equation}. Неверный формат уравнения.')
        raise ValueError("Неверный формат строки. Ожидается формат 'P(x) = Q(x) (mod m)'")
    module = int(match.group(3).strip())
    poly = Polynomial(match.group(1).strip()) - Polynomial(match.group(2).strip())
    logger.debug(f'parse_equation. equation={equation}, poly={poly}, module={module}')
    return poly, module

class TestComparisons(unittest.TestCase):
    
    def setUp(self):
        return super().setUp()
    
    def test_solve_comparisons(self):
        system = [
            (Polynomial("x^2 - 5"), 13),
            (Polynomial("x - 5"), 13),
        ]
        self.assertEqual(solve_comparisons(system), set(), "Ошибка solve_comparisons в случае отсуствия решения первого уравнения.")
        system = [
            (Polynomial("x - 5"), 13),
            (Polynomial("x^2 - 5"), 13),
        ]
        self.assertEqual(solve_comparisons(system), set(), "Ошибка solve_comparisons в случае отсуствия решения второго уравнения.")
        system = [
            (Polynomial("162x^10+135x^9+162x^8+56x^4+162x^3+162x^2+113x+188"), 189),
            (Polynomial("x^3+11x^2+2x+8"), 16),
            (Polynomial("x - 12"), 16)
        ]
        result = {
            (23692, 48384), (37516, 48384), (48316, 48384), (34492, 48384), (14620, 48384), (29740, 48384),
            (40540, 48384), (32764, 48384), (26716, 48384), (43564, 48384), (44860, 48384), (35788, 48384),
            (46588, 48384), (47884, 48384), (41836, 48384), (38812, 48384), (20668, 48384), (12028, 48384),
            (2956, 48384), (4252, 48384), (15052, 48384), (9004, 48384), (5980, 48384), (21100, 48384), (1228, 48384),
            (13324, 48384), (7276, 48384), (18076, 48384), (24124, 48384), (10300, 48384), (36220, 48384), (5548, 48384),
            (2524, 48384), (16348, 48384), (27148, 48384), (28444, 48384), (8572, 48384), (33196, 48384), (11596, 48384),
            (19372, 48384), (45292, 48384), (39244, 48384), (25420, 48384), (30172, 48384), (31468, 48384), (17644, 48384),
            (42268, 48384), (22396, 48384)
        }
        self.assertEqual(solve_comparisons(system), result, "Ошибка solve_comparisons в общем случае.")

if __name__ == '__main__':
    unittest.main()