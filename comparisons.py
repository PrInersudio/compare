""" Реализует функции:
solve_comparisons для решения систем полиномиальных сравнений;
parse_equation для для парсинга сравнения из строки.
"""

import re


import logging
from typing import List, Set
from polynomial import Polynomial
from integer_mod_ring import IntegerModRing


logger = logging.getLogger(__name__)
logger.setLevel(logging.DEBUG)
handler = logging.FileHandler(f'log/{__name__}.log')
handler.setLevel(logging.DEBUG)
handler.setFormatter(logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s'))
logger.addHandler(handler)


def solve_comparisons (equations: List[Polynomial.Element]) -> Set[IntegerModRing.Element]:
    """Решает системы полиномиальных сравнений."""
    logger.debug('solve_comparisons. equations = %s.', equations)
    result = equations[0].roots()
    for poly in equations[1:]:
        logger.debug('solve_comparisons. equations = %s, result = %s.', equations, result)
        new_result = set()
        for a in result:
            new_poly_ring = Polynomial(poly.ring.basic_ring)
            new_poly = new_poly_ring()
            for exp, coeff in poly.monomials():
                new_poly += coeff * new_poly_ring.binomial_theorem(a.ring.m,a.value,exp)
                logger.debug('solve_comparisons. result = %s, new_poly = %s.', result, new_poly)
            for a1 in new_poly.roots():
                new_result.add(IntegerModRing(a.ring.m*a1.ring.m)(a.ring.m*a1.value+a.value))
        result = new_result
    logger.debug('solve_comparisons. equations = %s, result = %s.', equations, result)
    return result

def parse_equation(equation: str) -> Polynomial.Element:
    """Парсит строку со сравним в многочлен."""
    logger.debug('parse_equation. equation = %s.', equation.strip())
    pattern = r'(.+?)\s*=\s*(.+?)\s*\(\s*mod\s*(\d+)\s*\)'
    match = re.match(pattern, equation)
    if not match:
        logger.critical('parse_equation. equation = %s. Неверный формат уравнения.', equation)
        raise ValueError("Неверный формат строки. Ожидается формат 'P(x) = Q(x) (mod m)'")
    ring = Polynomial(IntegerModRing(int(match.group(3).strip())))
    poly = ring(match.group(1).strip()) - ring(match.group(2).strip())
    logger.debug('parse_equation. equation = %s, poly = %s.', equation.strip(), poly)
    return poly
