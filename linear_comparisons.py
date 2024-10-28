""" Модуль реализует:
solve_linear_comparison - функцию, которая решает линейное сравнение.
crt - функцию, которая принимает список элементов колец вычетов по модулям m1,m2,...,mk и
возвращает эквивалентный им элемент по модулю m1*m2*...*mk.
"""

from math import gcd
import logging
from typing import List, Set
from integer_mod_ring import IntegerModRing


logger = logging.getLogger(__name__)
logger.setLevel(logging.DEBUG)
handler = logging.FileHandler(f'log/{__name__}.log')
handler.setLevel(logging.DEBUG)
handler.setFormatter(logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s'))
logger.addHandler(handler)


def solve_linear_comparison(a: int,b: int,m: int) -> Set[IntegerModRing.Element]:
    """ax = b (mod m)"""
    logger.debug('solve_linear_comparison. a = %s, b = %s, m = %s.', a, b, m)
    d = gcd(a,m)
    logger.debug('solve_linear_comparison. a = %s, b = %s, m = %s, d = %s.', a, b, m, d)
    if b % d != 0:
        logger.debug('solve_linear_comparison. a = %s, b = %s, m = %s. Нет решений.', a, b, m)
        return set()
    a0,b0,m0 = a//d,b//d,m//d
    logger.debug('solve_linear_comparison. a = %s, b = %s, m = %s, a0 = %s, b0 = %s, m0 = %s.',
                 a, b, m, a0, b0, m0)
    x0 = b0 * pow(a0,-1,m0)
    logger.debug('solve_linear_comparison.'
                 ' a = %s, b = %s, m = %s, d = %s, a0 = %s, b0 = %s, m0 = %s, x0 = %s.',
                 a, b, m, d, a0, b0, m0, x0)
    ring = IntegerModRing(m)
    result = {ring(x0 + i) for i in range(0,m,m0)}
    logger.debug('solve_linear_comparison. a = %s, b = %s, m = %s, result = %s.', a, b, m, result)
    return result

def crt (equations: List[IntegerModRing.Element]) -> IntegerModRing.Element:
    """x = a_i (mod m_i)"""
    logger.debug('crt. equations = %s.', equations)
    module = 1
    x = 0
    for element in equations:
        x += module*((int(element)-x)*pow(module,-1,element.ring.m))
        module *= element.ring.m
        logger.debug('crt. equations = %s, element = %s, element.ring.m = %s, x = %s, module = %s.',
                     equations, element, element.ring.m, x, module)
    result = IntegerModRing(module)(x)
    logger.debug('crt. equations = %s, result = %s.', equations, result)
    return result
