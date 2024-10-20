from math import gcd
import logging
from typing import List, Tuple
import unittest

logger = logging.getLogger(__name__)
logger.setLevel(logging.DEBUG)
handler = logging.FileHandler(f'log/{__name__}.log')
handler.setLevel(logging.DEBUG)
handler.setFormatter(logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s'))
logger.addHandler(handler)

def solve_linear_comparison(a: int,b: int,m: int) -> List[Tuple[int,int]]:
    """
    ax = b (mod m)
    """
    logger.debug(f'solve_linear_comparison. a={a}, b={b}, m={m}')
    d = gcd(a,m)
    logger.debug(f'solve_linear_comparison. a={a}, b={b}, m={m}, d={d}')
    if b % d != 0:
        logger.debug(f'solve_linear_comparison. a={a}, b={b}, m={m}. Нет решений.')
        return []
    a0,b0,m0 = a//d,b//d,m//d
    logger.debug(f'solve_linear_comparison. a={a}, b={b}, m={m}, a0={a0}, b0={b0}, m0={m0}')
    x0 = b0 * pow(a0,-1,m0)
    logger.debug(f'solve_linear_comparison. a={a}, b={b}, m={m}, d={d}, a0={a0}, b0={b0}, m0={m0}, x0={x0}')
    result = [((x0 + i) % m, m) for i in range(0,m,m0)]
    logger.debug(f'solve_linear_comparison. a={a}, b={b}, m={m}, result={result}')
    return result

def crt (equations: List[Tuple[int,int]]) -> Tuple[int,int]:
    """
    x = a_i (mod m_i)
    """
    logger.debug(f'crt. equations={equations}')
    M = 1
    x = 0
    for a,m in equations:
        x += M*((a-x)*pow(M,-1,m))
        M *= m
        logger.debug(f'crt. equations={equations}, a={a}, m={m}, x={x}, M={M}')
    result = (x % M,M)
    logger.debug(f'crt. equations={equations}, result={result}')
    return result

class TestLinearComparisons(unittest.TestCase):

    def test_solve_linear_comparison(self):
        self.assertCountEqual(solve_linear_comparison(2, 3, 5), [(4,5)], "Ошибка solve_linear_comparison в случае одного решения.")
        self.assertCountEqual(solve_linear_comparison(2, 3, 4), [], "Ошибка solve_linear_comparison в случае отсуствия решений.")
        self.assertCountEqual(solve_linear_comparison(2, 4, 4), [(0,4), (2,4)], "Ошибка solve_linear_comparison в случае нескольких решений.")
        self.assertCountEqual(solve_linear_comparison(5, 5, 5), [(0,5), (1,5), (2,5), (3,5), (4,5)], "Ошибка solve_linear_comparison в случае 0x=0.")
        self.assertCountEqual(solve_linear_comparison(2, 1, 1), [(0,1)], "Ошибка solve_linear_comparison в случае m=1.")

    def test_crt(self):
        self.assertEqual(crt([(4,5)]), (4,5), "Ошибка crt в случае одного сравнения.")
        self.assertEqual(crt([(2, 12), (16, 17), (0, 25)]), (50,5100), "Ошибка crt в случае нескольких сравнений.")

if __name__ == '__main__':
    unittest.main()