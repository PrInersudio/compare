from typing import Tuple, Set, Literal
import logging
from linear_comparisons import solve_linear_comparison
from sympy import factorint
import unittest
from math import gcd

logger = logging.getLogger(__name__)
logger.setLevel(logging.DEBUG)
handler = logging.FileHandler(f'log/{__name__}.log')
handler.setLevel(logging.DEBUG)
handler.setFormatter(logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s'))
logger.addHandler(handler)

def Jacobi_symbol(a: int, m: int) -> Literal[-1, 0, 1]:
    logger.debug(f'Jacobi_symbol. a={a}, m={m}.')
    result = 1
    a %= m
    while a:
        logger.debug(f'Шаг Jacobi_symbol. a={a}, m={m}, result={result}.')
        if not a&1: result *= 1 if m&7 in [1,7] else -1; a>>=1; continue
        result *= 1 if m&3 == 1 or a&3 == 1 else -1
        a,m = m%a,a
    if m != 1: result = 0
    logger.debug(f'Jacobi_symbol. a={a}, m={m}, result={result}.')
    return result

class QuadraticResidue:

    def __init__(self, a: int, p: int):
        self.a = a
        self.p = p
        
    def __repr__(self) -> str:
        return f'x^2 = {self.a} (mod {self.p})'
    
    def __str__(self) -> str:
        return self.__repr__()
    
    def __eq__(self, other: 'QuadraticResidue') -> bool:
        return self.p == other.p and self.a % self.p == other.a % other.p

    def __roots_p_modulo_4_equals_3(self) -> Set[int]:
        logger.debug(f'QuadraticResidue.__roots_p_modulo_4_equals_3. self={self}.')
        c  = pow(self.a, (self.p + 1) >> 2 , self.p)
        roots = {c, self.p-c}
        logger.debug(f'__roots_quadratic_residue_p_modulo_4_equals_3. self={self}, roots={roots}')
        return roots
    
    def __roots_p_modulo_8_equals_5(self) -> Set[int]:
        logger.debug(f'QuadraticResidue.__roots_p_modulo_8_equals_5. self={self}.')
        k = self.p >> 3
        t = 0 if pow(self.a, (self.p-1) >> 2, self.p) == 1 else 1
        c = (pow(self.a, k + 1, self.p) * pow(2, ((k<<1)+1) * t, self.p)) % self.p
        logger.debug(f'QuadraticResidue.__roots_p_modulo_8_equals_5. self={self}, k={k}, t={t}, c={c}.')
        roots = {c, self.p-c}
        logger.debug(f'QuadraticResidue.__roots_p_modulo_8_equals_5. self={self}. roots={roots}')
        return roots
    
    def __evaluate_s_t_p_modulo_8_equals_1(self) -> Tuple[int, int]:
        logger.debug(f'QuadraticResidue.__evaluate_s_t_p_modulo_8_equals_1. self={self}.')
        n = self.p - 1
        s = (n & -n).bit_length() - 1
        t = n >> s
        logger.debug(f'QuadraticResidue.__evaluate_s_t_p_modulo_8_equals_1. self={self}, s={s}, t={t}.')
        return s, t
    
    def __roots_a_exp_t_modulo_p_equals_1(self, t: int) -> Set[int]:
        logger.debug(f'QuadraticResidue.__roots_a_exp_t_modulo_p_equals_1. self={self}.')
        c = pow(self.a, (t+1) >> 1, self.p)
        roots = {c, self.p-c}
        logger.debug(f'QuadraticResidue.__roots_a_exp_t_modulo_p_equals_1. self={self}, roots={roots}.')
        return roots
    
    def __roots_a_exp_t_modulo_p_does_not_equal_1(self, s: int, t: int) -> Set[int]:
        logger.debug(f'QuadraticResidue.__roots_a_exp_t_modulo_p_does_not_equal_1. self={self}, s={s}, t={t}.')
        for n in range(2, self.p):
            if Jacobi_symbol(n,self.p) == - 1: break
        logger.debug(f'QuadraticResidue.__roots_a_exp_t_modulo_p_does_not_equal_1. self={self}, s={s}, t={t}, n={n}.')
        b = pow(n,t,self.p)
        inv_a = pow(self.a,-1,self.p)
        c = pow(self.a, (t+1) >> 1, self.p)
        exp_2 = 1 << (s - 2)
        logger.debug(f'QuadraticResidue.__roots_a_exp_t_modulo_p_does_not_equal_1. self={self}, s={s}, t={t}, b={b}, c={c}, exp_2={exp_2}.')
        for k in range(s-1):
            j = 0 if pow(pow(c,2,self.p) * inv_a, exp_2, self.p) == 1 else 1
            c *= pow(b, j<<k, self.p)
            exp_2 >>= 1
            logger.debug(f'QuadraticResidue.__roots_a_exp_t_modulo_p_does_not_equal_1. self={self}, s={s}, t={t}, k={k}, j={j}, c = {c}, exp_2={exp_2}.')
        c%=self.p
        roots = {c, self.p-c}
        logger.debug(f'QuadraticResidue.__roots_a_exp_t_modulo_p_does_not_equal_1. self={self}, roots={roots}.')
        return roots
    
    def __roots_p_modulo_8_equals_1(self) -> Set[int]:
        logger.debug(f'QuadraticResidue.__roots_p_modulo_8_equals_1. self={self}.')
        s, t = self.__evaluate_s_t_p_modulo_8_equals_1()
        if pow(self.a, t, self.p) == 1: return self.__roots_a_exp_t_modulo_p_equals_1(t)
        return self.__roots_a_exp_t_modulo_p_does_not_equal_1(s, t)
    
    def __roots_p_modulo_4_equals_1(self) -> Set[int]:
        logger.debug(f'QuadraticResidue.__roots_p_modulo_4_equals_1. self={self}.')
        if self.p&7 == 5: return self.__roots_p_modulo_8_equals_5()
        else: return self.__roots_p_modulo_8_equals_1()
    
    def __roots_Jacobi_symbol_equals_1(self) -> Set[int]:
        logger.debug(f'QuadraticResidue.__roots_Jacobi_symbol_equals_1. self={self}.')
        if self.p&3 == 3: return self.__roots_p_modulo_4_equals_3()
        else: return self.__roots_p_modulo_4_equals_1()
    
    def roots(self) -> Set[int]:
        logger.debug(f'QuadraticResidue.roots. self={self}.')
        if self.p == 2: return {self.a}
        if self.a % self.p == 1: return {1, self.p - 1}
        match Jacobi_symbol(self.a,self.p):
            case -1: return set()
            case 0: return {0}
        return self.__roots_Jacobi_symbol_equals_1()
    
def quadratic_comparison_to_QuadraticResidue(a: int, b: int, c: int, p: int) -> Tuple[QuadraticResidue, int]:
    """
    Получает из сравнения вида ax^2 + bx + c = 0 (mod p) экземпляр qr квадратного вычета вида
    y = x + b1 (mod qr.p)
    y^2 = qr.a (mod qr.p)
    Возвращает полученный экземпляр и b1
    """
    logger.debug(f'quadratic_comparison_to_QuadraticResidue. a={a}, b={b}, c={c}, p={p}.')
    inv_a = pow(a,-1,p)
    b1 = (inv_a * pow(2,-1, p) * b) % p
    qr = QuadraticResidue((pow(b1, 2, p) - inv_a * c) % p, p)
    logger.debug(f'quadratic_comparison_to_QuadraticResidue. a={a}, b={b}, c={c}, p={p}, qr={qr}, b1 = {b1}.')
    return qr, b1

def solve_quadratic_comparison(a: int, b: int, c: int, p: int) -> Set[int]:
    logger.debug(f'solve_quadratic_comparison. a={a}, b={b}, c={c}, p={p}.')
    qr, b1 = quadratic_comparison_to_QuadraticResidue(a, b, c, p)
    roots = set((y - b1) % p for y in qr.roots())
    logger.debug(f'solve_quadratic_comparison. a={a}, b={b}, c={c}, p={p} roots={roots}.')
    return roots


class TestQuadraticResidue(unittest.TestCase):
    
    def setUp(self):
        return super().setUp()
    
    def test_Jacobi_symbol(self):
        self.assertEqual(Jacobi_symbol(30, 25), 0, "Если a и m не взаимнопросты, то символ Якоби должен быть равен нулю.")
        self.assertEqual(Jacobi_symbol(5, 1), 1, "Символ Якоби (a,1) должен быть равен 1 для любого a.")
        self.assertEqual(Jacobi_symbol(-1, 7), -1, "Ошибка определения символа Якоби в случае (-1/m)=-1")
        self.assertEqual(Jacobi_symbol(-1, 13), 1, "Ошибка определения символа Якоби в случае (-1/m)=1")
        self.assertEqual(Jacobi_symbol(1, 5), 1, "Ошибка определения символа Якоби в случае (1/m)=1")
        self.assertEqual(Jacobi_symbol(4, 21), 1, "Ошибка определения символа Якоби в случае (2^(2k)/m)=1 при m % 8 = 5")
        self.assertEqual(Jacobi_symbol(4, 19), 1, "Ошибка определения символа Якоби в случае (2^(2k)/m)=1 при m % 8 = 3")
        self.assertEqual(Jacobi_symbol(4, 9), 1, "Ошибка определения символа Якоби в случае (2^(2k)/m)=1 при m % 8 = 1")
        self.assertEqual(Jacobi_symbol(4, 7), 1, "Ошибка определения символа Якоби в случае (2^(2k)/m)=1 при m % 8 = 7")
        self.assertEqual(Jacobi_symbol(8, 21), -1, "Ошибка определения символа Якоби в случае (2^(2k+1)/m)=1 при m % 8 = 5")
        self.assertEqual(Jacobi_symbol(8, 19), -1, "Ошибка определения символа Якоби в случае (2^(2k+1)/m)=1 при m % 8 = 3")
        self.assertEqual(Jacobi_symbol(8, 9), 1, "Ошибка определения символа Якоби в случае (2^(2k+1)/m)=1 при m % 8 = 1")
        self.assertEqual(Jacobi_symbol(8, 7), 1, "Ошибка определения символа Якоби в случае (2^(2k+1)/m)=1 при m % 8 = 7")
        self.assertEqual(Jacobi_symbol(3, 5), -1, "Ошибка определения символа Якоби в случае (-1)^(q-1)(m-1)//4=1")
        self.assertEqual(Jacobi_symbol(3, 7), -1, "Ошибка определения символа Якоби в случае (-1)^(q-1)(m-1)//4=-1")
        self.assertEqual(Jacobi_symbol(3, 7), -1, "Ошибка определения символа Якоби в случае (-1)^(q-1)(m-1)//4=-1")
        self.assertEqual(Jacobi_symbol(219, 233), 1, "Ошибка определения символа Якоби в общем случае при простом m.")
        self.assertEqual(Jacobi_symbol(219, 493), 1, "Ошибка определения символа Якоби в общем случае.")

    def test_quadratic_comparison_to_QuadraticResidue(self):
        self.assertEqual(quadratic_comparison_to_QuadraticResidue(7, 0, -8, 5), (QuadraticResidue(4,5), 0), "Ошибка quadratic_comparison_to_QuadraticResidue в случае b = 0.")
        self.assertEqual(quadratic_comparison_to_QuadraticResidue(6, 9, 8, 5), (QuadraticResidue(1,5), 2), "Ошибка quadratic_comparison_to_QuadraticResidue в случае b != 0.")

    def test_roots(self):
        self.assertEqual(QuadraticResidue(1, 2).roots(), {1}, "При p=2 корень квадратичного вычета должны быть равен a.")
        self.assertEqual(QuadraticResidue(4, 3).roots(), {1,2}, "При a = 1 (mod p) корень квадратичного вычета должны быть равен 1, -1.")
        self.assertEqual(QuadraticResidue(-1, 7).roots(), set(), "При символе Якоби равном -1, кореней нет.")
        self.assertEqual(QuadraticResidue(49, 7).roots(), {0}, "При символе Якоби равном 0, корень 0.")
        self.assertEqual(QuadraticResidue(4, 7).roots(), {2,5}, "Ошибка вычисления корней квадратичного вычета при p mod 4 = 3.")
        self.assertEqual(QuadraticResidue(4, 5).roots(), {2,3}, "Ошибка вычисления корней квадратичного вычета при p mod 8 = 5.")
        self.assertEqual(QuadraticResidue(10, 41).roots(), {16,25}, "Ошибка вычисления корней квадратичного вычета при p a^t % p = 1.")
        self.assertEqual(QuadraticResidue(9, 41).roots(), {3, 38}, "Ошибка вычисления корней квадратичного вычета при p a^t % p != 1.")
        self.assertEqual(QuadraticResidue(186, 401).roots(), {304, 97}, "Ошибка вычисления корней квадратичного вычета при p a^t % p != 1.")

    def test_solve_quadratic_comparison(self):
        self.assertEqual(solve_quadratic_comparison(2, -20, 32, 41), {2, 8}, "Ошибка вычисления корней квадратичного сравнения.")

if __name__ == '__main__':
    unittest.main()