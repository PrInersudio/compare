"""Модуль для поиска корней квадратичных сравнений."""

from typing import Tuple, Set, Literal
import logging
from field import GF


logger = logging.getLogger(__name__)
logger.setLevel(logging.DEBUG)
handler = logging.FileHandler(f'log/{__name__}.log')
handler.setLevel(logging.DEBUG)
handler.setFormatter(logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s'))
logger.addHandler(handler)


def jacobi_symbol(a: int, m: int) -> Literal[-1, 0, 1]:
    """Вычисление символа Якоби."""
    logger.debug('jacobi_symbol. a=%s, m=%s.', a, m)
    result = 1
    a %= m
    while a:
        logger.debug('Шаг jacobi_symbol. a=%s, m=%s, result=%s.', a, m, result)
        if not a&1:
            result *= 1 if m&7 in [1,7] else -1
            a>>=1
            continue
        result *= 1 if m&3 == 1 or a&3 == 1 else -1
        a,m = m%a,a
    if m != 1:
        result = 0
    logger.debug('jacobi_symbol. a=%s, m=%s, result=%s.', a, m, result)
    return result


class QuadraticResidue:
    """Класс квадратичных вычетов."""

    def __init__(self, a: int, p: int):
        self.a = a
        self.p = p

    def __repr__(self) -> str:
        return f'QuadraticResidue({self.a}, {self.p})'

    def __str__(self) -> str:
        return f'x^2 = {self.a} (mod {self.p})'

    def __eq__(self, other: 'QuadraticResidue') -> bool:
        if not isinstance(other, QuadraticResidue):
            return False
        return self.p == other.p and self.a % self.p == other.a % other.p

    def __roots_p_modulo_4_equals_3(self) -> Set[int]:
        logger.debug('QuadraticResidue.__roots_p_modulo_4_equals_3. self=%s.', self)
        c  = pow(self.a, (self.p + 1) >> 2 , self.p)
        roots = {c, self.p-c}
        logger.debug('__roots_quadratic_residue_p_modulo_4_equals_3. self=%s, roots=%s',
                     self, roots)
        return roots

    def __roots_p_modulo_8_equals_5(self) -> Set[int]:
        logger.debug('QuadraticResidue.__roots_p_modulo_8_equals_5. self=%s.', self)
        k = self.p >> 3
        t = 0 if pow(self.a, (self.p-1) >> 2, self.p) == 1 else 1
        c = (pow(self.a, k + 1, self.p) * pow(2, ((k<<1)+1) * t, self.p)) % self.p
        logger.debug('QuadraticResidue.__roots_p_modulo_8_equals_5.'
                     'self=%s, k=%s, t=%s, c=%s.', self, k, t, c)
        roots = {c, self.p-c}
        logger.debug('QuadraticResidue.__roots_p_modulo_8_equals_5.'
                     'self=%s. roots=%s', self, roots)
        return roots

    def __evaluate_s_t_p_modulo_8_equals_1(self) -> Tuple[int, int]:
        logger.debug('QuadraticResidue.__evaluate_s_t_p_modulo_8_equals_1. self=%s.', self)
        n = self.p - 1
        s = (n & -n).bit_length() - 1
        t = n >> s
        logger.debug('QuadraticResidue.__evaluate_s_t_p_modulo_8_equals_1.'
                     'self=%s, s=%s, t=%s.', self, s, t)
        return s, t

    def __roots_a_exp_t_modulo_p_equals_1(self, t: int) -> Set[int]:
        logger.debug('QuadraticResidue.__roots_a_exp_t_modulo_p_equals_1. self=%s.', self)
        c = pow(self.a, (t+1) >> 1, self.p)
        roots = {c, self.p-c}
        logger.debug('QuadraticResidue.__roots_a_exp_t_modulo_p_equals_1.'
                     'self=%s, roots=%s.', self, roots)
        return roots

    def __roots_a_exp_t_modulo_p_does_not_equal_1(self, s: int, t: int) -> Set[int]:
        logger.debug('QuadraticResidue.__roots_a_exp_t_modulo_p_does_not_equal_1.'
                     'self=%s, s=%s, t=%s.', self, s, t)
        for n in range(2, self.p):
            if jacobi_symbol(n,self.p) == - 1:
                break
        logger.debug('QuadraticResidue.__roots_a_exp_t_modulo_p_does_not_equal_1.'
                     'self=%s, s=%s, t=%s, n=%s.', self, s, t, n)
        b = pow(n,t,self.p)
        inv_a = pow(self.a,-1,self.p)
        c = pow(self.a, (t+1) >> 1, self.p)
        exp_2 = 1 << (s - 2)
        logger.debug('QuadraticResidue.__roots_a_exp_t_modulo_p_does_not_equal_1.'
                     'self=%s, s=%s, t=%s, b=%s, c=%s, exp_2=%s.',
                     self, s, t, b, c, exp_2)
        for k in range(s-1):
            j = 0 if pow(pow(c,2,self.p) * inv_a, exp_2, self.p) == 1 else 1
            c *= pow(b, j<<k, self.p)
            exp_2 >>= 1
            logger.debug('QuadraticResidue.__roots_a_exp_t_modulo_p_does_not_equal_1.'
                         'self=%s, s=%s, t=%s, k=%s, j=%s, c=%s, exp_2=%s.',
                         self, s, t, k, j, c, exp_2)
        c%=self.p
        roots = {c, self.p-c}
        logger.debug('QuadraticResidue.__roots_a_exp_t_modulo_p_does_not_equal_1.'
                     'self=%s, roots=%s.', self, roots)
        return roots

    def __roots_p_modulo_8_equals_1(self) -> Set[int]:
        logger.debug('QuadraticResidue.__roots_p_modulo_8_equals_1. self=%s.', self)
        s, t = self.__evaluate_s_t_p_modulo_8_equals_1()
        if pow(self.a, t, self.p) == 1:
            return self.__roots_a_exp_t_modulo_p_equals_1(t)
        return self.__roots_a_exp_t_modulo_p_does_not_equal_1(s, t)

    def __roots_p_modulo_4_equals_1(self) -> Set[int]:
        logger.debug('QuadraticResidue.__roots_p_modulo_4_equals_1. self=%s.', self)
        if self.p&7 == 5:
            return self.__roots_p_modulo_8_equals_5()
        else: return self.__roots_p_modulo_8_equals_1()

    def __roots_jacobi_symbol_equals_1(self) -> Set[int]:
        logger.debug('QuadraticResidue.__roots_jacobi_symbol_equals_1. self=%s.', self)
        if self.p&3 == 3:
            return self.__roots_p_modulo_4_equals_3()
        else: return self.__roots_p_modulo_4_equals_1()

    def roots(self) -> Set[int]:
        """Корни квадратичного вычета."""
        logger.debug('QuadraticResidue.roots. self=%s.', self)
        if self.p == 2:
            return {self.a}
        if self.a % self.p == 1:
            return {1, self.p - 1}
        match jacobi_symbol(self.a,self.p):
            case -1: return set()
            case 0: return {0}
        return self.__roots_jacobi_symbol_equals_1()


def quadratic_comparison_to_quadratic_residue(
        a: int, b: int, c: int, p: int,
) -> Tuple[QuadraticResidue, int]:
    """
    Получает из сравнения вида ax^2 + bx + c = 0 (mod p) экземпляр qr квадратного вычета вида
    y = x + b1 (mod qr.p)
    y^2 = qr.a (mod qr.p)
    Возвращает полученный экземпляр и b1
    """
    logger.debug('quadratic_comparison_to_QuadraticResidue. a=%s, b=%s, c=%s, p=%s.',
                 a, b, c, p)
    inv_a = pow(a,-1,p)
    b1 = (inv_a * pow(2,-1, p) * b) % p
    qr = QuadraticResidue((pow(b1, 2, p) - inv_a * c) % p, p)
    logger.debug('quadratic_comparison_to_QuadraticResidue.'
                 'a=%s, b=%s, c=%s, p=%s, qr=%s, b1=%s.',
                 a, b, c, p, qr, b1)
    return qr, b1


def solve_quadratic_comparison(a: int, b: int, c: int, p: int) -> Set[GF.Element]:
    """Поиск корней квадратичного сравнения."""
    logger.debug('solve_quadratic_comparison. a=%s, b=%s, c=%s, p=%s.',
                 a, b, c, p)
    qr, b1 = quadratic_comparison_to_quadratic_residue(a, b, c, p)
    gf = GF(p)
    roots = set(gf(y - b1) for y in qr.roots())
    logger.debug('solve_quadratic_comparison. a=%s, b=%s, c=%s, p=%s, roots=%s.',
                 a, b, c, p, roots)
    return roots
