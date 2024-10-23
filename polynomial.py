from __future__ import annotations
from math import comb
from sympy import factorint
from itertools import product
from typing import Tuple, Dict, Set, Generator
import logging
import re
from linear_comparisons import *
import unittest
from quadratic_comparisons import solve_quadratic_comparison

logger = logging.getLogger(__name__)
logger.setLevel(logging.DEBUG)
handler = logging.FileHandler(f'log/{__name__}.log')
handler.setLevel(logging.DEBUG)
handler.setFormatter(logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s'))
logger.addHandler(handler)

class Polynomial:
    """
    Полином с целыми коэффициентами.
    """

    def __init__(self, raw_poly: None | Dict[int,int] | str | int = None):
        logger.debug(f'Polynomial.__init__. raw_poly={raw_poly}')
        if not raw_poly:
            self.__coefficients = dict()
            return
        elif isinstance(raw_poly, int):
            self.__coefficients = {0:raw_poly} if raw_poly != 0 else dict()
            return
        elif isinstance(raw_poly, dict):
            self.__coefficients = {exp:coeff for exp, coeff in raw_poly.items() if isinstance(exp, int) and isinstance(coeff,int) and exp >= 0 and coeff != 0}
            return
        elif not isinstance(raw_poly, str):
            logger.critical(f'Polynomial.__init__. raw_poly={raw_poly}. Неправильный тип raw_poly.')
            raise TypeError("Полином создаётся только из строки или словаря")
        pattern = r'([+-]?[^-+]+)'
        monomials = re.findall(pattern, raw_poly)
        self.__coefficients = dict()
        logger.debug(f'Polynomial.__init__. self={self}, raw_poly={raw_poly}, monomials={monomials}')
        for monomial in monomials:
            splited = monomial.replace(' ', '').split('x')
            coeff = 1 if splited[0] in ['', '+'] else -1 if splited[0] == '-' else int(splited[0])
            if coeff == 0: continue
            exp = 0 if len(splited) != 2 else (1 if splited[1] == '' else int(splited[1][1:]))
            if exp < 0:
                logger.critical("Отрицательная степень.")
                raise ValueError("Отрицательные степени не поддерживаются.")
            self[exp] += coeff
        logger.debug(f'Polynomial.__init__ завершается. self={self}, raw_poly={raw_poly}')

    def __getitem__(self, exp: int) -> int:
        return self.__coefficients.get(exp, 0)
    
    def __setitem__(self, exp: int, coeff: int) -> None:
        if coeff == 0: self.__coefficients.pop(exp,0)
        else: self.__coefficients[exp] = coeff

    def __contains__(self, exp: int) -> bool:
        return exp in self.__coefficients
    
    def __iter__(self) -> Generator[int]:
        for exp, coeff in self.__coefficients.items():
            if coeff != 0: yield exp
    
    def __eq__(self, other: 'Polynomial') -> bool:
        if not isinstance(other, Polynomial): return False
        return all(self[exp] == other[exp] for exp in set(self) | set(other))
    
    def monomials(self)  -> Generator[Tuple[int,int]]:
        """
        Генератор наборов (exp, coeff).
        """
        for exp, coeff in self.__coefficients.items():
            if coeff != 0: yield exp, coeff

    def __hash__(self) -> int:
        return hash(tuple(sorted((exp, coeff) for exp, coeff in self.monomials())))

    def copy(self) -> 'Polynomial':
        return Polynomial(self.__coefficients.copy())

    def __add__(self, other: 'Polynomial' | int) -> 'Polynomial':
        logger.debug(f'Polynomial.__add__. self={self}, other={other}')
        if isinstance(other, int):
            other = Polynomial(other)
        if not isinstance(other, Polynomial):
            logger.critical(f'Polynomial.__add__. self={self}, other={other}. Неправильный тип other.')
            raise TypeError("Сложение возмозжно только с числом и другим полиномом")
        coeffs = dict()
        for i in set(self) | set(other):
            new_coeff = self[i] + other[i]
            if new_coeff != 0: coeffs[i] = new_coeff
        logger.debug(f'Polynomial.__add__. self={self}, other={other}, coeffs={coeffs}')
        return Polynomial(coeffs)
    
    def __iadd__(self, other: 'Polynomial' | int) -> 'Polynomial':
        logger.debug(f'Polynomial.__iadd__. self={self}, other={other}')
        return self + other
    
    def __radd__(self, other: 'Polynomial' | int) -> 'Polynomial':
        logger.debug(f'Polynomial.__radd__. self={self}, other={other}')
        return self + other
    
    def __neg__(self):
        return Polynomial({exp:-coeff for exp,coeff in self.monomials()})
    
    def __sub__(self, other: 'Polynomial' | int) -> 'Polynomial':
        logger.debug(f'Polynomial.__sub__. self={self}, other={other}')
        if isinstance(other, int):
            other = Polynomial(other)
        if not isinstance(other, Polynomial):
            logger.critical(f'Polynomial.__sub__. self={self}, other={other}. Неправильный тип other.')
            raise TypeError("Вычитание возмозжно только с числом и другим полиномом")
        coeffs = dict()
        for i in set(self) | set(other):
            new_coeff = self[i] - other[i]
            if new_coeff != 0: coeffs[i] = new_coeff
        logger.debug(f'Polynomial.__sub__. self={self}, other={other}, coeffs={coeffs}')
        return Polynomial(coeffs)
    
    def __isub__(self, other: 'Polynomial' | int) -> 'Polynomial':
        logger.debug(f'Polynomial.__isub__. self={self}, other={other}')
        return self - other
    
    def __rsub__(self, other: 'Polynomial' | int) -> 'Polynomial':
        logger.debug(f'Polynomial.__rsub__. self={self}, other={other}')
        return -self + other
    
    def __mul__(self, other: 'Polynomial' | int) -> 'Polynomial':
        logger.debug(f'Polynomial.__mul__. self={self}, other={other}')
        if isinstance(other, int):
            other = Polynomial(other)
        if not isinstance(other, Polynomial):
            logger.critical(f'Polynomial.__mul__. self={self}, other={other}. Неправильный тип other.')
            raise TypeError("Умножение многочлена возмозжно только на число или другой многочлен.")
        result = Polynomial()
        for exp1, coeff1 in self.monomials():
            for exp2, coeff2 in other.monomials():
                result[exp1+exp2] += coeff1 * coeff2
        result = Polynomial({exp:coeff for exp, coeff in result.monomials() if coeff != 0})
        logger.debug(f'Polynomial.__mul__. self={self}, other={other}, result={result}')
        return result
    
    def __imul__(self, other: 'Polynomial' | int) -> 'Polynomial':
        logger.debug(f'Polynomial.__imul__. self={self}, other={other}')
        return self * other
    
    def __rmul__(self, other: 'Polynomial' | int) -> 'Polynomial':
        logger.debug(f'Polynomial.__rmul__. self={self}, other={other}')
        return self * other
    
    def __bool__(self) -> bool:
        return any(coeff != 0 for coeff in self.__coefficients.values())
    
    def __repr__(self) -> str:
        if not self:
            return "0"
        monomials = []
        for exp in sorted(list(self), reverse=True):
            coeff = self[exp]
            if coeff == 0: continue
            if exp == 0: monomial = f"{coeff}"
            elif exp == 1: monomial = "x" if coeff == 1 else "-x" if coeff == -1 else f"{coeff}x"
            else: monomial = f"x^{exp}" if coeff == 1 else f"-x^{exp}" if coeff == -1 else f"{coeff}x^{exp}"
            monomials.append(monomial)
        polynomial_str = " + ".join(monomials)
        polynomial_str = polynomial_str.replace("+ -", "- ")
        return polynomial_str
    
    def __str__(self) -> str:
        return self.__repr__()
    
    def degree(self) -> int:
        logger.debug(f'Polynomial.degree. self={self}')
        if not self:
            logger.debug(f'Polynomial.degree. self={self}. Нулевой многочлен')
            return -1
        degree = max(exp for exp in self if self[exp] != 0)
        logger.debug(f'Polynomial.degree. self={self}, degree={degree}')
        return degree
    
    def __mod__(self,module: int) -> 'Polynomial':
        """
        Создаёт новый многочлен, у которого коэффициенты - остатки от деления на число.
        """
        logger.debug(f'Polynomial.__mod__. self={self}')
        poly = Polynomial({exp:coeff%module for exp,coeff in self.monomials() if coeff%module != 0})
        logger.debug(f'Polynomial.__mod__. self={self}, poly={poly}')
        return poly
    
    def value(self, x: int) -> int:
        """
        Значение многочлена в точке x
        """
        logger.debug(f'Polynomial.value. self={self}')
        value = sum(coeff * x ** exp for exp,coeff in self.monomials())
        logger.debug(f'Polynomial.value. self={self}, value={value}')
        return value
    
    def Fermat_little_theorem(self, p: int) -> 'Polynomial':
        """
        Сокращение многочлена по Малой теореме Ферма в GF(p)
        x^(p) = x (mod p) <==> x^(kp+t) = x^(k+t) (mod p)
        """
        logger.debug(f'Polynomial.Fermat_little_theorem. self={self}, p={p}')
        poly = Polynomial()
        for exp,coeff in self.monomials():
            while exp >= p: exp = exp // p + exp % p
            poly[exp] += coeff
        logger.debug(f'Polynomial.Fermat_little_theorem. self={self}, p={p}, poly={poly}')
        return poly % p
    
    def div(self, divisor: 'Polynomial' | int, p: int) -> Tuple['Polynomial', 'Polynomial']:
        """
        Деление многочленов в GF(p)[x]
        """
        logger.debug(f'Polynomial.div. self={self}, divisor={divisor}, p={p}')
        if isinstance(divisor, int):
            divisor = Polynomial(divisor)
        if not isinstance(divisor, Polynomial):
            logger.critical(f'Polynomial.div. self={self}, divisor={divisor}, p={p}. Неправильный тип divisor.')
            raise TypeError("Делитель должен быть экземпляром Polynomial или int.")
        divisor = divisor % p
        if not divisor:
            logger.critical(f'Polynomial.div. self={self}, divisor={divisor}, p={p}. Нулевой divisor.')
            raise ZeroDivisionError("Деление на нулевой многочлен.")
        remainder = self % p
        remainder_degree = remainder.degree()
        quotient = Polynomial()
        divisor_degree = divisor.degree()
        while remainder and remainder_degree >= divisor_degree:
            logger.debug(f'Polynomial.div. self={self}, divisor={divisor}, p={p}, quotient={quotient}, remainder={remainder}')
            degree_diff = remainder_degree - divisor_degree
            leading_coeff = solve_linear_comparison(divisor[divisor_degree], remainder[remainder_degree], p)[0][0]
            quotient[degree_diff] = leading_coeff
            for exp, coeff in divisor.monomials():
                remainder[exp + degree_diff] = (remainder[exp + degree_diff] - leading_coeff * coeff) % p
            remainder = Polynomial({exp: coeff for exp, coeff in remainder.monomials() if coeff != 0})
            remainder_degree = remainder.degree()
        logger.debug(f'Polynomial.div. self={self}, divisor={divisor}, p={p}, quotient={quotient}, remainder={remainder}')
        return quotient, remainder
    
    def __roots_BerlekampRabin(self, p: int) -> Set[int]:
        logger.debug(f'Polynomial.__roots_BerlekampRabin. self={self}, p={p}')
        roots = set()
        nonlinears = [self]
        for delta in range(p):
            if not nonlinears: break
            poly = nonlinears.pop()
            logger.debug(f'Polynomial.__roots_BerlekampRabin. self={self}, delta={delta}, roots={roots}, nonlinears = {nonlinears}, poly={poly}')
            d = gcd(poly, binomial_theorem(1,delta,(p-1) >> 1) - 1, p)
            if d.degree() <= 0: nonlinears.append(poly); continue
            for factor in [d, poly.div(d,p)[0]]:
                if factor.degree() == 1: roots |= {x for x,_ in solve_linear_comparison(factor[1], -factor[0], p)}
                else: nonlinears.append(factor)
        logger.debug(f'Polynomial.__roots_BerlekampRabin. self={self}, p={p}, roots={roots}')
        return roots
    
    def derivative(self) -> 'Polynomial':
        logger.debug(f'Polynomial.derivative. self={self}')
        derivative = Polynomial({exp-1:coeff*exp for exp,coeff in self.monomials() if exp > 0})
        logger.debug(f'Polynomial.derivative. self={self}, derivative={derivative}')
        return derivative

    def __roots_prime_module(self, p: int) -> Set[int]:
        logger.debug(f'Polynomial.__roots_prime_module. self={self}, p={p}')
        fermat = self.Fermat_little_theorem(p)
        if fermat.degree() <= 1:
            logger.debug(f'Polynomial.__roots_prime_module. self={self}, p={p} fermat={fermat}. Степень <= 1.')
            return set(a for a,_ in solve_linear_comparison(fermat[1], -fermat[0], p))
        reduced = gcd(fermat, Polynomial({p:1,1:-1}),p)
        if reduced.degree() == 2:
            logger.debug(f'Polynomial.__roots_prime_module. self={self}, p={p} reduced={reduced}. Степень 2.')
            return solve_quadratic_comparison(reduced[2], reduced[1], reduced[0], p)
        return reduced.__roots_BerlekampRabin(p)
    
    def __roots_primary_module(self, p: int, q: int = 1) -> Set[Tuple[int, int]]:
        logger.debug(f'Polynomial.__roots_primary_module. self={self}, p={p}, q={q}')
        result = self.__roots_prime_module(p)
        derivative = self.derivative().Fermat_little_theorem(p)
        for i in range(1,q):
            p_i = p**i
            result = {a:set(t for t,_ in solve_linear_comparison(derivative.value(a), - self.value(a) // p_i, p)) for a in result}
            logger.debug(f'Polynomial.__roots_primary_module. self={self}, p={p}, q={q}, result={result}')
            result = set(a + t * p_i for a, temp in result.items() for t  in temp)
        p_q = p**q
        logger.debug(f'Polynomial.__roots_primary_module. self={self}, p={p}, q={q}, result={result}')
        return set((a % p_q, p_q) for a in result)
    
    def roots(self, M: int) -> Set[Tuple[int, int]]:
        logger.debug(f'Polynomial.roots. self={self}, M={M}')
        if self.degree() <= 1:
            logger.debug(f'Polynomial.roots. self={self}, M={M}. Степень <= 1.')
            return set(solve_linear_comparison(self[1],-self[0],M))
        primary_results = []
        for p,q in factorint(M).items():
            primary_results.append(self.__roots_primary_module(p,q))
        if len(primary_results) == 1:
            logger.debug(f'Polynomial.roots. self={self}, M={M}, primary_results={primary_results}. M - примарное.')
            return primary_results[0]
        result = set()
        for element in product(*primary_results):
            result.add(crt(element))
        logger.debug(f'Polynomial.roots. self={self}, M={M}, result={result}')
        return result
    
    def monic(self, p) -> 'Polynomial':
        if not self: return Polynomial()
        return (self * pow(self[self.degree()], -1, p)) % p
    
def gcd(f: Polynomial, g: Polynomial, p:int) -> Polynomial:
    logger.debug(f'Начало gcd. f={f}, g={g}')
    while g:
        logger.debug(f'Шаг gcd. f={f}, g={g}')
        f, (_, g) = g, f.div(g,p)
    logger.debug(f'gcd, делаем унитарный. f={f}, g={g}')
    f = f.monic(p)
    logger.debug(f'Конец gcd. f={f}, g={g}')
    return f

def binomial_theorem(m: int, a: int, exp: int) -> Polynomial:
    """
    Разложение (mx+a)^exp
    """
    logger.debug(f'binomial_theorem. m={m}, a={a}, exp={exp}')
    poly = Polynomial()
    for k in range(exp + 1):
        coeff = comb(exp, k) * (m ** (exp - k)) * (a ** k)
        poly[exp - k] = poly[exp - k, 0] + coeff
        logger.debug(f'binomial_theorem. m={m}, a={a}, exp={exp}, poly={poly}')
    return poly

class TestPolynomial(unittest.TestCase):

    def setUp(self):
        return super().setUp()
    
    def test_init(self):
        # константы
        self.assertEqual(Polynomial("0"), Polynomial(), "Ошибка инициализации нулевого полинома.")
        self.assertEqual(Polynomial(""), Polynomial(), "Ошибка инициализации нулевого полинома (пустая строка).")
        self.assertEqual(Polynomial("5"), Polynomial(5), "Ошибка инициализации констатны.")
        self.assertEqual(Polynomial("+5"), Polynomial(5), "Ошибка инициализации положительной константы (плюс без пробела).")
        self.assertEqual(Polynomial("+ 5"), Polynomial(5), "Ошибка инициализации положительной константы (плюс с пробелом).")
        self.assertEqual(Polynomial("-5"), Polynomial(-5), "Ошибка инициализации отрицательной константы (без пробела).")
        self.assertEqual(Polynomial("- 5"), Polynomial(-5), "Ошибка инициализации отрицательной константы (с пробелом).")
        # линейный моном
        self.assertEqual(Polynomial("0x"), Polynomial(), "Ошибка инициализации линейного монома с нулевым коэффициентом.")
        self.assertEqual(Polynomial("x"), Polynomial({1:1}), "Ошибка инициализации линейного монома с коэффициентом, равным 1.")
        self.assertEqual(Polynomial("+x"), Polynomial({1:1}), "Ошибка инициализации линейного монома с коэффициентом, равным 1 (плюс без пробела).")
        self.assertEqual(Polynomial("+ x"), Polynomial({1:1}), "Ошибка инициализации линейного монома с коэффициентом, равным 1 (плюс с пробелом).")
        self.assertEqual(Polynomial("-x"), Polynomial({1:-1}), "Ошибка инициализации линейного монома с коэффициентом, равным -1 (без пробела).")
        self.assertEqual(Polynomial("- x"), Polynomial({1:-1}), "Ошибка инициализации линейного монома с коэффициентом, равным -1 (с пробелом).")
        self.assertEqual(Polynomial("5x"), Polynomial({1:5}), "Ошибка инициализации линейного монома с положительным коэффициентом.")
        self.assertEqual(Polynomial("+5x"), Polynomial({1:5}), "Ошибка инициализации линейного монома с положительным коэффициентом (плюс без пробела).")
        self.assertEqual(Polynomial("+ 5x"), Polynomial({1:5}), "Ошибка инициализации линейного монома с положительным коэффициентом (плюс с пробелом).")
        self.assertEqual(Polynomial("-5x"), Polynomial({1:-5}), "Ошибка инициализации линейного монома с отрицательным коэффициентом (без пробела).")
        self.assertEqual(Polynomial("- 5x"), Polynomial({1:-5}), "Ошибка инициализации линейного монома с отрицательным коэффициентом (с пробелом).")
        # произвольный моном
        self.assertEqual(Polynomial("0x^2"), Polynomial(), "Ошибка инициализации произвольного монома с нулевым коэффициентом.")
        self.assertEqual(Polynomial("5x^2"), Polynomial({2:5}), "Ошибка инициализации произвольного монома с положительным коэффициентом.")
        self.assertEqual(Polynomial("+5x^2"), Polynomial({2:5}), "Ошибка инициализации произвольного монома с положительным коэффициентом (плюс без пробела).")
        self.assertEqual(Polynomial("+ 5x^2"), Polynomial({2:5}), "Ошибка инициализации произвольного монома с положительным коэффициентом (плюс с пробелом).")
        self.assertEqual(Polynomial("-5x^2"), Polynomial({2:-5}), "Ошибка инициализации произвольного монома с отрицательным коэффициентом (без пробела).")
        self.assertEqual(Polynomial("- 5x^2"), Polynomial({2:-5}), "Ошибка инициализации произвольного монома с отрицательным коэффициентом (с пробелом).")
        # полиномы
        self.assertEqual(Polynomial("3x^3 - 4x^2 + 2x - 5"), Polynomial({0:-5, 1:2, 2:-4, 3:3}), "Ошибка инициализации в общем случае.")
        self.assertEqual(Polynomial("3x^3-4x^2+2x-5"), Polynomial({0:-5, 1:2, 2:-4, 3:3}), "Ошибка инициализации в общем случае (безпробельное написание).")
        self.assertEqual(Polynomial("2x^2 - 5"), Polynomial({0:-5, 2:2}), "Ошибка инициализации при нулевом коэффициенте.")
        self.assertEqual(Polynomial("2x^2 + 2x"), Polynomial({1:2, 2:2}), "Ошибка инициализации при нулевом свободном члене.")
        self.assertEqual(Polynomial("3x^2 + 5x^2 - 2x"), Polynomial({1:-2, 2:8}), "Ошибка инициализации в случае с одинаковми степенями.")
        # словарь
        self.assertEqual(Polynomial({1.4:1.4}), Polynomial(), "Ошибка проверки типов при инициализации словарём.")
        self.assertEqual(Polynomial({-5:5}), Polynomial(), "Ошибка проверки положительности степеней при инициализации словарём.")
        self.assertEqual(Polynomial({5:0}), Polynomial(), "Ошибка проверки коэффициентов при инициализации словарём.")
        # исключения
        with self.assertRaises(TypeError, msg="Ошибка проверки типа raw_poly при инициализации."): Polynomial(1.4)
        with self.assertRaises(ValueError, msg="init не поднял исключение при отрицательной степени в строке."): Polynomial("5x^-5")
    
    def test_str(self):
        self.assertEqual(str(Polynomial()), "0", "Ошибка преобразования в строку нулевого полинома.")
        self.assertEqual(str(Polynomial("5")), "5", "Ошибка преобразования в строку положительной констаты.")
        self.assertEqual(str(Polynomial("-5")), "-5", "Ошибка преобразования в строку отрицательной констаты.")
        self.assertEqual(str(Polynomial("5x")), "5x", "Ошибка преобразования в строку линейного монома с положительным коэффициентом.")
        self.assertEqual(str(Polynomial("-5x")), "-5x", "Ошибка преобразования в строку линейного монома с отрицательным коэффициентом.")
        self.assertEqual(str(Polynomial("5x^2")), "5x^2", "Ошибка преобразования в строку произвольного монома с положительным коэффициентом.")
        self.assertEqual(str(Polynomial("-5x^2")), "-5x^2", "Ошибка преобразования в строку произвольного монома с отрицательным коэффициентом.")
        self.assertEqual(str(Polynomial("3x^3 - 4x^2 + 2x - 5")), "3x^3 - 4x^2 + 2x - 5", "Ошибка преобразования в строку в общем случае.")
        self.assertEqual(str(Polynomial("2x^2 - 5")), "2x^2 - 5", "Ошибка преобразования в строку при нулевом коэффициенте.")
        self.assertEqual(str(Polynomial("2x^2 + 2x")), "2x^2 + 2x", "Ошибка преобразования в строку при нулевом свободном члене.")

    def test_getitem(self):
        self.assertEqual(Polynomial({2:5})[2], 5, "Ошибка при получении значения ненулевого коэффициента.")
        self.assertEqual(Polynomial()[2], 0, "Ошибка при получении значения нулевого коэффициента.")

    def test_setitem(self):
        p = Polynomial()
        p[2] = 5
        self.assertEqual(p, Polynomial({2:5}), "Ошибка при установке коэффициента.")
        p[2] = 2
        self.assertEqual(p, Polynomial({2:2}), "Ошибка при изменении коэффициента.")
        p[2] = 0
        self.assertEqual(p, Polynomial(), "Ошибка обнуления коэффициента.")
    
    def test_contains(self):
        p = Polynomial({2:5})
        self.assertTrue(2 in p, "Ошибка false negative при проверке содержания степени в полиноме.")
        self.assertFalse(3 in p, "Ошибка false positive при проверке содержания степени в полиноме.")

    def test_iter(self):
        self.assertEqual(list(Polynomial()), [], "Ошибка iter при нулевом многочлене.")
        self.assertEqual(list(Polynomial({0:1, 2:3, 4:5})), [0,2,4], "Ошибка iter при ненулевом многочлене.")
    
    def test_monomials(self):
        self.assertEqual(list(Polynomial().monomials()), [], "Ошибка monomials при нулевом многочлене.")
        self.assertEqual(list(Polynomial({0:1, 2:3, 4:5}).monomials()), [(0,1),(2,3),(4,5)], "Ошибка monomials при ненулевом многочлене.")

    def test_copy(self):
        p = Polynomial()
        self.assertEqual(p, p.copy(), "Ошибка copy при нулевом многочлене.")
        p = Polynomial({0:1, 2:3, 4:5})
        self.assertEqual(p, p.copy(), "Ошибка copy при ненулевом многочлене.")
    
    def test_add(self):
        self.assertEqual(Polynomial() + Polynomial(), Polynomial(), "Ошибка при суммировании нулевых многочленов.")
        self.assertEqual(Polynomial("3x^2 + 2x - 5") + Polynomial(), Polynomial("3x^2 + 2x - 5"), "Ошибка при суммировании нулевого и ненулевого многочлена.")
        self.assertEqual(Polynomial(5) + Polynomial(-3), Polynomial(2), "Ошибка при суммировании констант.")
        self.assertEqual(Polynomial("2x^2 + 4x + 5") + Polynomial("4x^2"), Polynomial("6x^2 + 4x + 5"), "Ошибка при суммировании многочленов одной степени.")
        self.assertEqual(Polynomial("2x^2 + 4x + 5") + Polynomial("- 4x^3"), Polynomial("- 4x^3 + 2x^2 + 4x + 5"), "Ошибка при суммировании многочленов разной степени.")
        self.assertEqual(Polynomial("3x^2 + 2x - 5") + Polynomial("-3x^2 - 2x + 5"), Polynomial(), "Ошибка при суммировании многочлена с его противополжным.")
        self.assertEqual(Polynomial("3x^2 + 2x - 5") + 5, Polynomial("3x^2 + 2x"), "Ошибка при суммировании многочлена с int.")
        self.assertEqual(5 + Polynomial("3x^2 + 2x - 5"), Polynomial("3x^2 + 2x"), "Ошибка при суммировании int с многочлена.")
        with self.assertRaises(TypeError, msg="add не поднял TypeError при неправильном типе правого слагаемого."): Polynomial() + 1.4
        with self.assertRaises(TypeError, msg="add не поднял TypeError при неправильном типе левого слагаемого."): 1.4 + Polynomial()

    def test_sub(self):
        self.assertEqual(Polynomial() - Polynomial(), Polynomial(), "Ошибка при вычитании нулевых многочленов.")
        self.assertEqual(Polynomial("3x^2 + 2x - 5") - Polynomial(), Polynomial("3x^2 + 2x - 5"), "Ошибка при вычитании нулевого и ненулевого многочлена.")
        self.assertEqual(Polynomial(5) - Polynomial(-3), Polynomial(8), "Ошибка при вычитании констант.")
        self.assertEqual(Polynomial("2x^2 + 4x + 5") - Polynomial("4x^2"), Polynomial("-2x^2 + 4x + 5"), "Ошибка при вычитании многочленов одной степени.")
        self.assertEqual(Polynomial("2x^2 + 4x + 5") - Polynomial("- 4x^3"), Polynomial("4x^3 + 2x^2 + 4x + 5"), "Ошибка при вычитании многочленов разной степени.")
        self.assertEqual(Polynomial("3x^2 + 2x - 5") - Polynomial("3x^2 + 2x - 5"), Polynomial(), "Ошибка при вычитании равных многочленов.")
        self.assertEqual(Polynomial("3x^2 + 2x + 5") - 5, Polynomial("3x^2 + 2x"), "Ошибка при левом вычитании int из многочлена.")
        self.assertEqual(5 - Polynomial("-3x^2 - 2x + 5"), Polynomial("3x^2 + 2x"), "Ошибка при правом вычитании многочлена из int.")
        with self.assertRaises(TypeError, msg="sub не поднял TypeError при неправильном типе правого слагаемого."): Polynomial() - 1.4
        with self.assertRaises(TypeError, msg="sub не поднял TypeError при неправильном типе левого слагаемого."): 1.4 - Polynomial()

    def test_mul(self):
        self.assertEqual(Polynomial() * Polynomial(), Polynomial(), "Ошибка при умножении нулевых многочленов.")
        self.assertEqual(Polynomial("3x^2 + 2x - 5") * Polynomial(), Polynomial(), "Ошибка при умножении ненулевого и нулевого многочлена.")
        self.assertEqual(Polynomial() * Polynomial("3x^2 + 2x - 5"), Polynomial(), "Ошибка при умножении нулевого и ненулевого многочлена.")
        self.assertEqual(Polynomial("5") * Polynomial("3x^2 + 2x - 5"), Polynomial("15x^2 + 10x - 25"), "Ошибка при умножении константы на многочлен.")
        self.assertEqual(Polynomial("x^2 + 2x + 1") * Polynomial("-3"), Polynomial("-3x^2 - 6x - 3"), "Ошибка при умножении многочлена на константу.")
        self.assertEqual(Polynomial(5) * Polynomial(-3), Polynomial(-15), "Ошибка при умножении констант.")
        self.assertEqual(Polynomial("3x^3 + 9x - 5") * Polynomial("-2x^2 + 4x + 6"), Polynomial("-6x^5 + 12x^4 + 46x^2 + 34x - 30"), "Ошибка умножения многочленов в общем случае.")
        self.assertEqual(Polynomial("3x^2 + 2x - 5") * 0, Polynomial(), "Ошибка при умножении многочлена на нулевой int.")
        self.assertEqual(0 * Polynomial("3x^2 + 2x - 5"), Polynomial(), "Ошибка при умножении нулевого int на многочлен.")
        self.assertEqual(5 * Polynomial("3x^2 + 2x - 5"), Polynomial("15x^2 + 10x - 25"), "Ошибка при умножении int на многочлен.")
        self.assertEqual(Polynomial("x^2 + 2x + 1") * -3, Polynomial("-3x^2 - 6x - 3"), "Ошибка при умножении многочлена на int.")
        with self.assertRaises(TypeError, msg="mul не поднял TypeError при неправильном типе правого множителя."): Polynomial() * 1.4
        with self.assertRaises(TypeError, msg="mul не поднял TypeError при неправильном типе левого множителя."): 1.4 * Polynomial()

    def test_bool(self):
        self.assertFalse(bool(Polynomial()), "Ненулевой многочлен даёт True, хотя должна быть False.")
        self.assertTrue(bool(Polynomial(5)), "Ненулевая константа даёт False, хотя должна быть True.")
        self.assertTrue(bool(Polynomial("5x^2")), "Моном даёт False, хотя должна быть True.")
        self.assertTrue(bool(Polynomial("5x^2 + 5x")), "Многочлен с нулевым свободным членом даёт False, хотя должна быть True.")
        self.assertTrue(bool(Polynomial("5x^2 + 5x + 5")), "Многочлен даёт False, хотя должна быть True.")

    def test_degree(self):
        self.assertEqual(Polynomial().degree(), -1, "Степень нулевого многочлена должна быть -1.")
        self.assertEqual(Polynomial(5).degree(), 0, "Степень константы должна быть нулём.")
        self.assertEqual(Polynomial("3x^2 + 2x + 1").degree(), 2, "Ошибка при определении степени в общем случае.")
        self.assertEqual(Polynomial("x^3 + 5").degree(), 3, "Ошибка при определении степени в многочлене с пропущенными степенями.")
        self.assertEqual(Polynomial("x^10 + 2x^9 + x^5").degree(), 10, "Ошибка при определении степени в многочлене без свободной и линейной части.")

    def test_mod(self):
        self.assertEqual(Polynomial() % 3, Polynomial(), "Ошибка при взятиии остатка нулевого многочлена.")
        self.assertEqual(Polynomial(5) % 3, Polynomial(2), "Ошибка при взятиии остатка положительной константы.")
        self.assertEqual(Polynomial(-4) % 3, Polynomial(2), "Ошибка при взятиии остатка отрицательной константы.")
        self.assertEqual(Polynomial("3x^2 + 2x + 1") % 3, Polynomial("2x + 1"), "Ошибка при взятиии остатка в общем случае.")
        self.assertEqual(Polynomial("3x^2 + 2x + 1") % 1, Polynomial(), "Ошибка при взятиии остатка деления на 1.")
        self.assertEqual(Polynomial("-1x^3 + 7x + 9") % -3, Polynomial("-x^3 - 2x"), "Ошибка при взятиии остатка деления на отрицательное число.")
        self.assertEqual(Polynomial("3x^4 + 2") % 2, Polynomial("x^4"), "Ошибка при взятиии остатка многочлена с пропущенными степенями.")

    def test_value(self):
        self.assertEqual(Polynomial().value(5), 0, "Значение нулевого многочлена должно быть 0 в любой точке.")
        self.assertEqual(Polynomial(5).value(10), 5, "Значение константы должно быть этой константой в любой точке.")
        self.assertEqual(Polynomial("3x^2 + 2x - 1").value(0), -1, "Значение многочлена в нуле должно быть равно его свободному члену.")
        self.assertEqual(Polynomial("3x^2 + 2x - 1").value(5), 84, "Ошибка вычисления значения полинома в общем случае.")
    
    def test_roots_BerlekampRabin(self):
        self.assertEqual(Polynomial("x + 1")._Polynomial__roots_BerlekampRabin(5), {4}, "Ошибка __roots_BerlekampRabin в случае линейного многочлена.")
        self.assertEqual(Polynomial("x^2 + 1")._Polynomial__roots_BerlekampRabin(3), set(), "Ошибка __roots_BerlekampRabin в случае неприводимого многочлена.")
        self.assertEqual(Polynomial("x^3 + x")._Polynomial__roots_BerlekampRabin(7), {0}, "Ошибка __roots_BerlekampRabin в случае, если многочлен расладывается на линейный и нелинейный.")
        self.assertEqual(Polynomial("3x^3 + 4x^2 + 2x + 1")._Polynomial__roots_BerlekampRabin(5), {1, 2, 4}, "Ошибка __roots_BerlekampRabin в общем случае.")

    def test_Fermat_little_theorem(self):
        self.assertEqual(Polynomial().Fermat_little_theorem(3), Polynomial(), "Нулевой многочлен после примениния Малой теоремы Ферма должен оставаться нулевым.")
        self.assertEqual(Polynomial(5).Fermat_little_theorem(3), Polynomial(2), "Константа после примениния Малой теоремы Ферма должена оставаться константой.")
        self.assertEqual(Polynomial("6x^8 - 3x^7 + 5x^6 + x^5 + x^3").Fermat_little_theorem(5), Polynomial("x^4 + 3x^3 + x"), "Ошибка примениния Малой теоремы Ферма в общем случае.")
    
    def test_div(self):
        self.assertEqual(Polynomial().div(Polynomial("x^2 + x + 1"), 3), (Polynomial(), Polynomial()), "Деление нулевого многочлена должно давать нулевой многочлен.")
        with self.assertRaises(ZeroDivisionError, msg="При попытке деления на нулевой многочлен должно подниматься исключение."): Polynomial("x^2 + x + 1").div(Polynomial(), 3)
        with self.assertRaises(TypeError, msg="div должен поднимать исключение при попытке деления на не полином или int."): Polynomial("x^2 + x + 1").div(1.4, 3)
        self.assertEqual(Polynomial("x^2 + x + 1").div(Polynomial("x^2 + x + 1"), 7), (Polynomial(1), Polynomial()), "При делении многочлена на самого себя должна получаться единица.")
        self.assertEqual(Polynomial("2x^3 + 5x + 1").div(2, 7), (Polynomial("x^3 + 6x + 4"), Polynomial()), "Ошибка при делении многочлена на константу.")
        self.assertEqual(Polynomial("x^2 + 2x + 1").div(Polynomial("x^3 + x + 1"), 5), (Polynomial(), Polynomial("x^2 + 2x + 1")), "Если степень делителя выше степени делимого, результатом должно быть частное 0 и остаток, равный самому делимому.")
        self.assertEqual(Polynomial("10x^2 + 2").div(Polynomial("3x"), 5), (Polynomial(), Polynomial("2")), "Если в делимом есть коэффициенты, кратные модулю, они должны быть занулены перед делением.")
        with self.assertRaises(ZeroDivisionError, msg="Если все коэффициенты делителя кратны модулю, должно подыматься исключение деления на нуль"): Polynomial("x^5 + x + 1").div(Polynomial("7x^2"), 7)
        self.assertEqual(Polynomial("3x^3 + x^2 + x + 5").div(Polynomial("x^2 + x + 1"),7), (Polynomial("3x + 5"), Polynomial()), "Ошибка при делении нацело.")
        self.assertEqual(Polynomial("3x^3 + x^2 + 2x + 6").div(Polynomial("x^2 + x + 1"),7), (Polynomial("3x + 5"), Polynomial("x + 1")), "Ошибка при делении в общем случае.")
    
    def test_derivative(self):
        self.assertEqual(Polynomial().derivative(), Polynomial(), "Производная нулевого многочлена должна быть нулевой.")
        self.assertEqual(Polynomial(5).derivative(), Polynomial(), "Производная константы должна быть нулевой.")
        self.assertEqual(Polynomial("3x^3 + x^2 + x + 5").derivative(), Polynomial("9x^2 + 2x + 1"), "Ошибка при определнии производной в общем случае.")

    def test_roots_prime_module(self):
        self.assertEqual(Polynomial()._Polynomial__roots_prime_module(5), {0, 1, 2 , 3, 4}, "Корнями нулевого многочлена должны являться все элементы поля.")
        self.assertEqual(Polynomial("5x^2 + 5x + 5")._Polynomial__roots_prime_module(5), {0, 1, 2 , 3, 4}, "Корнями нулевого многочлена должны являться все элементы поля (случай нулевого по модулю).")
        self.assertEqual(Polynomial("x - 3")._Polynomial__roots_prime_module(5), {3}, "Ошибка __roots_prime_module в случае линейного с одним корнем.")
        self.assertEqual(Polynomial(3)._Polynomial__roots_prime_module(5), set(), "Ошибка __roots_prime_module в случае константного без корней.")
        self.assertEqual(Polynomial("x^2 - 1")._Polynomial__roots_prime_module(7), {1, 6}, "Ошибка __roots_prime_module в случае квадратного с двумя корнями.")
        self.assertEqual(Polynomial("x^2 - 5")._Polynomial__roots_prime_module(13), set(), "Ошибка __roots_prime_module в случае квадратного без корней.")
        self.assertEqual(Polynomial("x^2 - 26")._Polynomial__roots_prime_module(13), {0}, "Ошибка __roots_prime_module в случае квадратного с одниим корнем.")
        self.assertEqual(Polynomial("3x^3 + 4x^2 + 2x + 1")._Polynomial__roots_prime_module(5), {1, 2, 4}, "Ошибка __roots_prime_module в общем случае.")

    def test_roots_primary_module(self):
        self.assertEqual(Polynomial("3x^3 + 4x^2 + 2x + 1")._Polynomial__roots_primary_module(5, 4), {(72, 625), (136, 625), (624, 625)}, "Ошибка _roots_primary_module в общем случае.")
        self.assertEqual(Polynomial("3x^3 + 4x^2 + 2x + 1")._Polynomial__roots_primary_module(5), {(1, 5), (2, 5), (4, 5)}, "Ошибка _roots_primary_module случае q = 1.")
        self.assertEqual(Polynomial("x^2 - 5")._Polynomial__roots_primary_module(13, 4), set(), "Ошибка __roots_primary_module в случае отсуствия корней по простому модулю.")

    def test_roots(self):
        self.assertEqual(Polynomial("162x^10+135x^9+162x^8+56x^4+162x^3+162x^2+113x+188").roots(189), {(94, 189), (67, 189), (121, 189)}, "Ошибка roots в общем случае.")
        self.assertEqual(Polynomial("x^3+11x^2+2x+8").roots(16), {(2, 16), (4, 16), (10, 16), (12, 16), (15, 16)}, "Ошибка roots в случае примарного модуля.")

    def test_gcd(self):
        self.assertEqual(gcd(Polynomial("2x^3 - 3x + 1"), Polynomial(), 5), Polynomial("x^3 + x + 3"), "Ошибка gcd в случае, если один из многочленов нулевой.")
        self.assertEqual(gcd(Polynomial("3x^4 + 6x^2 - 1"), Polynomial("x^2 - 7"), 13), Polynomial(1), "Ошибка gcd в случае взаимнопростых.")
        self.assertEqual(gcd(Polynomial(), Polynomial(), 13), Polynomial(), "Ошибка gcd в случае, если оба многочлена нулевые.")
        self.assertEqual(gcd(Polynomial("3x^4 + 6x^2 - 1"), Polynomial(5), 13), Polynomial(1), "Ошибка gcd в случае, если один многочлен констаната.")
        self.assertEqual(gcd(Polynomial("2x^3 - 3x + 1"), Polynomial("2x^3 - 3x + 1"), 5), Polynomial("x^3 + x + 3"), "Ошибка gcd в случае одинаковых многочленов.")
        self.assertEqual(gcd(Polynomial("2x^3 - 3x + 1"), Polynomial("x^3 + x + 3"), 5), Polynomial("x^3 + x + 3"), "Ошибка gcd в случае одинаковых унитарного и не унитарного многочленов.")
        self.assertEqual(gcd(Polynomial("2x - 2"), Polynomial("x^2 - 1"), 5), Polynomial("x + 4"), "Ошибка gcd в случае если один делитель друого.")
        self.assertEqual(gcd(Polynomial("x^3 + x"), Polynomial("x^3 + 3x^2 + 3x"), 5), Polynomial("x"), "Ошибка gcd в общем случае.")

    def test_binomial_theorem(self):
        self.assertEqual(binomial_theorem(2, 3, 0), Polynomial({0: 1}), "Ошибка binomial_theorem в случае exp=0.")
        self.assertEqual(binomial_theorem(2, 3, 1), Polynomial({1: 2, 0: 3}), "Ошибка binomial_theorem в случае exp=1.")
        self.assertEqual(binomial_theorem(0, 3, 2), Polynomial({0: 9}), "Ошибка binomial_theorem в случае m=0.")
        self.assertEqual(binomial_theorem(2, 0, 3), Polynomial({3: 8}), "Ошибка binomial_theorem в случае a=0.")
        self.assertEqual(binomial_theorem(1, 1, 5), Polynomial({5: 1, 4: 5, 3: 10, 2: 10, 1: 5, 0: 1}), "Ошибка binomial_theorem в общем случае.")
        self.assertEqual(binomial_theorem(-1, 2, 3), Polynomial({3: -1, 2: 6, 1: -12, 0: 8}), "Ошибка binomial_theorem случае m<0.")
        self.assertEqual(binomial_theorem(2, -3, 2), Polynomial({2: 4, 1: -12, 0: 9}), "Ошибка binomial_theorem случае a<0.")
        self.assertEqual(binomial_theorem(0, 0, 5), Polynomial(), "Ошибка binomial_theorem случае m=0, a=0." )

    def test_monic(self):
        self.assertEqual(Polynomial().monic(5), Polynomial(), "Ошибка monic в случае нулевого многочлена.")
        self.assertEqual(Polynomial("5x^3 - 3x^2 + 4x - 1").monic(7), Polynomial("x^3 + 5x^2 + 5x + 4"), "Ошибка monic в общем случае.")

if __name__ == '__main__':
    unittest.main()