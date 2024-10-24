from __future__ import annotations
from math import comb
from sympy import factorint
from itertools import product
from typing import Tuple, Dict, Set, Generator, Type
import logging
import re
from linear_comparisons import *
import unittest
from quadratic_comparisons import solve_quadratic_comparison
from ring import *
from field import *
import numpy as np

logger = logging.getLogger(__name__)
logger.setLevel(logging.DEBUG)
handler = logging.FileHandler(f'log/{__name__}.log')
handler.setLevel(logging.DEBUG)
handler.setFormatter(logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s'))
logger.addHandler(handler)

class Polynomial(Ring):

    def __init__(self, basic_ring: Ring):
        super().__init__()
        self.basic_ring = basic_ring
        self.has_multiplicative_identity = self.basic_ring.has_multiplicative_identity
        self.is_euclidean = isinstance(self.basic_ring, Field)
        self.is_commutative = self.basic_ring.is_commutative
        logger.debug(f'''
            Polynomial.__init__. self.basic_ring={self.basic_ring},
            self.has_multiplicative_identity={self.has_multiplicative_identity},
            self.is_euclidean={self.is_euclidean},
            self.is_commutative={self.is_commutative}
        ''')

    def __repr__(self) -> str:
        return f'Кольцо многочленов над {self.basic_ring}'
    
    def __str__(self) -> str:
        return self.__repr__()
    
    def __call__(self, raw_poly: None | Dict[int, Ring.Element] | Ring.Element | str = None):
        return self.Element(self, raw_poly)
    
    def __eq__(self, other: 'Polynomial') -> bool:
        if not isinstance(other, Polynomial): return False
        return self.basic_ring == other.basic_ring

    @staticmethod
    def gcd(f: Polynomial.Element, g: Polynomial.Element) -> Polynomial.Element:
        logger.debug(f'Начало Polynomial.gcd. f={f}, g={g}')
        while g:
            logger.debug(f'Шаг Polynomial.gcd. f={f}, g={g}')
            f, g = g, f % g
        logger.debug(f'Polynomial.gcd, делаем унитарный. f={f}, g={g}')
        f = f.monic()
        logger.debug(f'Конец Polynomial.gcd. f={f}, g={g}')
        return f
    
    def binomial_theorem(self, m: Ring.Element, a: Ring.Element, exp: int) -> Polynomial.Element:
        """
        Разложение (mx+a)^exp
        """
        logger.debug(f'Polynomial.binomial_theorem. self={self}, m={m}, a={a}, exp={exp}')
        poly = self()
        for k in range(exp + 1):
            coeff = comb(exp, k) * (m ** (exp - k)) * (a ** k)
            poly[exp - k] = coeff
            logger.debug(f'Polynomial.binomial_theorem. self={self}, m={m}, a={a}, exp={exp}, poly={poly}')
        return poly

    class Element(Ring.Element):
        def __init__(self, ring: Polynomial, raw_poly: None | Dict[int, Ring.Element] | Ring.Element | str = None):
            super().__init__(ring)
            if isinstance(self.ring.basic_ring, GF):
                self.Fermat_little_theorem = self.__Fermat_little_theorem
            if self.ring.basic_ring.has_multiplicative_inverses:
                self.monic = self.__monic
            logger.debug(f'Polynomial.Element.__init__. ring={self.ring}, raw_poly={raw_poly}')
            if not raw_poly:
                self.value = dict()
                logger.debug(f'Polynomial.Element.__init__ завершается. self={self}, raw_poly={raw_poly}')
                return
            elif isinstance(raw_poly, Polynomial.Element):
                self.value = {exp:self.ring.basic_ring(coeff) for exp, coeff in raw_poly.monomials()}
                logger.debug(f'Polynomial.Element.__init__ завершается. self={self}, raw_poly={raw_poly}')
                return 
            elif isinstance(raw_poly, dict):
                logger.debug(f'Polynomial.Element.__init__. type(list(raw_poly.values())[0]) = {type(list(raw_poly.values())[0])}')
                self.value = {exp:coeff for exp, coeff in raw_poly.items() if isinstance(exp, int) and exp >= 0 and coeff and isinstance(coeff, Ring.Element) and coeff.ring == self.ring.basic_ring}
                logger.debug(f'Polynomial.Element.__init__ завершается. self={self}, raw_poly={raw_poly}')
                return
            elif not isinstance(raw_poly, str):
                try:
                    self.value = {0:self.ring.basic_ring(raw_poly)} if raw_poly else dict()
                    logger.debug(f'Polynomial.Element.__init__ завершается. self={self}, raw_poly={raw_poly}')
                    return
                except Exception:
                    logger.critical(f'Polynomial.Element.__init__. ring={ring}, raw_poly={raw_poly}. Неправильный тип raw_poly.')
                    raise TypeError("Полином создаётся только из строки, словаря, элемента кольца (при этом тип добавляемых элементов должен соответствовать типу коэффициентов многочлена).")
            pattern = r'([+-]?[^-+]+)'
            monomials = re.findall(pattern, raw_poly)
            self.value = dict()
            logger.debug(f'Polynomial.Element.__init__. self={self}, raw_poly={raw_poly}, monomials={monomials}')
            for monomial in monomials:
                splited = monomial.replace(' ', '').split('x')
                coeff = self.ring.basic_ring(1) if splited[0] in ['', '+'] else self.ring.basic_ring(-1) if splited[0] == '-' else self.ring.basic_ring(splited[0])
                if not coeff: continue
                exp = 0 if len(splited) != 2 else (1 if splited[1] == '' else int(splited[1][1:]))
                if exp < 0:
                    logger.critical("Polynomial.Element.__init__. Отрицательная степень.")
                    raise ValueError("Отрицательные степени не поддерживаются.")
                self[exp] += coeff
            logger.debug(f'Polynomial.Element.__init__ завершается. self={self}, raw_poly={raw_poly}')

        def __getitem__(self, exp: int) -> Ring.Element:
            return self.value.get(exp, self.ring.basic_ring())
        
        def __setitem__(self, exp: int, coeff: Ring.Element) -> None:
            coeff = self.ring.basic_ring(coeff)
            if not coeff: self.value.pop(exp,0)
            else: self.value[exp] = coeff

        def __contains__(self, exp: int) -> bool:
            return exp in self.value
        
        def __iter__(self) -> Generator[int]:
            for exp, coeff in self.value.items():
                if coeff: yield exp
    
        def __eq__(self, other: 'Polynomial.Element') -> bool:
            if not isinstance(other, Polynomial.Element): return False
            if self.ring != other.ring: return False
            return all(self[exp] == other[exp] for exp in set(self) | set(other))
    
        def monomials(self)  -> Generator[Tuple[int,Ring.Element]]:
            """
            Генератор наборов (exp, coeff).
            """
            for exp, coeff in self.value.items():
                if coeff: yield exp, coeff

        def __hash__(self) -> int:
            return hash(tuple(sorted((exp, coeff) for exp, coeff in self.monomials())))

        def copy(self) -> 'Polynomial.Element':
            return self.ring(self.value)

        def __add__(self, other: 'Polynomial.Element' | Ring.Element) -> 'Polynomial.Element':
            logger.debug(f'Polynomial.Element.__add__. self={self}, other={other}')
            if not isinstance(other, Polynomial.Element):
                other = self.ring(other)
            coeffs = dict()
            for i in set(self) | set(other):
                new_coeff = self[i] + other[i]
                if new_coeff: coeffs[i] = new_coeff
            logger.debug(f'Polynomial.Element.__add__. self={self}, other={other}, coeffs={coeffs}')
            return self.ring(coeffs)
    
        def __neg__(self):
            return self.ring({exp:-coeff for exp,coeff in self.monomials()})
    
        def __sub__(self, other: 'Polynomial.Element' | Ring.Element) -> 'Polynomial.Element':
            logger.debug(f'Polynomial.__sub__. self={self}, other={other}')
            if not isinstance(other, Polynomial.Element):
                other = self.ring(other)
            coeffs = dict()
            for i in set(self) | set(other):
                new_coeff = self[i] - other[i]
                if new_coeff: coeffs[i] = new_coeff
            logger.debug(f'Polynomial.Element.__sub__. self={self}, other={other}, coeffs={coeffs}')
            return self.ring(coeffs)
    
        def __mul__(self, other: 'Polynomial.Element' | Ring.Element) -> 'Polynomial.Element':
            logger.debug(f'Polynomial.Element.__mul__. self={self}, other={other}')
            if not isinstance(other, Polynomial.Element):
                other = self.ring(other)
            result = self.ring()
            for exp1, coeff1 in self.monomials():
                for exp2, coeff2 in other.monomials():
                    result[exp1+exp2] += coeff1 * coeff2
            result = self.ring({exp:coeff for exp, coeff in result.monomials() if coeff})
            logger.debug(f'Polynomial.__mul__. self={self}, other={other}, result={result}')
            return result
        
        def __rmul__(self, other) -> 'Ring.Element':
            logger.debug(f'Polynomial.Element.__rmul__. self={self}, other={other}')
            if not isinstance(other, Polynomial.Element):
                other = self.ring(other)
            return other * self
    
        def __bool__(self) -> bool:
            return any(coeff for coeff in self.value.values())
    
        def __repr__(self) -> str:
            if not self: return str(self.ring.basic_ring())
            monomials = []
            for exp in sorted(list(self), reverse=True):
                coeff = self[exp]
                if not coeff: continue
                if exp == 0: monomial = f"{coeff}"
                elif exp == 1: monomial = f"{coeff}x"
                else: monomial = f"{coeff}x^{exp}"
                monomials.append(monomial)
            polynomial_str = " + ".join(monomials)
            polynomial_str = polynomial_str.replace("+ -", "- ")
            return polynomial_str
        
        def __gt__(self, other: 'Polynomial.Element' | Ring.Element) -> bool:
            if not isinstance(other, Polynomial.Element): other = self.ring(other)
            return self.degree() > other.degree()
        
        def __ge__(self, other: 'Polynomial.Element' | Ring.Element) -> bool:
            if not isinstance(other, Polynomial.Element): other = self.ring(other)
            return self.degree() >+ other.degree()
        
        def __lt__(self, other: 'Polynomial.Element' | Ring.Element) -> bool:
            if not isinstance(other, Polynomial.Element): other = self.ring(other)
            return self.degree() < other.degree()
        
        def __le__(self, other: 'Polynomial.Element' | Ring.Element) -> bool:
            if not isinstance(other, Polynomial.Element): other = self.ring(other)
            return self.degree() <= other.degree()
    
        def degree(self) -> int:
            logger.debug(f'Polynomial.Element.degree. self={self}')
            if not self:
                logger.debug(f'Polynomial.Element.degree. self={self}. Нулевой многочлен')
                return -1
            degree = max(exp for exp in self if self[exp])
            logger.debug(f'Polynomial.Element.degree. self={self}, degree={degree}')
            return degree
    
        def __call__(self, x: Ring.Element) -> Ring.Element:
            """
            Значение многочлена в точке x
            """
            logger.debug(f'Polynomial.Element.__call__. self={self}')
            if not isinstance(x, Ring.Element) or not x.ring == self.ring.basic_ring:
                x =  self.ring.basic_ring(x)
            value = sum(coeff * x ** exp for exp,coeff in self.monomials())
            logger.debug(f'Polynomial.Element.__call__. self={self}, value={value}')
            return value
    
        def __Fermat_little_theorem(self) -> 'Polynomial.Element':
            """
            Сокращение многочлена по Малой теореме Ферма в GF(p)
            x^(p) = x (mod p) <==> x^(kp+t) = x^(k+t) (mod p)
            """
            p = self.ring.basic_ring.p
            logger.debug(f'Polynomial.Element.Fermat_little_theorem. self={self}')
            poly = self.ring()
            for exp,coeff in self.monomials():
                while exp >= p: exp = exp // p + exp % p
                poly[exp] += coeff
            logger.debug(f'Polynomial.Element.Fermat_little_theorem. self={self}, poly={poly}')
            return poly
    
        def __divmod__(self, divisor: 'Polynomial.Element' | Ring.Element) -> Tuple['Polynomial.Element', 'Polynomial.Element']:
            logger.debug(f'Polynomial.Element.__divmod__. self={self}, divisor={divisor}')
            if not self.ring.is_euclidean:
                raise TypeError(f'Деление не определено в {self.ring}')
            if not isinstance(divisor, Polynomial.Element):
                divisor = self.ring(divisor)
            if not divisor:
                logger.critical(f'Polynomial.Element.__divmod__. self={self}, divisor={divisor}. Нулевой divisor.')
                raise ZeroDivisionError("Деление на нулевой многочлен.")
            remainder = self.copy()
            remainder_degree = remainder.degree()
            quotient = self.ring()
            divisor_degree = divisor.degree()
            while remainder and remainder_degree >= divisor_degree:
                logger.debug(f'Polynomial.Element.__divmod__. self={self}, divisor={divisor}, quotient={quotient}, remainder={remainder}')
                degree_diff = remainder_degree - divisor_degree
                leading_coeff = remainder[remainder_degree] / divisor[divisor_degree]
                quotient[degree_diff] = leading_coeff
                for exp, coeff in divisor.monomials():
                    remainder[exp + degree_diff] = remainder[exp + degree_diff] - leading_coeff * coeff
                remainder = self.ring({exp: coeff for exp, coeff in remainder.monomials() if coeff})
                remainder_degree = remainder.degree()
            logger.debug(f'Polynomial.Element.__divmod__. self={self}, divisor={divisor}, quotient={quotient}, remainder={remainder}')
            return quotient, remainder
        
        def __truediv__(self, divisor: 'Polynomial.Element' | Ring.Element) -> 'Polynomial.Element':
            return divmod(self, divisor)[0]
        
        def __mod__(self, divisor: 'Polynomial.Element' | Ring.Element) -> 'Polynomial.Element':
            return divmod(self, divisor)[1]
    
        def __roots_BerlekampRabin(self) -> Set[GF.Element]:
            p = self.ring.basic_ring.p
            logger.debug(f'Polynomial.Element.__roots_BerlekampRabin. self={self}, p={p}')
            roots = set()
            nonlinears = [self]
            for delta in range(p):
                if not nonlinears: break
                poly = nonlinears.pop()
                logger.debug(f'Polynomial.Element.__roots_BerlekampRabin. self={self}, delta={delta}, roots={roots}, nonlinears = {nonlinears}, poly={poly}')
                d = Polynomial.gcd(poly, self.ring.binomial_theorem(1,delta,(p-1) >> 1) - 1)
                if d.degree() <= 0: nonlinears.append(poly); continue
                for factor in [d, poly / d]:
                    if factor.degree() == 1: roots |= solve_linear_comparison(int(factor[1]), int(-factor[0]), p)
                    else: nonlinears.append(factor)
            logger.debug(f'Polynomial.Element.__roots_BerlekampRabin. self={self}, p={p}, roots={roots}')
            return roots
    
        def derivative(self) -> 'Polynomial.Element':
            logger.debug(f'Polynomial.Element.derivative. self={self}')
            derivative = self.ring({exp-1:coeff*exp for exp,coeff in self.monomials() if exp > 0})
            logger.debug(f'Polynomial.Element.derivative. self={self}, derivative={derivative}')
            return derivative

        def __roots_prime_module(self) -> Set[GF.Element]:
            p = self.ring.basic_ring.p
            logger.debug(f'Polynomial.Element.__roots_prime_module. self={self}, p={p}')
            fermat = self.Fermat_little_theorem()
            if fermat.degree() <= 1:
                logger.debug(f'Polynomial.Element.__roots_prime_module. self={self}, p={p} fermat={fermat}. Степень <= 1.')
                return solve_linear_comparison(int(fermat[1]), int(-fermat[0]), p)
            reduced = Polynomial.gcd(fermat, self.ring({p:1,1:-1}))
            if reduced.degree() == 2:
                logger.debug(f'Polynomial.Element.__roots_prime_module. self={self}, p={p} reduced={reduced}. Степень 2.')
                return solve_quadratic_comparison(int(reduced[2]), int(reduced[1]), int(reduced[0]), p)
            return reduced.__roots_BerlekampRabin()
    
        def __roots_primary_module(self, p: int, q: int = 1) -> Set[IntegerModRing.Element]:
            logger.debug(f'Polynomial.Element.__roots_primary_module. self={self}, p={p}, q={q}')
            gfp = GF(p)
            result = {int(a) for a in Polynomial(gfp)(self).__roots_prime_module()}
            derivative = Polynomial(gfp)(self.derivative()).Fermat_little_theorem()
            p_i = 1
            for _ in range(1,q):
                p_i *= p
                result = {a:solve_linear_comparison(int(derivative(a)), - int(self(a)) // p_i, p) for a in result}
                logger.debug(f'Polynomial.Element.__roots_primary_module. self={self}, p={p}, q={q}, result={result}')
                result = set(a + int(t) * p_i for a, temp in result.items() for t in temp)
            ring = IntegerModRing(p_i * p)
            logger.debug(f'Polynomial.Element.__roots_primary_module. self={self}, p={p}, q={q}, result={result}')
            return set(ring(a) for a in result)
    
        def roots(self):
            logger.debug(f'Polynomial.Element.roots. self={self}')
            if isinstance(self.ring.basic_ring, RealField):
                return set(np.roots([float(self[i]) for i in range(self.degree() + 1)]))
            if isinstance(self.ring.basic_ring, ComplexField):
                return set(np.roots([complex(self[i]) for i in range(self.degree() + 1)]))
            if isinstance(self.ring.basic_ring, IntegerModRing):
                M = self.ring.basic_ring.m
                if self.degree() <= 1:
                    logger.debug(f'Polynomial.Element.roots. self={self}, M={M}. Степень <= 1.')
                    return set(solve_linear_comparison(int(self[1]),int(-self[0]),M))
                primary_results = []
                for p,q in factorint(M).items():
                    primary_results.append(Polynomial(IntegerModRing(p**q))(self).__roots_primary_module(p,q))
                if len(primary_results) == 1:
                    logger.debug(f'Polynomial.Element.roots. self={self}, M={M}, primary_results={primary_results}. M - примарное.')
                    return primary_results[0]
                result = set()
                for element in product(*primary_results):
                    result.add(crt(element))
                logger.debug(f'Polynomial.Element.roots. self={self}, M={M}, result={result}')
                return result
            raise TypeError(f'Поиск корней в {self.ring.basic_ring} не поддерживается.')
    
        def __monic(self) -> 'Polynomial.Element':
            if not self: return self.ring()
            return self * (self[self.degree()] ** (-1))
class TestPolynomial(unittest.TestCase):

    def setUp(self):
        return super().setUp()
    
    def test_init(self):
        ring = Polynomial(IntegerRing())
        # константы
        self.assertEqual(ring("0"), ring(), "Ошибка инициализации нулевого полинома.")
        self.assertEqual(ring(""), ring(), "Ошибка инициализации нулевого полинома (пустая строка).")
        self.assertEqual(ring("5"), ring(5), "Ошибка инициализации констатны.")
        self.assertEqual(ring("+5"), ring(5), "Ошибка инициализации положительной константы (плюс без пробела).")
        self.assertEqual(ring("+ 5"), ring(5), "Ошибка инициализации положительной константы (плюс с пробелом).")
        self.assertEqual(ring("-5"), ring(-5), "Ошибка инициализации отрицательной константы (без пробела).")
        self.assertEqual(ring("- 5"), ring(-5), "Ошибка инициализации отрицательной константы (с пробелом).")
        # линейный моном
        self.assertEqual(ring("0x"), ring(), "Ошибка инициализации линейного монома с нулевым коэффициентом.")
        self.assertEqual(ring("x"), ring({1:IntegerRing()(1)}), "Ошибка инициализации линейного монома с коэффициентом, равным 1.")
        self.assertEqual(ring("+x"), ring({1:IntegerRing()(1)}), "Ошибка инициализации линейного монома с коэффициентом, равным 1 (плюс без пробела).")
        self.assertEqual(ring("+ x"), ring({1:IntegerRing()(1)}), "Ошибка инициализации линейного монома с коэффициентом, равным 1 (плюс с пробелом).")
        self.assertEqual(ring("-x"), ring({1:IntegerRing()(-1)}), "Ошибка инициализации линейного монома с коэффициентом, равным -1 (без пробела).")
        self.assertEqual(ring("- x"), ring({1:IntegerRing()(-1)}), "Ошибка инициализации линейного монома с коэффициентом, равным -1 (с пробелом).")
        self.assertEqual(ring("5x"), ring({1:IntegerRing()(5)}), "Ошибка инициализации линейного монома с положительным коэффициентом.")
        self.assertEqual(ring("+5x"), ring({1:IntegerRing()(5)}), "Ошибка инициализации линейного монома с положительным коэффициентом (плюс без пробела).")
        self.assertEqual(ring("+ 5x"), ring({1:IntegerRing()(5)}), "Ошибка инициализации линейного монома с положительным коэффициентом (плюс с пробелом).")
        self.assertEqual(ring("-5x"), ring({1:IntegerRing()(-5)}), "Ошибка инициализации линейного монома с отрицательным коэффициентом (без пробела).")
        self.assertEqual(ring("- 5x"), ring({1:IntegerRing()(-5)}), "Ошибка инициализации линейного монома с отрицательным коэффициентом (с пробелом).")
        # произвольный моном
        self.assertEqual(ring("0x^2"), ring(), "Ошибка инициализации произвольного монома с нулевым коэффициентом.")
        self.assertEqual(ring("5x^2"), ring({2:IntegerRing()(5)}), "Ошибка инициализации произвольного монома с положительным коэффициентом.")
        self.assertEqual(ring("+5x^2"), ring({2:IntegerRing()(5)}), "Ошибка инициализации произвольного монома с положительным коэффициентом (плюс без пробела).")
        self.assertEqual(ring("+ 5x^2"), ring({2:IntegerRing()(5)}), "Ошибка инициализации произвольного монома с положительным коэффициентом (плюс с пробелом).")
        self.assertEqual(ring("-5x^2"), ring({2:IntegerRing()(-5)}), "Ошибка инициализации произвольного монома с отрицательным коэффициентом (без пробела).")
        self.assertEqual(ring("- 5x^2"), ring({2:IntegerRing()(-5)}), "Ошибка инициализации произвольного монома с отрицательным коэффициентом (с пробелом).")
        # полиномы
        self.assertEqual(ring("3x^3 - 4x^2 + 2x - 5"), ring({0:IntegerRing()(-5), 1:IntegerRing()(2), 2:IntegerRing()(-4), 3:IntegerRing()(3)}), "Ошибка инициализации в общем случае.")
        self.assertEqual(ring("3x^3-4x^2+2x-5"), ring({0:IntegerRing()(-5), 1:IntegerRing()(2), 2:IntegerRing()(-4), 3:IntegerRing()(3)}), "Ошибка инициализации в общем случае (безпробельное написание).")
        self.assertEqual(ring("2x^2 - 5"), ring({0:IntegerRing()(-5), 2:IntegerRing()(2)}), "Ошибка инициализации при нулевом коэффициенте.")
        self.assertEqual(ring("2x^2 + 2x"), ring({1:IntegerRing()(2), 2:IntegerRing()(2)}), "Ошибка инициализации при нулевом свободном члене.")
        self.assertEqual(ring("3x^2 + 5x^2 - 2x"), ring({1:IntegerRing()(-2), 2:IntegerRing()(8)}), "Ошибка инициализации в случае с одинаковми степенями.")
        # словарь
        self.assertEqual(ring({1.4:1.4}), ring(), "Ошибка проверки типов при инициализации словарём.")
        self.assertEqual(ring({-5:IntegerRing()(5)}), ring(), "Ошибка проверки положительности степеней при инициализации словарём.")
        self.assertEqual(ring({5:IntegerRing()(0)}), ring(), "Ошибка проверки коэффициентов при инициализации словарём.")
        # исключенияY
        with self.assertRaises(ValueError, msg="init не поднял исключение при отрицательной степени в строке."): ring("5x^-5")
        # многочлены не над int
        ring = Polynomial(RealField())
        self.assertEqual(ring("-2.4x^3 + 5.33x^2 - x + 0.00"), ring({1:RealField()(-1), 2:RealField()(5.33), 3:RealField()(-2.4)}), "Ошибка инициализации вещественного многочлена.")
        ring = Polynomial(ComplexField())
        self.assertEqual(ring({1:5, 2:ComplexField()(complex(1,2))}), ring({2:ComplexField()(complex(1,2))}), "Ошибка инициализации комплексного многочлена.")
        gf = GF(7)
        ring = Polynomial(gf)
        self.assertEqual(ring("-2x^3 + 9x^2 - 5x + 1"), ring({0:gf(1), 1:gf(2), 2:gf(2), 3:gf(5)}), "Ошибка инициалицации многочлена над кольцом.")
    
    def test_str(self):
        ring = Polynomial(IntegerRing())
        self.assertEqual(str(ring()), "0", "Ошибка преобразования в строку нулевого полинома.")
        self.assertEqual(str(ring("5")), "5", "Ошибка преобразования в строку положительной констаты.")
        self.assertEqual(str(ring("-5")), "-5", "Ошибка преобразования в строку отрицательной констаты.")
        self.assertEqual(str(ring("5x")), "5x", "Ошибка преобразования в строку линейного монома с положительным коэффициентом.")
        self.assertEqual(str(ring("-5x")), "-5x", "Ошибка преобразования в строку линейного монома с отрицательным коэффициентом.")
        self.assertEqual(str(ring("5x^2")), "5x^2", "Ошибка преобразования в строку произвольного монома с положительным коэффициентом.")
        self.assertEqual(str(ring("-5x^2")), "-5x^2", "Ошибка преобразования в строку произвольного монома с отрицательным коэффициентом.")
        self.assertEqual(str(ring("3x^3 - 4x^2 + 2x - 5")), "3x^3 - 4x^2 + 2x - 5", "Ошибка преобразования в строку в общем случае.")
        self.assertEqual(str(ring("2x^2 - 5")), "2x^2 - 5", "Ошибка преобразования в строку при нулевом коэффициенте.")
        self.assertEqual(str(ring("2x^2 + 2x")), "2x^2 + 2x", "Ошибка преобразования в строку при нулевом свободном члене.")
        self.assertEqual(str(Polynomial(RealField())("-2.4x^3 + 5.33x^2 - x + 0.00")), "-2.4x^3 + 5.33x^2 - 1.0x", "Ошибка преобразования в строку вещественного многочлена.")
        self.assertEqual(str(Polynomial(ComplexField())({0:ComplexField()(complex(1,2)), 1:ComplexField()(complex(7,8)), 2:ComplexField()(complex(3,4)), 3:ComplexField()(complex(5,6))})), "(5+6j)x^3 + (3+4j)x^2 + (7+8j)x + (1+2j)", "Ошибка преобразования в строку комплексного многочлена.")
        self.assertEqual(str(Polynomial(GF(7))("-2x^3 + 9x^2 - 5x + 1")), "[5]_7x^3 + [2]_7x^2 + [2]_7x + [1]_7", "Ошибка преобразования в строку многочлена над кольцом.")

    def test_getitem(self):
        ring = Polynomial(IntegerRing())
        self.assertEqual(ring('5x^2')[2], 5, "Ошибка при получении значения ненулевого коэффициента.")
        self.assertEqual(ring()[2], 0, "Ошибка при получении значения нулевого коэффициента.")

    def test_setitem(self):
        ring = Polynomial(IntegerRing())
        p = ring()
        p[2] = 5
        self.assertEqual(p, ring({2:IntegerRing()(5)}), "Ошибка при установке коэффициента.")
        p[2] = 2
        self.assertEqual(p, ring({2:IntegerRing()(2)}), "Ошибка при изменении коэффициента.")
        p[2] = 0
        self.assertEqual(p, ring(), "Ошибка обнуления коэффициента.")
    
    def test_contains(self):
        p = Polynomial(IntegerRing())('5x^2')
        self.assertTrue(2 in p, "Ошибка false negative при проверке содержания степени в полиноме.")
        self.assertFalse(3 in p, "Ошибка false positive при проверке содержания степени в полиноме.")

    def test_iter(self):
        ring = Polynomial(IntegerRing())
        self.assertEqual(list(ring()), [], "Ошибка iter при нулевом многочлене.")
        self.assertEqual(list(ring({0:IntegerRing()(1), 2:IntegerRing()(3), 4:IntegerRing()(5)})), [0,2,4], "Ошибка iter при ненулевом многочлене.")
    
    def test_monomials(self):
        ring = Polynomial(IntegerRing())
        self.assertEqual(list(ring().monomials()), [], "Ошибка monomials при нулевом многочлене.")
        self.assertEqual(list(ring({0:IntegerRing()(1), 2:IntegerRing()(3), 4:IntegerRing()(5)}).monomials()), [(0,1),(2,3),(4,5)], "Ошибка monomials при ненулевом многочлене.")

    def test_copy(self):
        ring = Polynomial(IntegerRing())
        p = ring()
        self.assertEqual(p, p.copy(), "Ошибка copy при нулевом многочлене.")
        p = ring({0:1, 2:3, 4:5})
        self.assertEqual(p, p.copy(), "Ошибка copy при ненулевом многочлене.")
    
    def test_add(self):
        ring = Polynomial(IntegerRing())
        self.assertEqual(ring() + ring(), ring(), "Ошибка при суммировании нулевых многочленов.")
        self.assertEqual(ring("3x^2 + 2x - 5") + ring(), ring("3x^2 + 2x - 5"), "Ошибка при суммировании нулевого и ненулевого многочлена.")
        self.assertEqual(ring(5) + ring(-3), ring(2), "Ошибка при суммировании констант.")
        self.assertEqual(ring("2x^2 + 4x + 5") + ring("4x^2"), ring("6x^2 + 4x + 5"), "Ошибка при суммировании многочленов одной степени.")
        self.assertEqual(ring("2x^2 + 4x + 5") + ring("- 4x^3"), ring("- 4x^3 + 2x^2 + 4x + 5"), "Ошибка при суммировании многочленов разной степени.")
        self.assertEqual(ring("3x^2 + 2x - 5") + ring("-3x^2 - 2x + 5"), ring(), "Ошибка при суммировании многочлена с его противополжным.")
        self.assertEqual(ring("3x^2 + 2x - 5") + 5, ring("3x^2 + 2x"), "Ошибка при суммировании многочлена с int.")
        self.assertEqual(5 + ring("3x^2 + 2x - 5"), ring("3x^2 + 2x"), "Ошибка при суммировании int с многочлена.")

    def test_sub(self):
        ring = Polynomial(IntegerRing())
        self.assertEqual(ring() - ring(), ring(), "Ошибка при вычитании нулевых многочленов.")
        self.assertEqual(ring("3x^2 + 2x - 5") - ring(), ring("3x^2 + 2x - 5"), "Ошибка при вычитании нулевого и ненулевого многочлена.")
        self.assertEqual(ring(5) - ring(-3), ring(8), "Ошибка при вычитании констант.")
        self.assertEqual(ring("2x^2 + 4x + 5") - ring("4x^2"), ring("-2x^2 + 4x + 5"), "Ошибка при вычитании многочленов одной степени.")
        self.assertEqual(ring("2x^2 + 4x + 5") - ring("- 4x^3"), ring("4x^3 + 2x^2 + 4x + 5"), "Ошибка при вычитании многочленов разной степени.")
        self.assertEqual(ring("3x^2 + 2x - 5") - ring("3x^2 + 2x - 5"), ring(), "Ошибка при вычитании равных многочленов.")
        self.assertEqual(ring("3x^2 + 2x + 5") - 5, ring("3x^2 + 2x"), "Ошибка при левом вычитании int из многочлена.")
        self.assertEqual(5 - ring("-3x^2 - 2x + 5"), ring("3x^2 + 2x"), "Ошибка при правом вычитании многочлена из int.")

    def test_mul(self):
        ring = Polynomial(IntegerRing())
        self.assertEqual(ring() * ring(), ring(), "Ошибка при умножении нулевых многочленов.")
        self.assertEqual(ring("3x^2 + 2x - 5") * ring(), ring(), "Ошибка при умножении ненулевого и нулевого многочлена.")
        self.assertEqual(ring() * ring("3x^2 + 2x - 5"), ring(), "Ошибка при умножении нулевого и ненулевого многочлена.")
        self.assertEqual(ring("5") * ring("3x^2 + 2x - 5"), ring("15x^2 + 10x - 25"), "Ошибка при умножении константы на многочлен.")
        self.assertEqual(ring("x^2 + 2x + 1") * ring("-3"), ring("-3x^2 - 6x - 3"), "Ошибка при умножении многочлена на константу.")
        self.assertEqual(ring(5) * ring(-3), ring(-15), "Ошибка при умножении констант.")
        self.assertEqual(ring("3x^3 + 9x - 5") * ring("-2x^2 + 4x + 6"), ring("-6x^5 + 12x^4 + 46x^2 + 34x - 30"), "Ошибка умножения многочленов в общем случае.")
        self.assertEqual(ring("3x^2 + 2x - 5") * 0, ring(), "Ошибка при умножении многочлена на нулевой int.")
        self.assertEqual(0 * ring("3x^2 + 2x - 5"), ring(), "Ошибка при умножении нулевого int на многочлен.")
        self.assertEqual(5 * ring("3x^2 + 2x - 5"), ring("15x^2 + 10x - 25"), "Ошибка при умножении int на многочлен.")
        self.assertEqual(ring("x^2 + 2x + 1") * -3, ring("-3x^2 - 6x - 3"), "Ошибка при умножении многочлена на int.")

    def test_bool(self):
        ring = Polynomial(IntegerRing())
        self.assertFalse(bool(ring()), "Ненулевой многочлен даёт True, хотя должна быть False.")
        self.assertTrue(bool(ring(5)), "Ненулевая константа даёт False, хотя должна быть True.")
        self.assertTrue(bool(ring("5x^2")), "Моном даёт False, хотя должна быть True.")
        self.assertTrue(bool(ring("5x^2 + 5x")), "Многочлен с нулевым свободным членом даёт False, хотя должна быть True.")
        self.assertTrue(bool(ring("5x^2 + 5x + 5")), "Многочлен даёт False, хотя должна быть True.")

    def test_degree(self):
        ring = Polynomial(IntegerRing())
        self.assertEqual(ring().degree(), -1, "Степень нулевого многочлена должна быть -1.")
        self.assertEqual(ring(5).degree(), 0, "Степень константы должна быть нулём.")
        self.assertEqual(ring("3x^2 + 2x + 1").degree(), 2, "Ошибка при определении степени в общем случае.")
        self.assertEqual(ring("x^3 + 5").degree(), 3, "Ошибка при определении степени в многочлене с пропущенными степенями.")
        self.assertEqual(ring("x^10 + 2x^9 + x^5").degree(), 10, "Ошибка при определении степени в многочлене без свободной и линейной части.")

    def test_call(self):
        ring = Polynomial(IntegerRing())
        self.assertEqual(ring()(5), 0, "Значение нулевого многочлена должно быть 0 в любой точке.")
        self.assertEqual(ring(5)(10), 5, "Значение константы должно быть этой константой в любой точке.")
        self.assertEqual(ring("3x^2 + 2x - 1")(0), -1, "Значение многочлена в нуле должно быть равно его свободному члену.")
        self.assertEqual(ring("3x^2 + 2x - 1")(5), 84, "Ошибка вычисления значения полинома в общем случае.")
        self.assertEqual(Polynomial(GF(7))("5x^3 + 2x^2 + 3x + 1")(4), 1, "Ошибка вычисления значения полинома в поле.")
    
    def test_roots_BerlekampRabin(self):
        self.assertEqual(Polynomial(GF(5))("x + 1")._Element__roots_BerlekampRabin(), {GF(5)(4)}, "Ошибка __roots_BerlekampRabin в случае линейного многочлена.")
        self.assertEqual(Polynomial(GF(3))("x^2 + 1")._Element__roots_BerlekampRabin(), set(), "Ошибка __roots_BerlekampRabin в случае неприводимого многочлена.")
        self.assertEqual(Polynomial(GF(7))("x^3 + x")._Element__roots_BerlekampRabin(), {GF(7)()}, "Ошибка __roots_BerlekampRabin в случае, если многочлен расладывается на линейный и нелинейный.")
        self.assertEqual(Polynomial(GF(5))("3x^3 + 4x^2 + 2x + 1")._Element__roots_BerlekampRabin(), {GF(5)(1), GF(5)(2), GF(5)(4)}, "Ошибка __roots_BerlekampRabin в общем случае.")

    def test_Fermat_little_theorem(self):
        self.assertEqual(Polynomial(GF(3))().Fermat_little_theorem(), Polynomial(GF(3))(), "Нулевой многочлен после примениния Малой теоремы Ферма должен оставаться нулевым.")
        self.assertEqual(Polynomial(GF(3))(5).Fermat_little_theorem(), Polynomial(GF(3))(2), "Константа после примениния Малой теоремы Ферма должена оставаться константой.")
        self.assertEqual(Polynomial(GF(5))("6x^8 - 3x^7 + 5x^6 + x^5 + x^3").Fermat_little_theorem(), Polynomial(GF(5))("x^4 + 3x^3 + x"), "Ошибка примениния Малой теоремы Ферма в общем случае.")
    
    def test_divmod(self):
        self.assertEqual(divmod(Polynomial(GF(3))(), Polynomial(GF(3))("x^2 + x + 1")), (Polynomial(GF(3))(), Polynomial(GF(3))()), "Деление нулевого многочлена должно давать нулевой многочлен.")
        with self.assertRaises(ZeroDivisionError, msg="При попытке деления на нулевой многочлен должно подниматься исключение."): divmod(Polynomial(GF(3))("x^2 + x + 1"), Polynomial(GF(3))())
        with self.assertRaises(ZeroDivisionError, msg="При попытке деления на нулевой многочлен должно подниматься исключение (случай преобразования элемента в многочлен)."): divmod(Polynomial(GF(3))("x^2 + x + 1"), 3)
        self.assertEqual(divmod(Polynomial(GF(7))("x^2 + x + 1"), Polynomial(GF(7))("x^2 + x + 1")), (Polynomial(GF(7))(1), Polynomial(GF(7))()), "При делении многочлена на самого себя должна получаться единица.")
        self.assertEqual(divmod(Polynomial(GF(7))("2x^3 + 5x + 1"), 2), (Polynomial(GF(7))("x^3 + 6x + 4"), Polynomial(GF(7))()), "Ошибка при делении многочлена на константу.")
        self.assertEqual(divmod(Polynomial(GF(5))("x^2 + 2x + 1"), Polynomial(GF(5))("x^3 + x + 1")), (Polynomial(GF(5))(), Polynomial(GF(5))("x^2 + 2x + 1")), "Если степень делителя выше степени делимого, результатом должно быть частное 0 и остаток, равный самому делимому.")
        self.assertEqual(divmod(Polynomial(GF(5))("10x^2 + 2"), Polynomial(GF(5))("3x")), (Polynomial(GF(5))(), Polynomial(GF(5))("2")), "Если в делимом есть коэффициенты, кратные модулю, они должны быть занулены перед делением.")
        with self.assertRaises(ZeroDivisionError, msg="Если все коэффициенты делителя кратны модулю, должно подыматься исключение деления на нуль"): divmod(Polynomial(GF(7))("x^5 + x + 1"), Polynomial(GF(7))("7x^2"))
        self.assertEqual(divmod(Polynomial(GF(7))("3x^3 + x^2 + x + 5"), Polynomial(GF(7))("x^2 + x + 1")), (Polynomial(GF(7))("3x + 5"), Polynomial(GF(7))()), "Ошибка при делении нацело.")
        self.assertEqual(divmod(Polynomial(GF(7))("3x^3 + x^2 + 2x + 6"), Polynomial(GF(7))("x^2 + x + 1")), (Polynomial(GF(7))("3x + 5"), Polynomial(GF(7))("x + 1")), "Ошибка при делении в общем случае.")
        with self.assertRaises(TypeError, msg="divmod не поднимает исключение в случае многочлена над над Z."): divmod(Polynomial(IntegerRing())("x^2 - 1"), Polynomial(IntegerRing())("x - 1"))
    
    def test_derivative(self):
        ring =  Polynomial(IntegerRing())
        self.assertEqual(ring().derivative(), ring(), "Производная нулевого многочлена должна быть нулевой.")
        self.assertEqual(ring(5).derivative(), ring(), "Производная константы должна быть нулевой.")
        self.assertEqual(ring("3x^3 + x^2 + x + 5").derivative(), ring("9x^2 + 2x + 1"), "Ошибка при определнии производной в общем случае.")

    def test_roots_prime_module(self):
        gf = GF(5)
        self.assertEqual(Polynomial(gf)()._Element__roots_prime_module(), set(gf), "Корнями нулевого многочлена должны являться все элементы поля.")
        self.assertEqual(Polynomial(gf)("5x^2 + 5x + 5")._Element__roots_prime_module(), set(gf), "Корнями нулевого многочлена должны являться все элементы поля (случай нулевого по модулю).")
        self.assertEqual(Polynomial(gf)("x - 3")._Element__roots_prime_module(), {gf(3)}, "Ошибка __roots_prime_module в случае линейного с одним корнем.")
        self.assertEqual(Polynomial(gf)(3)._Element__roots_prime_module(), set(), "Ошибка __roots_prime_module в случае константного без корней.")
        self.assertEqual(Polynomial(GF(7))("x^2 - 1")._Element__roots_prime_module(), {GF(7)(1), GF(7)(6)}, "Ошибка __roots_prime_module в случае квадратного с двумя корнями.")
        self.assertEqual(Polynomial(GF(13))("x^2 - 5")._Element__roots_prime_module(), set(), "Ошибка __roots_prime_module в случае квадратного без корней.")
        self.assertEqual(Polynomial(GF(13))("x^2 - 26")._Element__roots_prime_module(), {GF(13)(0)}, "Ошибка __roots_prime_module в случае квадратного с одниим корнем.")
        self.assertEqual(Polynomial(gf)("3x^3 + 4x^2 + 2x + 1")._Element__roots_prime_module(), {gf(1), gf(2), gf(4)}, "Ошибка __roots_prime_module в общем случае.")

    def test_roots_primary_module(self):
        ring = IntegerModRing(625)
        self.assertEqual(Polynomial(ring)("3x^3 + 4x^2 + 2x + 1")._Element__roots_primary_module(5, 4), {ring(72), ring(136), ring(624)}, "Ошибка _roots_primary_module в общем случае.")
        gf = GF(5)
        self.assertEqual(Polynomial(gf)("3x^3 + 4x^2 + 2x + 1")._Element__roots_primary_module(5), {gf(1), gf(2), gf(4)}, "Ошибка _roots_primary_module случае q = 1.")
        ring = IntegerModRing(28561)
        self.assertEqual(Polynomial(ring)("x^2 - 5")._Element__roots_primary_module(13, 4), set(), "Ошибка __roots_primary_module в случае отсуствия корней по простому модулю.")

    def test_roots(self):
        ring = IntegerModRing(189)
        self.assertEqual(Polynomial(ring)("162x^10+135x^9+162x^8+56x^4+162x^3+162x^2+113x+188").roots(), {ring(94), ring(67), ring(121)}, "Ошибка roots в кольце.")
        ring = IntegerModRing(16)
        self.assertEqual(Polynomial(ring)("x^3+11x^2+2x+8").roots(), {ring(2), ring(4), ring(10), ring(12), ring(15)}, "Ошибка roots в кольце в случае примарного модуля.")
        self.assertAlmostEqual(Polynomial(RealField())("x^2 - 1").roots(), {np.complex64(-1), np.complex64(1)}, "Ошибка roots в вешественном поле.")
        with self.assertRaises(TypeError, msg="roots не поднимает исключение в случае многочлена над над Z."): Polynomial(IntegerRing())("x^2 - 1").roots()

    def test_gcd(self):
        gf = GF(5)
        self.assertEqual(Polynomial.gcd(Polynomial(gf)("2x^3 - 3x + 1"), Polynomial(gf)()), Polynomial(gf)("x^3 + x + 3"), "Ошибка gcd в случае, если один из многочленов нулевой.")
        self.assertEqual(Polynomial.gcd(Polynomial(gf)("2x^3 - 3x + 1"), Polynomial(gf)("2x^3 - 3x + 1")), Polynomial(gf)("x^3 + x + 3"), "Ошибка gcd в случае одинаковых многочленов.")
        self.assertEqual(Polynomial.gcd(Polynomial(gf)("2x^3 - 3x + 1"), Polynomial(gf)("x^3 + x + 3")), Polynomial(gf)("x^3 + x + 3"), "Ошибка gcd в случае одинаковых унитарного и не унитарного многочленов.")
        self.assertEqual(Polynomial.gcd(Polynomial(gf)("2x - 2"), Polynomial(gf)("x^2 - 1")), Polynomial(gf)("x + 4"), "Ошибка gcd в случае если один делитель друого.")
        self.assertEqual(Polynomial.gcd(Polynomial(gf)("x^3 + x"), Polynomial(gf)("x^3 + 3x^2 + 3x")), Polynomial(gf)("x"), "Ошибка gcd в общем случае.")
        gf = GF(13)
        self.assertEqual(Polynomial.gcd(Polynomial(gf)("3x^4 + 6x^2 - 1"), Polynomial(gf)("x^2 - 7")), Polynomial(gf)(1), "Ошибка gcd в случае взаимнопростых.")
        self.assertEqual(Polynomial.gcd(Polynomial(gf)(), Polynomial(gf)()), Polynomial(gf)(), "Ошибка gcd в случае, если оба многочлена нулевые.")
        self.assertEqual(Polynomial.gcd(Polynomial(gf)("3x^4 + 6x^2 - 1"), Polynomial(gf)(5)), Polynomial(gf)(1), "Ошибка gcd в случае, если один многочлен констаната.")

    def test_binomial_theorem(self):
        ring = Polynomial(IntegerRing())
        self.assertEqual(ring.binomial_theorem(2, 3, 0), ring(1), "Ошибка binomial_theorem в случае exp=0.")
        self.assertEqual(ring.binomial_theorem(2, 3, 1), ring('2x+3'), "Ошибка binomial_theorem в случае exp=1.")
        self.assertEqual(ring.binomial_theorem(0, 3, 2), ring(9), "Ошибка binomial_theorem в случае m=0.")
        self.assertEqual(ring.binomial_theorem(2, 0, 3), ring('8x^3'), "Ошибка binomial_theorem в случае a=0.")
        self.assertEqual(ring.binomial_theorem(1, 1, 5), ring('x^5 + 5x^4 + 10x^3 + 10x^2 + 5x + 1'), "Ошибка binomial_theorem в общем случае.")
        self.assertEqual(ring.binomial_theorem(-1, 2, 3), ring('-x^3 + 6x^2 - 12x + 8'), "Ошибка binomial_theorem случае m<0.")
        self.assertEqual(ring.binomial_theorem(2, -3, 2), ring('4x^2 - 12x + 9'), "Ошибка binomial_theorem случае a<0.")
        self.assertEqual(ring.binomial_theorem(0, 0, 5), ring(), "Ошибка binomial_theorem случае m=0, a=0." )

    def test_monic(self):
        self.assertEqual(Polynomial(GF(5))().monic(), Polynomial(GF(5))(), "Ошибка monic в случае нулевого многочлена.")
        self.assertEqual(Polynomial(GF(7))("5x^3 - 3x^2 + 4x - 1").monic(), Polynomial(GF(7))("x^3 + 5x^2 + 5x + 4"), "Ошибка monic над конечным полем.")

if __name__ == '__main__':
    unittest.main()