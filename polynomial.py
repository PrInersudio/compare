"""Реализует кольцо многочленов."""

from __future__ import annotations

from itertools import product
from typing import Tuple, Dict, Set, Generator
import logging
import re
from math import comb
from sympy import factorint
import numpy as np

from linear_comparisons import solve_linear_comparison, crt
from quadratic_comparisons import solve_quadratic_comparison
from ring import Ring
from field import Field, GF, RealField, ComplexField
from integer_mod_ring import IntegerModRing


logger = logging.getLogger(__name__)
logger.setLevel(logging.DEBUG)
handler = logging.FileHandler(f'log/{__name__}.log')
handler.setLevel(logging.DEBUG)
handler.setFormatter(logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s'))
logger.addHandler(handler)


class Polynomial(Ring):
    """Кольцо многочленов над другим кольцом."""

    def __init__(self, basic_ring: Ring):
        super().__init__()
        self.basic_ring = basic_ring
        self.has_multiplicative_identity = self.basic_ring.has_multiplicative_identity
        self.is_euclidean = isinstance(self.basic_ring, Field)
        self.is_commutative = self.basic_ring.is_commutative
        logger.debug(
            'Polynomial.__init__. self.basic_ring = %s,'
            ' self.has_multiplicative_identity = %s,'
            ' self.is_euclidean = %s,'
            ' self.is_commutative = %s.',
            self.basic_ring, self.has_multiplicative_identity,
            self.is_euclidean, self.is_commutative,
        )

    def __repr__(self) -> str:
        return f'Polynomial({repr(self.basic_ring)})'

    def __str__(self) -> str:
        return f'Кольцо многочленов над {self.basic_ring}'

    def __call__(self, raw_poly: None | Dict[int, Ring.Element] | Ring.Element | str = None):
        return self.Element(self, raw_poly)

    def __eq__(self, other: 'Polynomial') -> bool:
        if not isinstance(other, Polynomial):
            return False
        return self.basic_ring == other.basic_ring

    @staticmethod
    def gcd(f: Polynomial.Element, g: Polynomial.Element) -> Polynomial.Element:
        """НОД многочленов."""
        logger.debug('Начало Polynomial.gcd. f = %s, g = %s.', f, g)
        while g:
            logger.debug('Шаг Polynomial.gcd. f = %s, g = %s.', f, g)
            f, g = g, f % g
        logger.debug('Polynomial.gcd, делаем унитарный. f = %s, g = %s.', f, g)
        f = f.monic()
        logger.debug('Конец Polynomial.gcd. f = %s, g = %s.', f, g)
        return f

    def binomial_theorem(self, m: Ring.Element, a: Ring.Element, exp: int) -> Polynomial.Element:
        """
        Разложение (mx+a)^exp
        """
        logger.debug('Polynomial.binomial_theorem.'
                     ' self = %s, m = %s, a = %s, exp = %s.', self, m, a, exp)
        poly = self()
        for k in range(exp + 1):
            coeff = comb(exp, k) * (m ** (exp - k)) * (a ** k)
            poly[exp - k] = coeff
            logger.debug('Polynomial.binomial_theorem.'
                         ' self = %s, m = %s, a = %s, exp = %s, poly = %s.', self, m, a, exp, poly)
        return poly


    class Element(Ring.Element):
        """Многочлен из кольца многочленов."""

        def __init__(
                self, ring: Polynomial,
                raw_poly: None | Dict[int, Ring.Element] | Ring.Element | str = None
        ):
            super().__init__(ring)
            if isinstance(self.ring.basic_ring, GF):
                self.fermat_little_theorem = self.__fermat_little_theorem
            if self.ring.basic_ring.has_multiplicative_inverses:
                self.monic = self.__monic
            logger.debug('Polynomial.Element.__init__.'
                         ' ring = %s, raw_poly = %s.', self.ring, raw_poly)
            if not raw_poly:
                self.value = dict()
                logger.debug('Polynomial.Element.__init__ завершается.'
                             ' self = %s, raw_poly = %s.', self, raw_poly)
                return
            elif isinstance(raw_poly, Polynomial.Element):
                self.value = {
                    exp:self.ring.basic_ring(coeff)
                    for exp, coeff in raw_poly.monomials()
                }
                logger.debug('Polynomial.Element.__init__ завершается. self = %s, raw_poly = %s.',
                             self, raw_poly)
                return
            elif isinstance(raw_poly, dict):
                logger.debug('Polynomial.Element.__init__. type(list(raw_poly.values())[0]) = %s.',
                             type(list(raw_poly.values())[0]))
                self.value = {
                    exp:coeff
                    for exp, coeff in raw_poly.items()
                    if isinstance(exp, int) and exp >= 0
                    and coeff and isinstance(coeff, Ring.Element)
                    and coeff.ring == self.ring.basic_ring
                }
                logger.debug('Polynomial.Element.__init__ завершается. self = %s, raw_poly = %s.',
                             self, raw_poly)
                return
            elif not isinstance(raw_poly, str):
                try:
                    self.value = {0:self.ring.basic_ring(raw_poly)} if raw_poly else dict()
                    logger.debug('Polynomial.Element.__init__ завершается.'
                                 ' self = %s, raw_poly = %s.', self, raw_poly)
                    return
                except Exception as exc:
                    logger.critical(
                        'Polynomial.Element.__init__. elf = %s, raw_poly = %s.'
                        ' Неправильный тип raw_poly.',
                        self, raw_poly
                    )
                    raise TypeError("Полином создаётся только из строки,"
                                    "словаря, элемента кольца"
                                    "(при этом тип добавляемых элементов должен соответствовать"
                                    "типу коэффициентов многочлена).") from exc
            pattern = r'([+-]?[^-+]+)'
            monomials = re.findall(pattern, raw_poly)
            self.value = dict()
            logger.debug('Polynomial.Element.__init__. self = %s, raw_poly = %s, monomials = %s.',
                         self, raw_poly, monomials)
            for monomial in monomials:
                splited = monomial.replace(' ', '').split('x')
                coeff = (
                    self.ring.basic_ring(1) if splited[0] in ['', '+']
                    else self.ring.basic_ring(-1) if splited[0] == '-'
                    else self.ring.basic_ring(splited[0])
                )
                if not coeff:
                    continue
                exp = 0 if len(splited) != 2 else (1 if splited[1] == '' else int(splited[1][1:]))
                if exp < 0:
                    logger.critical("Polynomial.Element.__init__. Отрицательная степень.")
                    raise ValueError("Отрицательные степени не поддерживаются.")
                self[exp] += coeff
            logger.debug('Polynomial.Element.__init__ завершается. self = %s, raw_poly = %s.',
                         self, raw_poly)

        def __getitem__(self, exp: int) -> Ring.Element:
            return self.value.get(exp, self.ring.basic_ring())

        def __setitem__(self, exp: int, coeff: Ring.Element) -> None:
            coeff = self.ring.basic_ring(coeff)
            if not coeff:
                self.value.pop(exp,0)
            else: self.value[exp] = coeff

        def __contains__(self, exp: int) -> bool:
            return exp in self.value

        def __iter__(self) -> Generator[int]:
            for exp, coeff in self.value.items():
                if coeff:
                    yield exp

        def __eq__(self, other: 'Polynomial.Element') -> bool:
            if not isinstance(other, Polynomial.Element):
                return False
            if self.ring != other.ring:
                return False
            return all(self[exp] == other[exp] for exp in set(self) | set(other))

        def monomials(self)  -> Generator[Tuple[int,Ring.Element]]:
            """Генератор наборов (exp, coeff)."""
            for exp, coeff in self.value.items():
                if coeff:
                    yield exp, coeff

        def __hash__(self) -> int:
            return hash(tuple(sorted((exp, coeff) for exp, coeff in self.monomials())))

        def copy(self) -> 'Polynomial.Element':
            """Неглубокая копия многочлена."""
            return self.ring(self.value)

        def __add__(self, other: 'Polynomial.Element' | Ring.Element) -> 'Polynomial.Element':
            logger.debug('Polynomial.Element.__add__. self = %s, other = %s.', self, other)
            if not isinstance(other, Polynomial.Element):
                other = self.ring(other)
            coeffs = dict()
            for i in set(self) | set(other):
                new_coeff = self[i] + other[i]
                if new_coeff:
                    coeffs[i] = new_coeff
            logger.debug('Polynomial.Element.__add__. self = %s, other = %s, coeffs = %s.',
                         self, other, coeffs)
            return self.ring(coeffs)

        def __neg__(self):
            return self.ring({exp:-coeff for exp,coeff in self.monomials()})

        def __sub__(self, other: 'Polynomial.Element' | Ring.Element) -> 'Polynomial.Element':
            logger.debug('Polynomial.__sub__. self = %s, other = %s.', self, other)
            if not isinstance(other, Polynomial.Element):
                other = self.ring(other)
            coeffs = dict()
            for i in set(self) | set(other):
                new_coeff = self[i] - other[i]
                if new_coeff:
                    coeffs[i] = new_coeff
            logger.debug('Polynomial.Element.__sub__. self = %s, other = %s, coeffs = %s.',
                         self, other, coeffs)
            return self.ring(coeffs)

        def __mul__(self, other: 'Polynomial.Element' | Ring.Element) -> 'Polynomial.Element':
            logger.debug('Polynomial.Element.__mul__. self = %s, other = %s.', self, other)
            if not isinstance(other, Polynomial.Element):
                other = self.ring(other)
            result = self.ring()
            for exp1, coeff1 in self.monomials():
                for exp2, coeff2 in other.monomials():
                    result[exp1+exp2] += coeff1 * coeff2
            result = self.ring({exp:coeff for exp, coeff in result.monomials() if coeff})
            logger.debug('Polynomial.__mul__.'
                         ' self = %s, other = %s, result = %s.', self, other, result)
            return result

        def __rmul__(self, other) -> 'Ring.Element':
            logger.debug('Polynomial.Element.__rmul__. self = %s, other = %s.', self, other)
            if not isinstance(other, Polynomial.Element):
                other = self.ring(other)
            return other * self

        def __bool__(self) -> bool:
            return any(coeff for coeff in self.value.values())

        def __repr__(self) -> str:
            return f'Polynomial.Element({repr(self.ring)}, {self.value})'

        def __str__(self) -> str:
            if not self:
                return str(self.ring.basic_ring())
            monomials = []
            for exp in sorted(list(self), reverse=True):
                coeff = self[exp]
                if not coeff:
                    continue
                if exp == 0:
                    monomial = f"{coeff}"
                elif exp == 1:
                    monomial = f"{coeff}x"
                else: monomial = f"{coeff}x^{exp}"
                monomials.append(monomial)
            polynomial_str = " + ".join(monomials)
            polynomial_str = polynomial_str.replace("+ -", "- ")
            return polynomial_str

        def __gt__(self, other: 'Polynomial.Element' | Ring.Element) -> bool:
            if not isinstance(other, Polynomial.Element):
                other = self.ring(other)
            return self.degree() > other.degree()

        def __ge__(self, other: 'Polynomial.Element' | Ring.Element) -> bool:
            if not isinstance(other, Polynomial.Element):
                other = self.ring(other)
            return self.degree() >+ other.degree()

        def __lt__(self, other: 'Polynomial.Element' | Ring.Element) -> bool:
            if not isinstance(other, Polynomial.Element):
                other = self.ring(other)
            return self.degree() < other.degree()

        def __le__(self, other: 'Polynomial.Element' | Ring.Element) -> bool:
            if not isinstance(other, Polynomial.Element):
                other = self.ring(other)
            return self.degree() <= other.degree()

        def degree(self) -> int:
            """Степень многочлена."""
            logger.debug('Polynomial.Element.degree. self = %s.', self)
            if not self:
                logger.debug('Polynomial.Element.degree. self = %s. Нулевой многочлен.', self)
                return -1
            degree = max(exp for exp in self if self[exp])
            logger.debug('Polynomial.Element.degree. self = %s, degree = %s.', self, degree)
            return degree

        def __call__(self, x: Ring.Element) -> Ring.Element:
            """
            Значение многочлена в точке x
            """
            logger.debug('Polynomial.Element.__call__. self = %s.', self)
            if not isinstance(x, Ring.Element) or not x.ring == self.ring.basic_ring:
                x =  self.ring.basic_ring(x)
            value = sum(coeff * x ** exp for exp,coeff in self.monomials())
            logger.debug('Polynomial.Element.__call__. self = %s, value = %s.', self, value)
            return value

        def __fermat_little_theorem(self) -> 'Polynomial.Element':
            """
            Сокращение многочлена по Малой теореме Ферма в GF(p)
            x^(p) = x (mod p) <==> x^(kp+t) = x^(k+t) (mod p)
            """
            p = self.ring.basic_ring.p
            logger.debug('Polynomial.Element.fermat_little_theorem. self = %s.', self)
            poly = self.ring()
            for exp,coeff in self.monomials():
                while exp >= p:
                    exp = exp // p + exp % p
                poly[exp] += coeff
            logger.debug('Polynomial.Element.fermat_little_theorem.'
                         ' self = %s, poly = %s.', self, poly)
            return poly

        def __divmod__(
                self, divisor: 'Polynomial.Element' | Ring.Element
        ) -> Tuple['Polynomial.Element', 'Polynomial.Element']:
            logger.debug('Polynomial.Element.__divmod__. self = %s, divisor = %s.', self, divisor)
            if not self.ring.is_euclidean:
                raise TypeError(f'Деление не определено в {self.ring}')
            if not isinstance(divisor, Polynomial.Element):
                divisor = self.ring(divisor)
            if not divisor:
                logger.critical('Polynomial.Element.__divmod__. self = %s, divisor = %s.'
                                ' Нулевой divisor.', self, divisor)
                raise ZeroDivisionError("Деление на нулевой многочлен.")
            remainder = self.copy()
            remainder_degree = remainder.degree()
            quotient = self.ring()
            divisor_degree = divisor.degree()
            while remainder and remainder_degree >= divisor_degree:
                logger.debug('Polynomial.Element.__divmod__. self = %s, divisor = %s,'
                             ' quotient = %s, remainder = %s.',
                             self, divisor, quotient, remainder)
                degree_diff = remainder_degree - divisor_degree
                leading_coeff = remainder[remainder_degree] / divisor[divisor_degree]
                quotient[degree_diff] = leading_coeff
                for exp, coeff in divisor.monomials():
                    remainder[exp + degree_diff] = (remainder[exp + degree_diff]
                                                    - leading_coeff * coeff)
                remainder = self.ring({exp: coeff for exp, coeff in remainder.monomials() if coeff})
                remainder_degree = remainder.degree()
            logger.debug('Polynomial.Element.__divmod__. self = %s, divisor = %s, quotient = %s,'
                         ' remainder = %s.', self, divisor, quotient, remainder)
            return quotient, remainder

        def __truediv__(self, divisor: 'Polynomial.Element' | Ring.Element) -> 'Polynomial.Element':
            return divmod(self, divisor)[0]

        def __mod__(self, divisor: 'Polynomial.Element' | Ring.Element) -> 'Polynomial.Element':
            return divmod(self, divisor)[1]

        def _roots_berlekamp_rabin(self) -> Set[GF.Element]:
            p = self.ring.basic_ring.p
            logger.debug('Polynomial.Element._roots_berlekamp_rabin. self = %s, p = %s.', self, p)
            roots = set()
            nonlinears = [self]
            for delta in range(p):
                if not nonlinears:
                    break
                poly = nonlinears.pop()
                logger.debug('Polynomial.Element._roots_berlekamp_rabin. self = %s, delta = %s,'
                             ' roots = %s, nonlinears = %s, poly = %s.',
                             self, delta, roots, nonlinears, poly)
                d = Polynomial.gcd(poly, self.ring.binomial_theorem(1,delta,(p-1) >> 1) - 1)
                if d.degree() <= 0:
                    nonlinears.append(poly)
                    continue
                for factor in [d, poly / d]:
                    if factor.degree() == 1:
                        roots |= solve_linear_comparison(int(factor[1]), int(-factor[0]), p)
                    else: nonlinears.append(factor)
            logger.debug('Polynomial.Element._roots_berlekamp_rabin.'
                         ' self = %s, p = %s, roots = %s.', self, p, roots)
            return roots

        def derivative(self) -> 'Polynomial.Element':
            """Производная многочлена."""
            logger.debug('Polynomial.Element.derivative. self = %s.', self)
            derivative = self.ring({exp-1:coeff*exp for exp,coeff in self.monomials() if exp > 0})
            logger.debug('Polynomial.Element.derivative.'
                         ' self = %s, derivative = %s.', self, derivative)
            return derivative

        def _roots_prime_module(self) -> Set[GF.Element]:
            p = self.ring.basic_ring.p
            logger.debug('Polynomial.Element._roots_prime_module. self = %s, p = %s.', self, p)
            fermat = self.fermat_little_theorem()
            if fermat.degree() <= 1:
                logger.debug('Polynomial.Element._roots_prime_module.'
                             ' self = %s, p = %s, fermat = %s. Степень <= 1.', self, p, fermat)
                return solve_linear_comparison(int(fermat[1]), int(-fermat[0]), p)
            reduced = Polynomial.gcd(fermat, self.ring({p:1,1:-1}))
            if reduced.degree() == 2:
                logger.debug('Polynomial.Element._roots_prime_module.'
                             ' self = %s, p = %s, reduced = %s. Степень 2.', self, p, reduced)
                return solve_quadratic_comparison(
                    int(reduced[2]),
                    int(reduced[1]),
                    int(reduced[0]), p,
                )
            return reduced._roots_berlekamp_rabin() # pylint: disable=protected-access

        def _roots_primary_module(self, p: int, q: int = 1) -> Set[IntegerModRing.Element]:
            logger.debug('Polynomial.Element._roots_primary_module.'
                         ' self = %s, p = %s, q = %s.', self, p, q)
            gfp = GF(p)
            result = {int(a) for a in Polynomial(gfp)(self)._roots_prime_module()} # pylint: disable=protected-access
            derivative = Polynomial(gfp)(self.derivative()).fermat_little_theorem()
            p_i = 1
            for _ in range(1,q):
                p_i *= p
                result = {
                    a: solve_linear_comparison(int(derivative(a)), - int(self(a)) // p_i, p)
                    for a in result
                }
                logger.debug('Polynomial.Element._roots_primary_module.'
                             ' self = %s, p = %s, q = %s, result = %s.', self, p, q, result)
                result = set(a + int(t) * p_i for a, temp in result.items() for t in temp)
            ring = IntegerModRing(p_i * p)
            logger.debug('Polynomial.Element._roots_primary_module.'
                         ' self = %s, p = %s, q = %s, result = %s.', self, p, q, result)
            return set(ring(a) for a in result)

        def roots(self):
            """Корни многочлена в соответствующем кольце."""
            logger.debug('Polynomial.Element.roots. self = %s.', self)
            if isinstance(self.ring.basic_ring, RealField):
                return set(np.roots([float(self[i]) for i in range(self.degree() + 1)]))
            if isinstance(self.ring.basic_ring, ComplexField):
                return set(np.roots([complex(self[i]) for i in range(self.degree() + 1)]))
            if isinstance(self.ring.basic_ring, IntegerModRing):
                module = self.ring.basic_ring.m
                if self.degree() <= 1:
                    logger.debug('Polynomial.Element.roots. self = %s, module = %s. Степень <= 1.',
                                 self, module)
                    return set(solve_linear_comparison(int(self[1]),int(-self[0]), module))
                primary_results = []
                for p,q in factorint(module).items():
                    primary_results.append(
                        Polynomial(IntegerModRing(p**q))(self)._roots_primary_module(p,q) # pylint: disable=protected-access
                    )
                if len(primary_results) == 1:
                    logger.debug('Polynomial.Element.roots.'
                                 ' self = %s, module = %s, primary_results = %s.'
                                 ' module - примарный.', self, module, primary_results)
                    return primary_results[0]
                result = set()
                for element in product(*primary_results):
                    result.add(crt(element))
                logger.debug('Polynomial.Element.roots. self = %s, module = %s, result = %s.',
                             self, module, result)
                return result
            raise TypeError(f'Поиск корней в {self.ring.basic_ring} не поддерживается.')

        def __monic(self) -> 'Polynomial.Element':
            if not self:
                return self.ring()
            return self * (self[self.degree()] ** (-1))
