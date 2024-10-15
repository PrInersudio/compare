from __future__ import annotations
import sys
import re
from math import gcd as int_gcd
from math import comb
from sympy import factorint
from itertools import product
from typing import List, Tuple, Dict, Set
import logging

logging.basicConfig(filename='compare.log', level=logging.DEBUG, format=' %(asctime)s - %(levelname)s - %(message)s')

def solve_linear_comparison(a: int,b: int,m: int) -> List[Tuple[int,int]]:
    """
    ax = b (mod m)
    """
    logging.debug(f'solve_linear_comparison. a={a}, b={b}, m={m}')
    d = int_gcd(a,m)
    logging.debug(f'solve_linear_comparison. a={a}, b={b}, m={m}, d={d}')
    if b % d != 0:
        logging.debug(f'solve_linear_comparison. a={a}, b={b}, m={m}. Нет решений.')
        return []
    a0,b0,m0 = a//d,b//d,m//d
    logging.debug(f'solve_linear_comparison. a={a}, b={b}, m={m}, a0={a0}, b0={b0}, m0={m0}')
    x0 = b0 * pow(a0,-1,m0)
    logging.debug(f'solve_linear_comparison. a={a}, b={b}, m={m}, d={d}, a0={a0}, b0={b0}, m0={m0}, x0={x0}')
    result = [((x0 + i) % m, m) for i in range(0,m,m0)]
    logging.debug(f'solve_linear_comparison. a={a}, b={b}, m={m}, result={result}')
    return result

def crt (equations: List[Tuple[int,int]]) -> Tuple[int,int]:
    """
    x = a_i (mod m_i)
    """
    logging.debug(f'crt. equations={equations}')
    M = 1
    x = 0
    for a,m in equations:
        x += M*((a-x)*pow(M,-1,m))
        M *= m
        logging.debug(f'crt. equations={equations}, a={a}, m={m}, x={x}, M={M}')
    result = (x % M,M)
    logging.debug(f'crt. equations={equations}, result={result}')
    return result
class Polynomial:
    """
    Полином с целыми коэффициентами.
    """
    def __getitem__(self, exp: int) -> int:
        coeff = self.__coefficients.get(exp, 0)
        return coeff
    
    def __setitem__(self, exp: int, coeff: int) -> None:
        self.__coefficients[exp] = coeff

    def __contains__(self, exp: int) -> bool:
        contains = exp in self.__coefficients
        return contains
    
    def __iter__(self):
        it = iter(self.__coefficients)
        return it
    
    def monomials(self):
        """
        Возвращает наборы (exp, coeff).
        """
        monomials = self.__coefficients.items()
        return monomials

    def __init__(self, raw_poly: None | Dict[int,int] | str = None):
        logging.debug(f'Polynomial.__init__. raw_poly={raw_poly}')
        if not raw_poly:
            self.__coefficients = dict()
            return
        elif isinstance(raw_poly, dict):
            self.__coefficients = raw_poly
            return
        elif not isinstance(raw_poly, str):
            logging.critical(f'Polynomial.__init__. raw_poly={raw_poly}. Неправильный тип raw_poly.')
            raise TypeError("Полином создаётся только из строки или словаря")
        pattern = r'([+-]?[^-+]+)'
        monomials = re.findall(pattern, raw_poly)
        self.__coefficients = dict()
        logging.debug(f'Polynomial.__init__. self.__coefficients={self.__coefficients}, raw_poly={raw_poly}, monomials={monomials}')
        for monomial in monomials:
            if monomial == '0': continue
            splited = monomial.replace(' ', '').split('x')
            exp = 0 if len(splited) != 2 else (1 if splited[1] == '' else int(splited[1][1:]))
            self[exp] = 1 if splited[0] == '' else int(splited[0])
        logging.debug(f'Polynomial.__init__ завершается. self.__coefficients={self.__coefficients}, raw_poly={raw_poly}')

    def copy(self) -> 'Polynomial':
        return Polynomial(self.__coefficients.copy())

    def __add__(self, other: 'Polynomial' | int) -> 'Polynomial':
        logging.debug(f'Polynomial.__add__. self.__coefficients={self.__coefficients}, other={other}')
        if isinstance(other, int):
            poly = self.copy()
            poly[0] = poly[0] + other
            logging.debug(f'Polynomial.__add__. self.__coefficients={self.__coefficients}, other={other}, poly={poly}')
            return poly
        if not isinstance(other, Polynomial):
            logging.critical(f'Polynomial.__add__. self.__coefficients={self.__coefficients}, other={other}. Неправильный тип other.')
            raise TypeError("Сложение возмозжно только с числом и другим полиномом")
        coeffs = dict()
        for i in set(self) | set(other):
            coeffs[i] = self[i] + other[i]
        logging.debug(f'Polynomial.__add__. self.__coefficients={self.__coefficients}, other={other}, coeffs={coeffs}')
        return Polynomial(coeffs)
    
    def __iadd__(self, other: 'Polynomial' | int) -> 'Polynomial':
        logging.debug(f'Polynomial.__iadd__. self.__coefficients={self.__coefficients}, other={other}')
        return self + other
    
    def __sub__(self, other: 'Polynomial' | int) -> 'Polynomial':
        logging.debug(f'Polynomial.__sub__. self.__coefficients={self.__coefficients}, other={other}')
        if isinstance(other, int):
            poly = self.copy()
            poly[0] = poly[0] - other
            logging.debug(f'Polynomial.__sub__. self.__coefficients={self.__coefficients}, other={other}, poly={poly}')
            return poly
        if not isinstance(other, Polynomial):
            logging.critical(f'Polynomial.__sub__. self.__coefficients={self.__coefficients}, other={other}. Неправильный тип other.')
            raise TypeError("Вычитание возмозжно только с числом и другим полиномом")
        coeffs = dict()
        for i in set(self) | set(other):
            coeffs[i] = self[i] - other[i]
        logging.debug(f'Polynomial.__sub__. self.__coefficients={self.__coefficients}, other={other}, coeffs={coeffs}')
        return Polynomial(coeffs)
    
    def __isub__(self, other: 'Polynomial' | int) -> 'Polynomial':
        logging.debug(f'Polynomial.__isub__. self.__coefficients={self.__coefficients}, other={other}')
        return self - other
    
    def __mul__(self, scalar: int) -> 'Polynomial':
        logging.debug(f'Polynomial.__mul__. self.__coefficients={self.__coefficients}, scalar={scalar}')
        if not isinstance(scalar, int):
            logging.critical(f'Polynomial.__mul__. self.__coefficients={self.__coefficients}, scalar={scalar}. Неправильный тип scalar.')
            raise TypeError("Умножение возмозжно только с числом")
        result = Polynomial({exp:coeff*scalar for exp,coeff in self.monomials()})
        logging.debug(f'Polynomial.__mul__. self.__coefficients={self.__coefficients}, scalar={scalar}, result={result}')
        return result
    
    def __imul__(self, scalar: int) -> 'Polynomial':
        logging.debug(f'Polynomial.__imul__. self.__coefficients={self.__coefficients}, scalar={scalar}')
        return self * scalar
    
    def __rmul__(self, scalar: int) -> 'Polynomial':
        logging.debug(f'Polynomial.__rmul__. self.__coefficients={self.__coefficients}, scalar={scalar}')
        return self * scalar
    
    def __bool__(self) -> bool:
        result = any(coeff != 0 for coeff in self.__coefficients.values())
        return result
    
    def exponents(self):
        """
        Возвращает иттерабельный объкет с экспонентами многочлена.
        """
        exponents = self.__coefficients.keys()
        return exponents
    
    def __repr__(self) -> str:
        if not self:
            return "0"
        monomials = []
        for exp in sorted(self.exponents(), reverse=True):
            coeff = self[exp]
            if coeff == 0:
                continue
            if exp == 0:
                monomial = f"{coeff}"
            elif exp == 1:
                monomial = f"{coeff}x"
            else:
                monomial = f"{coeff}x^{exp}"
            monomials.append(monomial)
        polynomial_str = " + ".join(monomials)
        polynomial_str = polynomial_str.replace("+ -", "- ")
        return polynomial_str
    
    def __str__(self) -> str:
        return self.__repr__()
    
    def degree(self) -> int:
        logging.debug(f'Polynomial.degree. self.__coefficients={self.__coefficients}')
        if not self:
            logging.debug(f'Polynomial.degree. self.__coefficients={self.__coefficients}. Нулевой многочлен')
            return 0
        degree = max(exp for exp in self if self[exp] != 0)
        logging.debug(f'Polynomial.degree. self.__coefficients={self.__coefficients}, degree={degree}')
        return degree
    
    def __mod__(self,module: int) -> 'Polynomial':
        logging.debug(f'Polynomial.__mod__. self.__coefficients={self.__coefficients}')
        poly = Polynomial({exp:coeff%module for exp,coeff in self.monomials() if coeff%module != 0})
        logging.debug(f'Polynomial.__mod__. self.__coefficients={self.__coefficients}, poly={poly}')
        return poly
    
    def value(self, x: int) -> int:
        """
        Значение многочлена в точке x
        """
        logging.debug(f'Polynomial.value. self.__coefficients={self.__coefficients}')
        value = sum(coeff * x ** exp for exp,coeff in self.monomials())
        logging.debug(f'Polynomial.value. self.__coefficients={self.__coefficients}, value={value}')
        return value
    
    def Fermat_little_theorem(self, p: int) -> 'Polynomial':
        """
        Сокращение многочлена по Малой теореме Ферма в GF(p)
        x^(p) = x (mod p) <==> x^(kp+t) = x^(k+t) (mod p)
        """
        logging.debug(f'Polynomial.Fermat_little_theorem. self.__coefficients={self.__coefficients}, p={p}')
        poly = Polynomial()
        for exp,coeff in self.monomials():
            while exp >= p: exp = exp // p + exp % p
            poly[exp] += coeff
        logging.debug(f'Polynomial.Fermat_little_theorem. self.__coefficients={self.__coefficients}, p={p}, poly={poly}')
        return poly % p
    
    def div(self, divisor, p: int) -> 'Polynomial':
        """
        Деление многочленов в GF(p)
        """
        logging.debug(f'Polynomial.div. self.__coefficients={self.__coefficients}, divisor={divisor}, p={p}')
        if not isinstance(divisor, Polynomial):
            logging.critical(f'Polynomial.div. self.__coefficients={self.__coefficients}, divisor={divisor}, p={p}. Неправильный тип divisor.')
            raise TypeError("Делитель должен быть экземпляром Polynomial.")
        if not divisor:
            logging.critical(f'Polynomial.div. self.__coefficients={self.__coefficients}, divisor={divisor}, p={p}. Нулевой divisor.')
            raise ZeroDivisionError("Деление на нулевой многочлен.")
        remainder = self.copy()
        remainder_degree = remainder.degree()
        quotient = Polynomial()
        divisor_degree = divisor.degree()
        while remainder and remainder_degree >= divisor_degree:
            logging.debug(f'Polynomial.div. self.__coefficients={self.__coefficients}, divisor={divisor}, p={p}, quotient={quotient}, remainder={remainder}')
            degree_diff = remainder_degree - divisor_degree
            leading_coeff = solve_linear_comparison(divisor[divisor_degree], remainder[remainder_degree], p)[0][0]
            quotient[degree_diff] = leading_coeff
            for exp, coeff in divisor.monomials():
                remainder[exp + degree_diff] = (remainder[exp + degree_diff] - leading_coeff * coeff) % p
            remainder = Polynomial({exp: coeff for exp, coeff in remainder.monomials() if coeff != 0})
            remainder_degree = remainder.degree()
        logging.debug(f'Polynomial.div. self.__coefficients={self.__coefficients}, divisor={divisor}, p={p}, quotient={quotient}, remainder={remainder}')
        return quotient, remainder
    
    def factorize(self, p: int, start_delta: int = 0) -> Set['Polynomial']:
        """
        Вероятностная факторизация многочлена
        """
        logging.debug(f'Polynomial.factorize. self.__coefficients={self.__coefficients}, p={p}, start_delta={start_delta}')
        result = set()
        poly = self
        for delta in range(start_delta, p):
            logging.debug(f'Polynomial.factorize. self.__coefficients={self.__coefficients}, delta={delta}, poly={poly}, result={result}')
            if poly.degree() <= 1: break
            d = gcd(poly, binomial_theorem(1,delta,(p-1) // 2) - 1, p)
            logging.debug(f'Polynomial.factorize. self.__coefficients={self.__coefficients}, delta={delta}, poly={poly}, d={d}')
            if d.degree() == 0: continue
            quotient, remainder = poly.div(d,p)
            logging.debug(f'Polynomial.factorize. self.__coefficients={self.__coefficients}, delta={delta}, poly={poly}, d={d}, quotient={quotient}, remainder={remainder}')
            result |= quotient.factorize(p, delta + 1)
            poly = d
        result.add(poly)
        logging.debug(f'Polynomial.factorize. self.__coefficients={self.__coefficients}, result={result}')
        return result
    
    def derivative(self) -> 'Polynomial':
        logging.debug(f'Polynomial.derivative. self.__coefficients={self.__coefficients}')
        derivative = Polynomial({exp-1:coeff*exp for exp,coeff in self.monomials() if exp > 0})
        logging.debug(f'Polynomial.derivative. self.__coefficients={self.__coefficients}, derivative={derivative}')
        return derivative

    def __roots_prime_module(self, p: int) -> Set[int]:
        logging.debug(f'Polynomial.__roots_prime_module. self.__coefficients={self.__coefficients}, p={p}')
        fermat = self.Fermat_little_theorem(p)
        if fermat.degree() <= 1:
            logging.debug(f'Polynomial.__roots_prime_module. self.__coefficients={self.__coefficients}, p={p}. Степень <= 1.')
            return set(a for a,_ in solve_linear_comparison(fermat[1], -fermat[0], p))
        reduced = gcd(fermat, Polynomial({p:1,1:-1}),p)
        result = []
        for factor in reduced.factorize(p): result += set(solve_linear_comparison(factor[1],-factor[0],p))
        logging.debug(f'Polynomial.__roots_prime_module. self.__coefficients={self.__coefficients}, p={p}, result={result}')
        return set(a for a,_ in result)
    
    def __roots_primary_module(self, p: int, q: int = 1) -> Set[Tuple[int, int]]:
        logging.debug(f'Polynomial.__roots_primary_module. self.__coefficients={self.__coefficients}, p={p}, q={q}')
        result = self.__roots_prime_module(p)
        derivative = self.derivative().Fermat_little_theorem(p)
        for i in range(1,q):
            p_i = p**i
            result = {a:set(t for t,_ in solve_linear_comparison(derivative.value(a), - self.value(a) // p_i, p)) for a in result}
            logging.debug(f'Polynomial.__roots_primary_module. self.__coefficients={self.__coefficients}, p={p}, q={q}, result={result}')
            result = set(a + t * p_i for a, temp in result.items() for t  in temp)
        p_q = p**q
        logging.debug(f'Polynomial.__roots_primary_module. self.__coefficients={self.__coefficients}, p={p}, q={q}, result={result}')
        return set((a % p_q, p_q) for a in result)
    
    def roots(self, M: int) -> Set[Tuple[int, int]]:
        logging.debug(f'Polynomial.roots. self.__coefficients={self.__coefficients}, M={M}')
        if self.degree() <= 1:
            logging.debug(f'Polynomial.roots. self.__coefficients={self.__coefficients}, M={M}. Степень <= 1.')
            return set(solve_linear_comparison(self[1],-self[0],M))
        primary_results = []
        for p,q in factorint(M).items():
            primary_results.append(self.__roots_primary_module(p,q))
        if len(primary_results) == 1:
            logging.debug(f'Polynomial.roots. self.__coefficients={self.__coefficients}, M={M}, primary_results={primary_results}. M - примарное.')
            return primary_results[0]
        result = set()
        for element in product(*primary_results):
            result.add(crt(element))
        logging.debug(f'Polynomial.roots. self.__coefficients={self.__coefficients}, M={M}, result={result}')
        return result
    
def gcd(f: Polynomial, g: Polynomial, p:int) -> Polynomial:
    logging.debug(f'Начало gcd. f={f}, g={g}')
    while g:
        logging.debug(f'Шаг gcd. f={f}, g={g}')
        f, (_, g) = g, f.div(g,p)
    logging.debug(f'Конец gcd. f={f}, g={g}')
    return f

def binomial_theorem(m: int, a: int, exp: int) -> Polynomial:
    """
    Разложение (mx+a)^exp
    """
    logging.debug(f'binomial_theorem. m={m}, a={a}, exp={exp}')
    poly = Polynomial()
    for k in range(exp + 1):
        coeff = comb(exp, k) * (m ** (exp - k)) * (a ** k)
        poly[exp - k] = poly[exp - k, 0] + coeff
        logging.debug(f'binomial_theorem. m={m}, a={a}, exp={exp}, poly={poly}')
    return poly

def solve_comparisons (equations: List[Tuple[Polynomial, int]]) -> Set[Tuple[int,int]]:
    logging.debug(f'solve_comparisons. equations={equations}')
    result = equations[0][0].roots(equations[0][1])
    for poly, module in equations[1:]:
        logging.debug(f'solve_comparisons. equations={equations}, result={result}')
        new_result = set()
        for a,m in result:
            new_poly = Polynomial()
            for exp, coeff in poly.monomials():
                new_poly += coeff * binomial_theorem(m,a,exp)
                logging.debug(f'solve_comparisons. result={result}, new_poly={new_poly}')
            for a1,m1 in new_poly.roots(module):
                new_result.add((m*a1+a,m*m1))
        result = new_result
    logging.debug(f'solve_comparisons. equations={equations}, result={result}')
    return result    

def parse_equation(equation: str) -> Tuple[Polynomial, int]:
    logging.debug(f'parse_equation. equation={equation}')
    pattern = r'(.+?)\s*=\s*(.+?)\s*\(\s*mod\s*(\d+)\s*\)'
    match = re.match(pattern, equation)
    if not match:
        logging.critical(f'parse_equation. equation={equation}. Неверный формат уравнения.')
        raise ValueError("Неверный формат строки. Ожидается формат 'P(x) = Q(x) (mod m)'")
    module = int(match.group(3).strip())
    poly = Polynomial(match.group(1).strip()) - Polynomial(match.group(2).strip())
    logging.debug(f'parse_equation. equation={equation}, poly={poly}, module={module}')
    return poly, module

def main():
    logging.info("Начало программы.")
    if len(sys.argv) < 3:
        logging.error(f'Неправильное количество аргументов.')
        print("Запускать: python3", sys.argv[0], "<файл_с_сравнениями> <файл для вывода результата>")
        return
    logging.info(f'Считывание с файла {sys.argv[1]}.')
    with open(sys.argv[1]) as fp: equations_lines = fp.readlines()
    equations = []
    for equation in equations_lines:
        logging.info(f'Парсинг {equation}')
        if '#' in equation: continue
        equations.append(parse_equation(equation))
    if not equations:
        logging.info(f'Нет незакомментированных уравнений в {sys.argv[1]}.')
        print(f'Нет незакомментированных уравнений в {sys.argv[1]}.')
        return
    logging.info("Начинает решать.")
    result = solve_comparisons(equations)
    logging.info("Закончила решать. Начинает записывать ответ.")
    with open(sys.argv[2],"w") as fp:
        for x,m in result:
            print(f'[{x}]_{m}',end=' ', file=fp)
    logging.info("Записала ответ.")

if __name__ == '__main__':
    main()