"""Тесты quadratic_comparisons."""

import unittest

from quadratic_comparisons import (
    jacobi_symbol,
    QuadraticResidue,
    quadratic_comparison_to_quadratic_residue,
    solve_quadratic_comparison,
)
from field import GF


class TestQuadraticComparisons(unittest.TestCase):
    """Тесты quadratic_comparisons."""

    def test_jacobi_symbol(self):
        """Тест вычисления символа Якоби."""
        self.assertEqual(jacobi_symbol(30, 25), 0,
                         "Если a и m не взаимно просты, то символ Якоби должен быть равен нулю.")
        self.assertEqual(jacobi_symbol(5, 1), 1,
                         "Символ Якоби (a,1) должен быть равен 1 для любого a.")
        self.assertEqual(jacobi_symbol(-1, 7), -1,
                         "Ошибка определения символа Якоби в случае (-1/m)=-1")
        self.assertEqual(jacobi_symbol(-1, 13), 1,
                         "Ошибка определения символа Якоби в случае (-1/m)=1")
        self.assertEqual(jacobi_symbol(1, 5), 1,
                         "Ошибка определения символа Якоби в случае (1/m)=1")
        self.assertEqual(jacobi_symbol(4, 21), 1,
                         "Ошибка определения символа Якоби в случае (2^(2k)/m)=1 при m % 8 = 5")
        self.assertEqual(jacobi_symbol(4, 19), 1,
                         "Ошибка определения символа Якоби в случае (2^(2k)/m)=1 при m % 8 = 3")
        self.assertEqual(jacobi_symbol(4, 9), 1,
                         "Ошибка определения символа Якоби в случае (2^(2k)/m)=1 при m % 8 = 1")
        self.assertEqual(jacobi_symbol(4, 7), 1,
                         "Ошибка определения символа Якоби в случае (2^(2k)/m)=1 при m % 8 = 7")
        self.assertEqual(jacobi_symbol(8, 21), -1,
                         "Ошибка определения символа Якоби в случае (2^(2k+1)/m)=1 при m % 8 = 5")
        self.assertEqual(jacobi_symbol(8, 19), -1,
                         "Ошибка определения символа Якоби в случае (2^(2k+1)/m)=1 при m % 8 = 3")
        self.assertEqual(jacobi_symbol(8, 9), 1,
                         "Ошибка определения символа Якоби в случае (2^(2k+1)/m)=1 при m % 8 = 1")
        self.assertEqual(jacobi_symbol(8, 7), 1,
                         "Ошибка определения символа Якоби в случае (2^(2k+1)/m)=1 при m % 8 = 7")
        self.assertEqual(jacobi_symbol(3, 5), -1,
                         "Ошибка определения символа Якоби в случае (-1)^(q-1)(m-1)//4=1")
        self.assertEqual(jacobi_symbol(3, 7), -1,
                         "Ошибка определения символа Якоби в случае (-1)^(q-1)(m-1)//4=-1")
        self.assertEqual(jacobi_symbol(3, 7), -1,
                         "Ошибка определения символа Якоби в случае (-1)^(q-1)(m-1)//4=-1")
        self.assertEqual(jacobi_symbol(219, 233), 1,
                         "Ошибка определения символа Якоби в общем случае при простом m.")
        self.assertEqual(jacobi_symbol(219, 493), 1,
                         "Ошибка определения символа Якоби в общем случае.")

    def test_quadratic_comparison_to_quadratic_residue(self):
        """Тесты преобразования квадратичного сравнения в квадратичный вычет."""
        self.assertEqual(quadratic_comparison_to_quadratic_residue(7, 0, -8, 5),
                         (QuadraticResidue(4,5), 0),
                         "Ошибка quadratic_comparison_to_QuadraticResidue в случае b = 0.")
        self.assertEqual(quadratic_comparison_to_quadratic_residue(6, 9, 8, 5),
                         (QuadraticResidue(1,5), 2),
                         "Ошибка quadratic_comparison_to_QuadraticResidue в случае b != 0.")

    def test_roots(self):
        """Тесты вычисления корней квадратичного вычета."""
        self.assertEqual(QuadraticResidue(1, 2).roots(), {1},
                         "При p=2 корень квадратичного вычета должны быть равен a.")
        self.assertEqual(QuadraticResidue(4, 3).roots(), {1,2},
                         "При a = 1 (mod p) корень квадратичного вычета должны быть равен 1, -1.")
        self.assertEqual(QuadraticResidue(-1, 7).roots(), set(),
                         "При символе Якоби равном -1, корней нет.")
        self.assertEqual(QuadraticResidue(49, 7).roots(), {0},
                         "При символе Якоби равном 0, корень 0.")
        self.assertEqual(QuadraticResidue(4, 7).roots(), {2,5},
                         "Ошибка вычисления корней квадратичного вычета при p mod 4 = 3.")
        self.assertEqual(QuadraticResidue(4, 5).roots(), {2,3},
                         "Ошибка вычисления корней квадратичного вычета при p mod 8 = 5.")
        self.assertEqual(QuadraticResidue(10, 41).roots(), {16,25},
                         "Ошибка вычисления корней квадратичного вычета при p a^t % p = 1.")
        self.assertEqual(QuadraticResidue(9, 41).roots(), {3, 38},
                         "Ошибка вычисления корней квадратичного вычета при p a^t % p != 1.")
        self.assertEqual(QuadraticResidue(186, 401).roots(), {304, 97},
                         "Ошибка вычисления корней квадратичного вычета при p a^t % p != 1.")

    def test_solve_quadratic_comparison(self):
        """Тесты вычисления корней квадратичного сравнения."""
        self.assertEqual(solve_quadratic_comparison(2, -20, 32, 41),
                         {GF(41)(2), GF(41)(8)},
                         "Ошибка вычисления корней квадратичного сравнения.")


if __name__ == '__main__':
    unittest.main()
