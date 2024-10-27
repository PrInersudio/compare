"""Тестирование Polynomial."""

import unittest
import numpy as np

from polynomial import Polynomial
from integer_mod_ring import IntegerModRing
from ring import IntegerRing
from field import RealField, ComplexField, GF


class TestPolynomial(unittest.TestCase):
    """Тестирование Polynomial."""

    def test_init(self):
        """Тестирование инициализации Polynomial и Polynomial.Element."""
        ring = Polynomial(IntegerRing())
        # константы
        self.assertEqual(ring("0"), ring(), "Ошибка инициализации нулевого полинома.")
        self.assertEqual(ring(""), ring(),
                         "Ошибка инициализации нулевого полинома (пустая строка).")
        self.assertEqual(ring("5"), ring(5),
                         "Ошибка инициализации константы.")
        self.assertEqual(ring("+5"), ring(5),
                         "Ошибка инициализации положительной константы (плюс без пробела).")
        self.assertEqual(ring("+ 5"), ring(5),
                         "Ошибка инициализации положительной константы (плюс с пробелом).")
        self.assertEqual(ring("-5"), ring(-5),
                         "Ошибка инициализации отрицательной константы (без пробела).")
        self.assertEqual(ring("- 5"), ring(-5),
                         "Ошибка инициализации отрицательной константы (с пробелом).")
        # линейный моном
        self.assertEqual(ring("0x"), ring(),
                         "Ошибка инициализации линейного монома с нулевым коэффициентом.")
        self.assertEqual(ring("x"), ring({1:IntegerRing()(1)}),
                         "Ошибка инициализации линейного монома с коэффициентом, равным 1.")
        self.assertEqual(ring("+x"), ring({1:IntegerRing()(1)}),
                         "Ошибка инициализации линейного монома с коэффициентом,"
                         "равным 1 (плюс без пробела).")
        self.assertEqual(ring("+ x"), ring({1:IntegerRing()(1)}),
                         "Ошибка инициализации линейного монома с коэффициентом,"
                         "равным 1 (плюс с пробелом).")
        self.assertEqual(ring("-x"), ring({1:IntegerRing()(-1)}),
                         "Ошибка инициализации линейного монома с коэффициентом,"
                         "равным -1 (без пробела).")
        self.assertEqual(ring("- x"), ring({1:IntegerRing()(-1)}),
                         "Ошибка инициализации линейного монома с коэффициентом,"
                         "равным -1 (с пробелом).")
        self.assertEqual(ring("5x"), ring({1:IntegerRing()(5)}),
                         "Ошибка инициализации линейного монома с положительным коэффициентом.")
        self.assertEqual(ring("+5x"), ring({1:IntegerRing()(5)}),
                         "Ошибка инициализации линейного монома с положительным коэффициентом"
                         "(плюс без пробела).")
        self.assertEqual(ring("+ 5x"), ring({1:IntegerRing()(5)}),
                         "Ошибка инициализации линейного монома с положительным коэффициентом"
                         "(плюс с пробелом).")
        self.assertEqual(ring("-5x"), ring({1:IntegerRing()(-5)}),
                         "Ошибка инициализации линейного монома с отрицательным коэффициентом"
                         "(без пробела).")
        self.assertEqual(ring("- 5x"), ring({1:IntegerRing()(-5)}),
                         "Ошибка инициализации линейного монома с отрицательным коэффициентом"
                         "(с пробелом).")
        # произвольный моном
        self.assertEqual(ring("0x^2"), ring(),
                         "Ошибка инициализации произвольного монома с нулевым коэффициентом.")
        self.assertEqual(ring("5x^2"), ring({2:IntegerRing()(5)}),
                         "Ошибка инициализации произвольного монома с положительным коэффициентом.")
        self.assertEqual(ring("+5x^2"), ring({2:IntegerRing()(5)}),
                         "Ошибка инициализации произвольного монома с положительным коэффициентом"
                         "(плюс без пробела).")
        self.assertEqual(ring("+ 5x^2"), ring({2:IntegerRing()(5)}),
                         "Ошибка инициализации произвольного монома с положительным коэффициентом"
                         "(плюс с пробелом).")
        self.assertEqual(ring("-5x^2"), ring({2:IntegerRing()(-5)}),
                         "Ошибка инициализации произвольного монома с отрицательным коэффициентом"
                         "(без пробела).")
        self.assertEqual(ring("- 5x^2"), ring({2:IntegerRing()(-5)}),
                         "Ошибка инициализации произвольного монома с отрицательным коэффициентом"
                         "(с пробелом).")
        # полиномы
        self.assertEqual(ring("3x^3 - 4x^2 + 2x - 5"),
                         ring({
                             0:IntegerRing()(-5), 1:IntegerRing()(2),
                             2:IntegerRing()(-4), 3:IntegerRing()(3),
                         }), "Ошибка инициализации в общем случае.")
        self.assertEqual(ring("3x^3-4x^2+2x-5"),
                         ring({
                             0:IntegerRing()(-5), 1:IntegerRing()(2),
                             2:IntegerRing()(-4), 3:IntegerRing()(3),
                         }), "Ошибка инициализации в общем случае (написание без пробела).")
        self.assertEqual(ring("2x^2 - 5"),
                         ring({
                             0:IntegerRing()(-5), 2:IntegerRing()(2),
                         }), "Ошибка инициализации при нулевом коэффициенте.")
        self.assertEqual(ring("2x^2 + 2x"),
                         ring({
                             1:IntegerRing()(2), 2:IntegerRing()(2),
                         }), "Ошибка инициализации при нулевом свободном члене.")
        self.assertEqual(ring("3x^2 + 5x^2 - 2x"),
                         ring({
                             1:IntegerRing()(-2), 2:IntegerRing()(8),
                         }), "Ошибка инициализации в случае с одинаковыми степенями.")
        # словарь
        self.assertEqual(ring({1.4:1.4}), ring(),
                         "Ошибка проверки типов при инициализации словарём.")
        self.assertEqual(ring({-5:IntegerRing()(5)}), ring(),
                         "Ошибка проверки положительности степеней при инициализации словарём.")
        self.assertEqual(ring({5:IntegerRing()(0)}), ring(),
                         "Ошибка проверки коэффициентов при инициализации словарём.")
        # исключенияY
        with self.assertRaises(ValueError,
                               msg="init не поднял исключение при отрицательной степени в строке."):
            ring("5x^-5")
        # многочлены не над int
        ring = Polynomial(RealField())
        self.assertEqual(ring("-2.4x^3 + 5.33x^2 - x + 0.00"),
                         ring({
                             1:RealField()(-1), 2:RealField()(5.33), 3:RealField()(-2.4),
                         }), "Ошибка инициализации вещественного многочлена.")
        ring = Polynomial(ComplexField())
        self.assertEqual(ring({1:5, 2:ComplexField()(complex(1,2))}),
                         ring({2:ComplexField()(complex(1,2))}),
                         "Ошибка инициализации комплексного многочлена.")
        gf = GF(7)
        ring = Polynomial(gf)
        self.assertEqual(ring("-2x^3 + 9x^2 - 5x + 1"),
                         ring({0:gf(1), 1:gf(2), 2:gf(2), 3:gf(5)}),
                         "Ошибка инициализации многочлена над кольцом.")

    def test_str(self):
        """Тестирование преобразования многочлена в строку."""
        ring = Polynomial(IntegerRing())
        self.assertEqual(str(ring()), "0",
                         "Ошибка преобразования в строку нулевого полинома.")
        self.assertEqual(str(ring("5")), "5",
                         "Ошибка преобразования в строку положительной константы.")
        self.assertEqual(str(ring("-5")), "-5",
                         "Ошибка преобразования в строку отрицательной константы.")
        self.assertEqual(str(ring("5x")), "5x",
                         "Ошибка преобразования в строку линейного монома"
                         "с положительным коэффициентом.")
        self.assertEqual(str(ring("-5x")), "-5x",
                         "Ошибка преобразования в строку линейного монома"
                         "с отрицательным коэффициентом.")
        self.assertEqual(str(ring("5x^2")), "5x^2",
                         "Ошибка преобразования в строку произвольного монома"
                         "с положительным коэффициентом.")
        self.assertEqual(str(ring("-5x^2")), "-5x^2",
                         "Ошибка преобразования в строку произвольного монома"
                         "с отрицательным коэффициентом.")
        self.assertEqual(str(ring("3x^3 - 4x^2 + 2x - 5")), "3x^3 - 4x^2 + 2x - 5",
                         "Ошибка преобразования в строку в общем случае.")
        self.assertEqual(str(ring("2x^2 - 5")), "2x^2 - 5",
                         "Ошибка преобразования в строку при нулевом коэффициенте.")
        self.assertEqual(str(ring("2x^2 + 2x")), "2x^2 + 2x",
                         "Ошибка преобразования в строку при нулевом свободном члене.")
        self.assertEqual(str(Polynomial(RealField())("-2.4x^3 + 5.33x^2 - x + 0.00")),
                         "-2.4x^3 + 5.33x^2 - 1.0x",
                         "Ошибка преобразования в строку вещественного многочлена.")
        self.assertEqual(str(Polynomial(ComplexField())({
            0:ComplexField()(complex(1,2)), 1:ComplexField()(complex(7,8)),
            2:ComplexField()(complex(3,4)), 3:ComplexField()(complex(5,6)),
            })), "(5+6j)x^3 + (3+4j)x^2 + (7+8j)x + (1+2j)",
            "Ошибка преобразования в строку комплексного многочлена.")
        self.assertEqual(str(Polynomial(GF(7))("-2x^3 + 9x^2 - 5x + 1")),
                         "[5]_7x^3 + [2]_7x^2 + [2]_7x + [1]_7",
                         "Ошибка преобразования в строку многочлена над кольцом.")

    def test_getitem(self):
        """Тестирование получения коэффициентов многочлена."""
        ring = Polynomial(IntegerRing())
        self.assertEqual(ring('5x^2')[2], 5,
                         "Ошибка при получении значения ненулевого коэффициента.")
        self.assertEqual(ring()[2], 0, "Ошибка при получении значения нулевого коэффициента.")

    def test_setitem(self):
        """Тестирование установки коэффициентов многочлена."""
        ring = Polynomial(IntegerRing())
        p = ring()
        p[2] = 5
        self.assertEqual(p, ring({2:IntegerRing()(5)}), "Ошибка при установке коэффициента.")
        p[2] = 2
        self.assertEqual(p, ring({2:IntegerRing()(2)}), "Ошибка при изменении коэффициента.")
        p[2] = 0
        self.assertEqual(p, ring(), "Ошибка обнуления коэффициента.")

    def test_contains(self):
        """Тестирование проверки, является ли коэффициент ненулевым."""
        p = Polynomial(IntegerRing())('5x^2')
        self.assertTrue(2 in p,
                        "Ошибка false negative при проверке содержания степени в полиноме.")
        self.assertFalse(3 in p,
                         "Ошибка false positive при проверке содержания степени в полиноме.")

    def test_iter(self):
        """Тестирование итерации по степеням, соответствующим ненулевым коэффициентам."""
        ring = Polynomial(IntegerRing())
        self.assertEqual(list(ring()), [], "Ошибка iter при нулевом многочлене.")
        self.assertEqual(list(ring({
            0:IntegerRing()(1), 2:IntegerRing()(3), 4:IntegerRing()(5),
            })), [0,2,4], "Ошибка iter при ненулевом многочлене.")

    def test_monomials(self):
        """Тестирование итерации по парам (степень, коэффициент)."""
        ring = Polynomial(IntegerRing())
        self.assertEqual(list(ring().monomials()), [], "Ошибка monomials при нулевом многочлене.")
        self.assertEqual(list(ring({
            0:IntegerRing()(1), 2:IntegerRing()(3), 4:IntegerRing()(5),
            }).monomials()), [(0,1),(2,3),(4,5)], "Ошибка monomials при ненулевом многочлене.")

    def test_copy(self):
        """Тестирование копирования многочлена."""
        ring = Polynomial(IntegerRing())
        p = ring()
        self.assertEqual(p, p.copy(), "Ошибка copy при нулевом многочлене.")
        p = ring({0:1, 2:3, 4:5})
        self.assertEqual(p, p.copy(), "Ошибка copy при ненулевом многочлене.")

    def test_add(self):
        """Тестирование сложения многочленов."""
        ring = Polynomial(IntegerRing())
        self.assertEqual(ring() + ring(), ring(), "Ошибка при суммировании нулевых многочленов.")
        self.assertEqual(ring("3x^2 + 2x - 5") + ring(), ring("3x^2 + 2x - 5"),
                         "Ошибка при суммировании нулевого и ненулевого многочлена.")
        self.assertEqual(ring(5) + ring(-3), ring(2), "Ошибка при суммировании констант.")
        self.assertEqual(ring("2x^2 + 4x + 5") + ring("4x^2"), ring("6x^2 + 4x + 5"),
                         "Ошибка при суммировании многочленов одной степени.")
        self.assertEqual(ring("2x^2 + 4x + 5") + ring("- 4x^3"), ring("- 4x^3 + 2x^2 + 4x + 5"),
                         "Ошибка при суммировании многочленов разной степени.")
        self.assertEqual(ring("3x^2 + 2x - 5") + ring("-3x^2 - 2x + 5"), ring(),
                         "Ошибка при суммировании многочлена с его противоположным.")
        self.assertEqual(ring("3x^2 + 2x - 5") + 5, ring("3x^2 + 2x"),
                         "Ошибка при суммировании многочлена с int.")
        self.assertEqual(5 + ring("3x^2 + 2x - 5"), ring("3x^2 + 2x"),
                         "Ошибка при суммировании int с многочлена.")

    def test_sub(self):
        """Тестирование вычитания многочленов."""
        ring = Polynomial(IntegerRing())
        self.assertEqual(ring() - ring(), ring(),
                         "Ошибка при вычитании нулевых многочленов.")
        self.assertEqual(ring("3x^2 + 2x - 5") - ring(), ring("3x^2 + 2x - 5"),
                         "Ошибка при вычитании нулевого и ненулевого многочлена.")
        self.assertEqual(ring(5) - ring(-3), ring(8),
                         "Ошибка при вычитании констант.")
        self.assertEqual(ring("2x^2 + 4x + 5") - ring("4x^2"), ring("-2x^2 + 4x + 5"),
                         "Ошибка при вычитании многочленов одной степени.")
        self.assertEqual(ring("2x^2 + 4x + 5") - ring("- 4x^3"), ring("4x^3 + 2x^2 + 4x + 5"),
                         "Ошибка при вычитании многочленов разной степени.")
        self.assertEqual(ring("3x^2 + 2x - 5") - ring("3x^2 + 2x - 5"), ring(),
                         "Ошибка при вычитании равных многочленов.")
        self.assertEqual(ring("3x^2 + 2x + 5") - 5, ring("3x^2 + 2x"),
                         "Ошибка при левом вычитании int из многочлена.")
        self.assertEqual(5 - ring("-3x^2 - 2x + 5"), ring("3x^2 + 2x"),
                         "Ошибка при правом вычитании многочлена из int.")

    def test_mul(self):
        """Тестирование умножения многочленов."""
        ring = Polynomial(IntegerRing())
        self.assertEqual(ring() * ring(), ring(), "Ошибка при умножении нулевых многочленов.")
        self.assertEqual(ring("3x^2 + 2x - 5") * ring(), ring(),
                         "Ошибка при умножении ненулевого и нулевого многочлена.")
        self.assertEqual(ring() * ring("3x^2 + 2x - 5"), ring(),
                         "Ошибка при умножении нулевого и ненулевого многочлена.")
        self.assertEqual(ring("5") * ring("3x^2 + 2x - 5"), ring("15x^2 + 10x - 25"),
                         "Ошибка при умножении константы на многочлен.")
        self.assertEqual(ring("x^2 + 2x + 1") * ring("-3"), ring("-3x^2 - 6x - 3"),
                         "Ошибка при умножении многочлена на константу.")
        self.assertEqual(ring(5) * ring(-3), ring(-15), "Ошибка при умножении констант.")
        self.assertEqual(ring("3x^3 + 9x - 5") * ring("-2x^2 + 4x + 6"),
                         ring("-6x^5 + 12x^4 + 46x^2 + 34x - 30"),
                         "Ошибка умножения многочленов в общем случае.")
        self.assertEqual(ring("3x^2 + 2x - 5") * 0, ring(),
                         "Ошибка при умножении многочлена на нулевой int.")
        self.assertEqual(0 * ring("3x^2 + 2x - 5"), ring(),
                         "Ошибка при умножении нулевого int на многочлен.")
        self.assertEqual(5 * ring("3x^2 + 2x - 5"), ring("15x^2 + 10x - 25"),
                         "Ошибка при умножении int на многочлен.")
        self.assertEqual(ring("x^2 + 2x + 1") * -3, ring("-3x^2 - 6x - 3"),
                         "Ошибка при умножении многочлена на int.")

    def test_bool(self):
        """Тестирование преобразования многочлена в bool."""
        ring = Polynomial(IntegerRing())
        self.assertFalse(bool(ring()), "Ненулевой многочлен даёт True, хотя должна быть False.")
        self.assertTrue(bool(ring(5)), "Ненулевая константа даёт False, хотя должна быть True.")
        self.assertTrue(bool(ring("5x^2")), "Моном даёт False, хотя должна быть True.")
        self.assertTrue(bool(ring("5x^2 + 5x")),
                        "Многочлен с нулевым свободным членом даёт False, хотя должна быть True.")
        self.assertTrue(bool(ring("5x^2 + 5x + 5")), "Многочлен даёт False, хотя должна быть True.")

    def test_degree(self):
        """Тестирование получения степени многочлена."""
        ring = Polynomial(IntegerRing())
        self.assertEqual(ring().degree(), -1, "Степень нулевого многочлена должна быть -1.")
        self.assertEqual(ring(5).degree(), 0, "Степень константы должна быть нулём.")
        self.assertEqual(ring("3x^2 + 2x + 1").degree(), 2,
                         "Ошибка при определении степени в общем случае.")
        self.assertEqual(ring("x^3 + 5").degree(), 3,
                         "Ошибка при определении степени в многочлене с пропущенными степенями.")
        self.assertEqual(ring("x^10 + 2x^9 + x^5").degree(), 10,
                         "Ошибка при определении степени в многочлене"
                         "без свободной и линейной части.")

    def test_call(self):
        """Тестирование получения значения многочлена в точке."""
        ring = Polynomial(IntegerRing())
        self.assertEqual(ring()(5), 0, "Значение нулевого многочлена должно быть 0 в любой точке.")
        self.assertEqual(ring(5)(10), 5,
                         "Значение константы должно быть этой константой в любой точке.")
        self.assertEqual(ring("3x^2 + 2x - 1")(0), -1,
                         "Значение многочлена в нуле должно быть равно его свободному члену.")
        self.assertEqual(ring("3x^2 + 2x - 1")(5), 84,
                         "Ошибка вычисления значения полинома в общем случае.")
        self.assertEqual(Polynomial(GF(7))("5x^3 + 2x^2 + 3x + 1")(4), 1,
                         "Ошибка вычисления значения полинома в поле.")

    def test_roots_berlekamp_rabin(self):
        """Тестирование получения корней многочлена с помощью вероятностной факторизации."""
        # pylint: disable=protected-access
        self.assertEqual(Polynomial(GF(5))("x + 1")._roots_berlekamp_rabin(), {GF(5)(4)},
                         "Ошибка _roots_berlekamp_rabin в случае линейного многочлена.")
        self.assertEqual(Polynomial(GF(3))("x^2 + 1")._roots_berlekamp_rabin(), set(),
                         "Ошибка _roots_berlekamp_rabin в случае неприводимого многочлена.")
        self.assertEqual(Polynomial(GF(7))("x^3 + x")._roots_berlekamp_rabin(), {GF(7)()},
                         "Ошибка _roots_berlekamp_rabin в случае, если многочлен раскладывается"
                         "на линейный и нелинейный.")
        self.assertEqual(Polynomial(GF(5))("3x^3 + 4x^2 + 2x + 1")._roots_berlekamp_rabin(),
                         {GF(5)(1), GF(5)(2), GF(5)(4)},
                         "Ошибка _roots_berlekamp_rabin в общем случае.")
        # pylint: enable=protected-access

    def test_fermat_little_theorem(self):
        """Проверка редуцирования многочлена по малой теореме Ферма."""
        self.assertEqual(Polynomial(GF(3))().fermat_little_theorem(), Polynomial(GF(3))(),
                         "Нулевой многочлен после применения Малой теоремы Ферма"
                         "должен оставаться нулевым.")
        self.assertEqual(Polynomial(GF(3))(5).fermat_little_theorem(), Polynomial(GF(3))(2),
                         "Константа после применения Малой теоремы Ферма"
                         "должна оставаться константой.")
        self.assertEqual(Polynomial(
            GF(5))("6x^8 - 3x^7 + 5x^6 + x^5 + x^3").fermat_little_theorem(),
            Polynomial(GF(5))("x^4 + 3x^3 + x"),
            "Ошибка применения Малой теоремы Ферма в общем случае."
        )

    def test_divmod(self):
        """Тестирование деления с остатком многочлена."""
        self.assertEqual(divmod(Polynomial(GF(3))(), Polynomial(GF(3))("x^2 + x + 1")),
                         (Polynomial(GF(3))(), Polynomial(GF(3))()),
                         "Деление нулевого многочлена должно давать нулевой многочлен.")
        with self.assertRaises(ZeroDivisionError,
                               msg="При попытке деления на нулевой многочлен"
                               "должно подниматься исключение."):
            divmod(Polynomial(GF(3))("x^2 + x + 1"), Polynomial(GF(3))())
        with self.assertRaises(ZeroDivisionError,
                               msg="При попытке деления на нулевой многочлен"
                               "должно подниматься исключение"
                               "(случай преобразования элемента в многочлен)."):
            divmod(Polynomial(GF(3))("x^2 + x + 1"), 3)
        self.assertEqual(divmod(Polynomial(GF(7))("x^2 + x + 1"),
                                Polynomial(GF(7))("x^2 + x + 1")),
                                (Polynomial(GF(7))(1), Polynomial(GF(7))()),
                                "При делении многочлена на самого себя должна получаться единица.")
        self.assertEqual(divmod(Polynomial(GF(7))("2x^3 + 5x + 1"), 2),
                         (Polynomial(GF(7))("x^3 + 6x + 4"), Polynomial(GF(7))()),
                         "Ошибка при делении многочлена на константу.")
        self.assertEqual(divmod(Polynomial(GF(5))("x^2 + 2x + 1"),
                                Polynomial(GF(5))("x^3 + x + 1")),
                                (Polynomial(GF(5))(), Polynomial(GF(5))("x^2 + 2x + 1")),
                                "Если степень делителя выше степени делимого, результатом"
                                "должно быть частное 0 и остаток, равный самому делимому.")
        self.assertEqual(divmod(Polynomial(GF(5))("10x^2 + 2"), Polynomial(GF(5))("3x")),
                         (Polynomial(GF(5))(), Polynomial(GF(5))("2")),
                         "Если в делимом есть коэффициенты, кратные модулю, они должны"
                         "быть занулены перед делением.")
        with self.assertRaises(ZeroDivisionError,
                               msg="Если все коэффициенты делителя кратны модулю,"
                               "должно подыматься исключение деления на нуль"):
            divmod(Polynomial(GF(7))("x^5 + x + 1"), Polynomial(GF(7))("7x^2"))
        self.assertEqual(divmod(Polynomial(GF(7))("3x^3 + x^2 + x + 5"),
                                Polynomial(GF(7))("x^2 + x + 1")),
                                (Polynomial(GF(7))("3x + 5"), Polynomial(GF(7))()),
                                "Ошибка при делении нацело.")
        self.assertEqual(divmod(Polynomial(GF(7))("3x^3 + x^2 + 2x + 6"),
                                Polynomial(GF(7))("x^2 + x + 1")),
                                (Polynomial(GF(7))("3x + 5"), Polynomial(GF(7))("x + 1")),
                                "Ошибка при делении в общем случае.")
        with self.assertRaises(TypeError,
                               msg="divmod не поднимает исключение"
                               "в случае многочлена над над Z."):
            divmod(Polynomial(IntegerRing())("x^2 - 1"), Polynomial(IntegerRing())("x - 1"))

    def test_derivative(self):
        """Проверка взятия производной многочлена."""
        ring =  Polynomial(IntegerRing())
        self.assertEqual(ring().derivative(), ring(),
                         "Производная нулевого многочлена должна быть нулевой.")
        self.assertEqual(ring(5).derivative(), ring(),
                         "Производная константы должна быть нулевой.")
        self.assertEqual(ring("3x^3 + x^2 + x + 5").derivative(),
                         ring("9x^2 + 2x + 1"),
                         "Ошибка при определении производной в общем случае.")

    def test_roots_prime_module(self):
        """Тестирование вычисления корней по простому модулю."""
        gf = GF(5)
        # pylint: disable=protected-access
        self.assertEqual(Polynomial(gf)()._roots_prime_module(), set(gf),
                         "Корнями нулевого многочлена должны являться все элементы поля.")
        self.assertEqual(Polynomial(gf)("5x^2 + 5x + 5")._roots_prime_module(), set(gf),
                         "Корнями нулевого многочлена должны являться все элементы поля"
                         "(случай нулевого по модулю).")
        self.assertEqual(Polynomial(gf)("x - 3")._roots_prime_module(), {gf(3)},
                         "Ошибка _roots_prime_module в случае линейного с одним корнем.")
        self.assertEqual(Polynomial(gf)(3)._roots_prime_module(), set(),
                         "Ошибка _roots_prime_module в случае константного без корней.")
        self.assertEqual(Polynomial(GF(7))("x^2 - 1")._roots_prime_module(),
                         {GF(7)(1), GF(7)(6)},
                         "Ошибка _roots_prime_module в случае квадратного с двумя корнями.")
        self.assertEqual(Polynomial(GF(13))("x^2 - 5")._roots_prime_module(), set(),
                         "Ошибка _roots_prime_module в случае квадратного без корней.")
        self.assertEqual(Polynomial(GF(13))("x^2 - 26")._roots_prime_module(), {GF(13)(0)},
                         "Ошибка _roots_prime_module в случае квадратного с одним корнем.")
        self.assertEqual(Polynomial(gf)("3x^3 + 4x^2 + 2x + 1")._roots_prime_module(),
                         {gf(1), gf(2), gf(4)}, "Ошибка _roots_prime_module в общем случае.")
        # pylint: enable=protected-access

    def test_roots_primary_module(self):
        """Тестирование вычисления корней по примарному модулю."""
        # pylint: disable=protected-access
        ring = IntegerModRing(625)
        self.assertEqual(
            Polynomial(ring)("3x^3 + 4x^2 + 2x + 1")._roots_primary_module(5, 4),
            {ring(72), ring(136), ring(624)},
            "Ошибка _roots_primary_module в общем случае.",
        )
        gf = GF(5)
        self.assertEqual(Polynomial(gf)("3x^3 + 4x^2 + 2x + 1")._roots_primary_module(5),
                         {gf(1), gf(2), gf(4)}, "Ошибка _roots_primary_module случае q = 1.")
        ring = IntegerModRing(28561)
        self.assertEqual(Polynomial(ring)("x^2 - 5")._roots_primary_module(13, 4), set(),
                         "Ошибка _roots_primary_module в случае отсутствия"
                         "корней по простому модулю.")
        # pylint: enable=protected-access

    def test_roots(self):
        """Тестирование вычисления корней."""
        ring = IntegerModRing(189)
        self.assertEqual(
            Polynomial(ring)("162x^10+135x^9+162x^8+56x^4+162x^3+162x^2+113x+188").roots(),
            {ring(94), ring(67), ring(121)}, "Ошибка roots в кольце.",
        )
        ring = IntegerModRing(16)
        self.assertEqual(Polynomial(ring)("x^3+11x^2+2x+8").roots(),
                         {ring(2), ring(4), ring(10), ring(12), ring(15)},
                         "Ошибка roots в кольце в случае примарного модуля.")
        self.assertAlmostEqual(Polynomial(RealField())("x^2 - 1").roots(),
                               {np.complex64(-1), np.complex64(1)},
                               "Ошибка roots в вещественном поле.")
        with self.assertRaises(
            TypeError,
            msg="roots не поднимает исключение в случае многочлена над над Z.",
        ):
            Polynomial(IntegerRing())("x^2 - 1").roots()

    def test_gcd(self):
        """Тестирование получения НОД многочлена."""
        gf = GF(5)
        self.assertEqual(
            Polynomial.gcd(Polynomial(gf)("2x^3 - 3x + 1"), Polynomial(gf)()),
            Polynomial(gf)("x^3 + x + 3"),
            "Ошибка gcd в случае, если один из многочленов нулевой."
        )
        self.assertEqual(
            Polynomial.gcd(Polynomial(gf)("2x^3 - 3x + 1"), Polynomial(gf)("2x^3 - 3x + 1")),
            Polynomial(gf)("x^3 + x + 3"),
            "Ошибка gcd в случае одинаковых многочленов."
        )
        self.assertEqual(
            Polynomial.gcd(Polynomial(gf)("2x^3 - 3x + 1"), Polynomial(gf)("x^3 + x + 3")),
            Polynomial(gf)("x^3 + x + 3"),
            "Ошибка gcd в случае одинаковых унитарного и не унитарного многочленов."
        )
        self.assertEqual(
            Polynomial.gcd(Polynomial(gf)("2x - 2"), Polynomial(gf)("x^2 - 1")),
            Polynomial(gf)("x + 4"),
            "Ошибка gcd в случае если один делитель другого."
        )
        self.assertEqual(
            Polynomial.gcd(Polynomial(gf)("x^3 + x"), Polynomial(gf)("x^3 + 3x^2 + 3x")),
            Polynomial(gf)("x"),
            "Ошибка gcd в общем случае."
        )
        gf = GF(13)
        self.assertEqual(
            Polynomial.gcd(Polynomial(gf)("3x^4 + 6x^2 - 1"), Polynomial(gf)("x^2 - 7")),
            Polynomial(gf)(1),
            "Ошибка gcd в случае взаимно простых."
        )
        self.assertEqual(
            Polynomial.gcd(Polynomial(gf)(), Polynomial(gf)()),
            Polynomial(gf)(),
            "Ошибка gcd в случае, если оба многочлена нулевые."
        )
        self.assertEqual(
            Polynomial.gcd(Polynomial(gf)("3x^4 + 6x^2 - 1"), Polynomial(gf)(5)),
            Polynomial(gf)(1),
            "Ошибка gcd в случае, если один многочлен константа."
        )

    def test_binomial_theorem(self):
        """Тестирование разложения в многочлен по формуле бинома Ньютона."""
        ring = Polynomial(IntegerRing())
        self.assertEqual(
            ring.binomial_theorem(2, 3, 0), ring(1),
            "Ошибка binomial_theorem в случае exp=0.",
        )
        self.assertEqual(
            ring.binomial_theorem(2, 3, 1), ring('2x+3'),
            "Ошибка binomial_theorem в случае exp=1.",
        )
        self.assertEqual(
            ring.binomial_theorem(0, 3, 2), ring(9),
            "Ошибка binomial_theorem в случае m=0.",
        )
        self.assertEqual(
            ring.binomial_theorem(2, 0, 3), ring('8x^3'),
            "Ошибка binomial_theorem в случае a=0.",
        )
        self.assertEqual(
            ring.binomial_theorem(1, 1, 5),
            ring('x^5 + 5x^4 + 10x^3 + 10x^2 + 5x + 1'),
            "Ошибка binomial_theorem в общем случае.",
        )
        self.assertEqual(
            ring.binomial_theorem(-1, 2, 3),
            ring('-x^3 + 6x^2 - 12x + 8'),
            "Ошибка binomial_theorem случае m<0.",
        )
        self.assertEqual(
            ring.binomial_theorem(2, -3, 2),
            ring('4x^2 - 12x + 9'),
            "Ошибка binomial_theorem случае a<0.",
        )
        self.assertEqual(
            ring.binomial_theorem(0, 0, 5), ring(),
            "Ошибка binomial_theorem случае m=0, a=0.",
        )

    def test_monic(self):
        """Тестирование преобразование в унитарный многочлен."""
        self.assertEqual(Polynomial(GF(5))().monic(), Polynomial(GF(5))(),
                         "Ошибка monic в случае нулевого многочлена.")
        self.assertEqual(Polynomial(GF(7))("5x^3 - 3x^2 + 4x - 1").monic(),
                         Polynomial(GF(7))("x^3 + 5x^2 + 5x + 4"),
                         "Ошибка monic над конечным полем.")


if __name__ == '__main__':
    unittest.main()
