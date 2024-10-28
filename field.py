""" Модуль реализует:
Основной класс поля Field, от которого наследуются все остальные поля.
Класс конечного поле GF.
Класс RealField, который является обёрткой вокруг float.
Класс ComplexField, который является обёрткой вокруг complex.
"""

import logging
from sympy import factorint

from integer_mod_ring import IntegerModRing
from ring import Ring



logger = logging.getLogger(__name__)
logger.setLevel(logging.DEBUG)
handler = logging.FileHandler(f'log/{__name__}.log')
handler.setLevel(logging.DEBUG)
handler.setFormatter(logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s'))
logger.addHandler(handler)


class Field(Ring):
    """Основной класс поля, от которого наследуются все остальные поля."""

    def __init__(self):
        super().__init__()
        self.has_multiplicative_identity = True
        self.has_multiplicative_inverses = True
        self.is_euclidean = True
        self.is_commutative = True

    def __repr__(self) -> str:
        return "Field()"

    def __str__(self) -> str:
        return 'Абстрактное поле.'


    class Element(Ring.Element):
        """Элемент поля Field."""

        def __init__(self, field: 'Field', value = None):
            super().__init__(field, value)
            self.field = self.ring

        def __repr__(self) -> str:
            return f'Field.Element({repr(self.field)},{self.value})'

        def __str__(self) -> str:
            return f'Элемент {self.value} поля {self.field}'


class _GFPrime(IntegerModRing, Field):
    """Конечное поле по простому модулю."""

    def __init__(self, p: int):
        Field.__init__(self)
        IntegerModRing.__init__(self, p)
        self.p = p


    class Element(IntegerModRing.Element):
        """Элемент конечного поля по простому модулю."""

        def __init__(self, gf: 'GF', value = None):
            super().__init__(gf, value)
            self.field = self.ring

        def __repr__(self):
            return f'GF.Element({repr(self.field)}, {self.value})'


class _GFPrimary:
    """Конечное поле по примарному модулю."""

    def __init__(self, p: int, q: int):
        raise ValueError("Ещё не разработано.")


class GF(_GFPrime, _GFPrimary):
    """Конечное поле."""

    def __init__(self, P: int, dont_check_factors: bool = False, p: int = 0, q: int = 0):
        if dont_check_factors:
            if p and q:
                _GFPrimary.__init__(self, p, q)
            else:
                _GFPrime.__init__(self, P)
        factors = factorint(P)
        if len(factors) != 1:
            raise ValueError("Размер поля должен быть простым или примарным числом.")
        p, q = list(factors.items())[0]
        if q == 1:
            _GFPrime.__init__(self, p)
        else: _GFPrimary.__init__(self, p, q)

    def __repr__(self):
        return f'GF({self.m})'

    def __str__(self) -> str:
        return f'Конечное поле размера {self.p}.'


class RealField(Field):
    """Класс RealField является обёрткой вокруг float."""
    __instance = None

    def __new__(cls):
        if not cls.__instance:
            logger.debug("Создаётся RealField.")
            cls.__instance = super().__new__(cls)
        return cls.__instance

    def __init__(self):
        if hasattr(self, '_initialized'):
            return
        logger.debug("Инициализируется RealField.")
        super().__init__()
        self._initialized = True

    def __call__(self, value = 0) -> 'RealField.Element':
        return self.Element(self, float(value))

    def __repr__(self) -> str:
        return "RealField()"

    def __str__(self) -> str:
        return "Поле вещественных чисел."


    class Element(Field.Element):
        """Элемент RealField."""

        def __init__(self, field: 'RealField', value: float):
            super().__init__(field, value)

        def __repr__(self) -> str:
            return f'RealField.Element({repr(self.field)},{self.value})'

        def __str__(self) -> str:
            return str(self.value)


class ComplexField(Field):
    """Класс ComplexField является обёрткой вокруг complex."""
    __instance = None

    def __new__(cls):
        if not cls.__instance:
            logger.debug("Создаётся ComplexField.")
            cls.__instance = super().__new__(cls)
        return cls.__instance

    def __init__(self):
        if hasattr(self, '_initialized'):
            return
        logger.debug("Инициализируется ComplexField.")
        super().__init__()
        self._initialized = True

    def __call__(self, value = 0) -> 'ComplexField.Element':
        return self.Element(self, complex(value))

    def __repr__(self) -> str:
        return "ComplexField()"

    def __str__(self) -> str:
        return "Поле комплексных чисел."


    class Element(Field.Element):
        """Элемент ComplexField."""

        def __init__(self, field: 'ComplexField', value: complex):
            super().__init__(field, value)

        def __repr__(self) -> str:
            return f'ComplexField.Element({repr(self.field)},{self.value})'

        def __str__(self) -> str:
            return str(self.value)
