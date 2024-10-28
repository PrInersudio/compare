""" Модуль реализует:
Основной класс кольца Ring, от которого наследуются все остальные кольца.
Класс IntegerRing, который является обёрткой вокруг int.
"""

from typing import Tuple
import logging


logger = logging.getLogger(__name__)
logger.setLevel(logging.DEBUG)
handler = logging.FileHandler(f'log/{__name__}.log')
handler.setLevel(logging.DEBUG)
handler.setFormatter(logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s'))
logger.addHandler(handler)


class Ring:
    """Основной класс кольца, от которого наследуются все остальные кольца."""

    def __init__(self):
        self.has_multiplicative_identity = False
        self.has_multiplicative_inverses = False
        self.is_euclidean = False
        self.is_commutative = False

    def __call__(self, value = None) -> 'Ring.Element':
        return self.Element(self, value)

    def __repr__(self) -> str:
        return 'Ring()'

    def __str__(self) -> str:
        return 'Абстрактное кольцо.'

    def __eq__(self, other) -> bool:
        return isinstance(other, Ring)


    class Element:
        """Элемент кольца Ring."""

        def __init__(self, ring: 'Ring', value = None):
            if isinstance(value, Ring.Element):
                value = value.value
            self.value = value
            self.ring = ring

        def __str__(self):
            return f'Элемент {self.value} из {self.ring}'

        def __repr__(self) -> str:
            return f'Ring.Element({repr(self.ring)}, {self.value})'

        def _check(self, other) -> Tuple['Ring.Element', bool]:
            if not isinstance(other, Ring.Element):
                try:
                    other = self.ring(other)
                except Exception: # pylint: disable=broad-except
                    logger.info('Операция между %s из %s и %s не поддерживается.',
                                self, self.ring, other)
                    return other, False
            if self.ring != other.ring:
                logger.info('Операция между %s из %s и %s из %s не поддерживается.',
                            self, self.ring, other, other.ring)
                return other, False
            return other, True

        def __eq__(self, other: 'Ring.Element') -> bool:
            other, check_result = self._check(other)
            if check_result:
                return self.value == other.value
            return False

        def __int__(self):
            return int(self.value)

        def __float__(self):
            return float(self.value)

        def __complex__(self):
            return complex(self.value)

        def __bool__(self):
            return bool(self.value)

        def __add__(self, other) -> 'Ring.Element':
            other, check_result = self._check(other)
            if not check_result:
                return NotImplemented
            return self.ring(self.value + other.value)

        def __radd__(self, other) -> 'Ring.Element':
            return self + other

        def __iadd__(self, other) -> 'Ring.Element':
            return self + other

        def __sub__(self, other) -> 'Ring.Element':
            other, check_result = self._check(other)
            if not check_result:
                return NotImplemented
            return self.ring(self.value - other.value)

        def __neg__(self) -> 'Ring.Element':
            return self.ring(-self.value)

        def __rsub__(self, other) -> 'Ring.Element':
            return (-self) + other

        def __isub__(self, other) -> 'Ring.Element':
            return (-self) + other

        def __mul__(self, other) -> 'Ring.Element':
            other, check_result = self._check(other)
            if not check_result:
                return NotImplemented
            return self.ring(self.value * other.value)

        def __rmul__(self, other) -> 'Ring.Element':
            if self.ring.is_commutative:
                return self.__mul__(other)
            other, check_result = self._check(other)
            if not check_result:
                return NotImplemented
            return self.ring(other.value * self.value)

        def __imul__(self, other) -> 'Ring.Element':
            return self * other

        def __truediv__(self, other) -> 'Ring.Element':
            other, check_result = self._check(other)
            if not check_result:
                return NotImplemented
            return self.ring(self.value / other.value)

        def __mod__(self, other) -> 'Ring.Element':
            other, check_result = self._check(other)
            if not check_result:
                return NotImplemented
            return self.ring(self.value % other.value)

        def __divmod__(self, other) -> 'Ring.Element':
            other, check_result = self._check(other)
            if not check_result:
                return NotImplemented
            q, r = divmod(self.value, other.value)
            return self.ring(q), self.ring(r)

        def __floordiv__(self, other) -> 'Ring.Element':
            return self / other

        def __pow__(self, other: int) -> 'Ring.Element':
            return self.value ** other

        def __lt__(self, other) -> bool:
            other, check_result = self._check(other)
            if check_result:
                return self.value < other.value
            return False

        def __le__(self, other) -> bool:
            other, check_result = self._check(other)
            if check_result:
                return self.value <= other.value
            return False

        def __gt__(self, other) -> bool:
            other, check_result = self._check(other)
            if check_result:
                return self.value > other.value
            return False

        def __ge__(self, other) -> bool:
            other, check_result = self._check(other)
            if check_result:
                return self.value >= other.value
            return False

        def __hash__(self):
            return hash(self.value)


class IntegerRing(Ring):
    """Класс IntegerRing является обёрткой вокруг int."""
    __instance = None

    def __new__(cls):
        if not cls.__instance:
            logger.debug("Создаётся IntegerRing.")
            cls.__instance = super().__new__(cls)
        return cls.__instance

    def __init__(self):
        if hasattr(self, '_initialized'):
            return
        logger.debug("Инициализируется IntegerRing.")
        super().__init__()
        self._initialized = True
        self.has_multiplicative_identity = True
        self.is_euclidean = True
        self.is_commutative = True

    def __call__(self, value = 0) -> 'IntegerRing.Element':
        return self.Element(self, int(value))

    def __repr__(self) -> str:
        return "IntegerRing()"

    def __str__(self) -> str:
        return 'Кольцо целых чисел.'


    class Element(Ring.Element):
        """Элемент IntegerRing."""

        def __init__(self, ring: 'IntegerRing', value = None):
            super().__init__(ring, value)

        def __repr__(self) -> str:
            return f'IntegerRing.Element({repr(self.ring)}, {self.value})'

        def __str__(self) -> str:
            return str(self.value)
