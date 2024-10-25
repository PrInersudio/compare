from integer_mod_ring import IntegerModRing
from sympy import factorint
from ring import Ring
import logging

logger = logging.getLogger(__name__)
logger.setLevel(logging.DEBUG)
handler = logging.FileHandler(f'log/{__name__}.log')
handler.setLevel(logging.DEBUG)
handler.setFormatter(logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s'))
logger.addHandler(handler)

class Field(Ring):
    def __init__(self):
        super().__init__()
        self.has_multiplicative_identity = True
        self.has_multiplicative_inverses = True
        self.is_euclidean = True
        self.is_commutative = True

    def __repr__(self) -> str:
        return "Field()"

    def __str__(self) -> str:
        return f'Абстрактное поле.'

    class Element(Ring.Element):

        def __init__(self, field: 'Field', value = None):
            super().__init__(field, value)
            self.field = self.ring

        def __repr__(self) -> str:
            return f'Field.Element({repr(self.field)},{self.value})'

        def __str__(self) -> str:
            return f'Элемент {self.value} поля {self.field}.'
    
class GF_prime(IntegerModRing, Field):
    def __init__(self, p: int):
        IntegerModRing.__init__(self, p)
        self.p = p
        self.has_multiplicative_inverses = True
        self.is_euclidean = True
        self.is_commutative = True

    def __repr__(self) -> str:
        return f"GF_prime({self.p})"

    def __str__(self) -> str:
        return f'Конечное поле размера {self.p}.'
    
    class Element(IntegerModRing.Element):

        def __init__(self, gf: 'GF', value = None):
            super().__init__(gf, value)
            self.field = self.ring

        def __repr__(self) -> str:
            return f'GF_prime.Element({repr(self.field)},{self.value})'

    
class GF_primary:

    def __init__(self, p: int, q: int):
        raise ValueError("Ещё не разработано.")
    
    def __repr__(self) -> str:
        return f"GF_primary({self.p}, {self.q})."

    def __str__(self) -> str:
        return f'Конечное поле размера {self.p}^{self.q}.'
    
class GF(GF_prime, GF_primary):

    def __init__(self, P: int, dont_check_factors: bool = False, p: int = 0, q: int = 0):
        if dont_check_factors:
            if p and q:
                GF_primary.__init__(self, p, q)
            else:
                GF_prime.__init__(self, P)
        factors = factorint(P)
        if len(factors) != 1: raise ValueError("Размер поля должен быть простым или примарным числом.")
        p, q = list(factors.items())[0]
        if q == 1: GF_prime.__init__(self, p)
        else: GF_primary.__init__(self, p, q)
    

class RealField(Field):
    __instance = None

    def __new__(cls):
        if not cls.__instance:
            logger.debug("Созадаётся RealField.")
            cls.__instance = super().__new__(cls)
        return cls.__instance

    def __init__(self):
        if hasattr(self, '_initialized'): return
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

        def __init__(self, field: 'RealField', value: float):
            super().__init__(field, value)

        def __repr__(self) -> str:
            return f'RealField.Element({repr(self.field)},{self.value})'

        def __str__(self) -> str:
            return str(self.value)


class ComplexField(Field):
    __instance = None

    def __new__(cls):
        if not cls.__instance:
            logger.debug("Созадаётся ComplexField.")
            cls.__instance = super().__new__(cls)
        return cls.__instance

    def __init__(self):
        if hasattr(self, '_initialized'): return
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

        def __init__(self, field: 'ComplexField', value: complex):
            super().__init__(field, value)

        def __repr__(self) -> str:
            return str(self.value)
        
        def __repr__(self) -> str:
            return f'ComplexField.Element({repr(self.field)},{self.value})'

        def __str__(self) -> str:
            return str(self.value)