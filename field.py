from integer_mod_ring import IntegerModRing
from sympy import factorint
from ring import Ring

class Field(Ring):
    def __init__(self):
        super().__init__()
        self.has_multiplicative_identity = True
        self.has_multiplicative_inverses = True
        self.is_euclidean = True
        self.is_commutative = True

    def __repr__(self) -> str:
        return f'Абстрактное поле.'

    class Element(Ring.Element):

        def __init__(self, field: 'Field', value = None):
            super().__init__(field, value)
    
class GF_prime(IntegerModRing, Field):
    def __init__(self, p: int):
        IntegerModRing.__init__(self, p)
        self.p = p
        self.has_multiplicative_inverses = True
        self.is_euclidean = True
        self.is_commutative = True

    def __repr__(self) -> str:
        return f'Конечное поле размера {self.m}'
    
    class Element(IntegerModRing.Element):
        def __init__(self, gf: 'GF', value = None):
            super().__init__(gf, value)

    
class GF_primary:

    def __init__(self, p: int, q: int):
        raise ValueError("Ещё не разработано.")
    
class GF(GF_prime, GF_primary):

    def __init__(self, P: int, dont_check_consider_prime: bool = False):
        if dont_check_consider_prime: GF_prime.__init__(self, P)
        factors = factorint(P)
        if len(factors) != 1: raise ValueError("Размер поля должен быть простым или примарным числом.")
        p, q = list(factors.items())[0]
        if q == 1: GF_prime.__init__(self, p)
        else: GF_primary.__init__(self, p, q)
    

class RealField(Field):

    def __init__(self):
        super().__init__()

    def __call__(self, value = 0) -> 'RealField.Element':
        return self.Element(self, float(value))

    class Element(Field.Element):

        def __init__(self, field: 'RealField', value: float):
            super().__init__(field, value)

        def __repr__(self) -> str:
            return str(self.value)


class ComplexField(Field):

    def __init__(self):
        super().__init__()

    def __call__(self, value = 0) -> 'ComplexField.Element':
        return self.Element(self, complex(value))

    class Element(Field.Element):

        def __init__(self, field: 'ComplexField', value: complex):
            super().__init__(field, value)

        def __repr__(self) -> str:
            return str(self.value)