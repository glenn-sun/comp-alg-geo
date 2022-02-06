from typing import Any
import polynomial as poly

class Ideal:
    def __init__(self, generators: list[poly.Polynomial]):
        self._groebner_basis = poly.reduced_groebner_basis(generators)
    
    @property
    def groebner_basis(self):
        return self._groebner_basis

    def __eq__(self, other: Any):
        if isinstance(other, Ideal):
            return self.groebner_basis == other.groebner_basis
        else:
            return False

    def __contains__(self, f: poly.Polynomial) -> bool:
        _, r = poly.divmod_multiple(f, list(self.groebner_basis))
        return r == poly.ZERO

    def __repr__(self):
        return "Ideal("+ repr(self.groebner_basis) + ")"