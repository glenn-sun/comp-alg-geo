from __future__ import annotations
from collections import defaultdict
from functools import cache
from itertools import combinations
from typing import Any, Mapping
from fractions import Fraction

import config
from exceptions import NotDivisibleError

class Polynomial:
    def __init__(self, terms: Mapping[tuple[int, ...], float | Fraction]):
        if terms == dict():
            self._num_vars = 0
        else:
            self._num_vars = len(list(terms)[0])
        for degree, _ in terms.items():
            if len(degree) != self._num_vars:
                raise ValueError("All degrees must be same length tuples.")
        nonzero_terms = {d: Fraction(c) for d, c in terms.items() if c != 0}
        self._terms = defaultdict(Fraction, nonzero_terms)

    @property
    def terms(self) -> defaultdict[tuple[int, ...], Fraction]:
        return self._terms

    @property
    def num_vars(self) -> int:
        return self._num_vars

    @property
    @cache
    def degree(self) -> tuple[int, ...]:
        return max(self.terms, key=config.order)

    @property
    def lead_coeff(self) -> Fraction:
        return self.terms[self.degree]
    
    @property
    def lead_monomial(self) -> Monomial:
        return Monomial(self.degree)
    
    @property
    def lead_term(self) -> Polynomial:
        return Polynomial({self.degree: self.lead_coeff})

    def _del_if_zero(self, degree: tuple[int, ...]):
        if self.terms[degree] == 0:
            del self.terms[degree]

    def __hash__(self) -> int:
        return hash(tuple(sorted(self.terms.items())))

    def __eq__(self, other: Any) -> bool:
        if isinstance(other, Polynomial):
            return self.terms == other.terms
        else:
            return False

    def __add__(self, other: Polynomial) -> Polynomial:
        terms = defaultdict[tuple[int, ...], Fraction](Fraction)
        for degree, coeff in self.terms.items():
            terms[degree] += coeff
            self._del_if_zero(degree)
        for degree, coeff in other.terms.items():
            terms[degree] += coeff
            self._del_if_zero(degree)
        return Polynomial(terms)
    
    def __neg__(self) -> Polynomial:
        return Polynomial({d: -c for d, c in self.terms.items()})
    
    def __sub__(self, other: Polynomial) -> Polynomial:
        return self + (-other)

    def __mul__(self, other: Polynomial | Fraction) -> Polynomial:
        if isinstance(other, Polynomial):
            terms = defaultdict[tuple[int, ...], Fraction](Fraction)
            for d1, c1 in self.terms.items():
                for d2, c2 in other.terms.items():
                    degree = tuple(sum(d) for d in zip(d1, d2))
                    coeff = c1 * c2
                    terms[degree] += coeff
                    self._del_if_zero(degree)
            return Polynomial(terms)
        else:
            return Polynomial({d: c * other for d, c in self.terms.items()})
        
    def __rmul__(self, other: Fraction) -> Polynomial:
        return self * other

    def __divmod__(self, other: Polynomial) -> tuple[Polynomial, Polynomial]:
        q, r = divmod_multiple(self, [other])
        return q[0], r

    def __truediv__(self, other: Polynomial) -> Polynomial:
        q, r = divmod(self, other)
        if r == ZERO:
            return q
        else:
            raise NotDivisibleError

    def __floordiv__(self, other: Polynomial) -> Polynomial:
        q, _ = divmod(self, other)
        return q
        
    def __mod__(self, other: Polynomial) -> Polynomial:
        _, r = divmod(self, other)
        return r
    
    def __pow__(self, other: int) -> Polynomial:
        if not other >= 0:
            raise ValueError("Exponent must be a nonnegative integer.")
        output = ZERO
        for _ in range(other):
            output *= self
        return output

    def __repr__(self) -> str:
        if self == ZERO:
            return "0"

        term_strings = list[str]()
        for d in sorted(self.terms, key=config.order, reverse=True):
            coeff = self.terms[d]
            if coeff == int(coeff):
                coeff = int(coeff)
            
            if all(exp == 0 for exp in d):
                term_strings.append(str(coeff))
                continue

            factors = list[str]()
             
            if coeff != 1 and coeff != -1:
                factors.append(str(coeff))

            for i, exp in enumerate(d):
                if exp == 0:
                    continue
                elif exp == 1:
                    factors.append("x_" + str(i+1))
                else:
                    factors.append("x_" + str(i+1) + "^" + str(exp))
            
            times_string = "*".join(factors)
            if coeff == -1:
                times_string = "-" + times_string

            term_strings.append(times_string)
        plus_string = " + ".join(term_strings)
        return plus_string.replace("+ -", "- ")

class Monomial(Polynomial):
    def __init__(self, degree: tuple[int]):
        super().__init__({degree: Fraction(1)})

    def __mul__(self, other: Polynomial | Fraction) -> Monomial | Polynomial:
        if isinstance(other, Monomial):
            degree = tuple(i + j for i, j in zip(self.degree, other.degree))
            return Monomial(degree)
        elif isinstance(other, Fraction):
            return super().__mul__(other)
        else:
            return other * self
    
    def __divmod__(self, other: Polynomial) -> \
        tuple[Monomial, Polynomial] | tuple[Polynomial, Monomial] | \
        tuple[Polynomial, Polynomial]:
        if isinstance(other, Monomial):
            q_degree = tuple(i - j for i, j in zip(self.degree, other.degree))
            if all(i >= 0 for i in q_degree):
                return Monomial(q_degree), ZERO
            else:
                return ZERO, self
        else:
            return super().__divmod__(other)

    def __truediv__(self, other: Polynomial) -> Monomial | Polynomial:
        q, r = divmod(self, other)
        if r == ZERO:
            return q
        else:
            raise NotDivisibleError

    def __floordiv__(self, other: Polynomial) -> Monomial | Polynomial:
        q, _ = divmod(self, other)
        return q
        
    def __mod__(self, other: Polynomial) -> Monomial | Polynomial:
        _, r = divmod(self, other)
        return r

def lcm(f: Polynomial, g: Polynomial):
    if isinstance(f, Monomial) and isinstance(g, Monomial):
        degree = tuple(max(i, j) for i, j in zip(f.degree, g.degree))
        return Monomial(degree)
    else:
        raise NotImplementedError

def divmod_multiple(f: Polynomial, divisors: list[Polynomial]) -> \
    tuple[list[Polynomial], Polynomial]:
    def t_div(p: Polynomial, q: Polynomial) -> Polynomial:
        degree = tuple([i - j for i, j in zip(p.degree, q.degree)])
        if all(i >= 0 for i in degree):
            return Polynomial({degree: p.lead_coeff / q.lead_coeff})
        else:
            raise NotDivisibleError

    s = len(divisors)
    remainder = ZERO
    quotients = [ZERO] * s
    dividend = f
    while dividend != ZERO:
        i = 0
        division_occurred = False
        while i < s and not division_occurred:
            try: 
                q = t_div(dividend.lead_term, divisors[i].lead_term)
                quotients[i] += q
                dividend -= q * divisors[i]
                division_occurred = True
            except NotDivisibleError:
                i += 1
        if not division_occurred:
            remainder += dividend.lead_term
            dividend -= dividend.lead_term
    return quotients, remainder

def S_polynomial(f: Polynomial, g: Polynomial) -> Polynomial:
    x = lcm(f.lead_monomial, g.lead_monomial)
    return (x / f.lead_term) * f - x / g.lead_term * g

def is_groebner_basis(basis: list[Polynomial]) -> bool:
    for f, g in combinations(basis, 2):
        _, r = divmod_multiple(S_polynomial(f, g), basis)
        if r != ZERO:
            return False
    return True

def reduced_groebner_basis(generators: list[Polynomial]) -> set[Polynomial]:
    groebner_basis = generators.copy()
   
    while True:
        old_basis = groebner_basis.copy()
        for f, g in combinations(old_basis, 2):
            _, r = divmod_multiple(S_polynomial(f, g), old_basis)
            if r != ZERO:
                groebner_basis.append(r)
        if old_basis == groebner_basis:
            break
    
    minimal_basis = set(groebner_basis)
    for f in groebner_basis:
        for g in minimal_basis:
            if f == g:
                continue
            if f.lead_monomial % g.lead_monomial == ZERO:
                minimal_basis.discard(f)
                break
    
    minimal_basis = set(1 / f.lead_coeff * f for f in minimal_basis)

    reduced_basis = minimal_basis.copy()
    for f in minimal_basis:
        reduced_basis.discard(f)
        _, new_f = divmod_multiple(f, list(reduced_basis))
        reduced_basis.add(new_f)

    return reduced_basis

ZERO = Polynomial({})