from __future__ import annotations
from functools import total_ordering

class MonomialOrder(tuple[int, ...]):
    __lt__ = object.__lt__
    __le__ = object.__le__
    __ge__ = object.__ge__

@total_ordering
class Lex(MonomialOrder):
    def __gt__(self, other: tuple[int, ...]) -> bool:
        if isinstance(other, Lex):
            for xi, yi in zip(self, other):
                if xi > yi:
                    return True
                if xi < yi:
                    return False
            return False
        else:
            raise TypeError("Must compare with same ordering.")

@total_ordering
class GradedLex(MonomialOrder):
    def __gt__(self, other: tuple[int, ...]) -> bool:
        if isinstance(other, GradedLex):
            if sum(self) > sum(other):
                return True
            elif sum(self) < sum(other):
                return False
            else:
                for xi, yi in zip(self, other):
                    if xi > yi:
                        return True
                    if xi < yi:
                        return False
                return False
        else:
            raise TypeError("Must compare with same ordering.")

class GradedReverseLex(MonomialOrder):
    def __gt__(self, other: tuple[int, ...]) -> bool:
        if isinstance(other, GradedReverseLex):
            if sum(self) > sum(other):
                return True
            elif sum(self) < sum(other):
                return False
            else:
                for xi, yi in zip(reversed(self), reversed(other)):
                    if xi > yi:
                        return False
                    if xi < yi:
                        return True
                return False
        else:
            raise TypeError("Must compare with same ordering.")
