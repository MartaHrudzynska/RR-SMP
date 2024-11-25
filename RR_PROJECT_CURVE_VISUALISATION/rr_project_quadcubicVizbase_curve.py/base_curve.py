from abc import ABC, abstractmethod
import numpy as np
import matplotlib.pyplot as plt

class Curve(ABC):
    @abstractmethod
    def equation(self, x):
        """Define the curve equation."""
        pass

    def plot(self, x_range):
        x = np.linspace(*x_range, 500)
        y = self.equation(x)
        plt.plot(x, y)
        plt.xlabel('x')
        plt.ylabel('y')
        plt.title(self.__class__.__name__)
        plt.grid()
        plt.show()

class CurveOrder2(Curve):
    def __init__(self, a, b, c):
        self.a = a
        self.b = b
        self.c = c

    @abstractmethod
    def equation(self, x):
        """Define the second-order curve equation."""
        pass

    def __add__(self, other):
        if isinstance(other, CurveOrder2):
            return self.__class__(self.a + other.a, self.b + other.b, self.c + other.c)
        raise TypeError("Can only add another CurveOrder2 instance.")

    def __sub__(self, other):
        if isinstance(other, CurveOrder2):
            return self.__class__(self.a - other.a, self.b - other.b, self.c - other.c)
        raise TypeError("Can only subtract another CurveOrder2 instance.")

class CurveOrder3(Curve):
    def __init__(self, a, b, c, d):
        self.a = a
        self.b = b
        self.c = c
        self.d = d

    @abstractmethod
    def equation(self, x):
        """Define the third-order curve equation."""
        pass

    def __add__(self, other):
        if isinstance(other, CurveOrder3):
            return self.__class__(
                self.a + other.a, self.b + other.b, self.c + other.c, self.d + other.d
            )
        raise TypeError("Can only add another CurveOrder3 instance.")

    def __sub__(self, other):
        if isinstance(other, CurveOrder3):
            return self.__class__(
                self.a - other.a, self.b - other.b, self.c - other.c, self.d - other.d
            )
        raise TypeError("Can only subtract another CurveOrder3 instance.")
