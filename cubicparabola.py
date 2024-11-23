import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.optimize import fmin

class CubicParabola:
    def __init__(self, A, B, C, D):
        """
        Initializes a cubic parabola with the given coefficients.

        Args:
        - A, B, C, D: Coefficients for the equation y = Ax^3 + Bx^2 + Cx + D
        """
        self.A = A
        self.B = B
        self.C = C
        self.D = D

    def equation(self, x):
        """Computes the y-value for a given x based on the cubic parabola equation."""
        return self.A * x**3 + self.B * x**2 + self.C * x + self.D

    def derivative(self, x):
        """Computes the derivative (slope) of the cubic equation at x."""
        return 3 * self.A * x**2 + 2 * self.B * x + self.C

    def vertex(self):
        """Finds the vertex of the cubic parabola by minimizing the derivative (using numerical optimization)."""
        critical_x = fmin(lambda x: np.abs(self.derivative(x)), 0, disp=False)[0]   
        critical_y = self.equation(critical_x)
        return (critical_x, critical_y)

    
    def plot(self, x_range=(-10, 10)):
        """Plots the cubic parabola, its vertex, and the axis of symmetry."""
        x_vals = np.linspace(x_range[0], x_range[1], 500)
        y_vals = self.equation(x_vals)

        plt.figure(figsize=(8, 8))
        plt.plot(x_vals, y_vals, label="Cubic Parabola", color='purple', lw=2)

        x_vertex, y_vertex = self.vertex()
        plt.scatter(x_vertex, y_vertex, color='green', s=100, zorder=5, label="Vertex", marker='o')
        plt.axvline(x_vertex, color='gray', linestyle='--', label="Axis of Symmetry", lw=1)

        plt.title(f"Cubic Parabola: y = {self.A}x^3 + {self.B}x^2 + {self.C}x + {self.D}", fontsize=14)
        plt.xlabel("x", fontsize=12)
        plt.ylabel("y", fontsize=12)

        plt.legend(loc='upper left', fontsize=12)
        plt.grid(True)
        plt.show()

    def properties(self):
        """Displays key properties of the cubic parabola."""
        x_vertex, y_vertex = self.vertex()
        print(f"Vertex: ({x_vertex:.2f}, {y_vertex:.2f})")
