import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

class FoliumOfDescartes:
    def __init__(self, a, visible_asymptote=True):
        """
        Initializes the Folium of Descartes with parameter for customization.

        Args:
        - a: The parameter defining the shape of the curve.
        - visible_asymptote: Whether to display the asymptote (default: True).
        """
        self.a = a
        self.visible_asymptote = visible_asymptote
       
    def equation(self, x, y):
        """
        Computes the equation x^3 + y^3 = 3axy

        Args:
        - x: x-coordinate.
        - y: y-coordinate.
        
        Returns:
        - The value of the equation for the given x and y.
        """
        return x**3 + y**3-3*x*y*self.a


    def plot(self, x_range=(-10, 10), y_range=(-10, 10)):
        """
        Plots the Folium of Descartes and its asymptote.

        Args:
        - x_range: Range for x-axis (tuple of min and max values).
        - y_range: Range for y-axis (tuple of min and max values).
        """
        x_vals = np.linspace(x_range[0], x_range[1], 500)
        y_vals = np.linspace(y_range[0], y_range[1], 500)
        X, Y = np.meshgrid(x_vals, y_vals)
        Z = self.equation(X, Y)

        plt.figure(figsize=(8, 8))
        plt.contour(X, Y, Z, levels=[0], colors='purple', linewidths=2)  
        
        if self.visible_asymptote:
            plt.plot(x_vals, -x_vals-self.a, color='red', linestyle='--', label=f"Asymptote: y = - x - {self.a}")

        plt.title(f"Folium of Descartes: x^3 + y^3 = {self.a}xy", fontsize=14)
        plt.xlabel("x", fontsize=12)
        plt.ylabel("y", fontsize=12)
        plt.legend(loc='upper left', fontsize=12)
        plt.grid(True)
        plt.show()

