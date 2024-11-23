import numpy as np
import matplotlib.pyplot as plt

class CissoidOfDiocles:
    def __init__(self, a, visible_asymptote=True, visible_circle=True):
        """
        Initializes the Folium of Descartes with parameter for customization.

        Args:
        - a: The parameter defining the shape of the curve.
        - visible_asymptote: Whether to display the asymptote (default: True).
        - visible_circle: A cissoid intersects the auxiliary circle at points that belong to the diameter of this circle.(default: True).

        """
        self.a = a
        self.visible_asymptote = visible_asymptote
        self.visible_circle = visible_circle

    def equation(self, x, y):
        return (x**2 + y**2) * x - 2 * self.a * y**2

    def area(self):
        s = 3 * np.pi * (self.a)**2
        print(f"The area between the cissoid and the asymptote: {s}")
        return s

    def plot(self):
        """
        Plots the Cissoid of Diocles.
        """
        y_range=(-2 * self.a, 2 * self.a)
        # Цисоїда
        x_vals = np.linspace(0, 2 * self.a, 1000)
        y_vals = np.linspace(y_range[0], y_range[1], 1000)
        X, Y = np.meshgrid(x_vals, y_vals)
        Z = self.equation(X, Y)

        plt.figure(figsize=(5, 10))
        plt.contour(X, Y, Z, levels=[0], colors='purple', linewidths=2) 
        
        if self.visible_asymptote:
            plt.axvline(x=2*self.a, color='gray', linestyle='--', label=f'Asymptote x = {2 * self.a}') 
        if self.visible_circle:
            circle = plt.Circle((self.a, 0), self.a, color='yellow', fill=False, linewidth=2, label=f'Circle with radius a={self.a}')
            plt.gca().add_artist(circle)
        plt.xlim(min(x_vals)-0.1*self.a, max(x_vals) + 0.1*self.a)  
        plt.ylim(min(y_vals)-0.2*self.a, max(y_vals) + 0.2*self.a) 
        plt.title(f"Cissoid of Diocles with a = {self.a}")
        plt.xlabel("x")
        plt.ylabel("y")
        plt.legend()
        plt.grid(True)
        plt.show()
