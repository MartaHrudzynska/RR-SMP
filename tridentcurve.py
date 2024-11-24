import numpy as np
import matplotlib.pyplot as plt


class TridentCurve:
    """
    Class representing Trident curve.
    """
    def __init__(self, a=2, b=2, c=2, d=2):
        """
        Initializes the trident curve with parameters.
        
        Parameters:
        a, b, c, d (float): Parameters that influence the shape of the curve.
        """
        self.a = a
        self.b = b
        self.c = c
        self.d = d
        
    def calculate_points(self, x):
        """
        Calculate the value of y for a given x.
        
        :param x: The x-value.
        :return: The corresponding y-value.
        :raises ValueError: If x is zero to avoid division by zero.
        """
        if x == 0:
            raise ValueError("x cannot be zero (to avoid division by zero).")
        return self.d / x - self.a * x**2 - self.b * x - self.c
    
    def plot(self, x_range=(-3, 3), num_points=1000):
        """
        Plot the graph of the function.
        
        :param x_range: Tuple (min, max) for the x-axis.
        :param num_points: Number of points to plot the graph.
        """
        x = np.linspace(x_range[0], x_range[1], num_points)
        x = x[x != 0]  # Exclude x = 0 to avoid division by zero

        y = self.d / x - self.a * x**2 - self.b * x - self.c

        plt.figure(figsize=(8, 6))
        plt.plot(x, y, label=r"$y = \frac{d}{x} - ax^2 - bx - c$", color='purple')
        plt.axhline(0, color='black', linewidth=0.8, linestyle='--', label='y = 0')
        plt.axvline(0, color='black', linewidth=0.8, linestyle='--', label='x = 0 (asymptote)')
        plt.title("Graph of the Curve: $y = \\frac{d}{x} - ax^2 - bx - c$", fontsize=14)
        plt.xlabel("x", fontsize=12)
        plt.ylabel("y", fontsize=12)
        plt.legend()
        plt.grid()
        plt.ylim(-10, 10)  # y-axis limits for better visibility
        plt.show()





