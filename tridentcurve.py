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
    
    def plot(self, x_range=(-1, 1), y_range = (-10,10), num_points=1000, color='blue', 
            show_grid=True, show_axis=True, figsize=(10, 8)):
        """
        Plot the Trident curve within the specified ranges.

        Parameters:
        x_range (tuple): The range of x-values to plot as (min, max).
        y_range (tuple): The range of y-values to display on the plot as (min, max).
        num_points (int): The number of points to use for plotting the curve.
        color (str): The color of the curve line.
        show_grid (bool): Whether to display a grid on the plot.
        show_axis (bool): Whether to show axis lines (if False, axis lines are hidden).
        figsize (tuple): The size of the plot in inches as (width, height).
        
        Returns:
        None: Displays the curve plot.
        """
        x = np.linspace(x_range[0], x_range[1], num_points)
        x = x[x != 0]  # Exclude x = 0 to avoid division by zero

        y = self.d / x - self.a * x**2 - self.b * x - self.c

        plt.figure(figsize=figsize)
        plt.plot(x, y, color=color)
        
        if show_grid:
            plt.grid(True, linestyle='--', alpha=0.7)
        
        if not show_axis:
            plt.axis('off')
        
        plt.ylim(y_range)
        plt.title('Trident Curve')
        plt.xlabel('x')
        plt.ylabel('y')
        plt.axis('equal')
        plt.ylim(y_range)
        plt.show()





