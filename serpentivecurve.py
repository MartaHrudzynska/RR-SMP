import numpy as np
import matplotlib.pyplot as plt

class SerpentiveCurve:
    """
    Class representing Serpentine curve.
    """
    
    def __init__(self, a=1, b=1):
        """
        Initialize Serpentine curve with parameters.
        
        Args:
            a (float): First scale parameter
            b (float): Second scale parameter
        """
        self.a = a
        self.b = b
    
    def calculate_y(self, x):
        """
        Calculate y value for given x using the curve equation.
        
        Args:
            x (float or ndarray): x-coordinate(s)
            
        Returns:
            float or ndarray: Corresponding y-coordinate(s)
        """
        return (self.a * self.b * x) / (x**2 + self.a**2)
    
    def plot(self, x_range=(-1, 1), num_points=1000, color='blue', 
            show_grid=True, show_axis=True, figsize=(10, 8)):
        """
        Plot the Serpentine curve within the specified ranges.

        Parameters:
        x_range (tuple): The range of x-values to plot as (min, max).
        num_points (int): The number of points to use for plotting the curve.
        color (str): The color of the curve line.
        show_grid (bool): Whether to display a grid on the plot.
        show_axis (bool): Whether to show axis lines (if False, axis lines are hidden).
        figsize (tuple): The size of the plot in inches as (width, height).
        
        Returns:
        None: Displays the curve plot.
        """
        x = np.linspace(x_range[0], x_range[1], num_points)
        y = self.calculate_y(x)
        
        plt.figure(figsize=figsize)
        plt.plot(x, y, color=color)
        
        if show_grid:
            plt.grid(True, linestyle='--', alpha=0.7)
        
        if not show_axis:
            plt.axis('off')
            
        plt.title('Serpentive Curve')
        plt.xlabel('x')
        plt.ylabel('y')
        plt.show()