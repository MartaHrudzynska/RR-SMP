import numpy as np
import matplotlib.pyplot as plt

class TschirnhausenCubic:
    """
    Class representing Tschirnhausen Cubic curve.
    """
    def __init__(self, a=1):
        """
        Initialization of the Tschirnhausen cubic curve.
        :param a: scaling parameter of the curve (default is 1)
        """
        self.a = a
        
    def calculate_points(self, t_start=-1, t_end=1, num_points=1000):
        """
        Calculate points on the curve.
        :param t_start: starting value of the parameter t
        :param t_end: ending value of the parameter t
        :param num_points: number of points to calculate
        :return: arrays of x and y coordinates
        """
        t = np.linspace(t_start, t_end, num_points)
        x = 3 * self.a * (t**2 - 3)
        y = self.a * t * (3 - t**2)
        return x, y
    
    def plot(self, x_range=(-1, 1), num_points=1000, color='blue', 
            show_grid=True, show_axis=True, figsize=(10, 8)):
        """
        Plot the Tschirnhausen Cubic curve within the specified ranges.

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
        x, y = self.calculate_points(x_range[0], x_range[1], num_points)
        
        plt.figure(figsize=figsize)
        plt.plot(x, y, color=color)
        
        if show_grid:
            plt.grid(True, linestyle='--', alpha=0.7)
        
        if not show_axis:
            plt.axis('off')
            
        plt.title('Tschirnhausen Cubic Curve')
        plt.xlabel('x')
        plt.ylabel('y')
        plt.axis('equal')
        plt.show()
