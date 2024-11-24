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
    
    def plot(self, t_range=(-1, 1), num_points=1000, color='blue', 
            show_grid=True, show_axis=True, figsize=(10, 8)):
        """
        Visualize the curve.
        :param t_range: range of the parameter t as a tuple (start, end)
        :param num_points: number of points to calculate
        :param color: color of the curve
        :param show_grid: whether to show the grid
        :param show_axis: whether to show the axes
        :param figsize: size of the figure
        """
        x, y = self.calculate_points(t_range[0], t_range[1], num_points)
        
        plt.figure(figsize=figsize)
        plt.plot(x, y, color=color, label=f'Tschirnhausen Cubic (a={self.a})')
        
        if show_grid:
            plt.grid(True, linestyle='--', alpha=0.7)
        
        if not show_axis:
            plt.axis('off')
            
        plt.title('Tschirnhausen Cubic Curve')
        plt.xlabel('x')
        plt.ylabel('y')
        plt.legend()
        plt.axis('equal')  # Maintain aspect ratio
        plt.show()
