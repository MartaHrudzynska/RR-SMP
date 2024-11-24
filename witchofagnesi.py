import numpy as np
import matplotlib.pyplot as plt


class WitchOfAgnesi:
    """
    Class representing Witch of Agnesi curve.
    """
    
    def __init__(self, a=1):
        """
        Initialize Witch of Agnesi curve with parameter.
        
        Args:
            a (float): Scale parameter
        """
        self.a = a
    
    def calculate_y(self, x):
        """
        Calculate y value for given x using the curve equation.
        
        Args:
            x (float or ndarray): x-coordinate(s)
            
        Returns:
            float or ndarray: Corresponding y-coordinate(s)
        """
        return (8 * self.a**3) / (x**2 + 4 * self.a**2)
    
    def plot(self, x_range=(-10, 10), num_points=1000):
        """
        Plot the curve using matplotlib.
        
        Args:
            x_range (tuple): Range of x values (min, max)
            num_points (int): Number of points for plotting
        """
        x = np.linspace(x_range[0], x_range[1], num_points)
        y = self.calculate_y(x)
        
        plt.figure(figsize=(10, 6))
        plt.plot(x, y)
        plt.grid(True)
        plt.title("Witch of Agnesi")
        plt.xlabel('x')
        plt.ylabel('y')
        plt.axhline(y=0, color='k', linestyle='-', alpha=0.3)
        plt.axvline(x=0, color='k', linestyle='-', alpha=0.3)
        plt.show()
