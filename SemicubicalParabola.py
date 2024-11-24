import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import quad


class SemicubicParabola:
    """
    A class for constructing and visualizing a semicubic parabola in the form y^2 = ax^3.

    Attributes:
        a (float): The parameter `a` of the semicubic parabola, determining its shape.
        num_points (int): The number of points used to generate the plot.
        x_offset (float): Horizontal offset applied to the curve.
        y_offset (float): Vertical offset applied to the curve.

    Methods:
        __init__(a, num_points=1000, x_offset=0, y_offset=0):
            Initializes the semicubic parabola with the given parameters.

        arc_length(t0, t1):
            Calculates the arc length of the semicubic parabola between two parameter values t0 and t1.

        plot():
            Generates a plot of the semicubic parabola.
    """

    def __init__(self, a, num_points=1000, x_offset=0, y_offset=0):
        self.a = a  # Parameter a of the semicubic parabola
        self.num_points = num_points  # Number of points for the plot
        self.x_offset = x_offset  # Horizontal offset
        self.y_offset = y_offset  # Vertical offset

        # Generate the range for x (x >= 0 for semicubic parabola)
        self.x = np.linspace(0, 4 * a, num_points)
        self.y_positive = np.sqrt(self.a * self.x ** 3)  # Positive branch of y

    def arc_length(self, t0, t1):
        """
        Calculates the arc length of the semicubic parabola between t0 and t1.

        Args:
            t0 (float): The starting parameter value.
            t1 (float): The ending parameter value.

        Returns:
            float: The arc length of the parabola between t0 and t1.
        """

        def integrand(t):
            dxdt = 1
            dydt = 1.5 * self.a * np.sqrt(t)
            return np.sqrt(dxdt ** 2 + dydt ** 2)

        arc_length, _ = quad(integrand, t0, t1)
        return arc_length

    def plot(self):
        """
        Generates a plot of the semicubic parabola (positive branch only).
        """
        # Apply offsets
        x_plot = self.x + self.x_offset
        y_plot = self.y_positive + self.y_offset

        # Create the plot
        plt.figure(figsize=(8, 6))
        plt.plot(x_plot, y_plot, label="y^2 = ax^3", color="blue")

        # Add title and axis labels
        plt.title(f'Graph of the Semicubic Parabola y^2 = ax^3 (a={self.a})')
        plt.xlabel('X-axis')
        plt.ylabel('Y-axis')

        # Add coordinate axes
        plt.axhline(0, color='black', linewidth=0.7)  # Horizontal axis
        plt.axvline(0, color='black', linewidth=0.7)  # Vertical axis

        # Enable grid
        plt.grid(True)

        # Set limits
        plt.xlim(0, self.a * 4)
        plt.ylim(0, self.a * 4)

        plt.legend()
        plt.show()


# Example usage
semicubic_parabola = SemicubicParabola(a=0.5, num_points=1000, x_offset=0, y_offset=0)
semicubic_parabola.plot()

# Example arc length calculation
arc_length = semicubic_parabola.arc_length(0, 2)
print("Arc length between t=0 and t=2:", arc_length)
