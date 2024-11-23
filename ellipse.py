import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

class Ellipse:
    """
    A class for working with ellipses.

    The ellipse is defined by the equation:
        ((x - x0)^2 / a^2) + ((y - y0)^2 / b^2) = 1
    
    Args:
        a (float): The semi-major axis of the ellipse. Must be a positive number.
        b (float): The semi-minor axis of the ellipse. Must be a positive number.
        x0 (float): The x-coordinate of the center of the ellipse. Default is 0.
        y0 (float): The y-coordinate of the center of the ellipse. Default is 0.

    Raises:
        ValueError: If a <= 0 or b <= 0.

    Attributes:
        a (float): The semi-major axis of the ellipse.
        b (float): The semi-minor axis of the ellipse.
        x0 (float): The x-coordinate of the center of the ellipse.
        y0 (float): The y-coordinate of the center of the ellipse.
        c (float): The focal distance (calculated automatically).
    """

    def __init__(self, a: float, b: float, x0: float = 0, y0: float = 0):
        if a <= 0 or b <= 0:
            raise ValueError("Parameters 'a' and 'b' must be positive numbers.")
        self.a = max(a, b)  
        self.b = min(a, b)
        self.x0 = x0 
        self.y0 = y0
        self.c = np.sqrt(self.a**2 - self.b**2) 

    def focal_distance(self):
        """Returns the focal distance."""
        return self.c

    def eccentricity(self):
        """Calculates the eccentricity of the ellipse."""
        return self.c / self.a

    def focal_radii(self, x: float, y: float):
        """
        Calculates the focal radii (distances to the foci) for a given point (x, y).
        
        Args:
            x (float): x-coordinate of the point.
            y (float): y-coordinate of the point.

        Returns:
            tuple: Distances to the foci.
        """
        foci = self.focus_points()
        r1 = np.sqrt((x - foci[0][0])**2 + (y - foci[0][1])**2)
        r2 = np.sqrt((x - foci[1][0])**2 + (y - foci[1][1])**2)
        return r1, r2

    def focal_parameter(self):
        """Calculates the focal parameter of the ellipse."""
        return self.b**2 / self.a

    def compression_ratio(self):
        """Calculates the compression ratio (ellipticity) of the ellipse."""
        return 1 - (self.b / self.a)

    def area(self):
        """Calculates the area of the ellipse."""
        return np.pi * self.a * self.b

    def perimeter(self):
        """
        Calculates the approximate perimeter of the ellipse using Ramanujan's formula.
        
        Returns:
            float: The perimeter of the ellipse.
        """
        h = ((self.a - self.b)**2) / ((self.a + self.b)**2)
        return np.pi * (self.a + self.b) * (1 + (3 * h) / (10 + np.sqrt(4 - 3 * h)))

    def focus_points(self):
        """Calculates the coordinates of the foci."""
        return [(self.x0 - self.c, self.y0), (self.x0 + self.c, self.y0)]

    def describe(self, x, y):
        print(f"Focal distance: {self.focal_distance():.2f}")
        print(f"Eccentricity: {self.eccentricity():.2f}")
        print(f"Focal parameter: {self.focal_parameter():.2f}")
        print(f"Compression ratio: {self.compression_ratio():.2f}")
        print(f"Area: {self.area():.2f}")
        print(f"Perimeter (approximate): {self.perimeter():.2f}")
        r1, r2 = self.focal_radii(x, y)
        print(f"Focal radii for point ({x}, {y}): r1 = {r1:.2f}, r2 = {r2:.2f}")


    def plot(self):
        """Visualizes the ellipse, its foci, and the directrices, including the center."""
        theta = np.linspace(0, 2 * np.pi, 500)
        x = self.a * np.cos(theta) + self.x0
        y = self.b * np.sin(theta) + self.y0

        plt.figure(figsize=(8, 8))
        plt.plot(x, y, label="Ellipse", color='purple', lw=2)

        p = self.focal_parameter()
        plt.axvline(self.x0 - self.a**2 / self.c, color='blue', linestyle='--', label="Directrix", lw=1)
        plt.axvline(self.x0 + self.a**2 / self.c, color='blue', linestyle='--', lw=1)

        foci = self.focus_points()
        plt.scatter([foci[0][0], foci[1][0]], [foci[0][1], foci[1][1]], 
                    color='red', s=100, label='Foci', zorder=5)

        plt.scatter(self.x0, self.y0, color='green', s=100, zorder=5, label="Center", marker='o')

        plt.gca().set_aspect('equal', adjustable='box')
        plt.title(f"Ellipse: a={self.a}, b={self.b}, center=({self.x0}, {self.y0})", fontsize=14)
        plt.xlabel("x", fontsize=12)
        plt.ylabel("y", fontsize=12)
        plt.legend(loc='upper right', bbox_to_anchor=(1, 1), fontsize=12)
        plt.grid(True)
        plt.show()

