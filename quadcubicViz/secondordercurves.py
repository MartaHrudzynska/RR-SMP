import numpy as np
import matplotlib.pyplot as plt
from typing import Tuple, List

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


class Parabola:
    """
    Class representing Parabola.
    """
    def __init__(self, a: float, b: float, c: float):
        """
        Initialize a parabola with coefficients a, b, c for the equation y = ax² + bx + c
        
        Args:
            a: coefficient of x²
            b: coefficient of x
            c: constant term
        
        Raises:
            ValueError: if coefficient 'a' equals zero
        """
        if a == 0:
            raise ValueError("Coefficient 'a' cannot be zero for a parabola")
        
        self.a = a
        self.b = b
        self.c = c
        
    def calculate_y(self, x: float) -> float:
        """
        Calculate y value for a given x
        
        Args:
            x: x-coordinate
        
        Returns:
            y value at the given x
        """
        return self.a * x**2 + self.b * x + self.c
    
    def get_vertex(self) -> Tuple[float, float]:
        """
        Find the vertex coordinates of the parabola
        
        Returns:
            Tuple (x, y) representing the vertex coordinates
        """
        x = -self.b / (2 * self.a)
        y = self.calculate_y(x)
        return (x, y)
    
    def find_roots(self) -> List[float] | None:
        """
        Find roots of the quadratic equation ax² + bx + c = 0
        
        Returns:
            None if no real roots exist
            [x] if one root exists (when discriminant = 0)
            [x1, x2] if two roots exist (sorted in ascending order)
        """
        discriminant = self.b**2 - 4*self.a*self.c
        
        if discriminant < 0:
            return None
        elif discriminant == 0:
            x = -self.b / (2*self.a)
            return [x]
        else:
            x1 = (-self.b + np.sqrt(discriminant)) / (2*self.a)
            x2 = (-self.b - np.sqrt(discriminant)) / (2*self.a)
            return sorted([x1, x2])
    
    def plot(self, x_range: Tuple[float, float] = None, show_features: bool = True, color='blue', 
            show_grid=True, show_axis=True, figsize=(10, 8)):
        """
        Visualize the parabola and its key characteristics.
    
        Args:
            x_range (Tuple[float, float], optional): The range of x-values for the plot as (x_min, x_max).
                If None, the range is automatically set to span 5 units on either side of the vertex.
            show_features (bool): Whether to highlight key features of the parabola, including:
                - Vertex
                - Axis of symmetry
                - Roots (if they exist)
            color (str): The color of the parabola curve. Defaults to 'blue'.
            show_grid (bool): Whether to display a grid on the plot. Defaults to True.
            show_axis (bool): Whether to display the axis lines. If False, axis lines are hidden. Defaults to True.
            figsize (Tuple[float, float]): The size of the plot in inches as (width, height). Defaults to (10, 8).
    
        Returns:
            None: Displays the plot of the parabola.
    
        Additional Features:
            - The vertex is marked with a red dot and labeled as 'Vertex'.
            - The axis of symmetry is shown as a green dashed line.
            - The roots (if any) are marked with black dots and labeled as 'Roots'.
        """
        if x_range is None:
            vertex_x = self.get_vertex()[0]
            x_range = (vertex_x - 5, vertex_x + 5)
            
        x = np.linspace(x_range[0], x_range[1], 1000)
        y = [self.calculate_y(xi) for xi in x]
        
        plt.figure(figsize=figsize)
        plt.grid(True, linestyle='--', alpha=0.7)
        
        # Main parabola plot
        plt.plot(x, y, label=f'y = {self.a}x² + {self.b}x + {self.c}', color=color)
        
        if show_features:
            vertex = self.get_vertex()
            plt.plot(vertex[0], vertex[1], 'ro', label='Vertex')

            plt.axvline(x=vertex[0], color='g', linestyle='--', alpha=0.5, 
                       label='Axis of symmetry')
            
            roots = self.find_roots()
            if roots:
                root_y = [0] * len(roots)
                plt.plot(roots, root_y, 'ko', label='Roots')
        
        if show_grid:
            plt.grid(True, linestyle='--', alpha=0.7)
        
        if not show_axis:
            plt.axis('off')
            

        plt.xlabel('x')
        plt.ylabel('y')
        plt.title('Parabola Plot')
        plt.legend()
        plt.show()

class Hyperbole:
    """
    A class for constructing and visualizing a hyperbola of the second order in a Cartesian coordinate system.
    """

    def __init__(self, a, b, num_points=1000, x_offset=0, y_offset=0, visible_asymptote=False):
        self.Y_rot = None
        self.X_rot = None
        self.a = a  # The parameter of the hyperbola
        self.b = b  # The parameter of the hyperbola
        self.num_points = num_points  # Number of points in the grid
        self.x_offset = x_offset  # Horizontal offset
        self.y_offset = y_offset  # Vertical offset
        self.visible_asymptote = visible_asymptote  # Flag to show asymptotes

        # Create coordinate grids
        self.x = np.linspace(- 3 * a , 3 * a, num_points)
        self.y = np.linspace(-4 * b, 4 * b, num_points)
        self.X, self.Y = np.meshgrid(self.x, self.y)

    def calculate_Z(self):
        """
        Computes the values of Z based on the equation for the hyperbola.

        Returns:
            np.ndarray: A 2D array of computed Z values for the coordinate grid.
        """
        Z = (self.X ** 2 / self.a**2) - (self.Y ** 2 / self.b**2) - 1  # Standard form of a hyperbola
        return Z

    def rotate(self, angle_deg):
        """
        Rotates the contour plot by a specified angle.

        Args:
            angle_deg (float): The angle in degrees to rotate the plot.
        """
        angle_rad = np.deg2rad(angle_deg)  # Convert angle to radians
        rotation_matrix = np.array([[np.cos(angle_rad), -np.sin(angle_rad)],
                                    [np.sin(angle_rad), np.cos(angle_rad)]])  # Rotation matrix

        # Transform the X and Y coordinates
        points = np.vstack([self.X.ravel(), self.Y.ravel()])
        rotated_points = rotation_matrix @ points
        self.X_rot = rotated_points[0].reshape(self.X.shape)
        self.Y_rot = rotated_points[1].reshape(self.Y.shape)

        # Rotate the asymptote slopes
        slope_pos, slope_neg = self.find_asymptotes()
        self.slope_pos_rot, self.slope_neg_rot = self.rotate_slope(slope_pos, slope_neg, angle_rad)

    def rotate_slope(self, slope_pos, slope_neg, angle_rad):
        """
        Rotates the slopes of the asymptotes by the given angle.

        Args:
            slope_pos (float): The original positive slope of the asymptote.
            slope_neg (float): The original negative slope of the asymptote.
            angle_rad (float): The angle in radians by which to rotate the slopes.

        Returns:
            tuple: The new rotated slopes (positive and negative).
        """
        # Convert the slopes to angle form (tan⁻¹(slope)) and add the rotation angle
        angle_pos = np.arctan(slope_pos)
        angle_neg = np.arctan(slope_neg)

        angle_pos_rot = angle_pos + angle_rad
        angle_neg_rot = angle_neg + angle_rad

        # Convert back to slope form (tan of the rotated angles)
        slope_pos_rot = np.tan(angle_pos_rot)
        slope_neg_rot = np.tan(angle_neg_rot)

        return slope_pos_rot, slope_neg_rot

    def find_asymptotes(self):
        """
        Finds the asymptotes for the hyperbola.

        Returns:
            tuple: The slopes of the asymptotes (positive and negative).
        """
        return self.b / self.a, -self.b / self.a  # Asymptote slopes for a standard hyperbola

    def plot(self):
        """
        Generates a contour plot of the hyperbola.
        Displays the curve, coordinate axes, and optional asymptotes.
        """
        Z = self.calculate_Z()  # Compute Z values

        # Use rotated coordinates if rotation is applied
        X_plot = self.X_rot if hasattr(self, 'X_rot') else self.X
        Y_plot = self.Y_rot if hasattr(self, 'Y_rot') else self.Y

        # Create the plot
        plt.figure(figsize=(12, 8))
        plt.contour(X_plot, Y_plot, Z, levels=[0], colors='blue')  # Plot the contour for Z = 0

        # Add title and axis labels
        plt.title(f'Contour plot of the hyperbola (a={self.a}, b={self.b})')
        plt.xlabel('X-axis')
        plt.ylabel('Y-axis')

        # Add coordinate axes
        plt.axhline(0, color='black', linewidth=0.7)  # Horizontal axis
        plt.axvline(0, color='black', linewidth=0.7)  # Vertical axis

        # Add asymptotes if specified
        if self.visible_asymptote:
            # Asymptotes are lines y = m * x for the positive and negative slopes
            x_vals = np.linspace(-self.a * 5, self.a * 5, self.num_points)
            y_pos = self.slope_pos_rot * x_vals
            y_neg = self.slope_neg_rot * x_vals
            plt.plot(x_vals, y_pos, 'r--', label=f'Asymptote y = {self.slope_pos_rot:.2f}x')
            plt.plot(x_vals, y_neg, 'r--', label=f'Asymptote y = {self.slope_neg_rot:.2f}x')

        # Enable grid
        plt.grid(True)

        # Add legend if asymptotes are visible
        if self.visible_asymptote:
            plt.legend()

        # Set axis limits to compress the appearance
        plt.xlim(-self.a * 3 , self.a * 3 )  # Compress x-axis
        plt.ylim(-self.b * 1.5 , self.b * 1.5 )  # Compress y-axis

        # Adjust layout for better fit
        plt.tight_layout()

        # Display the plot
        plt.show()
