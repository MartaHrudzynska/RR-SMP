import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import fmin
from scipy.integrate import quad


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

class ConchoidOfDeSluze:
    def __init__(self, a):
        """
        Initializes the Conchoid of de Sluze with parameter for customization.

        Args:
        - a: The parameter defining the shape of the curve.
        """
        self.a = a
       
    def equation(self, x, y):
        """
        Computes the equation (x - 1)(x^2 + y^2) = ax^2

        Args:
        - x: x-coordinate.
        - y: y-coordinate.
        
        Returns:
        - The value of the equation for the given x and y.
        """
        return (x - 1) * (x**2 + y**2) - self.a * x**2


    def plot(self, x_range=(-10, 10), y_range=(-10, 10)):
        """
        Plots the Conchoid of de Sluze.

        Args:
        - x_range: Range for x-axis (tuple of min and max values).
        - y_range: Range for y-axis (tuple of min and max values).
        """
        x_vals = np.linspace(x_range[0], x_range[1], 500)
        y_vals = np.linspace(y_range[0], y_range[1], 500)
        X, Y = np.meshgrid(x_vals, y_vals)
        Z = self.equation(X, Y)

        plt.figure(figsize=(5, 5))
        plt.contour(X, Y, Z, levels=[0], colors='purple', linewidths=2)  
        plt.title(f"Conchoid of de Sluze: a = {self.a}", fontsize=14)
        plt.xlabel("x", fontsize=12)
        plt.ylabel("y", fontsize=12)
        plt.grid(True)
        plt.show()

    def plots_with_diff_a(self, row: int, col: int):
        """
        Plots the Conchoid of de Sluze for many parametrs a to see how the conchoid changes with different a.

        Args:
        - row: number of rows.
        - col: number of columns.
        """
        a = self.a
        if (len(a) != row*col):
            print("Wrong number of parameters a.")
        else:
            fig, axes = plt.subplots(row, col, figsize=(15,10))
            x_vals = np.linspace(-10, 10, 500)
            y_vals = np.linspace(-10, 10, 500)
            X, Y = np.meshgrid(x_vals, y_vals)
            for i, a in enumerate(a):
                row, col = divmod(i, 3) 
                conchoid = ConchoidOfDeSluze(a)
                Z = conchoid.equation(X, Y)
                ax = axes[row, col]
                ax.contour(X, Y, Z, levels=[0], colors='purple', linewidths=2)
                ax.set_title(f"a = {a}", fontsize=12)
                ax.set_xlabel("x")
                ax.set_ylabel("y")
                ax.grid(True)
            plt.tight_layout()
            plt.show()

class CubicParabola:
    def __init__(self, A, B, C, D):
        """
        Initializes a cubic parabola with the given coefficients.

        Args:
        - A, B, C, D: Coefficients for the equation y = Ax^3 + Bx^2 + Cx + D
        """
        self.A = A
        self.B = B
        self.C = C
        self.D = D

    def equation(self, x):
        """Computes the y-value for a given x based on the cubic parabola equation."""
        return self.A * x**3 + self.B * x**2 + self.C * x + self.D

    def derivative(self, x):
        """Computes the derivative (slope) of the cubic equation at x."""
        return 3 * self.A * x**2 + 2 * self.B * x + self.C

    def vertex(self):
        """Finds the vertex of the cubic parabola by minimizing the derivative (using numerical optimization)."""
        critical_x = fmin(lambda x: np.abs(self.derivative(x)), 0, disp=False)[0]   
        critical_y = self.equation(critical_x)
        return (critical_x, critical_y)

    
    def plot(self, x_range=(-10, 10)):
        """Plots the cubic parabola, its vertex, and the axis of symmetry."""
        x_vals = np.linspace(x_range[0], x_range[1], 500)
        y_vals = self.equation(x_vals)

        plt.figure(figsize=(8, 8))
        plt.plot(x_vals, y_vals, label="Cubic Parabola", color='purple', lw=2)

        x_vertex, y_vertex = self.vertex()
        plt.scatter(x_vertex, y_vertex, color='green', s=100, zorder=5, label="Vertex", marker='o')
        plt.axvline(x_vertex, color='gray', linestyle='--', label="Axis of Symmetry", lw=1)

        plt.title(f"Cubic Parabola: y = {self.A}x^3 + {self.B}x^2 + {self.C}x + {self.D}", fontsize=14)
        plt.xlabel("x", fontsize=12)
        plt.ylabel("y", fontsize=12)

        plt.legend(loc='upper left', fontsize=12)
        plt.grid(True)
        plt.show()

    def properties(self):
        """Displays key properties of the cubic parabola."""
        x_vertex, y_vertex = self.vertex()
        print(f"Vertex: ({x_vertex:.2f}, {y_vertex:.2f})")


class FoliumOfDescartes:
    def __init__(self, a, visible_asymptote=True):
        """
        Initializes the Folium of Descartes with parameter for customization.

        Args:
        - a: The parameter defining the shape of the curve.
        - visible_asymptote: Whether to display the asymptote (default: True).
        """
        self.a = a
        self.visible_asymptote = visible_asymptote
       
    def equation(self, x, y):
        """
        Computes the equation x^3 + y^3 = 3axy

        Args:
        - x: x-coordinate.
        - y: y-coordinate.
        
        Returns:
        - The value of the equation for the given x and y.
        """
        return x**3 + y**3-3*x*y*self.a


    def plot(self, x_range=(-10, 10), y_range=(-10, 10)):
        """
        Plots the Folium of Descartes and its asymptote.

        Args:
        - x_range: Range for x-axis (tuple of min and max values).
        - y_range: Range for y-axis (tuple of min and max values).
        """
        x_vals = np.linspace(x_range[0], x_range[1], 500)
        y_vals = np.linspace(y_range[0], y_range[1], 500)
        X, Y = np.meshgrid(x_vals, y_vals)
        Z = self.equation(X, Y)

        plt.figure(figsize=(8, 8))
        plt.contour(X, Y, Z, levels=[0], colors='purple', linewidths=2)  
        
        if self.visible_asymptote:
            plt.plot(x_vals, -x_vals-self.a, color='red', linestyle='--', label=f"Asymptote: y = - x - {self.a}")

        plt.title(f"Folium of Descartes: x^3 + y^3 = {self.a}xy", fontsize=14)
        plt.xlabel("x", fontsize=12)
        plt.ylabel("y", fontsize=12)
        plt.legend(loc='upper left', fontsize=12)
        plt.grid(True)
        plt.show()

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

class TrisectrixOfMaclaurin:
    """
    A class for constructing and visualizing the Trisectrix of Maclaurin in a Cartesian coordinate system.

    Attributes:
        a (float): The parameter of the trisectrix, affecting the curve's shape.
        num_points (int): The number of points used to create the coordinate grid (default is 1000).
        x_scale (float): Scaling factor for the X-axis.
        y_scale (float): Scaling factor for the Y-axis.
        x_offset (float): Horizontal offset applied to the curve.
        y_offset (float): Vertical offset applied to the curve.
        visible_asymptote (bool): Whether to display the asymptote on the plot.

    Methods:
        __init__(a, num_points=1000, x_scale=1, y_scale=1, x_offset=0, y_offset=0, visible_asymptote=False):
            Initializes the class with the given parameters and sets up the coordinate grid.

        calculate_Z():
            Computes the values of Z based on the trisectrix equation for the current configuration.

        rotate(angle_deg):
            Rotates the coordinate grid and curve by a specified angle in degrees.

        find_asymptotes():
            Calculates the location of the asymptote for the trisectrix.

        plot():
            Generates a contour plot of the trisectrix, including axes, grid, and optional asymptote.
    """

    def __init__(self, a, num_points=1000, x_scale=1, y_scale=1, x_offset=0, y_offset=0, visible_asymptote=False):
        self.Y_rot = None
        self.X_rot = None
        self.a = a  # The parameter of the trisectrix
        self.num_points = num_points  # Number of points in the grid
        self.x_scale = x_scale  # Scaling factor for the X-axis
        self.y_scale = y_scale  # Scaling factor for the Y-axis
        self.x_offset = x_offset  # Horizontal offset
        self.y_offset = y_offset  # Vertical offset
        self.visible_asymptote = visible_asymptote  # Flag to show asymptote

        # Create coordinate grids
        self.x = np.linspace(-2 * a, 16 * a, num_points)
        self.y = np.linspace(-16 * a, 16 * a, num_points)
        self.X, self.Y = np.meshgrid(self.x, self.y)

    def calculate_Z(self):
        """
        Computes the values of Z based on the Trisectrix of Maclaurin equation.

        Returns:
            np.ndarray: A 2D array of computed Z values for the coordinate grid.
        """
        Z = 2 * (self.X * self.x_scale + self.x_offset) * \
            ((self.X * self.x_scale + self.x_offset) ** 2 + (self.Y * self.y_scale - self.y_offset) ** 2) - \
            self.a * (3 * (self.X * self.x_scale + self.x_offset) ** 2 - (self.Y * self.y_scale - self.y_offset) ** 2)
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

    def find_asymptotes(self):
        """
        Finds the asymptote for the Trisectrix of Maclaurin.

        Returns:
            float: The x-coordinate of the asymptote.
        """
        return -self.a / 2 + self.x_offset

    def plot(self):
        """
        Generates a contour plot of the Trisectrix of Maclaurin.

        Displays the curve, coordinate axes, and an optional asymptote.
        """
        Z = self.calculate_Z()  # Compute Z values

        # Use rotated coordinates if rotation is applied
        X_plot = self.X_rot if hasattr(self, 'X_rot') else self.X
        Y_plot = self.Y_rot if hasattr(self, 'Y_rot') else self.Y

        # Create the plot
        plt.figure(figsize=(12, 8))
        plt.contour(X_plot, Y_plot, Z, levels=[0], colors='blue')  # Plot the contour for Z = 0

        # Add title and axis labels
        plt.title(f'Contour plot of the Trisectrix of Maclaurin (a={self.a})')
        plt.xlabel('X-axis')
        plt.ylabel('Y-axis')

        # Add coordinate axes
        plt.axhline(0, color='black', linewidth=0.7)  # Horizontal axis
        plt.axvline(0, color='black', linewidth=0.7)  # Vertical axis

        # Add asymptote if specified
        if self.visible_asymptote:
            asymptote = self.find_asymptotes()
            plt.axvline(x=asymptote, color='red', linestyle='--', label=f'Asymptote x = {asymptote}')

        # Enable grid
        plt.grid(True)

        # Add legend if asymptote is visible
        if self.visible_asymptote:
            plt.legend()

        # Set axis limits to ensure the plot fits nicely
        plt.xlim(-self.a * 2, self.a * 2)
        plt.ylim(-self.a * 2, self.a * 2)

        # Adjust layout for better fit
        plt.tight_layout()

        # Display the plot
        plt.show()

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
    
    def plot(self, x_range=(-1, 1), num_points=1000, color='blue', 
            show_grid=True, show_axis=True, figsize=(10, 8)):
        """
        Plot the Witch of Agnesi curve within the specified ranges.

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
            
        plt.title('Witch of Agnesi Curve')
        plt.xlabel('x')
        plt.ylabel('y')
        plt.axis('equal')
        plt.show()

class RightStrophoid:
    """
    A class for constructing and visualizing the right strophoid curve in a Cartesian coordinate system.

    Attributes:
        a (float): The parameter affecting the curve's shape.
        num_points (int): The number of points used to create the coordinate grid (default is 1000).
        x_scale (float): Scaling factor for the X-axis.
        y_scale (float): Scaling factor for the Y-axis.
        x_offset (float): Horizontal offset applied to the curve.
        y_offset (float): Vertical offset applied to the curve.
        visible_asymptote (bool): Whether to display the asymptote on the plot.

    Methods:
        __init__(a, num_points=1000, x_scale=1, y_scale=1, x_offset=0, y_offset=0, visible_asymptote=False):
            Initializes the class with the given parameters and sets up the coordinate grid.

        calculate_Z():
            Computes the values of Z based on the equation for the current configuration.

        rotate(angle_deg):
            Rotates the coordinate grid and curve by a specified angle in degrees.

        find_asymptotes():
            Calculates the location of the asymptote for the curve.

        plot():
            Generates a contour plot of the curve, including axes, grid, and optional asymptote.
    """

    def __init__(self, a, num_points=2000, x_scale=1, y_scale=1, x_offset=0, y_offset=0, visible_asymptote=False):
        self.Y_rot = None
        self.X_rot = None
        self.a = a  # The parameter of the curve
        self.num_points = num_points  # Number of points in the grid
        self.x_scale = x_scale  # Scaling factor for the X-axis
        self.y_scale = y_scale  # Scaling factor for the Y-axis
        self.x_offset = x_offset  # Horizontal offset
        self.y_offset = y_offset  # Vertical offset
        self.visible_asymptote = visible_asymptote  # Flag to show asymptote

        # Create coordinate grids
        self.x = np.linspace(-4 * a, 16 * a, num_points)
        self.y = np.linspace(-16 * a, 16 * a, num_points)
        self.X, self.Y = np.meshgrid(self.x, self.y)

    def calculate_Z(self):
        """
        Computes the values of Z based on the equation for the right strophoid.

        Returns:
            np.ndarray: A 2D array of computed Z values for the coordinate grid.
        """
        # Right strophoid equation: Z = (Y^2 * (a + X) - X^2 * (a - X - X)) 
        Z = (self.Y * self.y_scale - self.y_offset) ** 2 * (self.a + self.X * self.x_scale) - (self.X * self.x_scale) ** 2 * (self.a - self.X - self.X * self.x_scale)
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

    def find_asymptotes(self):
        """
        Finds the asymptote for the right strophoid.

        Returns:
            float: The x-coordinate of the asymptote.
        """
        return -self.a / 2 + self.x_offset

    def plot(self):
        """
        Generates a contour plot of the curve.

        Displays the right strophoid, coordinate axes, and an optional asymptote.
        """
        Z = self.calculate_Z()  # Compute Z values

        # Use rotated coordinates if rotation is applied
        X_plot = self.X_rot if hasattr(self, 'X_rot') else self.X
        Y_plot = self.Y_rot if hasattr(self, 'Y_rot') else self.Y

        # Create the plot
        plt.figure(figsize=(12, 8))
        plt.contour(X_plot, Y_plot, Z, levels=[0], colors='blue')  # Plot the contour for Z = 0

        # Add title and axis labels
        plt.title(f'Contour plot of the right strophoid curve (a={self.a})')
        plt.xlabel('X-axis')
        plt.ylabel('Y-axis')

        # Add coordinate axes
        plt.axhline(0, color='black', linewidth=0.7)  # Horizontal axis
        plt.axvline(0, color='black', linewidth=0.7)  # Vertical axis

        # Enable grid
        plt.grid(True)

        # Add legend if asymptote is visible
        if self.visible_asymptote:
            plt.legend()

        # Set axis limits to ensure the plot fits nicely
        plt.xlim(-self.a * 2, self.a * 2)
        plt.ylim(-self.a * 4, self.a * 4)

        # Adjust layout for better fit
        plt.tight_layout()

        # Display the plot
        plt.show()
