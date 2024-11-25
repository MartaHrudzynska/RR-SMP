import numpy as np
import matplotlib.pyplot as plt


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
