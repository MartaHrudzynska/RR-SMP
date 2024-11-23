import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

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
