import numpy as np
import matplotlib.pyplot as plt
from typing import Tuple, Optional, List

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
    
    def plot(self, x_range: Tuple[float, float] = None, show_features: bool = True):
        """
        Visualize the parabola and its characteristics
        
        Args:
            x_range: tuple (x_min, x_max) defining the range for plotting
                    If None, automatically sets range around vertex
            show_features: whether to show additional features (vertex, axis of symmetry, roots)
        """
        if x_range is None:
            vertex_x = self.get_vertex()[0]
            x_range = (vertex_x - 5, vertex_x + 5)
            
        x = np.linspace(x_range[0], x_range[1], 1000)
        y = [self.calculate_y(xi) for xi in x]
        
        plt.figure(figsize=(10, 6))
        plt.grid(True, linestyle='--', alpha=0.7)
        
        # Main parabola plot
        plt.plot(x, y, 'b-', label=f'y = {self.a}x² + {self.b}x + {self.c}')
        
        if show_features:
            vertex = self.get_vertex()
            plt.plot(vertex[0], vertex[1], 'ro', label='Vertex')

            plt.axvline(x=vertex[0], color='g', linestyle='--', alpha=0.5, 
                       label='Axis of symmetry')
            
            roots = self.find_roots()
            if roots:
                root_y = [0] * len(roots)
                plt.plot(roots, root_y, 'ko', label='Roots')
            
        plt.xlabel('x')
        plt.ylabel('y')
        plt.title('Parabola Plot')
        plt.legend()
        plt.show()