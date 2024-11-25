# Project: Visualization and Analysis of Curves in the Cartesian Coordinate System

## Subject:
Programming and Computational Methods

## Project Topic:
Visualization and Analysis of Curves in the Cartesian Coordinate System Using Python

## Participants:
- Dmytro Halaychuk
- Marta Hrudzynska
- Oleksandr Kondratiuk

---

## Description and Purpose of the Package:

This package is designed for the construction and analysis of various types of curves, specifically the **right strophoid**, in the Cartesian coordinate system. It allows for the generation of contour plots for curves, adjusting parameters such as axis scaling and offsets, and performing rotations of the curves within the coordinate system. The user can visualize the behavior of the curves depending on a given parameter and obtain precise data for further mathematical investigations.

## Main Functionality:
- **Curve Construction**: Generate contour plots for various types of curves, such as the right strophoid.
- **Rotation**: Rotate the curve by a specified angle.
- **Parameterization**: Modify the curve's parameters, such as axis scales and offsets.
- **Asymptote Calculation**: Compute and display the asymptotes for the curve (if applicable for the selected model).
- **Visualization**: Customize the appearance of the plot, including color, legend, and adding auxiliary lines.

## Example Usage:
```python
from thirdordercurves import RightStrophoid

# Create a RightStrophoid object with parameter a = 1
curve = RightStrophoid(a=1, visible_asymptote=True)

# Plot the curve
curve.plot()

# Rotate the curve by 45 degrees
curve.rotate(45)
curve.plot()
```
![image](https://github.com/user-attachments/assets/5cc39850-fb5a-4afc-a815-f5c0a63056ed)
```python
from quadcubicViz.secondordercurves import Ellipse

ellipse = Ellipse(5, 3, x0=2, y0=3)
ellipse.describe(3, 4)
ellipse.plot()
```
## Оптичні характеристики

Фокусна відстань: 4.00  
Ексцентриситет: 0.80  
Фокусний параметр: 1.80  
Коефіцієнт стиснення: 0.40  
Площа: 47.12  
Периметр (приблизно): 25.53  

## Фокусні радіуси для точки (3, 4)

r1: 5.10  
r2: 3.16


![image](https://github.com/user-attachments/assets/e77d48f0-916e-4b92-b0a2-4f8e625795ab)



