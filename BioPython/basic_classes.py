#!/usr/bin/env python3
# basic_classes.py
"""Write two classes with fields and methods and for each
class definition give two examples with at least one method"""

class Circle():
    """example of class and attributes + methods"""
    def __init__(self, color, radius):
        self.color = str(color)
        self.radius = float(radius)

    def diameter(self):
        """Calculates the diameter of the circle"""
        return self.radius * 2

    def circumference(self):
        """Calculates circumference"""
        return 2 * 3.14 * self.radius

    def isRed(self):
        """Boolean to determine if circle is Red"""
        if self.color.lower() == "red":
            return True

### examples for Circle:
C1 = Circle("Orange", 12)
C2 = Circle("red", 10)

# Calculate diameter for C1
print(C1.diameter())

# Get C2's color
print(C2.color)


class GraduateStudent():
    """another example of class and attributes + methods"""
    def __init__(self, first_name, last_name, year, major):
        self.first_name = str(first_name)
        self.last_name = str(last_name)
        self.year = int(year)
        self.major = str(major)

    def year_matriculated(self):
        """Assuming it is Jan 2020 per instructions, what year they started grad school"""
        return abs(2020 - self.year)

### examples for GraduateStudent:
G1 = GraduateStudent("dj", "khaled", 3, "remixing")
G2 = GraduateStudent("Tom", "Brady", 7, "champions")

# Calculate the year matriculated for G1
print(G1.year_matriculated())

# Get G2's major:
print(G2.major)
