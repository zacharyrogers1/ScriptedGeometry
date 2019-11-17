

import math

#-------------------------------------VARIABLES TO BE CHANGED-------------------------------------
inner_chamber_diameter = .05
inner_wall_thickness = .005
outer_wall_thickness = .005
number_of_beans_per_chamber = 10 
height_of_single_revolution = .2

inner_radius_bean = .002
distance_between_inner_outer_bean = .001
arc_length_in_degrees = 60.0
bean_area = .0000034034
#-------------------------------------------------------------Functions-----------------------------------------------------------------

def solve_quadratic(a, b, c):
    d = (b**2) + (-4.0*a*c)
    solution_plus = (-b + d**0.5)/(2.0*a)
    solution_minus = (-b - d**0.5)/(2.0*a)
    return (solution_plus, solution_minus)

# x = math.pi*((distance_between_inner_outer_bean/2)**2)
# y = math.pi*(distance_between_inner_outer_bean + inner_radius_bean)**2 - math.pi*inner_radius_bean**2
# z = y*arc_length_in_degrees/360
# bean_area = x+z
# print(bean_area)

a = ((arc_length_in_degrees/360) +0.25)
b = (2*inner_radius_bean*arc_length_in_degrees)/360
c = -bean_area/math.pi

distance_between_inner_outer_bean = solve_quadratic(a,b,c)
print(distance_between_inner_outer_bean[0])

