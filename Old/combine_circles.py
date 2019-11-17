delete_all = Selection.SelectAll()
ComponentHelper.DeleteEmptyComponents(delete_all)
Delete.Execute(delete_all)

import math
import time

initial_radius = MM(20)
tiers = 4
xyspread = MM(60)
distance_between_tiers = MM(60)
curve_resolution = 20
number_of_curves = 4**(tiers-1)

faux_origin = Point.Create(MM(0), MM(0), MM(-60))
origin_3d = Point.Create(MM(0),MM(0),MM(0))
z_axis = Direction.Create(0,0,1)
framed = Frame.Create(origin_3d, z_axis)


def create_circle(point, direction, radius):
#point: The center of the circle(Point Object)
#direction: The direction the circle is normal to (Direction Object)
#radius: The radius of the circle(float)
	circle = CircularSurface.Create(radius, direction, point)
	return circle.CreatedBody #returns IDOCobject


#Input a curve and how many points on that curve you want for example 20. Return a list of 20 point objects on that curve equally spaced
def dividing_curve(curve, num_divisions):
#curve: A SpaceClaim curve(Curve Object)
#num_divisions: How many points on the curve would you like returned.(Float)
	divided_curve_points = []
	for i in range(num_divisions+1):
		i = float(i)
		eval = curve.Evaluate(1-(i/num_divisions))
		divided_curve_points.append(eval.Point)
	return divided_curve_points #Returns a list of point objects on the curve.

#Input a parent point and have it split into 4 points in the cardinal directions. Note This function does not return anything. It only modifies the list labeled points.
def split_point(branch_pt, xysep, zsep, tier):
#branch_pt: The parent point which the 4 new points will branch off from (Single list of X,Y,Z coordinates)
#xysep: How far the next layer of points will branch out in the x and y direction. (Float)
#zsep: How far the next layer of points will be above or below the branch_pt in the z direction. (Float)
#tier: What tier/layer of points is the new layer going to be on. The origin is on tier 0 and its first children are on tier 1(Float)
	x0 = branch_pt[0]
	y0 = branch_pt[1]
	z0 = branch_pt[2]
	points[tier].append([x0+xysep, y0+xysep, z0+zsep])
	points[tier].append([x0-xysep, y0+xysep, z0+zsep])
	points[tier].append([x0-xysep, y0-xysep, z0+zsep])
	points[tier].append([x0+xysep, y0-xysep, z0+zsep])

#Input a curve and a number between 0-1. This equation will give you a point on that curve and the tangent direction of that point.
def tangent_direction(incoming_curve, evaluation_point):
#incoming_curve: A Spaceclaim curve(Curve Object)
#evaluation_point: A number between 0-1(Float)
	interim_1 = incoming_curve.EvalProportion(evaluation_point) #Find a point on the curve
	bottom_point = interim_1.Point
	interim_2 = incoming_curve.EvalProportion(evaluation_point +.001)#Find another point on the curve that is extremely close to the first point
	top_point = interim_2.Point
	X1 = bottom_point.X - top_point.X #Find the tangent direction by finding the vector between the two very close points
	Y1 = bottom_point.Y - top_point.Y
	Z1 = bottom_point.Z - top_point.Z
	direction = Direction.Create(X1,Y1,Z1)
	return (bottom_point, direction)#bottom_point: point object ||| direction: direction object

#Input a series of points and the radius of its parent circle. Return the radius which achieves the desired cross-sectional area and the shape.
#PoA 1: Create a geometry point from each list(point) that is given as an input
#PoA 2: Failsafe, Check to see if the circles are overlapping, if they are not overlapping then print an error.
#PoA 3: Run a while loop that creates overlapping circles, check their combined surface area against a reference and increase or decrease the radius until the desired cross-sectional area is achieved.
#PoA 3.1: Failsafe, If the radius has not converged after a certain number of iterations then exit the function and return nothing
#PoA 3.2: Create circles at all of the input points and select all of them
#PoA 3.3: Combine circles and find area
#PoA 3.4: Compare the combined circle area with the reference area. Evaluate how large to make the next radius or exit if correct radius found
#PoA 3.5: Delete all circles
def combine_circles(estimated_radius, parent_radius, tolerance, *points):
#estimated_radius: Where the solver will initially guess how large to make the overlapping circles(Float)
#parent_radius: How large is the parent circle radius that these branch out from? This is the area that we are trying to achieve.(Float)
#tolerance: What margin of error is acceptable when comparing to the desired cross-sectional area. (Float)
#*points: This is a variable arguement. You can add as many points as an input as you would like such as 2,3,4, or 5 points.(Single list of X,Y,Z coordinates)
	circle_centers = []
	circle_centers_check = []
	distances = []
	shifting_adder = float(2)
	new = 0 
	prev = 0
	iteration = 0
	target_area = math.pi * (parent_radius**2)
	radius_of_circle = estimated_radius	

#PoA 1: Create a geometry point from each list that is given as an input
	for p in points: #Take the variable input *points and create the centers of the circles within SpaceClaim
		circle_centers.append(Point.Create(MM(p[0]), MM(p[1]), MM(p[2])))
		circle_centers_check.append(p)
	
#PoA 2: Check to see if the circles have split, if they have then print an error
	lower_limit_radius = parent_radius * (len(circle_centers)**(-.5))
	for i in range(1, len(circle_centers)):
		a = (circle_centers_check[0][0]-circle_centers_check[i][0])**2
		b = (circle_centers_check[0][1]-circle_centers_check[i][1])**2
		c = (circle_centers_check[0][2]-circle_centers_check[i][2])**2
		d = (a + b +c)**.5
		distances.append(d)
	if(min(distances)/2 >=  lower_limit_radius):
		print("Error: The circles no longer overlap")
		return

#PoA 3: Run a while loop that creates overlapping circles, check their combined surface area against a reference and increase or decrease the radius until the desired cross-sectional area is achieved.
	while True:	
		#PoA 3.1: Failsafe, If the radius has not converged after a certain number of iterations then exit the function and return nothing
		if (iteration == 29):
			print("Error: The function combine_circles failed to find the correct radius")
			return
		
		#PoA 3.2: Create circles at all of the input points and select all of them
		sel = []
		trying_circles = []
		print("||||||||||Iteration: " + str(iteration) + "||||||||||||||||")
		print("Radius of Circle: " + str(radius_of_circle))
		print("Shifting Adder: " + str(shifting_adder))
		for i in range(len(circle_centers)): #Take the points of all the circles centers and make circles
			trying_circles.append(create_circle(circle_centers[i], z_axis, MM(radius_of_circle)))
			sel.append(Selection.Create(trying_circles[i]))	#Make a list of selections of each circle: sel[0], sel[1], sel[2]...

		#PoA 3.3: Combine circles and find area
		for i in range(1, len(circle_centers)):#Merge all of the created circles into one shape
			Combine.Merge(sel[0], sel[i])
		trying_circles[0].Faces[0].Area*(10**6)		
		combined_area = trying_circles[0].Faces[0].Area*(10**6)#find the area of the combined shape in MM2
		print("Area of shape: " + str(combined_area))
		
		#PoA 3.4: Compare the combined circle area with the reference area. Evaluate how large to make the next radius or exit if correct radius found
		if (combined_area >= target_area*(1-tolerance) and combined_area <= target_area*(1+tolerance)): #If the area of the shape is within the tolerance window return the golden radius
			correct_radius = MM(radius_of_circle)
			return correct_radius
		if(combined_area <= target_area*(1-tolerance)): #If the given radius undershoots the target area add the shifting_adder to the radius
			new = -1
			if(prev != new): #If you went from overshooting to undershooting then reduce shifting adder by half to close in on the right value
				shifting_adder = shifting_adder/2
			radius_of_circle = radius_of_circle + shifting_adder
		elif (combined_area >= target_area*(1+tolerance)):
			new = 1
			if(prev != new):
				shifting_adder = shifting_adder/2
			radius_of_circle = radius_of_circle - shifting_adder
		prev = new #This stores the last overshoot or undershoot estimate.
		iteration = iteration + 1
		
		#PoA 3.5: Delete all circles
		for i in range(len(circle_centers)):
			try: 
				Delete.Execute(sel[i]) #Because we do not know which circles were combined and which ones were not we could get errors deleting something that is not there.
			except:
				continue


creation0 = [5,5,0]
creation1 = [-5,5,0]
creation2 = [-5,-5,0]
creation3 = [5,-5,0]

help = combine_circles(4.2, 11, .01, creation0, creation1, creation2, creation3)
print(help)