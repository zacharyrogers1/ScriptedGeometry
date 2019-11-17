sel = Selection.SelectAll()
ComponentHelper.DeleteEmptyComponents(sel)
Delete.Execute(sel)

import math
import time

#-------------------------------------VARIABLES TO BE CHANGED(Only change these numbers)-------------------------------------
initial_radius = MM(20) #Initial radius of inlet that will be split
tiers = 4 #How many layers of splits 
xyspread = MM(60) #This describes how far in the xy direction each set of points splits. The larger the xy spread the more area it will take up
distance_between_tiers = MM(60) #This describes how far in the z direction the curve goes before splitting at each layer
curve_resolution = 20 #Each curve is made up of lofted circles. This number describes the number of circles that each curve is compromised of. The more circles, the higher the resolution, the more computing time needed.
#-----------------------------------------------------------------------------------------------------------------------------


#-------------------------------------Methodology-------------------------------------
#1. Create planes of points each with 4 times as many as the previous
#2. Create a curve starting from the top inlet winding through each tier and stopping at the lowest tier endpoint. You will have a collection of winding curves.
#3. Take a curve and create circles along regular intervals. Loft all of the circles in a curve. Repeat for all curves




#---------------------------------FUNCTIONAL CODE BEYOND THIS POINT-------------------------------------------
number_of_curves = 4**(tiers-1)
faux_origin = Point.Create(MM(0), MM(0), MM(-60))
origin_3d = Point.Create(MM(0),MM(0),MM(0))
z_axis = Direction.Create(0,0,1)


#Creates a circle that can be stored as a usable object
def create_circle(point, direction, radius):
	circle = CircularSurface.Create(radius, direction, point)
	return circle.CreatedBody
	
def dividing_curve(curve, num_divisions):
	divided_curve_points = []
	for i in range(num_divisions+1):
		i = float(i)
		eval = curve.Evaluate(1-(i/num_divisions))
		divided_curve_points.append(eval.Point)
	return divided_curve_points

def split_point(branch_pt, xysep, zsep, tier):
	x0 = branch_pt[0]
	y0 = branch_pt[1]
	z0 = branch_pt[2]
	points[tier].append([x0+xysep, y0+xysep, z0+zsep])
	points[tier].append([x0-xysep, y0+xysep, z0+zsep])
	points[tier].append([x0-xysep, y0-xysep, z0+zsep])
	points[tier].append([x0+xysep, y0-xysep, z0+zsep])

def tangent_direction(incoming_curve, evaluation_point):
	interim_1 = incoming_curve.EvalProportion(evaluation_point)
	bottom_point = interim_1.Point
	interim_2 = incoming_curve.EvalProportion(evaluation_point +.001)
	top_point = interim_2.Point
	X1 = bottom_point.X - top_point.X
	Y1 = bottom_point.Y - top_point.Y
	Z1 = bottom_point.Z - top_point.Z
	direction = Direction.Create(X1,Y1,Z1)
	return (bottom_point, direction)

def distance_between_points(point_a, point_b):
	a = (point_a[0]-point_b[0])**2
	b = (point_a[1]-point_b[1])**2
	c = (point_a[2]-point_b[2])**2
	d = (a+b+c)**.5
	return d

class constant_cross_curve:
	def __init__(self):
		self.circle_centers = []
		self.radius = []
	
#Input a series of points and the radius of its parent circle. Return the radius which achieves the desired cross-sectional area and the shape.
#PoA 1: Create a geometry point from each list(point) that is given as an input
#PoA 2: Failsafe, Check to see if the circles are overlapping, if they are not overlapping then print an error.
#PoA 3: Run a while loop that creates overlapping circles, check their combined surface area against a reference and increase or decrease the radius until the desired cross-sectional area is achieved.
#PoA 3.1: Failsafe, If the radius has not converged after a certain number of iterations then exit the function and return nothing
#PoA 3.2: Create circles at all of the input points and select all of them
#PoA 3.3: Combine circles and find area
#PoA 3.4: Use a linear approximation for the overlapping circles radius  v. Area plot to make the first guess. The first guess happens at iteration = 2.
#PoA 3.5: Compare the combined circle area with the reference area. Evaluate how large to make the next radius or exit if correct radius found
#PoA 3.6: Delete all circles
def combine_circles(target_area, tolerance, evaluation_points):
#target_area: What area are the overlapping circles trying to converge to. (Float)
#tolerance: What margin of error is acceptable when comparing to the desired cross-sectional area. (Float)
#evaluation_points: You can add as many points as an input as you would like such as 2,3,4, or 5 points.(List of SpaceClaim point objects)
	circle_centers = []
	circle_centers_check = []
	distances = []
	new = 0 
	prev = 0
	iteration = 0
	lower_limit_radius = (target_area/(math.pi * len(evaluation_points)))**0.5

#PoA 1: Create a geometry point from each list that is given as an input
	for p in evaluation_points: #Take the variable input *points and create the centers of the circles within SpaceClaim
		circle_centers.append(p)
		circle_centers_check.append([p.X, p.Y, p.Z])
	
#PoA 2: Check to see if the circles have split, if they have then print an error
	for i in range(1, len(circle_centers)):
		d = distance_between_points(circle_centers_check[0], circle_centers_check[i])
		distances.append(d)
	if(min(distances)/2 >=  lower_limit_radius):
		print("Error: The circles no longer overlap")
		return

#PoA 3: Run a while loop that creates overlapping circles, check their combined surface area against a reference and increase or decrease the radius until the desired cross-sectional area is achieved.
	linear_approx_1_radius = (min(distances)/2)*1.05
	radius_of_circle = linear_approx_1_radius
	shifting_adder = 0
	while True:	
		#PoA 3.1: Failsafe, If the radius has not converged after a certain number of iterations then exit the function and return nothing
		if (iteration == 29):
			print("Error: The function combine_circles failed to converge to the correct radius")
			print("Target Area: " + str(target_area))
			return
		
		#PoA 3.2: Create circles at all of the input points and select all of them
		sel = []
		trying_circles = []
		print("||||||||||Iteration: " + str(iteration) + "||||||||||||||||")
		print("Radius of Circle: " + str(radius_of_circle))
		print("Shifting Adder: " + str(shifting_adder))
		for i in range(len(circle_centers)): #Take the points of all the circles centers and make circles
			trying_circles.append(create_circle(circle_centers[i], z_axis, radius_of_circle))
			sel.append(Selection.Create(trying_circles[i]))	#Make a list of selections of each circle: sel[0], sel[1], sel[2]...

		#PoA 3.3: Combine circles and find area
		for i in range(1, len(circle_centers)):#Merge all of the created circles into one shape
			Combine.Merge(sel[0], sel[i])	
		combined_area = trying_circles[0].Faces[0].Area #find the area of the combined shape
		print("Area of shape: " + str(combined_area))
		
		#PoA 3.4: Use a linear approximation for the overlapping circles radius  v. Area plot to make the first guess. The first guess happens at iteration = 2.
		if (iteration==0):
			linear_approx_1_area = combined_area
			linear_approx_2_radius = 1.5*linear_approx_1_radius
			radius_of_circle = linear_approx_2_radius
		elif (iteration==1):
			linear_approx_2_area = combined_area
			slope = (linear_approx_2_area - linear_approx_1_area)/(linear_approx_2_radius - linear_approx_1_radius)
			y_intercept = linear_approx_2_area-(slope*linear_approx_2_radius)
			first_radius_attempt = (target_area - y_intercept)/slope
			radius_of_circle = first_radius_attempt
			shifting_adder = radius_of_circle *.25


		else:
			#PoA 3.5: Compare the combined circle area with the reference area. Evaluate how large to make the next radius or exit if correct radius found
			if (combined_area >= target_area*(1-tolerance) and combined_area <= target_area*(1+tolerance)): #If the area of the shape is within the tolerance window return the golden radius
				correct_radius = radius_of_circle
				print("Target area: " + str(target_area))
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
			
		#PoA 3.6: Delete all circles
		for i in range(len(circle_centers)):
			try: 
				Delete.Execute(sel[i]) #Because we do not know which circles were combined and which ones were not we could get errors deleting something that is not there.
			except:
				continue
		iteration = iteration + 1
















#Creates a 2D list that is filled with empty lists
points = [] #points is a 2D array where all planes of points are stored. points[0] contains the origin, points[1] contains the 4 points in the first tier. 
for i in range(tiers):
	points.append([])


points[0].append([0,0,0]) #Adds the origin to tier 0.

#This loop creates planes of points which splits in all 4 quadrents. Each successive plane will have 4 times as many points as the previous. All points are stored in the 2d array points[][]
for i in range(tiers-1):
	for j in range(len(points[i])):
		split_point(points[i][j], xyspread/(2**(i)), distance_between_tiers*math.exp(0), i+1)

#This loop creates circles at each one of the points on the planes and reduces their area by 1/4 each successive tier.		
# for i in range(tiers):
	# for j in range(len(points[i])):
		# happy = points[i][j]
		# point = Point.Create(happy[0],happy[1],happy[2])
		# create_circle(point, z_axis, math.sqrt(initial_radius**2/(tiers**i)))


#This loop creates a curve for every point in the bottom tier. If there are 64 points then 64 curves will be created which wind-up through each tier.
#Each curve will be created by a list of points. This list will contain a point from each tier. 
#We start with the bottom tier and find its parent point one level above it. We repeat this until all points of the tier are defined.
divisible = 0
curves_0 = [] #curves_0 is meant to be used to find where it intersects a plane
curves_1 = [] #curves_1 is meant to be used to divide the curve many times and draw out the curves
for i in range(number_of_curves): #loop through every single curve.
	transfer = points[tiers-1][i] #tiers-1 is the lowest tier so we add our first point.
	temp_point = Point.Create(transfer[0], transfer[1], transfer[2]) # Get the x, y, and z values from the starting point
	point_list = List[Point]()
	point_list.Add(temp_point) 
	divisible = i #If we know the bottom tier point # then we can find its parent
	for j in range(tiers-2, -1, -1):#j is going to loop through all of the tiers finding parents in each one. We start at the tier above the bottom tier and work our way to tier 0 going down one every time.
		for k in range(4):# To find a point's parent we subtract until the number is divisible by 4. Once divisible by 4 we divide and find its parent in the tier above.
			if ((divisible-k)%4 == 0):
				divisible = (divisible-k)/4
				break
		transfer = points[j][divisible]	
		temp_point = Point.Create(transfer[0], transfer[1], transfer[2])
		point_list.Add(temp_point) #Add the parent point to the list of points for which we will make a curve.
	point_list.Add(faux_origin) # add faux origin to curve
	ncurve = NurbsCurve.CreateThroughPoints(False, point_list,1) #make a curve using the list of points. THe third argument is a number which describes how'strict' the curve must be. Setting this value to 0 means you force every curve to hit every point while setting it to 1 allows the curve to be softer.
	curves_0.append(CurveSegment.Create(ncurve))
	curves_1.append(DesignCurve.Create(GetRootPart(), curves_0[i]))

villian = []
for i in range(4):
	villian.append(constant_cross_curve())

for i in range(10):
	intersections = []
	plane_intersect = DatumPlaneCreator.Create(Point.Create(0,0,MM(60 + i*5)), z_axis, True)
	test_plane = plane_intersect.CreatedPlanes[0].Shape.Geometry

	intersections.append(test_plane.IntersectCurve(curves_0[0].Geometry))
	intersections.append(test_plane.IntersectCurve(curves_0[1].Geometry))
	intersections.append(test_plane.IntersectCurve(curves_0[2].Geometry))
	intersections.append(test_plane.IntersectCurve(curves_0[3].Geometry))

	villian[0].circle_centers.append(intersections[0][0].Point)
	villian[1].circle_centers.append(intersections[1][0].Point)
	villian[2].circle_centers.append(intersections[2][0].Point)
	villian[3].circle_centers.append(intersections[3][0].Point)


	gunning_for_it = (math.pi * initial_radius * initial_radius)/16
	dessert = [intersections[0][0].Point, intersections[1][0].Point, intersections[2][0].Point, intersections[3][0].Point]
	result = combine_circles( gunning_for_it,.001, dessert)

	villian[0].radius.append(result)
	villian[1].radius.append(result)
	villian[2].radius.append(result)
	villian[3].radius.append(result)


distinct = []
#for i in range(10):
#	distinct.append(create_circle(villian[0].circle_centers[i], z_axis, villian[0].radius[i]))

# sel = Selection.Create(distinct)
# options = LoftOptions()
# options.GeometryCommandOptions = GeometryCommandOptions()
# Loft.Create(sel, None, options) #Loft all of the circles in a given curve








# start=time.clock()
# #This loop creates circles along regular intervals of the curve and then once created lofts all of the curves.
# for i in range(5):
# 	circles = []
# 	trying = []
# 	start = time.clock()
# 	trying = dividing_curve(curves_1[i], curve_resolution) #Divide each curve into a series of points from which we will create circles with their centerpoints on the curve
# 	for j in range(curve_resolution+1): #This for loop tailors the radius of each circle smaller and smaller depending on its z height.
# 		a = trying[j].Z
# 		b = a/distance_between_tiers + 1
# 		c = float(1)/(4**b)
# 		d = math.sqrt(c)
# 		e = d*initial_radius
# 		circles.append(create_circle(trying[j], z_axis, e)) 
# 	print("creating circles took: " + str(time.clock()-start))
# 	start = time.clock()		
# 	sel = Selection.Create(circles)
# 	options = LoftOptions()
# 	options.GeometryCommandOptions = GeometryCommandOptions()
# 	result = Loft.Create(sel, None, options) #Loft all of the circles in a given curve
# 	print("creating loft took: " + str(time.clock()-start))