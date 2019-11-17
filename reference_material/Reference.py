#-------------------------------------------------------------Functions-----------------------------------------------------------------
#Creates a circle that can be stored as a usable object
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

def average_children(parent_point, current_tier, current_place):
#parent_point: A point described in x,y,z coordinates(List)
#current_tier: What tier of points is the parent point on (int)
#current_place: What place in the current tier is the parent point. Is it the first point in the list or the fifteenth?(int)
	total_children = fractal_splits**(tiers - 1 - current_tier)
	starting_place = current_place*total_children
	x_avg = 0
	y_avg = 0
	z_avg = 0
	for i in range(starting_place, starting_place+total_children):
		x_avg = x_avg + points[tiers-1][i][0]
		y_avg = y_avg + points[tiers-1][i][1]
		z_avg = z_avg + points[tiers-1][i][2]
	x_avg = x_avg/total_children
	y_avg = y_avg/total_children
	z_avg = z_avg/total_children
	vector = [x_avg - parent_point[0], y_avg - parent_point[1], z_avg - parent_point[2]]
	print(vector)
	return vector

def distance_between_points(point_a, point_b):
		a = (point_a[0]-point_b[0])**2
		b = (point_a[1]-point_b[1])**2
		c = (point_a[2]-point_b[2])**2
		d = (a+b+c)**.5
		return d

def plane_intersect_curve(z_height, single_curve):
#z_height: At what z height will you create the plane to see where it intersects. (Float)
#single_curve: This is the curve object to be evaluated. It should be curves_0(SpaceClaim Curve Object)
	transfer_plane = DatumPlaneCreator.Create(Point.Create(0,0, z_height), z_axis, True)
	plane_intersect = transfer_plane.CreatedPlanes[0].Shape.Geometry
	intersections = plane_intersect.IntersectCurve(single_curve.Geometry)
	if (not intersections): #If intersections list is empty this will make the function return false
		return False
	intersection_point = intersections[0].Point
	return intersection_point #(SpaceClaim Point Object)


def hermite_curve(point_a, point_b, tangent_a, tangent_b):
    point_list = List[Point]()
    for i in range (21):
        t = i*.05
        x = (2*t**3 - 3*t**2 + 1)*point_a[0] + (3*t**2 - 2*t**3)*point_b[0] + (-2*t**2 + t**3 + t)*tangent_a[0] + (-t**2 + t**3)*tangent_b[0]
        y = (2*t**3 - 3*t**2 + 1)*point_a[1] + (3*t**2 - 2*t**3)*point_b[1] + (-2*t**2 + t**3 + t)*tangent_a[1] + (-t**2 + t**3)*tangent_b[1]
        z = (2*t**3 - 3*t**2 + 1)*point_a[2] + (3*t**2 - 2*t**3)*point_b[2] + (-2*t**2 + t**3 + t)*tangent_a[2] + (-t**2 + t**3)*tangent_b[2]
        temp_point = Point.Create(x, y, z)
        point_list.Add(temp_point)
    ncurve = NurbsCurve.CreateThroughPoints(False, point_list,.0001)
    curve_0 = CurveSegment.Create(ncurve)
    curve_1 = DesignCurve.Create(GetActivePart(), curve_0)
    return (curve_0, curve_1)

def generate_curve_from_points(point_list, reconnect_to_first_point_bool):
    point_list_holder = List[Point]()
    for tribot in point_list:
        point_list_holder.Add(tribot)
    ncurve = NurbsCurve.CreateThroughPoints(reconnect_to_first_point_bool, point_list_holder,.0001)
    curve_0 = CurveSegment.Create(ncurve)
    curve_1 = DesignCurve.Create(GetActivePart(), curve_0)
    return (curve_0, curve_1)  

def create_helix(center_point, radius, linear_height, single_rev_height):
    if (single_rev_height == 0):
        straight_start_point = Point.Create(center_point.X + radius, 0, 0)
        straight_end_point = Point.Create(center_point.X + radius, 0, linear_height)
        curveSegment = CurveSegment.Create(straight_start_point, straight_end_point)
        straight_curve = DesignCurve.Create(GetActivePart(), curveSegment)
        return (None, straight_curve)
    points_per_revolution = 30.0
    step_size = 2*math.pi/points_per_revolution
    final_t = 2*math.pi*linear_height/single_rev_height
    iterator_0=0
    helix_points = []
    while iterator_0 < final_t + 0.5*step_size:
        x= radius*math.cos(iterator_0) + center_point.X
        y= radius*math.sin(iterator_0) + center_point.Y
        z= (single_rev_height/(2*math.pi))*iterator_0 + center_point.Z
        temp_point = Point.Create(x,y,z)
        helix_points.append(temp_point)
        iterator_0 = iterator_0+step_size
    (curve_0, curve_1) = generate_curve_from_points(helix_points, False)
    return (curve_0, curve_1)

def create_circular_pattern(object_to_be_patterned, pattern_count, circular_axis, stop_count = None):
	if (isinstance(object_to_be_patterned, Selection)):
		result = Copy.ToClipboard(object_to_be_patterned)
	else:
		result = Copy.ToClipboard(Selection.Create(object_to_be_patterned))
	if (stop_count == None):
		repetitions = pattern_count
	else:
		repetitions = stop_count
	reference_trajectory = create_circle(origin_3d, circular_axis, .05)
	patterned_objects = []
	patterned_objects.append(object_to_be_patterned)
	for i in range(1,repetitions):
		copy_over = Paste.FromClipboard().CreatedObjects
		patterned_objects.append(copy_over)
		selection = Selection.Create(patterned_objects[i])
		trajectory = Selection.Create(reference_trajectory.Edges[0])
		localSystem = Frame.Create(origin_3d, z_axis)
		options = MoveOptions()
		result = Move.AlongTrajectory(selection, trajectory, i/float(pattern_count), options, localSystem)
	sel = Selection.Create(reference_trajectory)
	Delete.Execute(sel)
	return patterned_objects

def find_min_max_x_of_face(incoming_face):
	point_evaluation = []
	for single_edge in incoming_face.Edges:
		for i in range(21):
			interim_1 = single_edge.EvalProportion(i/20.0)
			x_value = interim_1.Point.X
			point_evaluation.append(x_value)
	maximum = max(point_evaluation)
	minimum = min(point_evaluation)
	return (maximum, minimum)

def points_along_circle(incoming_circle, num_points_wanted, offset):
    point_evaluation = []
    circle_edge = incoming_circle.Edges[0]
    if offset == True:
        for i in range(num_points_wanted):
            interim_1 = circle_edge.EvalProportion(0.5/num_points_wanted+i/float(num_points_wanted))
            point_evaluation.append(interim_1.Point)
    else:
         for i in range(num_points_wanted):
            interim_1 = circle_edge.EvalProportion(i/float(num_points_wanted))
            point_evaluation.append(interim_1.Point)       
    return point_evaluation

def orient_surfacenormal_to_direction(surface, direction):
	instance = DatumPlaneCreator.Create(origin_3d, direction, True)
	selection_1 = Selection.Create(surface)
	selection_2 = Selection.Create(instance.CreatedPlanes[0])
	mid_point = surface.Faces[0].MidPoint()
	mid_point = mid_point.Point
	localSystem = Frame.Create(mid_point, z_axis)
	options = MoveOptions()
	axis_type = HandleAxis.Z
	result = Move.OrientTo(selection_1,selection_2, localSystem, axis_type, options)
	instance.CreatedPlanes[0].Delete()

def move_surface_to_coordinate(surface,x,y,z):
	selection = Selection.Create(surface)
	mid_point = surface.Faces[0].MidPoint()
	mid_point = mid_point.Point
	localSystem = Frame.Create(mid_point, z_axis)
	options = MoveOptions()
	Move.MoveToCoordinate(selection,x,y,z, localSystem, options)

def create_point_between_two_points(point_a, point_b):
    x = (point_a.X + point_b.X)/2
    y = (point_a.Y + point_b.Y)/2
    z = (point_a.Z + point_b.Z)/2
    creation_point = Point.Create(x,y,z)
    return creation_point

def create_bean(inner_rad_bean, inner_outer_offset, arc_length_degrees):
	#Create Inner/Outer Bean Circle
    localSystem = Frame.Create(origin_3d, z_axis)
    outsideCircle = Circle.Create(localSystem, inner_rad_bean + inner_outer_offset)
    insideCircle = Circle.Create(localSystem, inner_rad_bean)
    arc_length_in_radians = (arc_length_degrees*math.pi)/180

    #Create Design Curves for inner/outer bean
    all_designcurves = []
    inner_curve_Segment = CurveSegment.Create(insideCircle, Interval.Create(2*math.pi-arc_length_in_radians/2, arc_length_in_radians/2))
    inner_curve = DesignCurve.Create(GetRootPart(), inner_curve_Segment)
    all_designcurves.append(inner_curve)
    outer_curve_Segment = CurveSegment.Create(outsideCircle, Interval.Create(2*math.pi-arc_length_in_radians/2, arc_length_in_radians/2))
    outer_curve = DesignCurve.Create(GetRootPart(), outer_curve_Segment)
    all_designcurves.append(outer_curve)

    #Create semi circles connecting inner and outer bean
    start_middle_point = create_point_between_two_points(inner_curve_Segment.StartPoint, outer_curve_Segment.StartPoint)
    end_middle_point = create_point_between_two_points(inner_curve_Segment.EndPoint, outer_curve_Segment.EndPoint)
    start_semi_circle = SketchArc.Create(start_middle_point, inner_curve_Segment.StartPoint, outer_curve_Segment.StartPoint).CreatedCurve[0]
    all_designcurves.append(start_semi_circle)
    end_semi_circle = SketchArc.Create(end_middle_point, outer_curve_Segment.EndPoint, inner_curve_Segment.EndPoint).CreatedCurve[0]
    all_designcurves.append(end_semi_circle)

    #Select all curves and then fill the space
    sel = Selection.Create2(all_designcurves)
    secondarySelection = Selection()
    options = FillOptions()
    result = Fill.Execute(sel, secondarySelection, options, FillMode.ThreeD)
    body_list = result.GetCreated[IDesignBody]()
    single_body = body_list[0]
    return single_body

def angle_between_vectors(vector_1, vector_2):
	numerator = vector_1[0]*vector_2[0] + vector_1[1]*vector_2[1] + vector_1[2]*vector_2[2]
	denominator = ((vector_1[0]**2 + vector_1[1]**2 +vector_1[2]**2) * (vector_2[0]**2 + vector_2[1]**2 +vector_2[2]**2))**0.5
	combine = numerator/denominator
	angle_in_radians = math.acos(combine)
	angle_in_degrees = angle_in_radians*180/math.pi
	return angle_in_degrees

def move_radially_about_axis(object_to_be_moved, axis_of_rotation, clockwise_degree_movement):
	reference_trajectory = create_circle(origin_3d, axis_of_rotation, .05)
	selection = Selection.Create(object_to_be_moved)
	trajectory = Selection.Create(reference_trajectory.Edges[0])
	localSystem = Frame.Create(origin_3d, z_axis)
	options = MoveOptions()
	result = Move.AlongTrajectory(selection, trajectory, clockwise_degree_movement/float(360), options, localSystem)
	sel = Selection.Create(reference_trajectory)
	Delete.Execute(sel)

def create_element(center_point, element_type, clockwise_rotation):
	if (element_type == 1):
		injector_point_o1 = Point.Create(center_point.X + element_separation, center_point.Y, center_point.Z)
		injector_point_o2 = Point.Create(center_point.X - element_separation, center_point.Y, center_point.Z)
		fuel_circle = create_circle(center_point, z_axis, fuel_injector_diameter/2)
		oxidizer_circle_1 = create_circle(injector_point_o1, z_axis, oxidizer_injector_diameter/2)
		oxidizer_circle_2 = create_circle(injector_point_o2, z_axis, oxidizer_injector_diameter/2)
		selection = Selection.Create(fuel_circle, oxidizer_circle_1, oxidizer_circle_2)
	localSystem = Frame.Create(center_point, z_axis)
	options = MoveOptions()
	moveType = TransformType.RotateZ
	result = Move.Execute(selection, localSystem, moveType, DEG(clockwise_rotation), options)
	return selection

class extremes_of_face_midpoint:
    def __init__(self, incoming_body):
        x_coordinates = []
        y_coordinates = []
        z_coordinates = []
        for try_out in incoming_body.Faces:
            x_coordinates.append(try_out.MidPoint().Point.X)
            y_coordinates.append(try_out.MidPoint().Point.Y)
            z_coordinates.append(try_out.MidPoint().Point.Z)
        self.highest_x_face = x_coordinates.index(max(x_coordinates))
        self.lowest_x_face = x_coordinates.index(min(x_coordinates))
        self.highest_y_face = y_coordinates.index(max(y_coordinates))
        self.lowest_y_face = y_coordinates.index(min(y_coordinates))
        self.highest_z_face = z_coordinates.index(max(z_coordinates))
        self.lowest_z_face = z_coordinates.index(min(z_coordinates))

class curve_class:
	def __init__(self, c0, c1):
		self.circle_centers = [] #(SpaceClaim point objects in an array)
		self.radius = []
		self.curves_0 = c0 #curves_0 is meant to be used to find where the curve intersects a plane
		self.curves_1 = c1 #curves_1 is meant to be used to divide the curve many times and draw out the curves
		interim_1 = c1.EvalProportion(0)
		interim_2 = c1.EvalProportion(1)
		self.start_point = interim_1.Point
		self.end_point = interim_2.Point

		iterator_0 = interim_1.Point.Z #Find where the curve intersects the z plane along regular intervals and add them to the circle_centers list. 
		while(True):
			temp_evaluate = plane_intersect_curve(iterator_0, c0)
			if (temp_evaluate == False):#Once the curve no longer intersects a plane break out of the loop.
				break
			self.circle_centers.append(temp_evaluate)
			iterator_0 = iterator_0 + constant_area_interval
		self.circle_centers.append(interim_2.Point)#Add the end_point of the curve to the circle_centers list



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

	if(evaluation_points[0]==evaluation_points[1]): # if the same points are input then just return the radius required for one circle to make the target area
		not_split_radius = (target_area/math.pi)**0.5
		return not_split_radius
	
#PoA 2: Check to see if the circles have split, if they have then print an error
	for i in range(1, len(circle_centers)):
		d = distance_between_points(circle_centers_check[0], circle_centers_check[i])
		distances.append(d)
	if(min(distances)/2 >=  lower_limit_radius):
		print("Error: The circles no longer overlap")
		return False

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
		print("Target Area: " + str(target_area))
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







#---------------------------------------------------------------TIPS----------------------------------------------------------------------------------------
#Creating named Selections. Creating groups
sel = Selection.Create(Face3, Face4)
sel.CreateAGroup("Top Face")

#Bringing in C# function into python Scripting
from SpaceClaim.Api.V16 import HeavyweightFacePattern

#Replacing the IDocObjects like Face1, Face 2. In this example circle[0] can be used in all areas Face1 can be used
circle = []
circle.append(Face1)

#Extracting coordinates from created points
point = Point.Create(1,2,3)
print(point.X)
>1

#Creating circles. First create 
point1 = Point.Create(MM(0), MM(0), MM(0)) #creates a point in3d Space
direction1 = Direction.Create(0,0,1) #Creates a direction vecotr
frame = Frame.Create(point1, direction1) #With both a point and vector you can describe a plane. A frame and plane are the same thing
plane = Plane.Create(frame) 
result = ViewHelper.SetSketchPlane(plane)#Create a sketch plane to draw circles on
origin = Point2D.Create(MM(0), MM(0)) #Creates a 2D point on the sketch plane
result = SketchCircle.Create(origin, MM(10)) #Sketches a circle on the 2D sketch plane
result = ViewHelper.SetViewMode(InteractionMode.Solid, Info1) #This just changes the mode from sketch mode to 3d mode. After creating the circle you MUST go into 3d mode or create a new sketch plane in order for the surface to be selectable.

#Setting a color
sel = Selection.Create(circle_inlet[i])
options = SetColorOptions()
options.FaceColorTarget = FaceColorTarget.Body
ColorHelper.SetColor(sel, options, Color.FromArgb(255, 255, 0, 0))

#Power Selection into extruding faces(Works)
selection = PowerSelection.Faces.ByArea(MM2(28274333.88), 
	PowerSelectOptions(True), 
	SearchCriteria.SizeComparison.Equal)
options = ExtrudeFaceOptions()
options.ExtrudeType = ExtrudeType.Add
result = ExtrudeFaces.Execute(selection, Direction.Create(0, 0, 1), MM(14585.43), options, None)

#Power Selection into named Faces(Works)
sel = PowerSelection.Faces.ByArea(MM2(28274333.88), 
	PowerSelectOptions(True), 
	SearchCriteria.SizeComparison.Equal)
sel.SetActive()
happy = Selection.GetActive()
happy.CreateAGroup("Outlet")

#Creating tiers of points all at once
for i in range(1,tiers):
	for j in range(len(points[i-1])):
		list0 = [points[i-1][j][0] + xyspread/(2**i), points[i-1][j][1] + xyspread/(2**i), distance_between_tiers*i] #take the x and y values of the parent point and add a spread that diminishes by one half every tier. The z distance is the same between all tiers.
		list1 = [points[i-1][j][0] - xyspread/(2**i), points[i-1][j][1] + xyspread/(2**i), distance_between_tiers*i]
		list2 = [points[i-1][j][0] - xyspread/(2**i), points[i-1][j][1] - xyspread/(2**i), distance_between_tiers*i]
		list3 = [points[i-1][j][0] + xyspread/(2**i), points[i-1][j][1] - xyspread/(2**i), distance_between_tiers*i]
		points[i].append(list0)
		points[i].append(list1)
		points[i].append(list2)
		points[i].append(list3)

#Finding a point on a curve
ratchet = curves[0].EvalProportion(0) #Eval Proportion takes an argument from 0-1 and evaluates the curve at that point. Think of 0-1 as the fraction of the curve. 1 is origin, 0 is bottom tier
joyful = ratchet.Point
print(joyful.X)
print(joyful.Y)
print(joyful.Z)

#Creating a plane
temp_point = Point.Create(3, 2, 1)
instance = DatumPlaneCreator.Create(temp_point, z_axis, True)

#Create a line
curveSegment = CurveSegment.Create(Point.Create(1,0,0), Point.Create(0,1,0))
DesignCurve.Create(GetRootPart(), curveSegment)

#Lofting Objects
circles = []
sel = Selection.Create(circles)
options = LoftOptions()
options.GeometryCommandOptions = GeometryCommandOptions()
result = Loft.Create(sel, None, options) #Note; Once a face has been lofted that face no longer exists. If you want to loft multiple things to the same face then you need to create duplicate faces

#Creating the subtraction of two selections
sel = Selection.Create(result.CreatedBodies[0].Faces)
sel.CreateAGroup("part_1")
sel1 = Selection.Create(result.CreatedBodies[0].Faces[0])
sel1.CreateAGroup("subpart")
difference = sel-sel1
difference.CreateAGroup("difference")

#Split a curve using a plane
targetCurves = Selection.Create(Curve1)
toolEntity = Selection.Create(DatumPlane1)
result = SplitSketchCurve.Execute(targetCurves, toolEntity)

#Copy and Paste
result = Copy.ToClipboard(Selection.Create(cirlce_1))
result = Paste.FromClipboard()

#Finding MidPoint of Face
plexus = stand_in_coolant_shape.Faces[0]
plexus_1 = plexus.MidPoint()
plexus_2 = plexus_1.Point
print(plexus_2.X)

# Sweep Face along a curve
shape = create_circle(origin_3d, z_axis, radius)
retry = shape.Faces[0] 
selection = Selection.Create(retry)#Select Face of circle
pathselection = Selection.Create(z_axis) #Select path
options = SweepFaceCommandOptions()
options.KeepMirror = True
options.SweepNormalTrajectory = False
options.KeepCompositeFaceRelationships = True
options.ExtrudeType = ExtrudeType.Add
options.KeepLayoutSurfaces = False
result = SweepFaces.Execute(selection, pathselection, MM(274.758810646276), options)

#Create new component and set as active
fluid_comp = ComponentHelper.CreateAtRoot("fluids")
ComponentHelper.SetActive(fluid_comp)
create_circle(origin_3d,z_axis, .3)

#rotate an object
selection = Selection.Create(circle_1)
localSystem = Frame.Create(midPoint, z_axis) #The loaction of the frame determines how the object will be rotated. Make sure to set the frame point equal to the center of the object.
options = MoveOptions()
moveType = TransformType.RotateZ
result = Move.Execute(selection, localSystem, moveType, DEG(62), options)

#Create a square surface
boundary = List[ITrimmedCurve]()
stand_in_coolant_radius = .002
curveSegment = CurveSegment.Create(Point.Create(0,0,0), Point.Create(stand_in_coolant_radius,0,0))
boundary.Add(curveSegment)
curveSegment = CurveSegment.Create(Point.Create(stand_in_coolant_radius,0,0), Point.Create(stand_in_coolant_radius, stand_in_coolant_radius,0))
boundary.Add(curveSegment)
curveSegment = CurveSegment.Create(Point.Create(stand_in_coolant_radius,stand_in_coolant_radius,0), Point.Create(0,stand_in_coolant_radius,0))
boundary.Add(curveSegment)
curveSegment = CurveSegment.Create(Point.Create(0,stand_in_coolant_radius,0), Point.Create(0,0,0))
boundary.Add(curveSegment)
designResult = PlanarBody.Create(Plane.PlaneXY, boundary)
stand_in_coolant_shape = designResult.CreatedBody

#GetCreated as a way to find all geometries of a particular type after completing a command
result = Fill.Execute(sel, secondarySelection, options, FillMode.ThreeD)
happy = result.GetCreated[IDesignBody]()#happy will be a IDesignBody list full of all created IDesignBodies. If only one was created just grab it at happy[0]

#Change Share Topology of part
ComponentHelper.SetRootActive()
only_part = DocumentHelper.GetActivePart()
only_part.ShareTopology = only_part.ShareTopology.Share

#Returns a list of all curves that collide with a given curve. Can print the length of the list to see if it has any neighbors.
connect_circ_start.GetNeighbors()

#Move bodies to component
fluid_comp = ComponentHelper.CreateAtRoot("Fluid_Component")
sel = Selection.Create2(patterned_coolant_bodies)
ComponentHelper.MoveBodiesToComponent(sel, fluid_comp)

#Loft with center line guiding the loft
selection = Selection.Create(Face5, Face6)
secondarySelection = Selection.Create(Curve4)
options = LoftOptions()
options.GeometryCommandOptions = GeometryCommandOptions()
result = Loft.CreateCenterLine(selection, secondarySelection, options, Info5)