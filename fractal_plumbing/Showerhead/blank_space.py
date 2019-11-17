ComponentHelper.SetRootActive()
root_part = DocumentHelper.GetActivePart()
root_part.ClearAllPartData()

import math
import random

#-------------------------------------VARIABLES TO BE CHANGED-------------------------------------
tiers = 2
tier_0_aspect_ratio = 2.0
final_tier_aspect_ratio = 1.5
aspect_ratio_exponent = 0.6
auto_strength_ratio = 1.0
centroid_initialization = "single_circle"
lofting_exponent = 2.0
stop_loft_taper_fraction = 1.0

circles_per_curve = 15
universal_direction_vector = [0,0,-1]
points_per_ring = [4, 8, 20, 32]
injector_diameter = 0.026797
outlet_radius = .0005
#-----------------------------------------------------------------------------------------------------------------------------

#-------------------------------------------------------------Universal Functions-----------------------------------------------------------------
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
	return vector

def distance_between_points(point_a, point_b):
		a = (point_a.X-point_b.X)**2
		b = (point_a.Y-point_b.Y)**2
		c = (point_a.Z-point_b.Z)**2
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

def hermite_curve(point_a, tangent_direction_a, magnitude_a, point_b, tangent_direction_b, magnitude_b):
    point_list = List[Point]()
    tangent_direction_a = [float(iterator_10) for iterator_10 in tangent_direction_a]
    vector_reduction_denominator_a = (tangent_direction_a[0]**2 + tangent_direction_a[1]**2 + tangent_direction_a[2]**2 )**0.5
    unit_vector_a = [iterator_10/vector_reduction_denominator_a for iterator_10 in tangent_direction_a]
    scaled_vector_a = [iterator_10 * magnitude_a for iterator_10 in unit_vector_a]

    tangent_direction_b = [float(iterator_10) for iterator_10 in tangent_direction_b]
    vector_reduction_denominator_b = (tangent_direction_b[0]**2 + tangent_direction_b[1]**2 + tangent_direction_b[2]**2 )**0.5
    unit_vector_b = [iterator_10/vector_reduction_denominator_b for iterator_10 in tangent_direction_b]
    scaled_vector_b = [iterator_10 * magnitude_b for iterator_10 in unit_vector_b]
    for i in range (21):
        t = i*.05
        x = (2*t**3 - 3*t**2 + 1)*point_a.X + (3*t**2 - 2*t**3)*point_b.X + (-2*t**2 + t**3 + t)*scaled_vector_a[0] + (-t**2 + t**3)*scaled_vector_b[0]
        y = (2*t**3 - 3*t**2 + 1)*point_a.Y + (3*t**2 - 2*t**3)*point_b.Y + (-2*t**2 + t**3 + t)*scaled_vector_a[1] + (-t**2 + t**3)*scaled_vector_b[1]
        z = (2*t**3 - 3*t**2 + 1)*point_a.Z + (3*t**2 - 2*t**3)*point_b.Z + (-2*t**2 + t**3 + t)*scaled_vector_a[2] + (-t**2 + t**3)*scaled_vector_b[2]
        temp_point = Point.Create(x, y, z)
        point_list.Add(temp_point)
    ncurve = NurbsCurve.CreateThroughPoints(False, point_list,.0001)
    curve_segment = CurveSegment.Create(ncurve)
    design_curve = DesignCurve.Create(GetActivePart(), curve_segment)
    return (curve_segment, design_curve)

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

def radius_of_circle_with_given_area(surface_area):
	return (surface_area/math.pi)**0.5

def create_average_point_from_list_of_points(point_list, weights = None):
	sum_of_weights = 0
    x_avg = 0
    y_avg = 0
    z_avg = 0	
	if (weights == None):
		weights = []
		for stand_in_object in point_list:
			weights.append(1.0)
	for stand_in_object in weights:
		sum_of_weights = sum_of_weights + stand_in_object
    for i in range(len(point_list)):
        x_avg = x_avg + point_list[i].X * (weights[i] / sum_of_weights)
        y_avg = y_avg + point_list[i].Y * (weights[i] / sum_of_weights)
        z_avg = z_avg + point_list[i].Z * (weights[i] / sum_of_weights)
    average_point = Point.Create(x_avg, y_avg, z_avg)
    return average_point

def compute_radius_of_circle_with_scalar_of_area(given_radius, area_scalar):
	adjusted_radius = ((given_radius**2) * area_scalar)**0.5
	return adjusted_radius
#----------------------------------------------------------------------------------------------------------------------------------------------------------














#------------------------------------------------------------------Showerhead Plumbing Specific Functions---------------------------------------------------------------
def hermite_curve_auto_strength(point_a, tangent_direction_a, point_b, tangent_direction_b, distance_to_strength_ratio):
    point_list = List[Point]()
    tangent_direction_a = [float(iterator_10) for iterator_10 in tangent_direction_a]
    vector_reduction_denominator_a = (tangent_direction_a[0]**2 + tangent_direction_a[1]**2 + tangent_direction_a[2]**2 )**0.5
    unit_vector_a = [iterator_10/vector_reduction_denominator_a for iterator_10 in tangent_direction_a]
	magnitude = distance_between_points(point_a, point_b)/distance_to_strength_ratio
    scaled_vector_a = [iterator_10 * magnitude for iterator_10 in unit_vector_a]

    tangent_direction_b = [float(iterator_10) for iterator_10 in tangent_direction_b]
    vector_reduction_denominator_b = (tangent_direction_b[0]**2 + tangent_direction_b[1]**2 + tangent_direction_b[2]**2 )**0.5
    unit_vector_b = [iterator_10/vector_reduction_denominator_b for iterator_10 in tangent_direction_b]
    scaled_vector_b = [iterator_10 * magnitude for iterator_10 in unit_vector_b]
    for i in range (21):
        t = i*.05
        x = (2*t**3 - 3*t**2 + 1)*point_a.X + (3*t**2 - 2*t**3)*point_b.X + (-2*t**2 + t**3 + t)*scaled_vector_a[0] + (-t**2 + t**3)*scaled_vector_b[0]
        y = (2*t**3 - 3*t**2 + 1)*point_a.Y + (3*t**2 - 2*t**3)*point_b.Y + (-2*t**2 + t**3 + t)*scaled_vector_a[1] + (-t**2 + t**3)*scaled_vector_b[1]
        z = (2*t**3 - 3*t**2 + 1)*point_a.Z + (3*t**2 - 2*t**3)*point_b.Z + (-2*t**2 + t**3 + t)*scaled_vector_a[2] + (-t**2 + t**3)*scaled_vector_b[2]
        temp_point = Point.Create(x, y, z)
        point_list.Add(temp_point)
    ncurve = NurbsCurve.CreateThroughPoints(False, point_list,.0001)
    curve_segment = CurveSegment.Create(ncurve)
    design_curve = DesignCurve.Create(GetActivePart(), curve_segment)
    return (curve_segment, design_curve)

def distance_between_pointsxy(point_a, point_b):
		a = (point_a.X-point_b.X)**2
		b = (point_a.Y-point_b.Y)**2
		d = (a+b)**.5
		return d

def create_random_injector_points_with_weights(injector_center_point, injector_radius_i, injector_normal_direction, number_points_to_create):
    temp_injector_circle = create_circle(injector_center_point, injector_normal_direction, injector_radius_i)
    injector_point_list = []
    weights_list = []
    number_succesfully_created_points = 0
    while (number_succesfully_created_points<number_points_to_create):
        float_1 = random.random()
        float_2 = random.random()
        dist_float_1 = (0.5-float_1)**2
        dist_float_2 = (0.5-float_2)**2
        total_distance = (dist_float_1 + dist_float_2)**0.5
        if (total_distance<= 0.5):
            my_eval = temp_injector_circle.Faces[0].EvalProportion(float_1,float_2)
            weights_list.append(1.0)
            injector_point_list.append(my_eval.Point)
            number_succesfully_created_points = number_succesfully_created_points + 1
    #sel = Selection.Create(temp_injector_circle)
    #Delete.Execute(sel)
    return (injector_point_list, weights_list)
    
def create_centroids(circle_center_point, circle_radius_i, plane_normal_direction, number_of_centroids, initialization_type = "single_circle"):
	if number_of_centroids == 1:
		initialization_type = "single_circle"
	if (initialization_type == "single_circle"):
		half_area_radius = compute_radius_of_circle_with_scalar_of_area(circle_radius_i, 0.5)
		temp_injector_circle = create_circle(circle_center_point, plane_normal_direction, half_area_radius)
		centroid_point_list = points_along_circle(temp_injector_circle, number_of_centroids, offset = False)
		sel = Selection.Create(temp_injector_circle)
		Delete.Execute(sel)
		return centroid_point_list
	elif (initialization_type == "double_circle"):
		centroid_number_inner = int(number_of_centroids/2.0)
		centroid_number_outer = int(round(number_of_centroids/2.0))
		quarter_area_radius = compute_radius_of_circle_with_scalar_of_area(circle_radius_i, 0.25)
		temp_injector_circle_1 = create_circle(circle_center_point, plane_normal_direction, quarter_area_radius)
		centroid_point_list = points_along_circle(temp_injector_circle_1, centroid_number_inner, offset = False) 

		half_area_radius = compute_radius_of_circle_with_scalar_of_area(circle_radius_i, 0.5)
		temp_injector_circle_2 = create_circle(circle_center_point, plane_normal_direction, half_area_radius)
		centroid_point_list = centroid_point_list + points_along_circle(temp_injector_circle_2, centroid_number_outer, offset = False)

		sel = Selection.Create(temp_injector_circle_1, temp_injector_circle_2)
		Delete.Execute(sel)
		return centroid_point_list		

def create_average_point_from_list_of_pointsxy(point_list, weights = None):
	sum_of_weights = 0
    x_avg = 0
    y_avg = 0
	if (weights == None):
		weights = []
		for stand_in_object in point_list:
			weights.append(1.0)
	for stand_in_object in weights:
		sum_of_weights = sum_of_weights + stand_in_object
    for i in range(len(point_list)):
        x_avg = x_avg + point_list[i].X * (weights[i] / sum_of_weights)
        y_avg = y_avg + point_list[i].Y * (weights[i] / sum_of_weights)
    average_point = Point.Create(x_avg, y_avg, 0)
    return average_point

class k_mean_point_class:
    def __init__(self, incoming_point, incoming_weight = 1.0, incoming_direction_vector = universal_direction_vector, incoming_radius = outlet_radius):
        self.point = incoming_point
		self.weight = incoming_weight
		self.direction_vector = incoming_direction_vector
		self.radius = incoming_radius
        self.centroid_point_designator = None
        self.centroid_point_object = None
        self.previous_centroid_point_designator = []
        self.previous_centroid_point_object = []
		self.xydistance_to_closest_centroid = 0
		self.is_current_designator_same_as_previous = False
	def find_and_set_closest_centroid(self, centroid_list):
		self.previous_centroid_point_designator.append(self.centroid_point_designator)
		self.previous_centroid_point_object.append(self.centroid_point_object)
		min_distance = 9000000000
		for single_centroid in centroid_list:
			distance_between_point_and_centroid = distance_between_pointsxy(self.point, single_centroid.current_point)
			if (distance_between_point_and_centroid < min_distance):
				min_distance = distance_between_point_and_centroid
				min_centroid_designator = single_centroid.designator
				min_centroid_point = single_centroid.current_point
		self.centroid_point_designator = min_centroid_designator
		self.centroid_point_object = min_centroid_point
		self.xydistance_to_closest_centroid = min_distance
		if (min_centroid_designator == self.previous_centroid_point_designator[-1]):
			self.is_current_designator_same_as_previous = True
		else :
			self.is_current_designator_same_as_previous = False

class centroid_class:
    def __init__(self, incoming_centroid_point, incoming_designator, incoming_direction_vector = universal_direction_vector):
        self.current_point = incoming_centroid_point
		self.direction_vector = incoming_direction_vector
        self.designator = incoming_designator
        self.previous_location_points = []
		self.points_in_group = []
		self.k_means_points_in_group = []
		self.weights_of_points_in_group = []
		self.radii_of_points_in_group = []
		self.combined_groups_weights = 0
		self.combined_groups_area = 0
	def move_centroid_to_average_point(self, k_mean_point_list):
		self.previous_location_points.append(self.current_point)
		self.points_in_group = []
		self.k_means_points_in_group = []
		self.weights_of_points_in_group = []
		self.radii_of_points_in_group = []
		for single_k_mean_point in k_mean_point_list:
			if(self.designator == single_k_mean_point.centroid_point_designator):
				self.points_in_group.append(single_k_mean_point.point)
				self.k_means_points_in_group.append(single_k_mean_point)
				self.weights_of_points_in_group.append(single_k_mean_point.weight)
				self.radii_of_points_in_group.append(single_k_mean_point.radius)
		if(len(self.points_in_group) >= 1):
			average_point = create_average_point_from_list_of_pointsxy(point_list = self.points_in_group, weights = self.weights_of_points_in_group)
			self.current_point = average_point
			self.combined_groups_weights = sum(self.weights_of_points_in_group)
			self.combined_groups_area = sum(radius_list_2_area_list(self.radii_of_points_in_group))
		else :
			print("Error this centroid has no close points: " + str(self.designator))
	def raise_centroid_to_correct_aspect_ratio(self, min_aspect_ratio):
		required_z_heights = []
		for iterator_10 in self.k_means_points_in_group:
			z_height = (min_aspect_ratio * iterator_10.xydistance_to_closest_centroid) + iterator_10.point.Z
			required_z_heights.append(z_height)
		maximum_z = max(required_z_heights)
		self.current_point = Point.Create(self.current_point.X, self.current_point.Y, maximum_z)
	def raise_tier_0_to_have_equal_aspect(self, min_aspect_ratio):
		for iterator_10 in self.k_means_points_in_group:
			new_z_height = self.current_point.Z - (min_aspect_ratio*iterator_10.xydistance_to_closest_centroid)
			iterator_10.point = Point.Create(iterator_10.point.X, iterator_10.point.Y, new_z_height)


		
        
def converge_k_mean(list_of_k_mean_point_objects, list_of_k_mean_centroid_objects):
	iterator_3 = 0
	while (True):
		continue_convergance = False
		iterator_3 = iterator_3 + 1
		for single_point in list_of_k_mean_point_objects: #Find closest centroid to each point object
			single_point.find_and_set_closest_centroid(list_of_k_mean_centroid_objects)

		for single_centroid in list_of_k_mean_centroid_objects: #Go through each centroid and move its point to the average of its collection
			single_centroid.move_centroid_to_average_point(list_of_k_mean_point_objects)
		
		for single_point_check in list_of_k_mean_point_objects:#Go through all points to see if the current centroid point is the same as the previous
			if (single_point_check.is_current_designator_same_as_previous == False):
				continue_convergance = True
				break
		print("Convergence iteration complete: " + str(iterator_3))

		if (continue_convergance == False):# Exit condition. If all points have the same grouping as the previous convergence then exit while loop
			break

def transfer_converged_centroids_2_points(previous_tier_centroid_objects):
	current_tier_k_mean_point_objects = []
	for old_centroid in previous_tier_centroid_objects:
		temp_point = old_centroid.current_point
		temp_weight = old_centroid.combined_groups_weights
		temp_radius = radius_of_circle_with_given_area(old_centroid.combined_groups_area)
		current_tier_k_mean_point_objects.append(k_mean_point_class(temp_point, temp_weight, incoming_radius = temp_radius))
	return current_tier_k_mean_point_objects

def point_list_2_centroid_object(incoming_point_list):
    centroid_object_list = []
    iterator_3 = 0
    for single_centroid in incoming_point_list:
        centroid_object_list.append(centroid_class(single_centroid, iterator_3))
        iterator_3 = iterator_3 + 1	
    return centroid_object_list

def point_list_2_k_mean_point_object(incoming_point_list, incoming_weights, incoming_direction_vectors = None, incoming_radius = None):
    k_mean_object_list = []
    for i in range(len(incoming_point_list)):
        k_mean_object_list.append(k_mean_point_class(incoming_point_list[i], incoming_weights[i], incoming_radius = incoming_radius[i]))
    return k_mean_object_list

#This is for modification of the tier 1 points. Once the k-means points have all been moved up the previous centroids are moved as well.
def centroid_backpedal(incoming_higher_tier_point_list, incoming_lower_tier_centroids):
	for i in range(len(incoming_lower_tier_centroids)):
		incoming_lower_tier_centroids[i].current_point = incoming_higher_tier_point_list[i].point

def delete_parents_without_children(incoming_centroid_objects):
	iterator_11 = 0
	while(True):
		if iterator_11 == len(incoming_centroid_objects):
			break
		if (not incoming_centroid_objects[iterator_11].points_in_group):
			incoming_centroid_objects.pop(iterator_11)
			print("Deleted Empty!")
			continue
		iterator_11 = iterator_11 + 1



class curve_class:
	def __init__(self, incoming_curve_segment, incoming_design_curve, incoming_k_means_point, incoming_centroid):
		self.curve_segment = incoming_curve_segment #curve_segment is meant to be used to find where the curve intersects a plane
		self.design_curve = incoming_design_curve #design_curve is meant to be used to divide the curve many times and draw out the curves
		self.k_means_point = incoming_k_means_point
		self.centroid = incoming_centroid
		self.circle_centers = []
		self.circle_center_normal_direction = []
		self.radii = []
		self.circle_surfaces = []
	def populate_curve_with_circle_info(self, incoming_lofting_exponent):
		#populate list with first parent circle
		self.circle_centers.append(self.centroid.current_point) 
		transfer_list = self.centroid.direction_vector
		transfer_direction = Direction.Create(transfer_list[0], transfer_list[1], transfer_list[2])
		self.circle_center_normal_direction.append(transfer_direction)
		target_start_surface_area = self.centroid.combined_groups_area
		centroid_radius = radius_of_circle_with_given_area(target_start_surface_area)
		self.radii.append(centroid_radius)

		#Find radius difference
		destination_radius = self.k_means_point.radius
		radius_difference = centroid_radius-destination_radius

		#Iterate through entire path length and fill in circles while tapering.
		circles_during_split = int(round(circles_per_curve * stop_loft_taper_fraction))
		zero_hit_after_steps = radius_difference/(circles_during_split**lofting_exponent)
		for iterator_11 in range(1, circles_per_curve):
			(eval_point, normal_direction) = tangent_direction(self.design_curve, iterator_11/float(circles_per_curve))
			self.circle_centers.append(eval_point)
			self.circle_center_normal_direction.append(normal_direction)
			if iterator_11 <= circles_during_split: #If the path is still tapering the radius
				transfer_radius = -zero_hit_after_steps*(iterator_11**incoming_lofting_exponent) + radius_difference + destination_radius
			else: #If the radius is no longer tapering and being held constant
				transfer_radius = destination_radius
			self.radii.append(transfer_radius)

		#The final circle on the child point is added to the lists
		self.circle_centers.append(self.k_means_point.point)
		transfer_list = self.k_means_point.direction_vector
		transfer_direction = Direction.Create(transfer_list[0], transfer_list[1], transfer_list[2])
		self.circle_center_normal_direction.append(transfer_direction)
		self.radii.append(destination_radius)
	def populate_curve_with_circles(self):
		for i in range(len(self.radii)):
			self.circle_surfaces.append(create_circle(self.circle_centers[i], self.circle_center_normal_direction[i], self.radii[i]))
	def populate_curve_with_circles_with_added_thickness(self, incoming_thickness):
		for i in range(len(self.radii)):
			self.circle_surfaces.append(create_circle(self.circle_centers[i], self.circle_center_normal_direction[i], self.radii[i] + incoming_thickness))		
	def loft_circle_surfaces(self):
		sel = Selection.Create(self.circle_surfaces)
		options = LoftOptions()
		options.GeometryCommandOptions = GeometryCommandOptions()
		options.ExtrudeType = ExtrudeType.ForceIndependent
		result = Loft.Create(sel, None, options)
		self.circle_surfaces = []		

def link_point_tiers_with_hermite(incoming_destination_points, incoming_centroid_points):
	curve_class_objects = [] 
	iterator_13 = 0
	for stand_in_centroid in incoming_centroid_points:
		for stand_in_destination in stand_in_centroid.k_means_points_in_group:
			destination_point = stand_in_destination.point
			destination_vector = stand_in_destination.direction_vector
			centroid_point = stand_in_centroid.current_point
			centroid_vector = stand_in_centroid.direction_vector
			if (centroid_point==destination_point): #A centroid point will be equal to one of its children when it only has one child. This will skip that hermite curve.
				continue
			(transfer_curve_segment, transfer_design_curve) = hermite_curve_auto_strength(centroid_point, centroid_vector, destination_point, destination_vector, auto_strength_ratio)
			curve_class_objects.append(curve_class(transfer_curve_segment, transfer_design_curve, stand_in_destination, stand_in_centroid))
	return curve_class_objects

def python_2d_list_2_spaceclaim_points_list(incoming_2d_array):
	spaceclaim_points = []
	for stand_in_object in incoming_2d_array:
		single_point = Point.Create(stand_in_object[0], stand_in_object[1], stand_in_object[2])
		spaceclaim_points.append(single_point)
	return spaceclaim_points

def radius_list_2_area_list(incoming_radius_list):
	area_list = []
	for single_radius in incoming_radius_list:
		area_list.append(single_radius*single_radius*math.pi)
	return area_list
#----------------------------------------------------------------------------------------------------------------------------------------------------------


#---------------------------------------------------------Constants-----------------------------------------------------------
origin_3d = Point.Create(0,0,0)
z_axis = Direction.Create(0,0,1)
injector_radius = injector_diameter/2.0
injector_points = sum(points_per_ring)
aspect_difference = tier_0_aspect_ratio-final_tier_aspect_ratio
centroid_difference = injector_points-1
#-----------------------------------------------------------------------------------------------------------------------------

#---------------------------------------------------Initialization of injector points-----------------------------------------------------------
point_2d_array = [[0, 0, 0], [0.0044661666666666669, 0.0, 0.0], [0.0022330833333333339, 0.0038678137908686005, 0.0], [-0.0022330833333333326, 0.003867813790868601, 0.0], [-0.0044661666666666669, 5.4692960772297073e-19, 0.0], [-0.0022330833333333356, -0.0038678137908685992, 0.0], [0.0022330833333333339, -0.0038678137908686005, 0.0], [0.0089323333333333338, 0.0, 0.0], [0.0079091883731430227, 0.0041510622804189559, 0.0], [0.005074143672718261, 0.0073511662181174374, 0.0], [0.0010766738069339637, 0.0088672065664018515, 0.0], [-0.003167449052692943, 0.0083518767517470217, 0.0], [-0.0066859475062470003, 0.0059232326242928636, 0.0], [-0.0086727759605219716, 0.0021376472853045611, 0.0], [-0.0086727759605219733, -0.002137647285304559, 0.0], [-0.0066859475062469977, -0.0059232326242928644, 0.0], [-0.0031674490526929452, -0.0083518767517470199, 0.0], [0.0010766738069339654, -0.0088672065664018515, 0.0], [0.0050741436727182593, -0.0073511662181174392, 0.0], [0.0079091883731430227, -0.0041510622804189542, 0.0]]
starting_injector_points = python_2d_list_2_spaceclaim_points_list(point_2d_array)
injector_points = len(starting_injector_points)
starting_weights = [1.0] * injector_points
starting_radii = [outlet_radius] * injector_points
combined_outlet_surface_area = math.pi*outlet_radius*outlet_radius*injector_points
#-----------------------------------------------------------------------------------------------------------------------------


#---------------------------------------------------Geometry Creation-----------------------------------------------------------
master_k_mean_points = []
master_k_mean_centroids = []
master_curves = []
centroids_per_tier = []
aspect_ratio_per_tier = []

#Create empty 2D arrays for points, centroids, and curves. Also create list of centroids and aspect ratios depending on tier.
centroids_constant = math.exp(math.log(injector_points)/tiers)
aspect_ratio_constant = aspect_difference/((tiers-1)**aspect_ratio_exponent)
for i in range(tiers):
	master_k_mean_points.append([])
	master_k_mean_centroids.append([])
	master_curves.append([])
	transfer_num = injector_points/(centroids_constant**(i+1))
	centroids_per_tier.append(int(round(transfer_num)))
	transfer_num = -aspect_ratio_constant * (i**aspect_ratio_exponent) + tier_0_aspect_ratio
	aspect_ratio_per_tier.append(transfer_num)

for iterator_4 in range(tiers):
	#Create next tier of centroids
	reduced_radius = compute_radius_of_circle_with_scalar_of_area(injector_radius, 1.0/(2.0**iterator_4))
	temp_centroid_points = create_centroids(origin_3d, reduced_radius, z_axis, centroids_per_tier[iterator_4], centroid_initialization)

	#convert all centroid points to k_mean_centroid objects
	master_k_mean_centroids[iterator_4] = point_list_2_centroid_object(temp_centroid_points)

    if (iterator_4 == 0):
        master_k_mean_points[iterator_4] = point_list_2_k_mean_point_object(starting_injector_points,  starting_weights, incoming_radius = starting_radii)
    else:
	    #Convert converged centroids from lower tier to k_mean_point objects
	    master_k_mean_points[iterator_4] = transfer_converged_centroids_2_points(master_k_mean_centroids[iterator_4-1])

	#K-mean Converge
	converge_k_mean(master_k_mean_points[iterator_4], master_k_mean_centroids[iterator_4])

	delete_parents_without_children(master_k_mean_centroids[iterator_4])

	#Raise centroids to z location which satisfies min aspect ratio
	for iterator_5 in master_k_mean_centroids[iterator_4]:
		iterator_5.raise_centroid_to_correct_aspect_ratio(aspect_ratio_per_tier[iterator_4])
	
	#Raise 1st tier centroids to highest limiting point
	if (iterator_4 == 1):
		for iterator_5 in master_k_mean_centroids[1]:
			iterator_5.raise_tier_0_to_have_equal_aspect(aspect_ratio_per_tier[iterator_4])
		centroid_backpedal(master_k_mean_points[1], master_k_mean_centroids[0])
		
#Iterate through all tiers and connect points with their centroid point.
iterator_5 = 0
for iterator_4 in range(tiers):
	print(iterator_5)
	master_curves[iterator_4] = link_point_tiers_with_hermite(master_k_mean_points[iterator_4], master_k_mean_centroids[iterator_4])
	iterator_5 = iterator_5 + 1

#Populating the curves with circles and lofting them
for iterator_5 in range(len(master_curves)):
	for iterator_4 in range(len(master_curves[iterator_5])):
		master_curves[iterator_5][iterator_4].populate_curve_with_circle_info(lofting_exponent)
		master_curves[iterator_5][iterator_4].populate_curve_with_circles()
		master_curves[iterator_5][iterator_4].loft_circle_surfaces()



# fluid_comp = ComponentHelper.CreateAtRoot("fluids")
# ComponentHelper.SetActive(fluid_comp)
# wall_thickness = 0.0005

# for iterator_5 in range(len(master_curves)):
# 	for iterator_4 in range(len(master_curves[iterator_5])):
#     	master_curves[iterator_5][iterator_4].populate_curve_with_circles_with_added_thickness(wall_thickness)
#     	master_curves[iterator_5][iterator_4].loft_circle_surfaces()
#-----------------------------------------------------------------------------------------------------------------------------




#---------------------------------------------------Naming Faces-----------------------------------------------------------
#Outlet
outlet_area = outlet_radius*outlet_radius*math.pi
sel = PowerSelection.Faces.ByArea(outlet_area, 
    PowerSelectOptions(True), 
    SearchCriteria.SizeComparison.Equal)
sel.SetActive()
power_to_selection = Selection.GetActive()
all_outlets = power_to_selection.GetItems[IDesignFace]()
for iterator_3 in range(len(all_outlets)):
	sel = Selection.Create(all_outlets[iterator_3])
	sel.CreateAGroup("Outlet " + str(iterator_3))