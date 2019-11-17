ComponentHelper.SetRootActive()
root_part = DocumentHelper.GetActivePart()
root_part.ClearAllPartData()

import math
import random

#-------------------------------------VARIABLES TO BE CHANGED-------------------------------------
injector_diameter = .050
distance_between_tiers = .0125
vector_strength = 0.02
tiers = 3
injector_points = 64
number_of_centroid_points = [16,4,1]
points_per_ring = [4, 8, 20, 32]
circles_per_curve = 7
outlet_radius = .0005
#-----------------------------------------------------------------------------------------------------------------------------

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
    
def create_centroids(circle_center_point, circle_radius_i, plane_normal_direction, number_of_centroids):
	half_area_radius = compute_radius_of_circle_with_scalar_of_area(circle_radius_i, 0.5)
    temp_injector_circle = create_circle(circle_center_point, plane_normal_direction, half_area_radius)
	centroid_point_list = points_along_circle(temp_injector_circle, number_of_centroids, offset = True)
	sel = Selection.Create(temp_injector_circle)
	Delete.Execute(sel)
	return centroid_point_list

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

class k_mean_point_class:
    def __init__(self, incoming_point, incoming_weight = 1.0, incoming_direction_vector = [0,0,-1.0], incoming_magnitude = vector_strength):
        self.point = incoming_point
		self.weight = incoming_weight
		self.direction_vector = incoming_direction_vector
		self.magnitude = incoming_magnitude
        self.centroid_point_designator = None
        self.centroid_point_object = None
        self.previous_centroid_point_designator = []
        self.previous_centroid_point_object = []
		self.is_current_designator_same_as_previous = False
	def find_and_set_closest_centroid(self, centroid_list):
		self.previous_centroid_point_designator.append(self.centroid_point_designator)
		self.previous_centroid_point_object.append(self.centroid_point_object)
		min_distance = 9000000000
		for single_centroid in centroid_list:
			distance_between_point_and_centroid = distance_between_points(self.point, single_centroid.current_point)
			if (distance_between_point_and_centroid < min_distance):
				min_distance = distance_between_point_and_centroid
				min_centroid_designator = single_centroid.designator
				min_centroid_point = single_centroid.current_point
		self.centroid_point_designator = min_centroid_designator
		self.centroid_point_object = min_centroid_point

		if (min_centroid_designator == self.previous_centroid_point_designator[-1]):
			self.is_current_designator_same_as_previous = True
		else :
			self.is_current_designator_same_as_previous = False

class centroid_class:
    def __init__(self, incoming_centroid_point, incoming_designator, incoming_direction_vector = [0,0,-1.0], incoming_magnitude = vector_strength):
        self.current_point = incoming_centroid_point
		self.elevated_current_point = Point.Create(incoming_centroid_point.X, incoming_centroid_point.Y, incoming_centroid_point.Z + distance_between_tiers)
		self.direction_vector = incoming_direction_vector
		self.magnitude = incoming_magnitude
        self.designator = incoming_designator
        self.previous_location_points = []
		self.points_in_group = []
		self.combined_groups_weights = 0
	def move_centroid_to_average_point(self, k_mean_point_list):
		self.previous_location_points.append(self.current_point)
		self.points_in_group = []
        
		for single_k_mean_point in k_mean_point_list:
			if(self.designator == single_k_mean_point.centroid_point_designator):
				self.points_in_group.append(single_k_mean_point.point)
		if(len(self.points_in_group) >= 1):
			average_point = create_average_point_from_list_of_points(point_list = self.points_in_group)
			self.current_point = average_point
			self.elevated_current_point = Point.Create(average_point.X, average_point.Y, average_point.Z + distance_between_tiers)
		else :
			print("Error this centroid has no close points: " + str(self.designator))
        
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
			for single_point in list_of_k_mean_point_objects: #Before exiting go through all point objects and add them to whichever centroid they belong to.
				my_centroid_designator = single_point.centroid_point_designator
				my_weight = single_point.weight
				list_of_k_mean_centroid_objects[my_centroid_designator].points_in_group.append(single_point.point)
				list_of_k_mean_centroid_objects[my_centroid_designator].combined_groups_weights = list_of_k_mean_centroid_objects[my_centroid_designator].combined_groups_weights + my_weight
			break

def transfer_converged_centroids_2_points(previous_tier_centroid_objects):
	current_tier_k_mean_point_objects = []
	for old_centroid in previous_tier_centroid_objects:
		temp_point = old_centroid.current_point
		temp_weight = old_centroid.combined_groups_weights
		temp_point = Point.Create(temp_point.X, temp_point.Y, temp_point.Z + distance_between_tiers)
		current_tier_k_mean_point_objects.append(k_mean_point_class(temp_point, temp_weight))
	return current_tier_k_mean_point_objects

def point_list_2_centroid_object(incoming_point_list):
    centroid_object_list = []
    iterator_3 = 0
    for single_centroid in incoming_point_list:
        centroid_object_list.append(centroid_class(single_centroid, iterator_3))
        iterator_3 = iterator_3 + 1	
    return centroid_object_list

def point_list_2_k_mean_point_object(incoming_point_list, incoming_weights):
    k_mean_object_list = []
    for i in range(len(incoming_point_list)):
        k_mean_object_list.append(k_mean_point_class(incoming_point_list[i], incoming_weights[i]))
    return k_mean_object_list

def radius_of_circle_with_given_area(surface_area):
	return (surface_area/math.pi)**0.5

class curve_class:
	def __init__(self, incoming_curve_segment, incoming_design_curve, incoming_k_means_point, incoming_centroid):
		self.curve_segment = incoming_curve_segment #curve_segment is meant to be used to find where the curve intersects a plane
		self.design_curve = incoming_design_curve #design_curve is meant to be used to divide the curve many times and draw out the curves
		self.k_means_point = incoming_k_means_point
		self.centroid = incoming_centroid
		#interim_1 = incoming_design_curve.EvalProportion(0)
		#interim_2 = incoming_design_curve.EvalProportion(1)
		#self.start_point = interim_1.Point
		#self.end_point = interim_2.Point
		self.circle_centers = []
		self.circle_center_normal_direction = []
		self.radii = []
		self.circle_surfaces = []
	def populate_curve_with_circle_info(self, incoming_total_weights, total_area):
		self.circle_centers.append(self.centroid.elevated_current_point) 
		transfer_list = self.centroid.direction_vector
		transfer_direction = Direction.Create(transfer_list[0], transfer_list[1], transfer_list[2])
		self.circle_center_normal_direction.append(transfer_direction)
		target_start_surface_area = (self.centroid.combined_groups_weights/incoming_total_weights) * total_area
		centroid_radius = radius_of_circle_with_given_area(target_start_surface_area)
		self.radii.append(centroid_radius)

		#Only used to create linear taper
		target_end_surface_area = (self.k_means_point.weight/incoming_total_weights) * total_area
		destination_radius = radius_of_circle_with_given_area(target_end_surface_area)
		radius_difference = centroid_radius-destination_radius


		for i in range(1, circles_per_curve):
			(eval_point, normal_direction) = tangent_direction(self.design_curve, i/float(circles_per_curve))
			self.circle_centers.append(eval_point)
			self.circle_center_normal_direction.append(normal_direction)
			linear_taper_radius = centroid_radius-(radius_difference*i/float(circles_per_curve))
			self.radii.append(linear_taper_radius)
		
		self.circle_centers.append(self.k_means_point.point)
		transfer_list = self.k_means_point.direction_vector
		transfer_direction = Direction.Create(transfer_list[0], transfer_list[1], transfer_list[2])
		self.circle_center_normal_direction.append(transfer_direction)
		self.radii.append(destination_radius)
	def populate_curve_with_circles(self):
		for i in range(len(self.radii)):
			self.circle_surfaces.append(create_circle(self.circle_centers[i], self.circle_center_normal_direction[i], self.radii[i]))
	def loft_circle_surfaces(self):
		sel = Selection.Create(self.circle_surfaces)
		options = LoftOptions()
		options.GeometryCommandOptions = GeometryCommandOptions()
		result = Loft.Create(sel, None, options)		
		



def link_point_tiers_with_hermite(incoming_destination_points, incoming_centroid_points):
	curve_class_objects = []
	for stand_in_object in incoming_destination_points:
		destination_point = stand_in_object.point
		destination_vector = stand_in_object.direction_vector
		destination_magnitude =  stand_in_object.magnitude
		point_designator = stand_in_object.centroid_point_designator
		centroid_point = incoming_centroid_points[point_designator].elevated_current_point
		centroid_vector = incoming_centroid_points[point_designator].direction_vector
		centroid_magnitude = incoming_centroid_points[point_designator].magnitude
		(transfer_curve_segment, transfer_design_curve) = hermite_curve(centroid_point, centroid_vector, centroid_magnitude, destination_point, destination_vector, destination_magnitude)
		curve_class_objects.append(curve_class(transfer_curve_segment, transfer_design_curve, stand_in_object, incoming_centroid_points[point_designator]))
	return curve_class_objects

#---------------------------------------------------------Constants-----------------------------------------------------------
origin_3d = Point.Create(0,0,0)
z_axis = Direction.Create(0,0,1)
injector_radius = injector_diameter/2.0
#-----------------------------------------------------------------------------------------------------------------------------

#---------------------------------------------------Functional Code-----------------------------------------------------------


this_point = Point.Create(1,1,0)
direction_1 = [1,0,0]
ratio_b = 1.0
hermite_curve_auto_strength(origin_3d, direction_1, this_point, direction_1, ratio_b)