sel = Selection.SelectAll()
ComponentHelper.DeleteEmptyComponents(sel)
Delete.Execute(sel)

import math

#-------------------------------------VARIABLES TO BE CHANGED-------------------------------------
inner_chamber_diameter = .05
inner_wall_thickness = .005
outer_wall_thickness = .005
number_of_beans_per_chamber = 5
height_of_single_revolution = .2
chamber_height = .2
shrinkwrap_recess = .002

#-----------------------------------------------------------------------------------------------------------------------------


#-------------------------------------------------------Functions/Classes----------------------------------------------------
#Creates a circle that can be stored as a usable object
def create_circle(point, direction, radius):
	circle = CircularSurface.Create(radius, direction, point)
	return circle.CreatedBody #returns IDOCobject


def dividing_curve(curve, num_divisions):
	divided_curve_points = []
	for i in range(num_divisions+1):
		i = float(i)
		eval = curve.Evaluate(1-(i/num_divisions))
		divided_curve_points.append(eval.Point)
	return divided_curve_points #Returns a list of point objects on the curve.

def split_point(branch_pt, xysep, zsep, tier):
	x0 = branch_pt[0]
	y0 = branch_pt[1]
	z0 = branch_pt[2]
	points[tier].append([x0+xysep, y0+xysep, z0+zsep])
	points[tier].append([x0-xysep, y0+xysep, z0+zsep])
	points[tier].append([x0-xysep, y0-xysep, z0+zsep])
	points[tier].append([x0+xysep, y0-xysep, z0+zsep])


def tangent_direction(incoming_curve, evaluation_point):
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
	curve_1 = DesignCurve.Create(GetRootPart(), curve_0)
	return (curve_0, curve_1)

def generate_curve_from_points(point_list, reconnect_to_first_point_bool):
    point_list_holder = List[Point]()
    for tribot in point_list:
        point_list_holder.Add(tribot)
    ncurve = NurbsCurve.CreateThroughPoints(reconnect_to_first_point_bool, point_list_holder,.0001)
    curve_0 = CurveSegment.Create(ncurve)
    curve_1 = DesignCurve.Create(GetRootPart(), curve_0)
    return (curve_0, curve_1)    

def create_helix(center_point, radius, linear_height, single_rev_height):
    revolution_count = linear_height/single_rev_height
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

def create_circular_pattern(object_to_be_patterned, pattern_count, circular_axis):
	if (isinstance(object_to_be_patterned, Selection)):
		result = Copy.ToClipboard(object_to_be_patterned)
	else :
		result = Copy.ToClipboard(Selection.Create(object_to_be_patterned))
	reference_trajectory = create_circle(origin_3d, circular_axis, .05)
	patterned_objects = []
	patterned_objects.append(object_to_be_patterned)
	for i in range(1,pattern_count):
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
#-----------------------------------------------------------------------------------------------------------------------------

#---------------------------------------------------------Constants----------------------------------------------------------------
origin_3d = Point.Create(0,0,0)
z_axis = Direction.Create(0,0,1)
inner_chamber_radius = inner_chamber_diameter/2
#-----------------------------------------------------------------------------------------------------------------------------



#---------------------------------------------------Functional Code-----------------------------------------------------------------
#Create inner chamber circle
inner_circle = create_circle(origin_3d, z_axis, inner_chamber_radius)

#Create coolant face
stand_in_coolant_radius = .002
coolant_center_circle_point = Point.Create(inner_chamber_radius+inner_wall_thickness, 0,0)
stand_in_coolant_shape = create_circle(coolant_center_circle_point, z_axis, stand_in_coolant_radius)

#Evaluate all edges of coolant face to find min/max X values
x_maximum, x_minimum = find_min_max_x_of_face(stand_in_coolant_shape.Faces[0])
move_step = inner_chamber_radius + inner_wall_thickness - x_minimum

#Move Coolant face to be precisely a wall thickness away from inner chamber
selection = Selection.Create(stand_in_coolant_shape)
localSystem = Frame.Create(origin_3d, Direction.Create(1, 0, 0), Direction.Create(0, 1, 0))
options = MoveOptions()
moveType = TransformType.TranslateX
result = Move.Execute(selection, localSystem, moveType, move_step, options)

#Find Mid Point of coolant Face
plexus = stand_in_coolant_shape.Faces[0].MidPoint()
coolant_face_mid_point = plexus.Point

"""
#Create Helix curve to sweep coolant Face
(helix_0, helix_1) = create_helix(center_point = origin_3d, radius=coolant_face_mid_point.X, linear_height = chamber_height, single_rev_height = height_of_single_revolution)

#Sweep coolant Face along helix curve
shape_face = stand_in_coolant_shape.Faces[0]
selection = Selection.Create(shape_face)
pathselection = Selection.Create(helix_1)
options = SweepFaceCommandOptions()
options.KeepMirror = True
options.SweepNormalTrajectory = False
options.KeepCompositeFaceRelationships = True
options.ExtrudeType = ExtrudeType.Add
options.KeepLayoutSurfaces = False
result = SweepFaces.Execute(selection, pathselection, 9, options)
base_body = GetRootPart().Bodies[1]

#Pattern swept curve around the chamber
patterned_coolant_bodies = create_circular_pattern(base_body, number_of_beans_per_chamber, z_axis)
"""

#Create Outer chamber circle
outer_wall_radius = x_maximum + move_step + outer_wall_thickness
outer_chamber_circle = create_circle(origin_3d,z_axis,outer_wall_radius)
holding_array_1 = points_along_circle(outer_chamber_circle, number_of_beans_per_chamber, False)

#Create Inner shrink wrap inner circle based on Midpoint of coolant shape
shrinkwrap_inner_circle = create_circle(origin_3d, z_axis, outer_wall_radius - shrinkwrap_recess)
holding_array_2 = points_along_circle(shrinkwrap_inner_circle, number_of_beans_per_chamber, True)

#Merge inner and outer point list into a single array
combined_array = []
for i in range(len(holding_array_1)):
	combined_array.append(holding_array_1[i])
	combined_array.append(holding_array_2[i])




generate_curve_from_points(combined_array, True)





"""
#Create new component for wall
wall_comp = ComponentHelper.CreateAtRoot("Wall")
ComponentHelper.SetActive(wall_comp)


#Create Outer Wall
outer_wall_radius = x_maximum + move_step + outer_wall_thickness
outer_chamber_circle = create_circle(origin_3d,z_axis,outer_wall_radius)

#Create pipe using inner chamber circle as the cutter
targets = Selection.Create(outer_chamber_circle)
tools = Selection.Create(inner_circle)
options = MakeSolidsOptions()
options.KeepCutter = False
result = Combine.Intersect(targets, tools, options)

#Extrude chamber wall
selection = Selection.Create(outer_chamber_circle.Faces[0])
options = ExtrudeFaceOptions()
options.ExtrudeType = ExtrudeType.Add
extrude_chamber_wall_result = ExtrudeFaces.Execute(selection, z_axis, chamber_height, options, None)
chamber_wall = extrude_chamber_wall_result.CreatedBodies[0]

#Cut coolant channels out of chamber wall
targets = Selection.Create(chamber_wall)
tools = Selection.Create(patterned_coolant_bodies)
options = MakeSolidsOptions()
options.KeepCutter = True
result = Combine.Intersect(targets, tools, options)

#Delete remenants of the cutting operation
result = Copy.ToClipboard(Selection.Create(chamber_wall))
sel = Selection.SelectAll()
ComponentHelper.DeleteEmptyComponents(sel)
Delete.Execute(sel)
Paste.FromClipboard().PastedObjects[0]

#Create Fluid Domain Component
fluid_comp = ComponentHelper.CreateAtRoot("Fluid")
ComponentHelper.MoveBodiesToComponent(patterned_coolant_bodies, fluid_comp)
"""


