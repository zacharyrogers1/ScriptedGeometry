ComponentHelper.SetRootActive()
root_part = DocumentHelper.GetActivePart()
root_part.ClearAllPartData()

import math

#-------------------------------------VARIABLES TO BE CHANGED-------------------------------------
inner_chamber_diameter = .05
inner_wall_thickness = .005
outer_wall_thickness = .005
number_of_beans_per_chamber = 10 
height_of_single_revolution = .2

inner_radius_bean = .002
distance_between_inner_outer_bean = .003
arc_length_in_degrees = 120.0

wall_height = .2
beans_to_sim = 5
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
    sketch_frame = Frame.Create(origin_3d, z_axis)
    sketch_plane = Plane.Create(sketch_frame)

    #Create Inner/Outer Bean Circle
    localSystem = Frame.Create(origin_3d, z_axis)
    outsideCircle = Circle.Create(localSystem, inner_rad_bean + inner_outer_offset)
    insideCircle = Circle.Create(localSystem, inner_rad_bean)
    arc_length_in_radians = (arc_length_degrees*math.pi)/180

    #Create Design Curves for inner/outer bean
    all_designcurves = []
    inner_curve_Segment = CurveSegment.Create(insideCircle, Interval.Create(2*math.pi-arc_length_in_radians/2, arc_length_in_radians/2))
    inner_curve = DesignCurve.Create(GetActivePart(), inner_curve_Segment)
    all_designcurves.append(inner_curve)
    outer_curve_Segment = CurveSegment.Create(outsideCircle, Interval.Create(2*math.pi-arc_length_in_radians/2, arc_length_in_radians/2))
    outer_curve = DesignCurve.Create(GetActivePart(), outer_curve_Segment)
    all_designcurves.append(outer_curve)

    #Find midpoint between inner/outerbean ends
    start_middle_point = create_point_between_two_points(inner_curve_Segment.StartPoint, outer_curve_Segment.StartPoint)
    end_middle_point = create_point_between_two_points(inner_curve_Segment.EndPoint, outer_curve_Segment.EndPoint)

    #Create circle and trim away inner half
    connect_circ_start = SketchCircle.Create(start_middle_point, inner_outer_offset/2, sketch_plane).CreatedCurve[0]
    curveSelPoint_start = SelectionPoint.Create(connect_circ_start, 0)
    start_semi_circle = TrimSketchCurve.Execute(curveSelPoint_start).CreatedCurve[0]
    all_designcurves.append(start_semi_circle)
    connect_circ_end = SketchCircle.Create(end_middle_point, inner_outer_offset/2, sketch_plane).CreatedCurve[0]
    curveSelPoint_end = SelectionPoint.Create(connect_circ_end, 0)
    end_semi_circle = TrimSketchCurve.Execute(curveSelPoint_end).CreatedCurve[0]
    all_designcurves.append(end_semi_circle)

    #Select all curves and then fill the space
    sel = Selection.Create2(all_designcurves)
    secondarySelection = Selection()
    options = FillOptions()
    result = Fill.Execute(sel, secondarySelection, options, FillMode.ThreeD)
    body_list = result.GetCreated[IDesignBody]()
    single_body = body_list[0]

    sel = Selection.Create(inner_curve, outer_curve, start_semi_circle, end_semi_circle)
    Delete.Execute(sel)
    return single_body


def create_partial_chamber_surface(surface, inner_radius, outer_radius, total_coolant_channels, channels_to_sim):
    resolving = surface.Faces[0].MidPoint()
    mid_point = resolving.Point
    pi_positioning = math.atan2(mid_point.Y, mid_point.X)

    raised_point = Point.Create(0,0,mid_point.Z)
    localSystem = Frame.Create(raised_point, z_axis)
    local_plane = Plane.Create(localSystem)
    outsideCircle = Circle.Create(localSystem, outer_radius)
    insideCircle = Circle.Create(localSystem, inner_radius)
    radial_move_step = (2*math.pi)/total_coolant_channels
    local_interval = Interval.Create(pi_positioning - radial_move_step/2, pi_positioning + radial_move_step*(channels_to_sim-1) + radial_move_step/2)

    #Create Design Curves for inner/outer circle
    boundary = List[ITrimmedCurve]()
    inner_curve_Segment = CurveSegment.Create(insideCircle, local_interval)
    boundary.Add(inner_curve_Segment)
    outer_curve_Segment = CurveSegment.Create(outsideCircle, local_interval)
    boundary.Add(outer_curve_Segment)


    #Attach Inner/Outer circles
    attach_1 = CurveSegment.Create(inner_curve_Segment.StartPoint, outer_curve_Segment.StartPoint)
    boundary.Add(attach_1)
    attach_2 = CurveSegment.Create(inner_curve_Segment.EndPoint, outer_curve_Segment.EndPoint)
    boundary.Add(attach_2)

    designResult = PlanarBody.Create(local_plane, boundary)
    partial_chamber_surface = designResult.CreatedBody
    return partial_chamber_surface

def project_surface_to_plane(incoming_surface, projection_plane):
    #Project Curves to plane
    sel = Selection.Create(incoming_surface.Edges)
    print(len(incoming_surface.Edges))
    projected_curves = ProjectToSketch.Create(sel, projection_plane).CreatedCurve
    print(type(projected_curves))
    print(len(projected_curves))

    #Fill in projected sketch
    sel = Selection.Create(projected_curves)
    secondarySelection = Selection()
    options = FillOptions()
    result = Fill.Execute(sel, secondarySelection, options, FillMode.ThreeD)
    body_list = result.GetCreated[IDesignBody]()
    projected_body = body_list[0]

    #Delete leftover curves
    sel = Selection.Create(projected_curves)
    Delete.Execute(sel)
    return projected_body

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

#-----------------------------------------------------------------------------------------------------------------------------

#---------------------------------------------------------Constants----------------------------------------------------------------
origin_3d = Point.Create(0,0,0)
z_axis = Direction.Create(0,0,1)
y_axis = Direction.Create(0,1,0)
inner_chamber_radius = inner_chamber_diameter/2
#-----------------------------------------------------------------------------------------------------------------------------



#---------------------------------------------------Functional Code-----------------------------------------------------------------
stand_in_coolant_shape = create_bean(inner_radius_bean, distance_between_inner_outer_bean, arc_length_in_degrees)



#Evaluate all edges of coolant face to find min/max X values
x_maximum, x_minimum = find_min_max_x_of_face(stand_in_coolant_shape.Faces[0])
plexus = stand_in_coolant_shape.Faces[0].MidPoint()
coolant_face_mid_point = plexus.Point
protrude_from_center_min = coolant_face_mid_point.X - x_minimum
protrude_from_center_max = x_maximum - coolant_face_mid_point.X 
coolant_center_x = inner_chamber_radius + inner_wall_thickness + protrude_from_center_min

#Move Coolant face to be precisely a wall thickness away from inner chamber
move_surface_to_coordinate(stand_in_coolant_shape,coolant_center_x,0,0)

#Find Mid Point of coolant Face
plexus = stand_in_coolant_shape.Faces[0].MidPoint()
coolant_face_mid_point = plexus.Point


#Create Helix curve to sweep coolant Face
(helix_0, helix_1) = create_helix(center_point = origin_3d, radius=coolant_face_mid_point.X, linear_height = wall_height, single_rev_height = height_of_single_revolution)

faces_per_curve = 30
faces_to_be_lofted = []
for i in range(1, faces_per_curve+1):
    #Copy/Paste original coolant shape and place it periodically around helix 
    sel = Selection.Create(stand_in_coolant_shape)
    result = Copy.ToClipboard(sel)
    copy_over = Paste.FromClipboard().CreatedObjects[0]
    faces_to_be_lofted.append(copy_over)
    fraction_of_curve = i/float(faces_per_curve)
    (helix_eval_point, tangent_direction_1) = tangent_direction(helix_1, fraction_of_curve)
    move_surface_to_coordinate(copy_over, helix_eval_point.X, helix_eval_point.Y, helix_eval_point.Z)

    #Rotate the coolant face depending on how far it is up the helix so it maintains the same orientation
    sel = Selection.Create(copy_over)
    patterned_coolant_face_midpoint = copy_over.Faces[0].MidPoint()
    patterned_coolant_face_midpoint = patterned_coolant_face_midpoint.Point
    localSystem = Frame.Create(patterned_coolant_face_midpoint, z_axis)
    options = MoveOptions()
    moveType = TransformType.RotateZ
    if(height_of_single_revolution == 0):
        continue
    revolutions_completed_by_helix = wall_height/height_of_single_revolution
    degree_of_rotation = 360 * revolutions_completed_by_helix * fraction_of_curve
    result = Move.Execute(sel, localSystem, moveType, DEG(degree_of_rotation), options)

    #Make the coolant shape normal to the helix path
    orient_surfacenormal_to_direction(copy_over, tangent_direction_1)

#Orient original coolant shape to the normal of helix curve
(helix_start_point, tangent_direction_1) = tangent_direction(helix_1, 0)
orient_surfacenormal_to_direction(stand_in_coolant_shape, tangent_direction_1)
faces_to_be_lofted = [stand_in_coolant_shape] + faces_to_be_lofted



#Create loft of all faces
sel = Selection.Create2(faces_to_be_lofted)
options = LoftOptions()
options.GeometryCommandOptions = GeometryCommandOptions()
result = Loft.Create(sel, None, options)
base_body = result.CreatedBodies[0]

#Create plane normal to z_axis
temp_point = Point.Create(0, 0, .040)
sketch_frame = Frame.Create(temp_point, z_axis)
sketch_plane = Plane.Create(sketch_frame)
instance = DatumPlaneCreator.Create(temp_point, z_axis, True)

# Slice Bodies by Plane
selection = Selection.Create(base_body)
cut_plane = Selection.Create(instance.CreatedPlanes[0])
result = SplitBody.Execute(selection, cut_plane)


#Copy Intersecting Face then delete lofted bodies
sel = Selection.Create(result.CreatedBodies[0].Faces[1])
result = Copy.ToClipboard(sel)
sel = Selection.SelectAll()
Delete.Execute(sel)
stand_in_coolant_shape = Paste.FromClipboard().CreatedObjects[0]

#Switch to sub Component
wall_comp = ComponentHelper.CreateAtRoot("Wall_Component")
ComponentHelper.SetActive(wall_comp)

#Create partial chamber surface and  extrude
outer_wall_radius = coolant_center_x + protrude_from_center_max + outer_wall_thickness
partial_chamber_surface = create_partial_chamber_surface(stand_in_coolant_shape, inner_chamber_radius, outer_wall_radius, number_of_beans_per_chamber, beans_to_sim)
sel = Selection.Create(partial_chamber_surface.Faces[0])
options = ExtrudeFaceOptions()
options.ExtrudeType = ExtrudeType.Add
extrude_chamber_wall_result = ExtrudeFaces.Execute(sel, z_axis, wall_height, options, None)
partial_chamber_wall = extrude_chamber_wall_result.CreatedBodies[0]

#Switch to root part
ComponentHelper.SetRootActive()

#Extrude Coolant Shape
sel = Selection.Create(stand_in_coolant_shape.Faces[0])
options = ExtrudeFaceOptions()
options.ExtrudeType = ExtrudeType.Add
extrude_wall_result = ExtrudeFaces.Execute(sel, z_axis, wall_height, options, None)
single_coolant_body = extrude_wall_result.CreatedBodies[0]

#pattern projected coolant shape
patterned_coolant_bodies = create_circular_pattern(single_coolant_body, number_of_beans_per_chamber, z_axis, stop_count = beans_to_sim)

#Change Active component to wall
ComponentHelper.SetActive(wall_comp)

#Cut coolant channels out of chamber wall
targets = Selection.Create(partial_chamber_wall)
tools = Selection.Create2(patterned_coolant_bodies)
options = MakeSolidsOptions()
options.KeepCutter = True
result = Combine.Intersect(targets, tools, options)

#Delete remenants of the cutting operation
result = Copy.ToClipboard(Selection.Create(partial_chamber_wall))
sel = Selection.SelectAll()
ComponentHelper.DeleteEmptyComponents(sel)
Delete.Execute(sel)
Paste.FromClipboard()

#Create Fluid Domain Component
fluid_comp = ComponentHelper.CreateAtRoot("Fluid_Component")
sel = Selection.Create2(patterned_coolant_bodies)
ComponentHelper.MoveBodiesToComponent(sel, fluid_comp)

#Set the root part as the active component and set it to Share Topolgy
ComponentHelper.SetRootActive()
only_part = DocumentHelper.GetActivePart()
only_part.ShareTopology = only_part.ShareTopology.Share




#----------------------------------------------------------------------Naming Bodies/Faces------------------------------------------------------
all_fluid_bodies = fluid_comp.GetAllBodies()
all_wall_bodies = wall_comp.GetAllBodies()

sel = sel.Create(all_fluid_bodies)
sel.CreateAGroup("fluid_volume")

sel = sel.Create(all_wall_bodies)
sel.CreateAGroup("wall_volume")

inner_area_estimate = 2 * math.pi*inner_chamber_radius*(beans_to_sim/float(number_of_beans_per_chamber))*wall_height
sel = PowerSelection.Faces.ByArea(inner_area_estimate, 
    PowerSelectOptions(True), 
    SearchCriteria.SizeComparison.Equal)
sel.SetActive()
sel = Selection.GetActive()
sel.CreateAGroup("inner_chamber_wall")


iterator_3 = 0
fluid_walls_list = []
for change in all_fluid_bodies:
    sel = sel.Create(change.Faces[4])
    sel.CreateAGroup("fluid_inlet_" + str(iterator_3))

    sel = sel.Create(change.Faces[5])
    sel.CreateAGroup("fluid_outlet_" + str(iterator_3))

    fluid_walls_list.append(change.Faces[0])
    fluid_walls_list.append(change.Faces[1])
    fluid_walls_list.append(change.Faces[2])
    fluid_walls_list.append(change.Faces[3])

    iterator_3 = iterator_3 + 1

sel = sel.Create2(fluid_walls_list)
sel.CreateAGroup("fluid_wall")