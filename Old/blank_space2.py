sel = Selection.SelectAll()
ComponentHelper.DeleteEmptyComponents(sel)
Delete.Execute(sel)

import math
import time

#-------------------------------------VARIABLES TO BE CHANGED(Only change these numbers)-------------------------------------
throat_radius = .0099949 #.787in Diamater. Constant
throat_z_mag = .02 #.002-.03
exit_throat_ratio = 5 #Area ratio from 2-5 exit_radius/throat_radius

exit_z_mag = .06 # .05 -.135
exit_x_mag = 0.01 #0-0.03
nozzle_length = .04 #.04-.2

throat_area = math.pi*(throat_radius**2)
exit_area = throat_area*exit_throat_ratio
exit_radius = (exit_area/math.pi)**0.5



#---------------------------------FUNCTIONAL CODE BEYOND THIS POINT-------------------------------------------
faux_origin = Point.Create(MM(0), MM(0), MM(-60))
origin_3d = Point.Create(MM(0),MM(0),MM(0))
z_axis = Direction.Create(0,0,1)
#cross_sectional_area = (math.pi * initial_radius * initial_radius)


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
    degree_shift = 0
    for i in range(fractal_splits):
        points[tier].append([x0+xysep*math.cos(math.radians(degree_shift)), y0+xysep*math.sin(math.radians(degree_shift)), z0+zsep])
        degree_shift = degree_shift + 360/fractal_splits


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





def plane_intersect_curve(z_height, single_curve):
    transfer_plane = DatumPlaneCreator.Create(Point.Create(0,0, z_height), z_axis, True)
    plane_intersect = transfer_plane.CreatedPlanes[0].Shape.Geometry
    intersections = plane_intersect.IntersectCurve(single_curve.Geometry)
    delete_plane = transfer_plane.CreatedPlanes[0]
    delete_plane.Delete()
    if (not intersections): #If intersections list is empty this will make the function return false
        return False
    intersection_point = intersections[0].Point
    return intersection_point





    
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
            sel.append(Selection.Create(trying_circles[i])) #Make a list of selections of each circle: sel[0], sel[1], sel[2]...

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
            shifting_adder = radius_of_circle *.3


        else:
            #PoA 3.5: Compare the combined circle area with the reference area. Evaluate how large to make the next radius or exit if correct radius found
            if (combined_area >= target_area*(1-tolerance) and combined_area <= target_area*(1+tolerance)): #If the area of the shape is within the tolerance window return the golden radius
                correct_radius = radius_of_circle
                for i in range(len(circle_centers)):
                    try: 
                        Delete.Execute(sel[i]) #Because we do not know which circles were combined and which ones were not we could get errors deleting something that is not there.
                    except:
                        continue
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
                if(shifting_adder>=radius_of_circle):
                    shifting_adder = radius_of_circle/2
                radius_of_circle = radius_of_circle - shifting_adder
            prev = new #This stores the last overshoot or undershoot estimate.
            
        #PoA 3.6: Delete all circles
        for i in range(len(circle_centers)):
            try: 
                Delete.Execute(sel[i]) #Because we do not know which circles were combined and which ones were not we could get errors deleting something that is not there.
            except:
                continue
        iteration = iteration + 1


#-----------------------------Create Contour-----------------------------------------------------
point_1 = [throat_radius,0, 0]
point_2 = [exit_radius, 0, nozzle_length]
point_1_t = [0,0, throat_z_mag]
point_2_t = [exit_x_mag , 0 , exit_z_mag]
(transfer_0, transfer_1) = hermite_curve(point_1, point_2, point_1_t, point_2_t)



point_1_m = [-throat_radius,0, 0]
point_2_m = [-exit_radius, 0, nozzle_length]
point_1_t_m = [0,0, throat_z_mag]
point_2_t_m = [-exit_x_mag , 0 , exit_z_mag]
(transfer_2, transfer_3) = hermite_curve(point_1_m, point_2_m, point_1_t_m, point_2_t_m)




#-----------------------------Connect Curves-----------------------------------------------------

#1. Create Exhaust box
curveSegment = CurveSegment.Create(Point.Create(exit_radius,0,nozzle_length), Point.Create(exit_radius+.05,0,nozzle_length))
connector_1 = DesignCurve.Create(GetRootPart(), curveSegment)

curveSegment = CurveSegment.Create(Point.Create(exit_radius+.05,0,nozzle_length), Point.Create(exit_radius+.05,0,nozzle_length+0.2))
connector_2 = DesignCurve.Create(GetRootPart(), curveSegment)

curveSegment = CurveSegment.Create(Point.Create(exit_radius+.05,0,nozzle_length+0.2), Point.Create(-exit_radius-.05,0,nozzle_length+0.2))
connector_3 = DesignCurve.Create(GetRootPart(), curveSegment)

curveSegment = CurveSegment.Create(Point.Create(-exit_radius,0,nozzle_length), Point.Create(-exit_radius-.05,0,nozzle_length))
connector_4 = DesignCurve.Create(GetRootPart(), curveSegment)

curveSegment = CurveSegment.Create(Point.Create(-exit_radius-.05,0,nozzle_length), Point.Create(-exit_radius-.05,0,nozzle_length+0.2))
connector_5 = DesignCurve.Create(GetRootPart(), curveSegment)


#2. Create converging throat section
curveSegment = CurveSegment.Create(Point.Create(throat_radius,0,0), Point.Create(throat_radius+.05,0,-.04))
connector_6 = DesignCurve.Create(GetRootPart(), curveSegment)

curveSegment = CurveSegment.Create(Point.Create(-throat_radius,0,0), Point.Create(-throat_radius-.05,0,-.04))
connector_7 = DesignCurve.Create(GetRootPart(), curveSegment)

curveSegment = CurveSegment.Create(Point.Create(throat_radius+.05,0,-.04), Point.Create(-throat_radius-.05,0,-.04))
connector_8 = DesignCurve.Create(GetRootPart(), curveSegment)


#3. Split into separate mesh areas
curveSegment = CurveSegment.Create(Point.Create(throat_radius,0,0), Point.Create(-throat_radius,0,0))
connector_9 = DesignCurve.Create(GetRootPart(), curveSegment)

curveSegment = CurveSegment.Create(Point.Create(exit_radius,0,nozzle_length), Point.Create(-exit_radius,0,nozzle_length))
connector_10 = DesignCurve.Create(GetRootPart(), curveSegment)

curveSegment = CurveSegment.Create(Point.Create(exit_radius+.05,0,nozzle_length+0.1), Point.Create(-exit_radius-.05,0,nozzle_length+0.1))
connector_11 = DesignCurve.Create(GetRootPart(), curveSegment)



#-----------------------------Separate Surfaces-----------------------------------------------------
#1. Converging Surface
selection = Selection.Create(connector_6, connector_7, connector_8, connector_9)
secondarySelection = Selection()
options = FillOptions()
converge_surface = Fill.Execute(selection, secondarySelection, options, FillMode.ThreeD)


#2. Nozzle Surface
selection = Selection.Create(transfer_1, transfer_3, connector_9, connector_10)
secondarySelection = Selection()
options = FillOptions()
nozzle_surface = Fill.Execute(selection, secondarySelection, options, FillMode.ThreeD)

#3. exhaust box Surface close
selection = Selection.Create(connector_10, connector_11, connector_1, connector_4, connector_2, connector_5)
secondarySelection = Selection()
options = FillOptions()
exhaust_1_surface = Fill.Execute(selection, secondarySelection, options, FillMode.ThreeD)


#4. Exhuast box surface Far
selection = Selection.Create(connector_11, connector_3, connector_2, connector_5)
secondarySelection = Selection()
options = FillOptions()
exhaust_2_surface = Fill.Execute(selection, secondarySelection, options, FillMode.ThreeD)


# #-----------------------------Named Faces-----------------------------------------------------
sel = Selection.Create(converge_surface.CreatedFaces[0].Edges[2])
sel.CreateAGroup("Inlet")


sel = Selection.Create(exhaust_1_surface.CreatedFaces[0].Edges[1], exhaust_1_surface.CreatedFaces[0].Edges[3], exhaust_2_surface.CreatedFaces[0].Edges[0], exhaust_2_surface.CreatedFaces[0].Edges[2], exhaust_2_surface.CreatedFaces[0].Edges[3])
sel.CreateAGroup("Outlet")

sel = Selection.Create(converge_surface.CreatedFaces[0])
sel.CreateAGroup("converge_surface")

sel = Selection.Create(nozzle_surface.CreatedFaces[0])
sel.CreateAGroup("nozzle_surface")

sel = Selection.Create(exhaust_1_surface.CreatedFaces[0])
sel.CreateAGroup("exhaust_1_surface")

sel = Selection.Create(exhaust_2_surface.CreatedFaces[0])
sel.CreateAGroup("exhaust_2_surface")

