sel = Selection.SelectAll()
ComponentHelper.DeleteEmptyComponents(sel)
Delete.Execute(sel)

point_1 =  [0.0,0.0,0.0]
point_2 = [0.06, 0.0, 0.06]

direction_1 = [0,0,0.1]
direction_2 = [0,0,0.1]

curves_0 = []
curves_1 = []
point_list = List[Point]()


for i in range (21):
	t = i*.05
	x = (2*t**3 - 3*t**2 + 1)*point_1[0] + (3*t**2 - 2*t**3)*point_2[0] + (-2*t**2 + t**3 + t)*direction_1[0] + (-t**2 + t**3)*direction_2[0]
	y = (2*t**3 - 3*t**2 + 1)*point_1[1] + (3*t**2 - 2*t**3)*point_2[1] + (-2*t**2 + t**3 + t)*direction_1[1] + (-t**2 + t**3)*direction_2[1]
	z = (2*t**3 - 3*t**2 + 1)*point_1[2] + (3*t**2 - 2*t**3)*point_2[2] + (-2*t**2 + t**3 + t)*direction_1[2] + (-t**2 + t**3)*direction_2[2]
	x_prime = (12*t -6)*point_1[0] + (6-12*t)*point_2[0] + (-4 + 6*t)*direction_1[0] + (-4 + 6*t)*direction_2[0]
	z_prime = (12*t -6)*point_1[2] + (6-12*t)*point_2[2] + (-4 + 6*t)*direction_1[2] + (-4 + 6*t)*direction_2[2]
	magnitude = (x_prime**2 + z_prime**2)**0.5

	trident = [x,y,z]
	print("This is i: " + str(i) + " : This is the point i produced: " + str(trident))
	print("This is z double prime evaluation: " +str(z_prime))
	print("This is magnitude: " + str(magnitude))
	print("\n\n\n")

	DatumPointCreator.Create(Point.Create(x,y,z))

	temp_point = Point.Create(x, y, z)
	point_list.Add(temp_point)


ncurve = NurbsCurve.CreateThroughPoints(False, point_list,1)
curves_0.append(CurveSegment.Create(ncurve))
curves_1.append(DesignCurve.Create(GetRootPart(), curves_0[0])) 




