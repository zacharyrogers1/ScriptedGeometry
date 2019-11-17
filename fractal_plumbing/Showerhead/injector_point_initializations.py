#----------------------------------------------------------Random Injection Points----------------------------------------------------
combined_outlet_surface_area = math.pi*outlet_radius*outlet_radius*injector_points

(starting_injector_points, starting_weights) = create_random_injector_points_with_weights(origin_3d, injector_radius, z_axis, injector_points)
total_weights = sum(starting_weights)
#-----------------------------------------------------------------------------------------------------------------------------







#---------------------------------------------------Circular Inection Points----------------------------------------------------------------
combined_outlet_surface_area = math.pi*outlet_radius*outlet_radius*injector_points

number_of_rings = len(points_per_ring)
starting_injector_points = []
for iterator_4 in range(number_of_rings):
    single_ring_radius = injector_radius/(number_of_rings+1) * (iterator_4+1)
    temp_circle = create_circle(origin_3d, z_axis, single_ring_radius)
    starting_injector_points = starting_injector_points + points_along_circle(temp_circle, points_per_ring[iterator_4], False)

create_circle(origin_3d,z_axis,injector_radius)
starting_weights = [1.0] * injector_points
total_weights = sum(starting_weights)
#-----------------------------------------------------------------------------------------------------------------------------

#---------------------------------------------------2D array Injection Points----------------------------------------------------------------
point_2d_array = [[0, 0, 0], [0.0044661666666666669, 0.0, 0.0], [0.0022330833333333339, 0.0038678137908686005, 0.0], [-0.0022330833333333326, 0.003867813790868601, 0.0], [-0.0044661666666666669, 5.4692960772297073e-19, 0.0], [-0.0022330833333333356, -0.0038678137908685992, 0.0], [0.0022330833333333339, -0.0038678137908686005, 0.0], [0.0089323333333333338, 0.0, 0.0], [0.0079091883731430227, 0.0041510622804189559, 0.0], [0.005074143672718261, 0.0073511662181174374, 0.0], [0.0010766738069339637, 0.0088672065664018515, 0.0], [-0.003167449052692943, 0.0083518767517470217, 0.0], [-0.0066859475062470003, 0.0059232326242928636, 0.0], [-0.0086727759605219716, 0.0021376472853045611, 0.0], [-0.0086727759605219733, -0.002137647285304559, 0.0], [-0.0066859475062469977, -0.0059232326242928644, 0.0], [-0.0031674490526929452, -0.0083518767517470199, 0.0], [0.0010766738069339654, -0.0088672065664018515, 0.0], [0.0050741436727182593, -0.0073511662181174392, 0.0], [0.0079091883731430227, -0.0041510622804189542, 0.0]]
starting_injector_points = python_2d_list_2_spaceclaim_points_list(point_2d_array)
injector_points = len(starting_injector_points)
starting_weights = [1.0] * injector_points
total_weights = sum(starting_weights)
combined_outlet_surface_area = math.pi*outlet_radius*outlet_radius*injector_points
