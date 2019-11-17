outlet_area = .0000013744
sel = PowerSelection.Faces.ByArea(outlet_area, 
    PowerSelectOptions(True), 
    SearchCriteria.SizeComparison.Equal)
sel.SetActive()
power_to_selection = Selection.GetActive()
all_outlets = power_to_selection.GetItems[IDesignFace]()
for iterator_3 in range(len(all_outlets)):
	sel = Selection.Create(all_outlets[iterator_3])
	sel.CreateAGroup("Outlet " + str(iterator_3))