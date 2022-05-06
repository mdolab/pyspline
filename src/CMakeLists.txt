set(pyspline_src ${SRC_DIR}/adtProjections.F90
				 ${SRC_DIR}/basis.f90
				 ${SRC_DIR}/compute_curve.f90
				 ${SRC_DIR}/compute_surface.f90
				 ${SRC_DIR}/compute_volume.f90
				 ${SRC_DIR}/eval_curve.f90
				 ${SRC_DIR}/eval_surface.f90
				 ${SRC_DIR}/eval_volume.f90
				 ${SRC_DIR}/evaluations.f90
				 ${SRC_DIR}/findSpan.f90
				 ${SRC_DIR}/getBasisPt.f90
				 ${SRC_DIR}/insertKnot.f90
				 ${SRC_DIR}/knots.f90
				 ${SRC_DIR}/parameterizations.f90
				 ${SRC_DIR}/precision.f90
				 ${SRC_DIR}/projections.F90
				 ${SRC_DIR}/tfi2d.f90
				 )

# Add a custom target for f2py so it is in the dependency graph
add_custom_target( gen_lib_c
	DEPENDS ${SRC_DIR}/f2py/pyspline.pyf
			${SRC_DIR}/f2py/f2py_f2cmap.ref)
add_custom_command(
	OUTPUT "${CMAKE_CURRENT_BINARY_DIR}/${f2py_module_c}"
		   "${CMAKE_CURRENT_BINARY_DIR}/libspline-f2pywrappers2.f90"
	COMMAND ${Python3_EXECUTABLE} -m "numpy.f2py"
						--f2cmap "${SRC_DIR}/f2py/f2py_f2cmap.ref"
						"${SRC_DIR}/f2py/pyspline.pyf"
	DEPENDS ${SRC_DIR}/f2py/pyspline.pyf 
			${SRC_DIR}/f2py/f2py_f2cmap.ref
)

# library target
add_library(spline SHARED
	"${CMAKE_CURRENT_BINARY_DIR}/${f2py_module_c}"
	"${CMAKE_CURRENT_BINARY_DIR}/libspline-f2pywrappers2.f90"
	"${F2PY_INCLUDE_DIR}/fortranobject.c"
	"${pyspline_src}"	
)

# make the library depend on the f2py target
add_dependencies(spline gen_lib_c)

# install it to top level
install(TARGETS spline DESTINATION ${PROJECT_SOURCE_DIR}/pyspline)

# add the test import target
add_test(NAME test_import
	COMMAND ${Python3_EXECUTABLE} "${SRC_DIR}/f2py/importTest.py"
	WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
)
add_custom_target(check
	COMMAND ${CMAKE_CTEST_COMMAND}
	DEPENDS spline)
