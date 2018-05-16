# Install python packages via pip and setup.py
if(BUILD_PYTHON)
  set(PYTHON_SETUP_PKGS "sirf" CACHE INTERNAL "list of provided python packages")

  # alias sirf.p* -> p* for backward-compatibility
  function(python_pkg_alias PY_PKG_NEW PY_PKG_OLD)
    list(APPEND PYTHON_SETUP_PKGS ${PY_PKG_NEW})
    set(PYTHON_SETUP_PKGS "${PYTHON_SETUP_PKGS}" PARENT_SCOPE)
    set(SETUP_PY_INIT_IN "${CMAKE_CURRENT_LIST_DIR}/__init__.py.in")
    set(SETUP_PY_INIT "${PYTHON_DEST}/${PY_PKG_NEW}/__init__.py")
    # configure_file("${SETUP_PY_INIT_IN}" "${SETUP_PY_INIT}")  # v1
    # install(CODE "configure_file(\"${SETUP_PY_INIT_IN}\" \"${SETUP_PY_INIT}\")")  # v2
    configure_file("${SETUP_PY_INIT_IN}"
      "${CMAKE_CURRENT_BINARY_DIR}/${PY_PKG_NEW}/__init__.py")
    install(FILES "${CMAKE_CURRENT_BINARY_DIR}/${PY_PKG_NEW}/__init__.py"
      DESTINATION "${PYTHON_DEST}/${PY_PKG_NEW}")
    # message(STATUS "setup.py:${SETUP_PY_INIT}")
    message(STATUS "setup.py:${PY_PKG_NEW}<-${PY_PKG_OLD}")
  endfunction(python_pkg_alias)
  python_pkg_alias(pGadgetron "sirf.Gadgetron")
  python_pkg_alias(pSTIR "sirf.STIR")
  python_pkg_alias(pUtilities "sirf.Utilities")
  python_pkg_alias(pygadgetron "sirf.pygadgetron")
  python_pkg_alias(pystir "sirf.pystir")
  python_pkg_alias(pyiutilities "sirf.pyiutilities")
  # convert to python CSV tuple for setup.py configure_file
  string(REPLACE ";" "', '" PYTHON_SETUP_PKGS_CSV "${PYTHON_SETUP_PKGS}")
  set(PYTHON_SETUP_PKGS_CSV "'${PYTHON_SETUP_PKGS_CSV}'")
  # message(STATUS "setup.py:pacakges:${PYTHON_SETUP_PKGS_CSV}")

  # Create setup.py
  set(SETUP_PY_IN "${CMAKE_CURRENT_LIST_DIR}/setup.py.in")
  set(SETUP_PY "${PYTHON_DEST}/setup.py")
  set(SETUP_PY_INIT "${PYTHON_DEST}/sirf/__init__.py")
  message(STATUS "setup.py:${SETUP_PY}")
  # configure_file("${SETUP_PY_IN}" "${SETUP_PY}")  # v1
  # install(CODE "configure_file(\"${SETUP_PY_IN}\" \"${SETUP_PY}\")")  # v2
  configure_file("${SETUP_PY_IN}" "${CMAKE_CURRENT_BINARY_DIR}/setup.py")
  install(FILES "${CMAKE_CURRENT_BINARY_DIR}/setup.py"
    DESTINATION "${PYTHON_DEST}")

  if(PYTHONINTERP_FOUND)
    # python setup.py build
    add_custom_command(OUTPUT "${SETUP_PY_INIT}"
      COMMAND "${CMAKE_COMMAND}" -E make_directory "${PYTHON_DEST}/sirf"
      COMMAND "${CMAKE_COMMAND}" -E touch "${SETUP_PY_INIT}"
      # COMMAND "${PYTHON_EXECUTABLE}" setup.py build
      DEPENDS "${SETUP_PY}"
      WORKING_DIRECTORY "${PYTHON_DEST}")
    add_custom_target(pybuild_sirf ALL DEPENDS "${SETUP_PY_INIT}")

    # python setup.py install
    if("${PYTHON_STRATEGY}" STREQUAL "SETUP_PY")
      install(CODE "execute_process(COMMAND\n\
        \"${PYTHON_EXECUTABLE}\" setup.py install\n\
        WORKING_DIRECTORY \"${PYTHON_DEST}\")")
    endif()
  endif(PYTHONINTERP_FOUND)
endif(BUILD_PYTHON)
