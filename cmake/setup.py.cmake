# Install python packages via pip and setup.py
if(BUILD_PYTHON)
  set(PYTHON_SETUP_PKGS "sirf" CACHE INTERNAL "list of provided python packages")

  # alias sirf.p* -> p* for backward-compatibility
  function(python_pkg_alias PY_PKG_NEW PY_PKG_OLD)
    list(APPEND PYTHON_SETUP_PKGS ${PY_PKG_NEW})
    set(PYTHON_SETUP_PKGS "${PYTHON_SETUP_PKGS}" PARENT_SCOPE)
    set(PY_PKG_INIT_IN "${CMAKE_CURRENT_LIST_DIR}/__init__.py.in")
    set(PY_PKG_INIT "${CMAKE_CURRENT_BINARY_DIR}/cmake/${PY_PKG_NEW}/__init__.py")
    configure_file("${PY_PKG_INIT_IN}" "${PY_PKG_INIT}")  # depends: PY_PKG_OLD
    install(FILES "${PY_PKG_INIT}"
      DESTINATION "${PYTHON_DEST}/${PY_PKG_NEW}")
    # message(STATUS "__init__.py:${PY_INIT}")
    message(STATUS "python alias:${PY_PKG_NEW}<-${PY_PKG_OLD}")
  endfunction(python_pkg_alias)
  python_pkg_alias(pGadgetron "sirf.Gadgetron")
  python_pkg_alias(pSTIR "sirf.STIR")
  python_pkg_alias(pUtilities "sirf.Utilities")
  python_pkg_alias(pygadgetron "sirf.pygadgetron")
  python_pkg_alias(pystir "sirf.pystir")
  python_pkg_alias(pyiutilities "sirf.pyiutilities")
  python_pkg_alias(pReg "sirf.Reg")
  python_pkg_alias(pyreg "sirf.pyreg")
  # convert to python CSV tuple for setup.py configure_file
  string(REPLACE ";" "', '" PYTHON_SETUP_PKGS_CSV "${PYTHON_SETUP_PKGS}")
  set(PYTHON_SETUP_PKGS_CSV "'${PYTHON_SETUP_PKGS_CSV}'")
  # message(STATUS "setup.py:pacakges:${PYTHON_SETUP_PKGS_CSV}")

  # Create setup.py
  set(SETUP_PY_IN "${CMAKE_CURRENT_LIST_DIR}/setup.py.in")
  set(SETUP_PY "${CMAKE_CURRENT_BINARY_DIR}/cmake/setup.py")
  configure_file("${SETUP_PY_IN}" "${SETUP_PY}")
  install(FILES "${SETUP_PY}" DESTINATION "${PYTHON_DEST}")
  message(STATUS "setup.py:${PYTHON_DEST}/setup.py")
  # create sirf/__init__.py
  set(PY_MOD_INIT_IN "${CMAKE_CURRENT_LIST_DIR}/sirf.__init__.py.in")
  set(PY_MOD_INIT "${CMAKE_CURRENT_BINARY_DIR}/cmake/__init__.py")
  configure_file("${PY_MOD_INIT_IN}" "${PY_MOD_INIT}")
  install(FILES "${PY_MOD_INIT}" DESTINATION "${PYTHON_DEST}/sirf")

  if(PYTHONINTERP_FOUND)
    # python setup.py install
    if("${PYTHON_STRATEGY}" STREQUAL "SETUP_PY")
      install(CODE "execute_process(COMMAND\n\
        \"${PYTHON_EXECUTABLE}\" setup.py build install\n\
        WORKING_DIRECTORY \"${PYTHON_DEST}\")")
    endif()
  endif(PYTHONINTERP_FOUND)
endif(BUILD_PYTHON)
