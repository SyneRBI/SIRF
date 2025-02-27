if (BUILD_PYTHON)
  set(PythonFiles torch.py)
  INSTALL(FILES ${PythonFiles} DESTINATION "${PYTHON_DEST}/sirf")
endif()