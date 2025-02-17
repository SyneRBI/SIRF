if (BUILD_PYTHON)
  set(PythonFiles SIRF_torch.py)
  INSTALL(FILES ${PythonFiles} DESTINATION "${PYTHON_DEST}/sirf")
endif()