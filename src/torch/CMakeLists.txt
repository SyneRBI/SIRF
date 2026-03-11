if (BUILD_PYTHON)
  INSTALL(FILES torch.py DESTINATION "${PYTHON_DEST}/sirf")
  add_test(NAME TORCH_TESTS_PYTHON
    COMMAND ${Python_EXECUTABLE} -m pytest -v --cov=sirf.torch --cov-config=.coveragerc-torch src/torch/tests
    WORKING_DIRECTORY ${CMAKE_SOURCE_DIR})
endif()
