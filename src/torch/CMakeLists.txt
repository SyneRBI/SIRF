if (BUILD_PYTHON)
  INSTALL(FILES torch.py DESTINATION "${PYTHON_DEST}/sirf")
  add_test(NAME TORCH_TESTS_PYTHON
    COMMAND ${Python_EXECUTABLE} -m pytest -v -k "not (PET and (test_objective_function_gradcheck or test_objective_function_gradient_gradcheck))" --cov=sirf.torch --cov-config=.coveragerc-torch src/torch/tests
    WORKING_DIRECTORY ${CMAKE_SOURCE_DIR})
endif()
