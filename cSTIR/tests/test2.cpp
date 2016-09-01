#include <string>

#include "cstir.h"
#include "data_handle.h"
#include "iutilities.h"

int execution_status(void* handle)
{
	int s = executionStatus(handle);
	if (s)
		std::cout << executionError(handle) << '\n';
	return s;
}

void test2()
{
	int status;
	void* handle;

	void* matrix = cSTIR_newObject("RayTracingMatrix");
	handle = cSTIR_setParameter
		(matrix, "RayTracingMatrix", "num_tangential_LORs", intDataHandle(2));
	status = execution_status(handle);
	deleteDataHandle(handle);
	deleteDataHandle(matrix);
}

