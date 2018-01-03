#pragma once

#include "data_handle.h"
#include "iutilities.h"

#define HANDLE(H, F) H = F; if (execution_status(H)) break
#define CALL(F) if (execution_status(F, 1)) break
//#define CALL(F) TMP_HANDLE = F; if (execution_status(TMP_HANDLE, 1)) break

inline int execution_status(void* handle, int clear = 0)
{
	int s = executionStatus(handle);
	if (s)
		std::cout << executionError(handle) << '\n';
	if (clear)
		deleteDataHandle(handle);
	return s;
}
