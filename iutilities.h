#ifndef INTERFACE_UTILITIES
#define INTERFACE_UTILITIES

#ifndef CGADGETRON_FOR_MATLAB
extern "C" {
#endif
	void* newDataHandle();
	void deleteDataHandle(void* ptr);
	void* copyOfObjectHandle(void* ptr);
	void deleteObject(void* ptr);
	void* charDataHandle(const char* s);
	void* intDataHandle(int i);
	void* floatDataHandle(float i);
	void* doubleDataHandle(double i);
	char* charDataFromHandle(const void* ptr);
	int intDataFromHandle(const void* ptr);
	float floatDataFromHandle(const void* ptr);
	double doubleDataFromHandle(const void* ptr);
	double doubleReDataFromHandle(const void* ptr);
	double doubleImDataFromHandle(const void* ptr);
	int executionStatus(const void* ptr);
	const char* executionError(const void* ptr);
	const char* executionErrorFile(const void* ptr);
	int executionErrorLine(const void* ptr);
#ifndef CGADGETRON_FOR_MATLAB
}
#endif

#endif