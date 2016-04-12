#ifndef CGADGETRON_PARAMETERS_HANDLERS
#define CGADGETRON_PARAMETERS_HANDLERS

extern "C"
void* cGT_acquisitionParameter(void* ptr_acq, const char* name);

extern "C"
void* cGT_setCSParameter(void* ptr, const char* par, const void* val);

#endif
