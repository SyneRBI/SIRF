//void* STIR_OSMAPOSLReconstructedImage(const void* par);
//void* STIR_OSMAPOSLReconstruction(const char* parFile);
//void* STIR_OSSPSReconstructedImage(const void* par);
//void* STIR_OSSPSReconstruction(const char* parFile);
//void* STIR_RayTracingMatrixParameters();
//void* STIR_ProjectorPairUsingMatrixParameters();
//void* STIR_MedianFilterParameters();
//void* STIR_FilterRootPriorParameters();
//void* STIR_QuadraticPriorParameters();
//void* STIR_TruncateToCylindricalFOVParameters();
//void* STIR_SeparableCartesianMetzFilterParameters();
//void* STIR_ChainedDataProcessorParameters();
//void* STIR_PoissonLogLikelihoodWithLinearModelForMeanAndProjDataParameters();
//void* STIR_OSMAPOSLReconstructionParameters();
//void* STIR_OSSPSReconstructionParameters();

//void writeText(void* ptr, char* text);

//#define SWIG_FILE_WITH_INIT

/* Let's just grab the original header file here */

/*
%include "carrays.i"
%array_class(int, intArray);
%array_class(double, doubleArray);
*/

/*
%include "numpy.i"
%init %{
import_array();
%}
%apply(int INPLACE_ARRAY1[ANY]){void* ptr, (int dim[3])}
*/

//void STIR_getImageDimensions(void* ptr, void* dim);
//void STIR_getImageDimensions(void* ptr, int dim[3]);
//void STIR_getImageData(void* ptr, double data[]);

//double* STIR_imageData(void* ptr);

/* File : pystir.i */

#include "ex.h"
#include "dh.h"
#include "tw.h"

%include "ex.h"
%include "dh.h"
%include "tw.h"

//void* STIR_reconstructedImage(const char* method, const void* par);
//void* STIR_imageFromFile(const char* filename);
//void* STIR_imagesDifference(void* first, void* second, int rimsize);
//void STIR_getImageDimensions(void* ptr, size_t pd);
//void STIR_getImageData(void* ptr, size_t pd);
//void STIR_deleteImage(void* ptr);
//const char* STIR_parameter(void* parameters, const char* name);
//void* STIR_parameters(const char* name);
////void* STIR_reconstructionParameters(const char* method);
//void* STIR_parametersFromFile(const char* name);
//int STIR_setParameterList(void* parameters, const char* name, const void* value);
//int STIR_setParameter(void* parameters, const char* name, const char* value);
//void STIR_printParameters(const void* par);
//void STIR_deleteParameters(void* par);

void* cSTIR_update(void* ptr_r, void* ptr_i);

char* cSTIR_newInt(int i);
void cSTIR_deleteInt(char* ptr);
void* cSTIR_charToVoid(char* ptr);
void printCharData(const void* ptr);

void* cSTIR_reconstructionSetup(void* ptr_r, void* ptr_i);

//void* cSTIR_setupObjective(const char* obj, void* ptr_obj);

//void* cSTIR_setupReconstruction(const char* method, void* ptr_r, void* ptr_i);
//void* cSTIR_update(const char* method, void* ptr_r, void* ptr_i);
//void* cSTIR_setSubiteration(const char* method, void* ptr_r, int i);
//void* cSTIR_getSubiteration(const char* method, void* ptr_r);

