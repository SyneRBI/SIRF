extern "C"
int doOSMAPOSLReconstruction(char* parFile) {

	USING_NAMESPACE_STIR

		try {
		OSMAPOSLReconstruction<DiscretisedDensity<3, float> >
			reconstruction_object(parFile);
		if (reconstruction_object.reconstruct() == Succeeded::yes) {
			//_getch();
			return EXIT_SUCCESS;
		}
		else {
			return EXIT_FAILURE;
		}
	}
	catch (string s) {
		cout << s << endl;
		return -1;
	}
	catch (StirException& se) {
		cout << se.what() << endl;
		cout << se.file() << endl;
		cout << se.line() << endl;
		return -1;
	}
	catch (...) {
		cout << "exception occured" << endl;
		return -1;
	}
}

//cout << se.what() << endl;
//cout << se.file() << endl;
//cout << se.line() << endl;

int doOSMAPOSLReconstruction(char* parFile);

//boost::shared_ptr< DiscretisedDensity<3, float> >* ptr_first;
//boost::shared_ptr< DiscretisedDensity<3, float> >* ptr_second;

//	ptr_first = &first_operand;
//}
//catch (StirException& se) {
//	ExecutionStatus status(se);
//	handle->set(0, &status);
//	return (void*)handle;
//}

//try {

//	ptr_second = &second_operand;
//}
//catch (StirException& se) {
//	ExecutionStatus status(se);
//	handle->set(0, &status);
//	return (void*)handle;
//}

//PostFiltering::PostFiltering()
//{
//	filter_ptr.reset();
//	parser.add_start_key("PostFilteringParameters");
//	parser.add_parsing_key("PostFilter type", &filter_ptr);
//	parser.add_stop_key("END PostFiltering Parameters");
//}

//boost::shared_ptr< DiscretisedDensity<3, float> >* ptr =
//	new boost::shared_ptr< DiscretisedDensity<3, float> >
//	(read_from_file<DiscretisedDensity<3, float> >(filename));

//delete (boost::shared_ptr< DiscretisedDensity<3, float> >*)handle->data();

//boost::shared_ptr< DiscretisedDensity<3, float> >
//	first_operand(read_from_file<DiscretisedDensity<3, float> >(image_a));
//if (is_null_ptr(first_operand))
//{
//	ExecutionStatus status("Could not read first file", __FILE__, __LINE__);
//	handle->set(0, &status);
//	return (void*)handle;
//}

//if (is_null_ptr(second_operand))
//{
//	ExecutionStatus status("Could not read second file", __FILE__, __LINE__);
//	handle->set(0, &status);
//	return (void*)handle;
//}

//void* ptr = STIR_imageFromFile(image_a);
//DataHandle* first_h = (DataHandle*)ptr;

//boost::shared_ptr< DiscretisedDensity<3, float> >* first_ptr =
//	(boost::shared_ptr< DiscretisedDensity<3, float> >*)first_h->data();
//boost::shared_ptr< DiscretisedDensity<3, float> >& first_operand = *first_ptr;

//boost::shared_ptr< DiscretisedDensity<3, float> >
//	second_operand(read_from_file<DiscretisedDensity<3, float> >(image_b));

//ptrImage3DF* first_ptr = (ptrImage3DF*)first_h->data();
//ptrImage3DF& first_operand = *first_ptr;
//ptrImage3DF* second_ptr = (ptrImage3DF*)second_h->data();
//ptrImage3DF& second_operand = *second_ptr;
//{
//	std::string explanation;
//	if (!first_operand->has_same_characteristics(*second_operand,
//		explanation))
//	{
//		warning("input images do not have the same characteristics.\n%s",
//			explanation.c_str());
//		ExecutionStatus status(
//			"input images do not have the same characteristics",
//			__FILE__, __LINE__);
//		handle->set(0, &status);
//		return (void*)handle;
//	}
//}
//if (rimsize >= 0)
//{
//	truncate_rim(*first_operand, rimsize);
//	truncate_rim(*second_operand, rimsize);
//}

//float reference_max = first_operand->find_max();
//float reference_min = first_operand->find_min();

//float amplitude = fabs(reference_max)>fabs(reference_min) ?
//	fabs(reference_max) : fabs(reference_min);

//*first_operand -= *second_operand;
//const float max_error = first_operand->find_max();
//const float min_error = first_operand->find_min();
//in_place_abs(*first_operand);
//const float max_abs_error = first_operand->find_max();

//const bool same = (max_abs_error / amplitude <= tolerance);

//cout << "\nMaximum absolute error = " << max_abs_error
//	<< "\nMaximum in (1st - 2nd) = " << max_error
//	<< "\nMinimum in (1st - 2nd) = " << min_error << endl;
//cout << "Error relative to sup-norm of first image = " 
//	<< (max_abs_error / amplitude) * 100 << " %" << endl;

//cout << "\nImage arrays ";
//if (same)
//{
//	cout << (max_abs_error == 0 ? "are " : "deemed ")
//		<< "identical\n";
//}
//else
//{
//	cout << "deemed different\n";
//}
//cout << "(tolerance used: " << tolerance * 100 << " %)\n\n";
//if (!same) {
//	int* result = (int*)malloc(sizeof(int));
//	*result = 1;
//	handle->set(result, 0, GRAB);
//}

#if 0
extern "C"
void* STIR_imageComparisonResult
(void* first, void* second, int rimsize, double tolerance) {

	DataHandle* handle = new DataHandle;
	try {

		DataHandle* first_h = (DataHandle*)first;
		DataHandle* second_h = (DataHandle*)second;

		ptrImage3DF* first_ptr = (ptrImage3DF*)first_h->data();
		ptrImage3DF& first_operand = *first_ptr;
		ptrImage3DF* second_ptr = (ptrImage3DF*)second_h->data();
		ptrImage3DF& second_operand = *second_ptr;

		{
			std::string explanation;
			if (!first_operand->has_same_characteristics(*second_operand,
				explanation))
			{
				warning("input images do not have the same characteristics.\n%s",
					explanation.c_str());
				ExecutionStatus status(
					"input images do not have the same characteristics",
					__FILE__, __LINE__);
				handle->set(0, &status);
				return (void*)handle;
			}
		}

		if (rimsize >= 0)
		{
			truncate_rim(*first_operand, rimsize);
			truncate_rim(*second_operand, rimsize);
		}

		float reference_max = first_operand->find_max();
		float reference_min = first_operand->find_min();

		float amplitude = fabs(reference_max) > fabs(reference_min) ?
			fabs(reference_max) : fabs(reference_min);

		*first_operand -= *second_operand;
		const float max_error = first_operand->find_max();
		const float min_error = first_operand->find_min();

		in_place_abs(*first_operand);
		const float max_abs_error = first_operand->find_max();

		const bool same = (max_abs_error / amplitude <= tolerance);

		cout << "\nMaximum absolute error = " << max_abs_error
			<< "\nMaximum in (1st - 2nd) = " << max_error
			<< "\nMinimum in (1st - 2nd) = " << min_error << endl;
		cout << "Error relative to sup-norm of first image = "
			<< (max_abs_error / amplitude) * 100 << " %" << endl;

		cout << "\nImage arrays ";
		if (same)
		{
			cout << (max_abs_error == 0 ? "are " : "deemed ")
				<< "identical\n";
		}
		else
		{
			cout << "deemed different\n";
		}
		cout << "(tolerance used: " << tolerance * 100 << " %)\n\n";
		if (!same) {
			int* result = (int*)malloc(sizeof(int));
			*result = 1;
			handle->set(result, 0, GRAB);
		}
	}
	catch (StirException& se) {
		ExecutionStatus status(se);
		handle->set(0, &status);
		return (void*)handle;
	}
	catch (exception& e) {
		cout << "exception " << e.what() << endl;
	}
	catch (...) {
		cout << "unhandled exception occured" << endl;
		ExecutionStatus status("unhandled exception", __FILE__, __LINE__);
		handle->set(0, &status);
		return (void*)handle;
	}
	return (void*)handle;

}

extern "C"
void* STIR_compareImagesResult
(const char* image_a, const char* image_b, int rimsize, double tolerance) {

	DataHandle* handle = new DataHandle;
	try {
		DataHandle* first_h = (DataHandle*)STIR_imageFromFile(image_a);
		DataHandle* second_h = (DataHandle*)STIR_imageFromFile(image_b);
		handle = (DataHandle*)STIR_imageComparisonResult
			(first_h, second_h, rimsize, tolerance);
		STIR_deleteImage(first_h);
		STIR_deleteImage(second_h);
	}
	catch (StirException& se) {
		ExecutionStatus status(se);
		handle->set(0, &status);
		return (void*)handle;
	}
	return (void*)handle;
}

extern "C"
void* STIR_postFilterResult(char* in, char* out, char* par) {

	//	USING_NAMESPACE_STIR

	boost::shared_ptr<DiscretisedDensity<3, float> > input_image_ptr;
	PostFiltering post_filtering;
	std::string out_filename;
	bool verbose = false;

	DataHandle* handle = new DataHandle;

	try {
		input_image_ptr =
			read_from_file<DiscretisedDensity<3, float> >(in);
		if (is_null_ptr(input_image_ptr))
		{
			warning("postfilter: No input image. Not writing any output.\n");
			ExecutionStatus status("postfilter: no input image", __FILE__, __LINE__);
			handle->set(0, &status);
			return (void*)handle;
		}

		if (post_filtering.parser.parse(par) == false)
		{
			warning
				("postfilter aborting because error in parsing. Not writing any output");
			ExecutionStatus status("postfilter: error in parsing", __FILE__, __LINE__);
			handle->set(0, &status);
			return (void*)handle;
		}
		if (is_null_ptr(post_filtering.filter_ptr))
		{
			warning("postfilter: No filter set. Not writing any output.\n");
			ExecutionStatus status("postfilter: no filter set", __FILE__, __LINE__);
			handle->set(0, &status);
			return (void*)handle;
		}

		post_filtering.filter_ptr->apply(*input_image_ptr);

		if (!(OutputFileFormat<DiscretisedDensity<3, float> >::default_sptr()->
			write_to_file(out, *input_image_ptr) == Succeeded::yes)) {
			ExecutionStatus status("postfilter: error writing", __FILE__, __LINE__);
			handle->set(0, &status);
			return (void*)handle;
		}
	}
	catch (StirException& se) {
		ExecutionStatus status(se);
		handle->set(0, &status);
		return (void*)handle;
	}

	return (void*)handle;
}

extern "C"
void* STIR_generateImageResult(char* par) {

	//	USING_NAMESPACE_STIR

	DataHandle* handle = new DataHandle;

	try {
		GenerateImage application(par);
		Succeeded success = application.compute();
		if (success != Succeeded::yes) {
			ExecutionStatus status("generateImage failure", __FILE__, __LINE__);
			handle->set(0, &status);
		}
		return (void*)handle;
		//return success == Succeeded::yes ? EXIT_SUCCESS : EXIT_FAILURE;
	}
	catch (StirException& se) {
		ExecutionStatus status(se);
		handle->set(0, &status);
		return (void*)handle;
	}

}
#endif

#if 0
extern "C"
void* STIR_OSMAPOSLReconstructionResult(const char* parFile) {

	DataHandle* handle = new DataHandle;
	try {
		OSMAPOSLReconstruction<Image3DF> reconstruction_object(parFile);
		ptrImage3DF image(reconstruction_object.get_initial_data_ptr());
		if (
			reconstruction_object.set_up(image) != Succeeded::yes ||
			reconstruction_object.reconstruct(image) != Succeeded::yes) {
			ExecutionStatus status("OSMAPOSLReconstruction::reconstruct() failed",
				__FILE__, __LINE__);
			handle->set(0, &status);
		}
	}
	catch (StirException& se) {
		ExecutionStatus status(se);
		handle->set(0, &status);
	}
	return (void*)handle;
}
#endif

//void* STIR_generateImageResult(char* par);
//void* STIR_postFilterResult(char* in, char* out, char* par);
//void* STIR_OSMAPOSLReconstructionResult(const char* parFile);

//void* STIR_imageComparisonResult
//	(void* first, void* second, int rimsize, double tolerance);
//void* STIR_compareImagesResult
//	(const char* image_a, const char* image_b, int rimsize, double tolerance);

//if (
//	reconstruction_object.set_up(*image) != Succeeded::yes ||
//	reconstruction_object.reconstruct(*image) != Succeeded::yes) {
//	ExecutionStatus status("OSMAPOSLReconstruction::reconstruct() failed",
//		__FILE__, __LINE__);
//	handle->set(0, &status);
//	return (void*)handle;
//}

//if (
//	reconstruction_object.set_up(*image) != Succeeded::yes ||
//	reconstruction_object.reconstruct(*image) != Succeeded::yes) {
//	ExecutionStatus status("OSMAPOSLReconstruction::reconstruct() failed",
//		__FILE__, __LINE__);
//	handle->set(0, &status);
//	return (void*)handle;
//}

//int parse_parameters_file(ifstream& in, NestedList* list);

//par->set(0, "Matrix Type", (NestedList*)STIR_RayTracingMatrixParameters());

//par->set(1, "Filter Type", (NestedList*)STIR_MedianFilterParameters());

//par->set(0, "Data Processor to apply first",
//	(NestedList*)STIR_TruncateToCylindricalFOVParameters());
//par->set(1, "Data Processor to apply second",
//	(NestedList*)STIR_SeparableCartesianMetzFilterParameters());

//par->set(5, "projector pair type",
//	(NestedList*)STIR_ProjectorPairUsingMatrixParameters());
//par->set(6, "prior type",
//	(NestedList*)STIR_FilterRootPriorParameters());

//par->set(0, "objective function type",
//	(NestedList*)
//	STIR_PoissonLogLikelihoodWithLinearModelForMeanAndProjDataParameters());

//par->set(7, "inter-iteration filter type",
//	(NestedList*)STIR_ChainedDataProcessorParameters());

//if (!par->set("inter-iteration filter type.Data Processor to apply first",
//	(NestedList*)STIR_ChainedDataProcessorParameters()))
//	cout << "first item not found" << endl << endl;
//if (!par->set("inter-iteration filter type.Data Processor to apply second",
//	(NestedList*)STIR_TruncateToCylindricalFOVParameters()))
//	cout << "second item not found" << endl << endl;

//return (void*)par;

//NestedList* rp = (NestedList*)par;

//return ((NestedList*)parameters)->set(name, (NestedList*)value);

//return ((NestedList*)parameters)->set(name, value);

//NestedList* par = (NestedList*)STIR_OSMAPOSLReconstructionParameters();

//int status = parse_parameters_file(in, par);
//if (status) {

//ExecutionStatus status(
//	"parametersFromFile parsing failure",
//	__FILE__, __LINE__);

// if (!par->set("inter-iteration filter type.Data Processor to apply first",
//(NestedList*)hp->data()))
//cout << "first item not found" << endl << endl;

//if (!par->set("inter-iteration filter type.Data Processor to apply second",
//	(NestedList*)hp->data()))
//	cout << "second item not found" << endl << endl;

//if (!par->set(name,	"1"))
//	cout << "third item not found" << endl << endl;

// EO 14 Aug 15, 25 Aug 15
//char buff[128];
//sprintf_s(buff, 128, "Now processing subset #%d\n", subset_num);
//writeText(buff);

//char buff[128];
//sprintf_s(buff, 128, 
//	"Number of (cancelled) singularities in Sensitivity division: %d\n", 
//	count);
//writeText(buff);

//char buff[128];
//sprintf_s(buff, 128, 
//	"Update image old min,max: %f, %f, new min, max: %f, %f\n",
//	current_min, current_max, 
//	max(current_min, new_min), 
//	min(current_max, new_max));
//writeText(buff);

void* STIR_parameters0(char* name) {
	void* par;
	if (boost::iequals(name, "Ray Tracing"))
		par = STIR_RayTracingMatrixParameters();
	//else if (boost::iequals(name,
	//	"Forward Projector Using Ray Tracing Parameters"))
	//	par = STIR_RayTracingForwardProjectorParameters();
	else if (boost::iequals(name, "Matrix"))
		par = STIR_ProjectorPairUsingMatrixParameters();
	else if (boost::iequals(name,
		"Separate Projectors"))
		par = STIR_ProjectorPairUsingSeparateProjectorsParameters();
	else if (boost::iequals(name, "Median"))
		par = STIR_MedianFilterParameters();
	else if (boost::iequals(name, "FilterRootPrior"))
		par = STIR_FilterRootPriorParameters();
	else if (boost::iequals(name, "Quadratic"))
		par = STIR_QuadraticPriorParameters();
	else if (boost::iequals(name, "Truncate To Cylindrical FOV"))
		par = STIR_TruncateToCylindricalFOVParameters();
	else if (boost::iequals(name, "Separable Cartesian Metz"))
		par = STIR_SeparableCartesianMetzFilterParameters();
	else if (boost::iequals(name, "Chained Data Processor"))
		par = STIR_ChainedDataProcessorParameters();
	else if (boost::iequals(name, "Interfile Output File Format"))
		par = STIR_InterfileOutputFileFormatParameters();
	else if (boost::iequals(name,
		"PoissonLogLikelihoodWithLinearModelForMeanAndProjData"))
		par =
		STIR_PoissonLogLikelihoodWithLinearModelForMeanAndProjDataParameters();
	else if (boost::iequals(name, "OSMAPOSL"))
		par = STIR_OSMAPOSLReconstructionParameters();
	else if (boost::iequals(name, "OSSPS"))
		par = STIR_OSSPSReconstructionParameters();
	else {
		DataHandle* h = new DataHandle;
		string error = "unknown parameter '";
		error += name;
		error += "'";
		ExecutionStatus status(error.c_str(), __FILE__, __LINE__);
		h->set(0, &status);
		return (void*)h;
	}
	return par;
}

#ifdef _WIN32
#define NULLDEV  "NUL"
#define DELETE_FILE "del "
#else
#define NULLDEV "/dev/null"
#define DELETE_FILE "rm -f"
#endif
#define COMMAND_SIZE 1024

#if 0
void redirectOutput(char command[], char* log, char* errlog) {
	if (log == 0 || log && strcmp(log, "stdout")) {
		strcat_s(command, COMMAND_SIZE, " 1>");
		//if (log && strlen(log) > 0)
		if (log && strcmp(log, "NUL") != 0)
			strcat_s(command, COMMAND_SIZE, log);
		else
			strcat_s(command, COMMAND_SIZE, NULLDEV);
	}
	if (errlog == 0 || errlog && strcmp(errlog, "stderr")) {
		strcat_s(command, COMMAND_SIZE, " 2>");
		if (errlog && strcmp(errlog, "stdout") == 0)
			strcat_s(command, COMMAND_SIZE, "&1");
		//else if (strlen(errlog) > 0)
		else if (errlog && strcmp(errlog, "NUL"))
			strcat_s(command, COMMAND_SIZE, errlog);
		else if (log && strcmp(log, "NUL"))
			strcat_s(command, COMMAND_SIZE, NULLDEV);
		else
			strcat_s(command, COMMAND_SIZE, "&1");
	}
}

extern "C"
int generateImage(char stirPath[], char parFile[], char* log, char* errlog, int verb) {
	char command[COMMAND_SIZE] = { 0 };
	char code[] = "generate_image";
	strcat_s(command, COMMAND_SIZE, stirPath);
	strcat_s(command, COMMAND_SIZE, code);
	strcat_s(command, COMMAND_SIZE, " ");
	strcat_s(command, COMMAND_SIZE, parFile);
	redirectOutput(command, log, errlog);
	if (verb)
		printf("%s\n", command);
	return system(command);
}

extern "C"
int postFilter(
char stirPath[], char fin[], char fout[], char parFile[],
char* log, char* errlog, int verb
) {
	char command[COMMAND_SIZE] = { 0 };
	char code[] = "postfilter";
	strcat_s(command, COMMAND_SIZE, stirPath);
	strcat_s(command, COMMAND_SIZE, code);
	strcat_s(command, COMMAND_SIZE, " ");
	strcat_s(command, COMMAND_SIZE, fout);
	strcat_s(command, COMMAND_SIZE, " ");
	strcat_s(command, COMMAND_SIZE, fin);
	strcat_s(command, COMMAND_SIZE, " ");
	strcat_s(command, COMMAND_SIZE, parFile);
	redirectOutput(command, log, errlog);
	if (verb)
		printf("%s\n", command);
	return system(command);
}

extern "C"
int runOSMAPOSL
(char stirPath[], char parFile[], char* log, char* errlog, int verb) {
	char command[COMMAND_SIZE] = { 0 };
	char code[] = "OSMAPOSL";
	strcat_s(command, COMMAND_SIZE, stirPath);
	strcat_s(command, COMMAND_SIZE, code);
	strcat_s(command, COMMAND_SIZE, " ");
	strcat_s(command, COMMAND_SIZE, parFile);
	redirectOutput(command, log, errlog);
	if (verb)
		printf("%s\n", command);
	return system(command);
}

extern "C"
int runOSSPS
(char stirPath[], char parFile[], char* log, char* errlog, int verb) {
	char command[COMMAND_SIZE] = { 0 };
	char code[] = "OSSPS";
	strcat_s(command, COMMAND_SIZE, stirPath);
	strcat_s(command, COMMAND_SIZE, code);
	strcat_s(command, COMMAND_SIZE, " ");
	strcat_s(command, COMMAND_SIZE, parFile);
	redirectOutput(command, log, errlog);
	if (verb)
		printf("%s\n", command);
	return system(command);
}

extern "C"
int compareImage
(char stirPath[], char a[], char b[], char* log, char* errlog, int verb) {
	char command[COMMAND_SIZE] = { 0 };
	char code[] = "compare_image";
	strcat_s(command, COMMAND_SIZE, stirPath);
	strcat_s(command, COMMAND_SIZE, code);
	strcat_s(command, COMMAND_SIZE, " ");
	strcat_s(command, COMMAND_SIZE, a);
	strcat_s(command, COMMAND_SIZE, " ");
	strcat_s(command, COMMAND_SIZE, b);
	redirectOutput(command, log, errlog);
	if (verb)
		printf("%s\n", command);
	return system(command);
}
#endif

void* STIR_OSMAPOSLReconstruction(const char* parFile) {

	DataHandle* handle = new DataHandle;
	try {
		OSMAPOSLReconstruction<Image3DF> reconstruction_object(parFile);
		ptrImage3DF* image =
			new ptrImage3DF(reconstruction_object.get_initial_data_ptr());
		if (reconstruction_object.set_up(*image) != Succeeded::yes) {
			ExecutionStatus status("OSMAPOSLReconstruction::set_up failed",
				__FILE__, __LINE__);
			handle->set(0, &status);
			return (void*)handle;
		}
		if (reconstruction_object.reconstruct(*image) != Succeeded::yes) {
			ExecutionStatus status("OSMAPOSLReconstruction::reconstruct failed",
				__FILE__, __LINE__);
			handle->set(0, &status);
			return (void*)handle;
		}
		handle->set((void*)image, 0, GRAB);
	}
	catch (StirException& se) {
		ExecutionStatus status(se);
		handle->set(0, &status);
	}
	return (void*)handle;
}

//extern "C"
//void* STIR_OSMAPOSLReconstructedImage(const void* par) {
//	char filename[64];
//	make_par_file(par, filename, 64);
//	void* image = STIR_OSMAPOSLReconstruction(filename);
//	string cmd = DELETE_FILE;
//	cmd += filename;
//	//cout << cmd << endl;
//	system(cmd.c_str());
//	return image;
//}

//extern "C"
void* STIR_OSSPSReconstruction(const char* parFile) {

	DataHandle* handle = new DataHandle;
	try {
		OSSPSReconstruction<Image3DF> reconstruction_object(parFile);
		ptrImage3DF* image =
			new ptrImage3DF(reconstruction_object.get_initial_data_ptr());
		if (reconstruction_object.set_up(*image) != Succeeded::yes) {
			ExecutionStatus status("OSSPSReconstruction::set_up failed",
				__FILE__, __LINE__);
			handle->set(0, &status);
			return (void*)handle;
		}
		if (reconstruction_object.reconstruct(*image) != Succeeded::yes) {
			ExecutionStatus status("OSSPSReconstruction::reconstruct failed",
				__FILE__, __LINE__);
			handle->set(0, &status);
			return (void*)handle;
		}
		handle->set((void*)image, 0, GRAB);
	}
	catch (StirException& se) {
		ExecutionStatus status(se);
		handle->set(0, &status);
	}
	return (void*)handle;
}

//extern "C"
//void* STIR_OSSPSReconstructedImage(const void* par) {
//	char filename[64];
//	make_par_file(par, filename, 64);
//	void* image = STIR_OSSPSReconstruction(filename);
//	string cmd = DELETE_FILE;
//	cmd += filename;
//	system(cmd.c_str());
//	return image;
//}

//image = STIR_OSMAPOSLReconstruction(filename);

//image = STIR_OSSPSReconstruction(filename);

//string cmd = DELETE_FILE;
//cmd += filename;
//system(cmd.c_str());

//int generateImage
//	(char stirPath[], char parFile[], char* log, char* errlog, int verb);
//int postFilter
//	(char stirPath[], char fin[], char fout[], char parFile[],
//	char* log, char* errlog, int verb);
//int runOSMAPOSL
//	(char stirPath[], char parFile[], char* log, char* errlog, int verb);
//int runOSSPS
//	(char stirPath[], char parFile[], char* log, char* errlog, int verb);
//int compareImage
//	(char stirPath[], char a[], char b[], char* log, char* errlog, int verb);

//char* keyValue(char parFile[], char key[]);
//void deleteString(char* s);

//void* STIR_OSMAPOSLReconstruction(const char* parFile);
//void* STIR_OSMAPOSLReconstructedImage(const void* par);
//void* STIR_OSSPSReconstruction(const char* parFile);
//void* STIR_OSSPSReconstructedImage(const void* par);

//void STIR_displayImageData(void* ptr);

//void STIR_getImageDimensions(void* ptr, void* dim);
//void STIR_getImageDimensions(void* ptr, int dim[3]);
//void STIR_getImageData(void* ptr, double data[]);

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

int cSTIR_setPoissonLogLikelihoodWithLinearModelForMeanAndProjDataParameter
(PoissonLogLikelihoodWithLinearModelForMeanAndProjData<Image3DF>& obj,
string name, string value)
{
	if (boost::iequals(name, "sensitivity_filename")) {
		obj.set_sensitivity_filename(value);
		return 0;
	}
	return 1;
}

	else if (boost::iequals(set, "PoissonLogLikelihoodWithLinearModelForMeanAndProjData")) {
		PoissonLogLikelihoodWithLinearModelForMeanAndProjData<Image3DF>* obj =
			(PoissonLogLikelihoodWithLinearModelForMeanAndProjData<Image3DF>*)h->data();
		status = cSTIR_setPoissonLogLikelihoodWithLinearModelForMeanAndProjDataParameter(*obj, name, value);
	}

	//size_t len;
	//const char* w;
	//w = ex.what();
	//len = strlen(w) + 1;
	//_error = new char[len];
	//memcpy(_error, w, len);
	//w = ex.file();
	//len = strlen(w) + 1;
	//_file = new char[len];
	//memcpy(_file, w, len);
	//_line = ex.line();

	//#include "nlist.h"
	//#include "stir_util.h"

#if 0
	// old stuff
	template<class Method>
	void* STIR_reconstruction(const char* parFile) {

		DataHandle* handle = new DataHandle;
		try {
			Method reconstruction_object(parFile);
			//cout << reconstruction_object.parameter_info() << endl;
			ptrImage3DF* image =
				new ptrImage3DF(reconstruction_object.get_initial_data_ptr());
			if (reconstruction_object.set_up(*image) != Succeeded::yes) {
				ExecutionStatus status("STIR_reconstruction::set_up failed",
					__FILE__, __LINE__);
				handle->set(0, &status);
				return (void*)handle;
			}
			if (reconstruction_object.reconstruct(*image) != Succeeded::yes) {
				ExecutionStatus status("STIR_reconstruction::reconstruct failed",
					__FILE__, __LINE__);
				handle->set(0, &status);
				return (void*)handle;
			}
			handle->set((void*)image, 0, GRAB);
		}
		catch (StirException& se) {
			ExecutionStatus status(se);
			handle->set(0, &status);
		}
		return (void*)handle;
	}

	static void make_par_file(const void* par, char* filename, int len) {
		DataHandle* hp = (DataHandle*)par;
		const NestedList* rp = (const NestedList*)hp->data();
		time_t seconds;
		seconds = time(NULL);
		sprintf_s(filename, len, "tmp%d.par", seconds);
		//cout << filename << endl;
		ofstream out;
		out.open(filename, ios::out);
		TextWriter tw(&out);
		//TextWriterHandle twh;
		aTextWriter* ptr_tw = TextWriterHandle::writer;
		//aTextWriter* ptr_tw = twh.writer;
		//twh.writer = &tw;
		TextWriterHandle::writer = &tw;
		printNestedList(*rp);
		out.close();
		TextWriterHandle::writer = ptr_tw;
		//twh.writer = ptr_tw;
	}

	extern "C"
		void* STIR_reconstructedImage(const char* method, const void* par) {
		char filename[64];
		make_par_file(par, filename, 64);
		void* image;
		if (boost::iequals(method, "OSMAPOSL"))
			image = STIR_reconstruction<OSMAPOSLReconstruction<Image3DF>>(filename);
		else if (boost::iequals(method, "OSSPS"))
			image = STIR_reconstruction<OSSPSReconstruction<Image3DF>>(filename);
		else {
			DataHandle* h = new DataHandle;
			string error = "unknown method ";
			error += method;
			ExecutionStatus status(error.c_str(), __FILE__, __LINE__);
			h->set(0, &status);
			image = (void*)h;
		}
		remove(filename);
		return image;
	}

	extern "C"
		void* STIR_imageFromFile(const char* filename) {
		DataHandle* handle = new DataHandle;
		try {
			ptrImage3DF* ptr = new ptrImage3DF(read_from_file<Image3DF>(filename));
			handle->set((void*)ptr, 0, GRAB);
			return (void*)handle;
		}
		catch (StirException& se) {
			ExecutionStatus status(se);
			handle->set(0, &status);
			return (void*)handle;
		}
		catch (...) {
			cout << "unhandled exception occured" << endl;
			ExecutionStatus status("unhandled exception", __FILE__, __LINE__);
			handle->set(0, &status);
			return (void*)handle;
		}
	}

	extern "C"
		void STIR_displayImageData(void* ptr) {
		DataHandle* h = (DataHandle*)ptr;
		ptrImage3DF* ptrImage = (ptrImage3DF*)h->data();
		ptrImage3DF& image = *ptrImage;
		Coordinate3D<int> min_indices;
		Coordinate3D<int> max_indices;
		if (!image->get_regular_range(min_indices, max_indices)) {
			cout << "image is not regular" << endl;
		}
		else {
			cout << "voxels indices ranges:" << endl;
			for (int i = 1; i <= 3; i++)
				cout << min_indices[i] << " to " << max_indices[i] << endl;
		}
	}

	extern "C"
		//void STIR_getImageDimensions(void* ptr, void* pd) {
		//void STIR_getImageDimensions(void* ptr, int dim[3]) {
		void STIR_getImageDimensions(void* ptr, size_t pd) {
		DataHandle* h = (DataHandle*)ptr;
		ptrImage3DF* ptrImage = (ptrImage3DF*)h->data();
		ptrImage3DF& image = *ptrImage;
		Coordinate3D<int> min_indices;
		Coordinate3D<int> max_indices;
		int* dim = (int*)pd;
		if (!image->get_regular_range(min_indices, max_indices)) {
			cout << "image is not regular" << endl;
			dim[0] = 0;
			dim[1] = 0;
			dim[2] = 0;
		}
		else {
			for (int i = 0; i < 3; i++)
				dim[i] = max_indices[i + 1] - min_indices[i + 1] + 1;
		}
	}

	extern "C"
		//void STIR_getImageData(void* ptr, double data[]) {
		void STIR_getImageData(void* ptr, size_t pd) {
		DataHandle* h = (DataHandle*)ptr;
		ptrImage3DF* ptrImage = (ptrImage3DF*)h->data();
		ptrImage3DF& image = *ptrImage;
		Coordinate3D<int> min_indices;
		Coordinate3D<int> max_indices;
		double* data = (double*)pd;
		if (!image->get_regular_range(min_indices, max_indices)) {
			cout << "image is not regular" << endl;
			return;
		}
		else {
			for (int z = min_indices[1], i = 0; z <= max_indices[1]; z++) {
				for (int y = min_indices[2]; y <= max_indices[2]; y++) {
					for (int x = min_indices[3]; x <= max_indices[3]; x++, i++) {
						data[i] = (*image)[z][y][x];
					}
				}
			}
		}
	}

	//extern "C"
	//double* STIR_imageData(void* ptr) {
	//	int dim[3];
	//	STIR_getImageDimensions(ptr, (size_t)dim);
	//	double* data = new double[dim[0] * dim[1] * dim[2]];
	//	STIR_getImageData(ptr, (size_t)data);
	//	return data;
	//}

	extern "C"
		void STIR_deleteImage(void* ptr) {
		DataHandle* handle = (DataHandle*)ptr;
		//delete (ptrImage3DF*)handle->data();
		delete handle;
	}

	extern "C"
		void* STIR_imagesDifference
		(void* first, void* second, int rimsize) {

		DataHandle* handle = new DataHandle;
		try {

			DataHandle* first_h = (DataHandle*)first;
			DataHandle* second_h = (DataHandle*)second;

			ptrImage3DF* first_ptr = (ptrImage3DF*)first_h->data();
			ptrImage3DF& first_operand = *first_ptr;
			ptrImage3DF* second_ptr = (ptrImage3DF*)second_h->data();
			ptrImage3DF& second_operand = *second_ptr;

			{
				std::string explanation;
				if (!first_operand->has_same_characteristics(*second_operand,
					explanation))
				{
					warning("input images do not have the same characteristics.\n%s",
						explanation.c_str());
					ExecutionStatus status(
						"input images do not have the same characteristics",
						__FILE__, __LINE__);
					handle->set(0, &status);
					return (void*)handle;
				}
			}

			if (rimsize >= 0)
			{
				truncate_rim(*first_operand, rimsize);
				truncate_rim(*second_operand, rimsize);
			}

			float reference_max = first_operand->find_max();
			float reference_min = first_operand->find_min();

			float amplitude = fabs(reference_max) > fabs(reference_min) ?
				fabs(reference_max) : fabs(reference_min);

			*first_operand -= *second_operand;
			const float max_error = first_operand->find_max();
			const float min_error = first_operand->find_min();

			float max_abs_error = fabs(min_error);
			if (max_error > max_abs_error)
				max_abs_error = max_error;
			*first_operand += *second_operand;

			double* result = (double*)malloc(sizeof(double));
			*result = max_abs_error / amplitude;
			handle->set(result, 0, GRAB);
		}
		catch (StirException& se) {
			ExecutionStatus status(se);
			handle->set(0, &status);
			return (void*)handle;
		}
		catch (exception& e) {
			cout << "exception " << e.what() << endl;
		}
		catch (...) {
			cout << "unhandled exception occured" << endl;
			ExecutionStatus status("unhandled exception", __FILE__, __LINE__);
			handle->set(0, &status);
			return (void*)handle;
		}

		return (void*)handle;
	}

#endif

//Reconstruction<Image3DF>* recon = (Reconstruction<Image3DF>*)hr->data();

//s = setup<OSMAPOSLReconstruction<Image3DF>>((void*)recon, *image);
//s = ((STIR_Reconstruction<OSMAPOSLReconstruction<Image3DF>>*)recon)->setup(*image);

//#include <conio.h>
//#include <stdio.h>
//#include <stdlib.h>
//#include <string.h>
//#include <time.h>
//
//#include <boost\algorithm\string.hpp>
//
//#include "stir/DiscretisedDensity.h"
//#include "stir/ArrayFunction.h"
//#include "stir/CartesianCoordinate3D.h"
//#include "stir/CPUTimer.h"
//#include "stir/DataProcessor.h"
//#include "stir/HighResWallClockTimer.h"
//#include "stir/IndexRange3D.h"
//#include "stir/is_null_ptr.h"
//#include "stir/KeyParser.h"
//#include "stir/recon_array_functions.h"
//#include "stir/Succeeded.h"
//#include "stir/utilities.h"
//#include "stir/VoxelsOnCartesianGrid.h"
//#include "stir/IO/OutputFileFormat.h"
//#include "stir/IO/read_from_file.h"
//#include "stir/OSMAPOSL/OSMAPOSLReconstruction.h"
//#include "stir/OSSPS/OSSPSReconstruction.h"
//#include "stir/recon_buildblock/distributable_main.h"
//#include "stir/recon_buildblock/PoissonLogLikelihoodWithLinearModelForMeanAndProjData.h"
//#include "stir/Shape/Shape3D.h"
//
//#include "StirException.h"
//#include "TextWriter.h"
//
//#include <fstream>
//#include <iostream>
//#include <numeric>
//#include <string>
////using namespace std;

//#define GRAB 1
//
//USING_NAMESPACE_STIR

//typedef DiscretisedDensity<3, float> Image3DF;
//typedef boost::shared_ptr<Image3DF> ptrImage3DF;

	//string error = "parameter ";
	//error += name;
	//error += " not found";
	//ExecutionStatus status(error.c_str(), __FILE__, __LINE__);
	//h->set(0, &status);
	//return (void*)h;

	//string error = "unknown method ";
	//error += set;
	//ExecutionStatus status(error.c_str(), __FILE__, __LINE__);
	//h->set(0, &status);
	//return (void*)h;

	//recon = cSTIR_newReconstructionMethod<OSMAPOSLReconstruction<Image3DF>>

	//DataHandle* h = new DataHandle;
	//string error = "unknown method ";
	//error += method;
	//ExecutionStatus status(error.c_str(), __FILE__, __LINE__);
	//h->set(0, &status);
	//recon = (void*)h;

	//if (recon->set_up(*image) != Succeeded::yes) {

	//string error = "unknown method ";
	//error += method;
	//ExecutionStatus status(error.c_str(), __FILE__, __LINE__);
	//handle->set(0, &status);
	//return (void*)handle;

	//string error = "unknown method ";
	//error += method;
	//ExecutionStatus status(error.c_str(), __FILE__, __LINE__);
	//handle->set(0, &status);
	//return (void*)handle;

	//string error = "unknown method ";
	//error += method;
	//ExecutionStatus status(error.c_str(), __FILE__, __LINE__);
	//handle->set(0, &status);
	//return (void*)handle;

	//string error = "unknown method ";
	//error += method;
	//ExecutionStatus status(error.c_str(), __FILE__, __LINE__);
	//handle->set(0, &status);
	//return (void*)handle;

//template class xSTIR_Reconstruction < OSMAPOSLReconstruction<Image3DF> > ;

//int* ptr = new int;
//*ptr = recon.get_num_subsets();
//h->set((void*)ptr, 0, GRAB);
//return (void*)h;

//int* ptr = new int;
//*ptr = recon.get_num_subiterations();
//h->set((void*)ptr, 0, GRAB);
//return (void*)h;

//int* ptr = new int;
//*ptr = recon.get_start_subiteration_num();
//h->set((void*)ptr, 0, GRAB);
//return (void*)h;

//int n;

////char* value = (char*)ptr;
////sscanf_s(value, "%d", &n);
//n = atoi((char*)ptr);
//recon.set_num_subsets(n);

//return 0;

//char* value = (char*)ptr;
//sscanf_s(value, "%d", &n);
//recon.set_num_subiterations(n);

//char* value = (char*)ptr;
//sscanf_s(value, "%d", &n);
//recon.set_start_subiteration_num(n);

//char* value = (char*)ptr;
//obj.set_sensitivity_filename(value);

//char* value = (char*)ptr;
//if (boost::iequals(value, "true"))
//	obj.set_recompute_sensitivity(true);
//else
//	obj.set_recompute_sensitivity(false);
//return 0;

//(OSMAPOSLReconstruction<Image3DF>& recon, const char* name, const void* ptr)

//(OSMAPOSLReconstruction<Image3DF>& recon, string name)

//OSMAPOSLReconstruction<Image3DF>* recon =
//	(OSMAPOSLReconstruction<Image3DF>*)h->data();
//status = cSTIR_setOSMAPOSLParameter(*recon, name, value);

//OSMAPOSLReconstruction<Image3DF>* recon =
//	(OSMAPOSLReconstruction<Image3DF>*)h->data();
//return cSTIR_OSMAPOSLParameter(*recon, name);

void* cSTIR_intHandle(int i) {
	DataHandle* h = new DataHandle;
	int* ptr = new int;
	*ptr = i;
	h->set((void*)ptr, 0, GRAB);
	return (void*)h;
}

extern "C"
char* cSTIR_newInt(int i) {
	int* ptr = new int;
	*ptr = i;
	return (char*)ptr;
}

extern "C"
void cSTIR_deleteInt(char* ptr) {
	delete (int*)ptr;
}

extern "C"
void* cSTIR_charToVoid(char* ptr) {
	return (void*)ptr;
}

//recon->set_num_subsets(atoi((char*)hv->data()));
//recon->set_num_subsets(atoi((char*)ptr_v));

//recon->set_num_subiterations(atoi((char*)ptr_v));

//	recon->set_start_subiteration_num(atoi((char*)ptr_v));

//if (boost::iequals(name, "num_subsets"))
//	recon->set_num_subsets(atoi((char*)ptr_v));
//else if (boost::iequals(name, "num_subiterations"))
//	recon->set_num_subiterations(atoi((char*)ptr_v));
//else if (boost::iequals(name, "start_subiteration_num"))
//	recon->set_start_subiteration_num(atoi((char*)ptr_v));
//else
//	not_found = 1;

//if (boost::iequals(name, "num_subsets"))
//	return cSTIR_intHandle(recon->get_num_subsets());
//if (boost::iequals(name, "num_subiterations"))
//	return cSTIR_intHandle(recon->get_num_subiterations());
//if (boost::iequals(name, "start_subiteration_num"))
//	return cSTIR_intHandle(recon->get_start_subiteration_num());
//if (boost::iequals(name, "objective_function")) {
//	DataHandle* h = new DataHandle;
//	void *ptr = (void*)&recon->get_objective_function();
//	h->set((void*)ptr, 0);
//	return (void*)h;
//}

//if (boost::iequals(name, "num_subsets"))
//	recon->set_num_subsets(atoi((char*)ptr_v));
//else if (boost::iequals(name, "num_subiterations"))
//	recon->set_num_subiterations(atoi((char*)ptr_v));
//else if (boost::iequals(name, "start_subiteration_num"))
//	recon->set_start_subiteration_num(atoi((char*)ptr_v));
//else
//	not_found = 1;

//if (boost::iequals(name, "num_subsets"))
//	return cSTIR_intHandle(recon->get_num_subsets());
//if (boost::iequals(name, "num_subiterations"))
//	return cSTIR_intHandle(recon->get_num_subiterations());
//if (boost::iequals(name, "start_subiteration_num"))
//	return cSTIR_intHandle(recon->get_start_subiteration_num());
//if (boost::iequals(name, "objective_function")) {
//	DataHandle* h = new DataHandle;
//	void *ptr = (void*)&recon->get_objective_function();
//	h->set((void*)ptr, 0);
//	return (void*)h;
//}

//cout << (char*)h->data() << endl;

//cout << s << endl;

//obj.set_recompute_sensitivity(boost::iequals((char*)ptr, "true"));

//DataHandle* hv = (DataHandle*)value;
//cout << (char*)hv->data() << endl;

extern "C"
void* cSTIR_reconstructionSetup(void* ptr_r, void* ptr_i) {
	DataHandle* hr = (DataHandle*)ptr_r;
	DataHandle* hi = (DataHandle*)ptr_i;
	Reconstruction<Image3DF>* recon = (Reconstruction<Image3DF>*)hr->data();
	ptrImage3DF* image = (ptrImage3DF*)hi->data();
	DataHandle* handle = new DataHandle;
	try {
		Succeeded s = recon->set_up(*image);
		if (s != Succeeded::yes) {
			ExecutionStatus status("cSTIR_reconstructionSetup failed",
				__FILE__, __LINE__);
			handle->set(0, &status);
			return (void*)handle;
		}
	}
	catch (StirException& se) {
		ExecutionStatus status(se);
		handle->set(0, &status);
	}
	return (void*)handle;
}

//cout << (char*)h->data() << endl;

void printCharData(const void* ptr) {
	DataHandle* hv = (DataHandle*)ptr;
	cout << (char*)hv->data() << endl;
}

void* cSTIR_intHandle(int i) {
	DataHandle* h = new DataHandle;
	int* ptr = new int;
	*ptr = i;
	h->set((void*)ptr, 0, GRAB);
	return (void*)h;
}

//recon->update_estimate(**image);
//recon->end_of_iteration_processing(**image);

//xSTIR_subiteration<OSMAPOSLReconstruction<Image3DF>>(hr->data())++;

//xSTIR_subiteration<OSSPSReconstruction<Image3DF>>(hr->data())++;

//(PoissonLogLikelihoodWithLinearModelForMean<Image3DF>& obj,

//PoissonLogLikelihoodWithLinearModelForMean<Image3DF>* obj =
//	(PoissonLogLikelihoodWithLinearModelForMean<Image3DF>*)h->data();
//status = cSTIR_setPoissonLogLikelihoodWithLinearModelForMeanParameter(*obj, name, value);

//(void* ptr_r, const char* name, const void* ptr_v)

//DataHandle* hr = (DataHandle*)ptr_r;
//DataHandle* hv = (DataHandle*)ptr_v;

//(void* ptr_r, const char* name, const void* ptr_v)

//DataHandle* hr = (DataHandle*)ptr_r;
//DataHandle* hv = (DataHandle*)ptr_v;

//cSTIR_IterativeReconstructionParameter(void* ptr_r, string name)

//DataHandle* h = (DataHandle*)ptr_r;

//int cSTIR_setOSMAPOSLParameter(void* ptr_r, const char* name, const void* ptr_v)

//DataHandle* h = (DataHandle*)ptr_r;

//void* cSTIR_OSMAPOSLParameter(void* ptr_r, string name)

//(void* ptr_r, const char* name, const void* ptr_v)

//(void* ptr_r, string name)

//(void* ptr_of, const char* name, const void* ptr)

//DataHandle* ho = (DataHandle*)ptr_of;

//DataHandle* h = (DataHandle*)ptr;

//OSMAPOSLReconstruction<Image3DF>* recon =
//	(OSMAPOSLReconstruction<Image3DF>*)hr->data();

//OSMAPOSLReconstruction<Image3DF>* recon =
//	(OSMAPOSLReconstruction<Image3DF>*)hr->data();

//OSMAPOSLReconstruction<Image3DF>* recon =
//	(OSMAPOSLReconstruction<Image3DF>*)handle->data();

//DataHandle* h = new DataHandle;
//int* ptr = new int;
//*ptr = i;
//h->set((void*)ptr, 0, GRAB);
//return (void*)h;

//if (!ptr)
//	return 0;
//DataHandle* ptr_h = (DataHandle*)ptr;
//void* ptr_d = ptr_h->data();
//if (!ptr_d)
//	return 0;
//else
//	return *((int*)ptr_d);

//if (!ptr)
//	return 0;
//DataHandle* ptr_h = (DataHandle*)ptr;
//void* ptr_d = ptr_h->data();
//if (!ptr_d)
//	return 0;
//else
//	return *((double*)ptr_d);

//static void*
//cSTIR_unknownMethod(const char* set, const char* file, int line)
//{
//	DataHandle* handle = new DataHandle;
//	string error = "unknown method ";
//	error += set;
//	ExecutionStatus status(error.c_str(), file, line);
//	handle->set(0, &status);
//	return (void*)handle;
//}
//
//static void*
//unknownParameterSet(const char* set, const char* file, int line)
//{
//	DataHandle* h = new DataHandle;
//	string error = "unknown set of parameters ";
//	error += set;
//	ExecutionStatus status(error.c_str(), file, line);
//	h->set(0, &status);
//	return (void*)h;
//}

//return unknownParameterSet(set, __FILE__, __LINE__);

//recon = cSTIR_unknownMethod(method, __FILE__, __LINE__);

//return cSTIR_unknownMethod(method, __FILE__, __LINE__);

//return cSTIR_unknownMethod(method, __FILE__, __LINE__);

//QuadraticPrior < float >* p 

//*pp = 
//	(boost::shared_ptr<GeneralisedPrior<Image3DF> >)new QuadraticPrior < float >();

//else if (boost::iequals(name, "prior")) {
//	void* ptr = hv->data();
//	obj_fun->set_prior_sptr((boost::shared_ptr<stir::GeneralisedPrior<Image3DF>> &)ptr);
//	//obj_fun->set_prior_sptr((GeneralisedPrior<Image3DF>*)hv->data());
//}

//GeneralisedPrior<Image3DF>* prior = (GeneralisedPrior<Image3DF>*)hp->data();

//GeneralisedPrior<Image3DF>* prior = (GeneralisedPrior<Image3DF>*)handle->data();

//boost::shared_ptr<ProjMatrixByBinUsingRayTracing>* pp =
//	new boost::shared_ptr<ProjMatrixByBinUsingRayTracing>
//	(new ProjMatrixByBinUsingRayTracing());
//handle->set((void*)pp, 0, GRAB);
//return handle;

//boost::shared_ptr<GeneralisedPrior<Image3DF> >* pp = 
//	new boost::shared_ptr<GeneralisedPrior<Image3DF> >
//	(new QuadraticPrior < float >());
//handle->set((void*)pp, 0, GRAB);
//return handle;

extern "C"
void* cSTIR_newMatrixProjector(const char* matrix_type)
{
	DataHandle* handle = new DataHandle;
	if (boost::iequals(matrix_type, "RayTracingMatrix"))
		return newObjectHandle<ProjMatrixByBinUsingRayTracing>();
	else
		return unknownObject("matrix projector", matrix_type, __FILE__, __LINE__);
}

extern "C"
void* cSTIR_newPrior(const char* prior_type)
{
	DataHandle* handle = new DataHandle;
	if (boost::iequals(prior_type, "QuadraticPrior"))
		return newObjectHandle<GeneralisedPrior<Image3DF>, QuadraticPrior<float>>();
	else
		return unknownObject("prior", prior_type, __FILE__, __LINE__);
}

void* cSTIR_newMatrixProjector(const char* matrix_type);
void* cSTIR_newPrior(const char* prior_type);

//boost::shared_ptr<ProjMatrixByBin>* pp =
//	(boost::shared_ptr<ProjMatrixByBin>*)hv->data();
//ProjMatrixByBinUsingRayTracing* matrix = (ProjMatrixByBinUsingRayTracing*)pp->get();

ProjMatrixByBinUsingRayTracing* matrix = (ProjMatrixByBinUsingRayTracing*)sptr.get();
printf("%p %p %p\n", hv, hv->data(), matrix);
cout << matrix->get_num_tangential_LORs() << endl;

//projectors->set_proj_matrix_ptr(*pp);

printf("%p %p %p\n", handle, handle->data(), matrix);

//boost::shared_ptr<ProjMatrixByBinUsingRayTracing> matrix =
//	*(boost::shared_ptr<ProjMatrixByBinUsingRayTracing>*)hp->data();

//prior->set_penalisation_factor(*(float*)hv->data());

//string s((char*)hv->data());
//obj_fun->set_sensitivity_filename(s);

//obj_fun->set_sensitivity_filename(charDataFromHandle((void*)hv));

//string s((char*)hv->data());

//string s((char*)hv->data());

//cout << recon->get_objective_function().get_prior_ptr()->get_penalisation_factor() << endl;

//cout << sptr->get()->get_prior_ptr()->get_penalisation_factor() << endl;
//cout << recon->get_objective_function().get_prior_ptr()->get_penalisation_factor() << endl;
//TextWriterHandle h;
//h.write(recon->parameter_info().c_str());

//GeneralisedObjectiveFunction<Image3DF>* obj_fun =
//	(GeneralisedObjectiveFunction<Image3DF>*)ho->data();

//printf("%p\n", (void*)pp);

//cout << pp->get()->get_penalisation_factor() << endl;
//cout << obj_fun->get_prior_ptr()->get_penalisation_factor() << endl;

//PoissonLogLikelihoodWithLinearModelForMean<Image3DF>* obj_fun =
//	(PoissonLogLikelihoodWithLinearModelForMean<Image3DF>*)ho->data();

//cout << (char*)hv->data() << endl;

//cout << s << endl;

//void *ptr = (void*)&recon->get_objective_function();
//h->set((void*)ptr, 0);

//boost::shared_ptr<ProjectorByBinPairUsingProjMatrixByBin> sptr_p =
//	*(boost::shared_ptr<ProjectorByBinPairUsingProjMatrixByBin>*)hp->data();
//ProjectorByBinPairUsingProjMatrixByBin* projectors = 
//	(ProjectorByBinPairUsingProjMatrixByBin*)sptr_p.get();
//ProjectorByBinPairUsingProjMatrixByBin* projectors =
//	objectFromHandle<ProjectorByBinPairUsingProjMatrixByBin>(hp);

//boost::shared_ptr<ProjMatrixByBin> sptr_m =
//	*(boost::shared_ptr<ProjMatrixByBin>*)hv->data();

//sptr_p->set_proj_matrix_ptr(sptr_m);
//projectors->set_proj_matrix_ptr(sptr_m);

//boost::shared_ptr<GeneralisedObjectiveFunction<Image3DF>> sptr =
//	*(boost::shared_ptr<GeneralisedObjectiveFunction<Image3DF>>*)ho->data();
//GeneralisedObjectiveFunction<Image3DF>* obj_fun = sptr.get();

//boost::shared_ptr<GeneralisedPrior<Image3DF> >* pp =
//	(boost::shared_ptr<GeneralisedPrior<Image3DF> >*)hv->data();

//obj_fun->set_prior_sptr(*pp);

//boost::shared_ptr<GeneralisedObjectiveFunction<Image3DF>> sptr =
//	*(boost::shared_ptr<GeneralisedObjectiveFunction<Image3DF>>*)ho->data();
//xSTIR_ObjectiveFunction
//	<PoissonLogLikelihoodWithLinearModelForMean
//	<Image3DF >> *obj_fun =
//	(xSTIR_ObjectiveFunction
//	<PoissonLogLikelihoodWithLinearModelForMean
//	<Image3DF >> *) sptr.get();

//boost::shared_ptr<GeneralisedObjectiveFunction<Image3DF>> sptr =
//	*(boost::shared_ptr<GeneralisedObjectiveFunction<Image3DF>>*)ho->data();
//xSTIR_ObjectiveFunction
//	<PoissonLogLikelihoodWithLinearModelForMeanAndProjData
//	<Image3DF>>* obj_fun =
//	(xSTIR_ObjectiveFunction
//	<PoissonLogLikelihoodWithLinearModelForMeanAndProjData
//	<Image3DF>>*) sptr.get();

//boost::shared_ptr<ProjectorByBinPairUsingProjMatrixByBin> sptr_p =
//	*(boost::shared_ptr<ProjectorByBinPairUsingProjMatrixByBin>*)hv->data();

//obj_fun->set_projector_pair_sptr(sptr_p);

//xSTIR_ObjectiveFunction<
//> 

//<xSTIR_Reconstruction<OSMAPOSLReconstruction<Image3DF>>>

//<xSTIR_Reconstruction<OSSPSReconstruction<Image3DF>>>

//void* cSTIR_setupObjective(const char* obj, void* ptr_obj)

//if (boost::iequals
//	(obj, "PoissonLogLikelihoodWithLinearModelForMeanAndProjData"))
//	failure = xSTIR_setupObjective
//	<PoissonLogLikelihoodWithLinearModelForMeanAndProjData<Image3DF>>
//	(h->data());
//else
//	return unknownObject("objective function", obj, __FILE__, __LINE__);

//void* cSTIR_setupReconstruction(const char* method, void* ptr_r, void* ptr_i)

//if (boost::iequals(method, "OSMAPOSL"))
//	s = xSTIR_setupReconstruction<OSMAPOSLReconstruction<Image3DF>>
//	(hr->data(), *image);
//else if (boost::iequals(method, "OSSPS"))
//	s = xSTIR_setupReconstruction<OSSPSReconstruction<Image3DF>>
//	(hr->data(), *image);
//else
//	return unknownObject("method", method, __FILE__, __LINE__);

//void* cSTIR_setSubiteration(const char* method, void* ptr_r, int i)

//if (boost::iequals(method, "OSMAPOSL"))
//	xSTIR_subiteration<OSMAPOSLReconstruction<Image3DF>>(hr->data()) = i;
//else if (boost::iequals(method, "OSSPS"))
//	xSTIR_subiteration<OSSPSReconstruction<Image3DF>>(hr->data()) = i;
//else
//	return unknownObject("method", method, __FILE__, __LINE__);

//void* cSTIR_getSubiteration(const char* method, void* ptr_r)

//if (boost::iequals(method, "OSMAPOSL"))
//	iter = xSTIR_subiteration<OSMAPOSLReconstruction<Image3DF>>(hr->data());
//else if (boost::iequals(method, "OSSPS"))
//	iter = xSTIR_subiteration<OSSPSReconstruction<Image3DF>>(hr->data());
//else
//	return unknownObject("method", method, __FILE__, __LINE__);

//void* cSTIR_update(const char* method, void* ptr_r, void* ptr_i)

//if (boost::iequals(method, "OSMAPOSL"))
//	xSTIR_updateReconstruction<OSMAPOSLReconstruction<Image3DF>>
//	(hr->data(), **image);
//else if (boost::iequals(method, "OSSPS"))
//	xSTIR_updateReconstruction<OSSPSReconstruction<Image3DF>>
//	(hr->data(), **image);
//else
//	return unknownObject("method", method, __FILE__, __LINE__);

//boost::shared_ptr<ProjMatrixByBin> sptr =
//	*(boost::shared_ptr<ProjMatrixByBin>*)hp->data();
//ProjMatrixByBinUsingRayTracing* matrix = (ProjMatrixByBinUsingRayTracing*)sptr.get();

//boost::shared_ptr<ProjMatrixByBin> sptr =
//	*(boost::shared_ptr<ProjMatrixByBin>*)handle->data();
//ProjMatrixByBinUsingRayTracing* matrix = (ProjMatrixByBinUsingRayTracing*)sptr.get();

//void* cSTIR_setupObjective(const char* obj, void* ptr_obj);

//void* cSTIR_setupReconstruction(const char* method, void* ptr_r, void* ptr_i);

//void* cSTIR_update(const char* method, void* ptr_r, void* ptr_i);

//void* cSTIR_setupReconstruction(const char* method, void* ptr_r, void* ptr_i);

//xSTIR_Reconstruction<OSMAPOSLReconstruction<Image3DF>>* recon =
//	(xSTIR_Reconstruction<OSMAPOSLReconstruction<Image3DF>>*)hr->data();

//else if (boost::iequals(name, "initial_estimate"))
//	xSTIR_set_initial_estimate_file<IterativeReconstruction<Image3DF>>
//	(hr->data(), (char*)hv->data());
//	//recon->set_initial_estimate_file((char*)hv->data());

//try {
//	xSTIR_subiteration<IterativeReconstruction<Image3DF>>(hr->data()) = iter;
//}
//catch (StirException& se) {
//	ExecutionStatus status(se);
//	handle->set(0, &status);
//}

//DataHandle* handle = new DataHandle;
//int iter;
//try {
//	iter = xSTIR_subiteration<IterativeReconstruction<Image3DF>>(hr->data());
//	int* i = new int;
//	*i = iter;
//	handle->set((void*)i, 0, GRAB);
//}
//catch (StirException& se) {
//	ExecutionStatus status(se);
//	handle->set(0, &status);
//}
//return (void*)handle;

//template < class R >
//class xSTIR_IterativeReconstruction : public R {
//public:
//	xSTIR_IterativeReconstruction() : R() {}
//	xSTIR_IterativeReconstruction(const char* parFile) : R(parFile) {}
//	Succeeded setup(ptrImage3DF const& image) {
//		return set_up(image);
//	}
//	void update(Image3DF &image) {
//		update_estimate(image);
//		end_of_iteration_processing(image);
//		subiteration_num++;
//	}
//	int& subiteration() {
//		return subiteration_num;
//	}
//	int subiteration() const {
//		return subiteration_num;
//	}
//	void set_initial_estimate_file(const char* filename) {
//		initial_data_filename = filename;
//	}
//};

//xSTIR_IterativeReconstruction3DF() : IterativeReconstruction<Image3DF>() {}
//xSTIR_IterativeReconstruction3DF(const char* parFile) : 
//	IterativeReconstruction<Image3DF>(parFile) {}

//DataHandle* handle = (DataHandle*)ptr;

//template < class R >
//Succeeded xSTIR_setupReconstruction(void* ptr, ptrImage3DF const& image) {
//	xSTIR_IterativeReconstruction<R>* recon = (xSTIR_IterativeReconstruction<R>*)ptr;
//	Succeeded s = recon->setup(image);
//	recon->subiteration() = recon->get_start_subiteration_num();
//	return s;
//}

//template < class R >
//void xSTIR_updateReconstruction(void* ptr, Image3DF& image) {
//	xSTIR_IterativeReconstruction<R>* recon = (xSTIR_IterativeReconstruction<R>*)ptr;
//	recon->update(image);
//}

//template < class R >
//int& xSTIR_subiteration(void* ptr) {
//	xSTIR_IterativeReconstruction<R>* recon = (xSTIR_IterativeReconstruction<R>*)ptr;
//	//std::cout << recon->subiteration() << endl;
//	return recon->subiteration();
//}

//std::cout << recon->subiteration() << endl;

//template < class R >
//void xSTIR_set_initial_estimate_file(void* ptr, const char* filename) {
//	xSTIR_IterativeReconstruction<R>* recon = (xSTIR_IterativeReconstruction<R>*)ptr;
//	recon->set_initial_estimate_file(filename);
//}

//DataHandle* handle = new DataHandle;

//return (void*)handle;

//bool post_process() {
//	return post_processing();
//}

//template < class F >
//bool xSTIR_setupObjective(void* ptr) 
//{
//	boost::shared_ptr< GeneralisedObjectiveFunction<Image3DF> >* sptr_obj =
//		(boost::shared_ptr< GeneralisedObjectiveFunction<Image3DF> >*)ptr;
//	xSTIR_ObjectiveFunction<F>* obj = 
//		(xSTIR_ObjectiveFunction<F>*)sptr_obj->get();
//	bool status = obj->post_process();
//	return status;
//}

//#include "ci.h"

//bool failure = xSTIR_setupObjective<GeneralisedObjectiveFunction<Image3DF> >
//	(h->data());
//if (failure) {

//template < class F >
//class xSTIR_ObjectiveFunction : public F {
//public:
//	void set_input_file(const char* filename) {
//		this->input_filename = filename;
//	}
//};

//xSTIR_ObjectiveFunction<
//	PoissonLogLikelihoodWithLinearModelForMean<Image3DF> >
//	*obj_fun = objectFromHandle<GeneralisedObjectiveFunction<Image3DF>,
//	xSTIR_ObjectiveFunction<
//	PoissonLogLikelihoodWithLinearModelForMean<Image3DF> > >
//	(ho);

//xSTIR_ObjectiveFunction<
//	PoissonLogLikelihoodWithLinearModelForMeanAndProjData<Image3DF> >*
//	obj_fun1 =
//	objectFromHandle<GeneralisedObjectiveFunction<Image3DF>,
//	xSTIR_ObjectiveFunction<
//	PoissonLogLikelihoodWithLinearModelForMeanAndProjData<Image3DF> > >
//	(ho);

//boost::shared_ptr< GeneralisedObjectiveFunction<Image3DF> >* sptr_obj =
//	(boost::shared_ptr< GeneralisedObjectiveFunction<Image3DF> >*)ptr;
//xSTIR_GeneralisedObjectiveFunction3DF* obj =
//	(xSTIR_GeneralisedObjectiveFunction3DF*)sptr_obj->get();

//xSTIR_IterativeReconstruction3DF* recon =
//	(xSTIR_IterativeReconstruction3DF*)ptr;

//xSTIR_IterativeReconstruction3DF* recon = 
//	(xSTIR_IterativeReconstruction3DF*)ptr;

//xSTIR_IterativeReconstruction3DF* recon =
//	(xSTIR_IterativeReconstruction3DF*)ptr;

//xSTIR_IterativeReconstruction3DF* recon = 
//	(xSTIR_IterativeReconstruction3DF*)ptr;

//boost::shared_ptr<Base>* sptr = (boost::shared_ptr<Base>*)handle->data();

//boost::shared_ptr<ProjectorByBinPairUsingProjMatrixByBin>* sptr_p =
//	(boost::shared_ptr<ProjectorByBinPairUsingProjMatrixByBin>*)hp->data();

////boost::shared_ptr<TruncateToCylindricalFOVImageProcessor<float> > sptr_p =
////	*(boost::shared_ptr<TruncateToCylindricalFOVImageProcessor<float> >*)
////	hp->data();

//boost::shared_ptr<GeneralisedPrior<Image3DF> >* prior =
//	(boost::shared_ptr<GeneralisedPrior<Image3DF> >*)hp->data();

//boost::shared_ptr<GeneralisedPrior<Image3DF> >* prior =
//	(boost::shared_ptr<GeneralisedPrior<Image3DF> >*)handle->data();

//boost::shared_ptr<GeneralisedObjectiveFunction<Image3DF> >* sptr =
//	(boost::shared_ptr<GeneralisedObjectiveFunction<Image3DF> >*)ho->data();

//Reconstruction<Image3DF>* recon = (Reconstruction<Image3DF>*)hr->data();

//IterativeReconstruction<Image3DF>* recon =
//	(IterativeReconstruction<Image3DF>*)hr->data();

//IterativeReconstruction<Image3DF>* recon =
//	(IterativeReconstruction<Image3DF>*)handle->data();

//OSMAPOSLReconstruction<Image3DF>* recon =
//	(OSMAPOSLReconstruction<Image3DF>*)hr->data();

//OSMAPOSLReconstruction<Image3DF>* recon =
//	(OSMAPOSLReconstruction<Image3DF>*)handle->data();

//OSSPSReconstruction<Image3DF>* recon =
//	(OSSPSReconstruction<Image3DF>*)hr->data();

//OSSPSReconstruction<Image3DF>* recon =
//	(OSSPSReconstruction<Image3DF>*)handle->data();

//DataHandle* hs = (DataHandle*)ptr_s;
//DataHandle* hv = (DataHandle*)ptr_v;

//DataHandle* handle = (DataHandle*)ptr;

//DataHandle* h = (DataHandle*)ptr_obj;

//DataHandle* hr = (DataHandle*)ptr_r;
//DataHandle* hi = (DataHandle*)ptr_i;
//ptrImage3DF* image = (ptrImage3DF*)hi->data();

//DataHandle* hr = (DataHandle*)ptr_r;
//DataHandle* hi = (DataHandle*)ptr_i;
//Reconstruction<Image3DF>* recon = (Reconstruction<Image3DF>*)hr->data();
//ptrImage3DF* image = (ptrImage3DF*)hi->data();

//DataHandle* hr = (DataHandle*)ptr_r;

//DataHandle* hr = (DataHandle*)ptr_r;

//DataHandle* hr = (DataHandle*)ptr_r;
//DataHandle* hi = (DataHandle*)ptr_i;
//IterativeReconstruction<Image3DF>* recon = 
//	(IterativeReconstruction<Image3DF>*)hr->data();
//ptrImage3DF* image = (ptrImage3DF*)hi->data();

//DataHandle* h = (DataHandle*)ptr;
//ptrImage3DF* ptrImage = (ptrImage3DF*)h->data();

//DataHandle* h = (DataHandle*)ptr;
//ptrImage3DF* ptrImage = (ptrImage3DF*)h->data();

//DataHandle* first_h = (DataHandle*)first;
//DataHandle* second_h = (DataHandle*)second;

//ptrImage3DF* first_ptr = (ptrImage3DF*)first_h->data();
//ptrImage3DF* second_ptr = (ptrImage3DF*)second_h->data();

//string s((char*)hv->data());
//obj_fun->set_zero_seg0_end_planes(boost::iequals(s, "true"));

//string s((char*)hv->data());
//recon->set_output_filename_prefix(s);

//boost::shared_ptr<DataProcessor<Image3DF> >* sptr_p =
//	(boost::shared_ptr<DataProcessor<Image3DF> >*)hv->data();

//boost::shared_ptr<GeneralisedObjectiveFunction<Image3DF> >* sptr =
//	(boost::shared_ptr<GeneralisedObjectiveFunction<Image3DF> >*)hv->data();

//xSTIR_OSSPSReconstruction3DF* rec = (xSTIR_OSSPSReconstruction3DF*)recon;

//DataHandle* h = new DataHandle;
//boost::shared_ptr<GeneralisedObjectiveFunction<Image3DF> >* sptr =
//	new boost::shared_ptr<GeneralisedObjectiveFunction<Image3DF> > ;

//string s((char*)hv->data());
//recon->set_MAP_model(s);

//Object* object = (Object*)sptr->get();

//void* cSTIR_getSubiteration(void* ptr_r)

//return intDataHandle(xSTIR_subiteration(hr->data()));

//PoissonLogLikelihoodWithLinearModelForMeanAndProjData<Image3DF>*
//	obj_fun =
//	objectFromHandle<GeneralisedObjectiveFunction<Image3DF>,
//	PoissonLogLikelihoodWithLinearModelForMeanAndProjData<Image3DF> >
//	(hp);

//xSTIR_PoissonLogLikelihoodWithLinearModelForMeanAndProjData3DF*
//	obj_fun =
//	objectFromHandle<GeneralisedObjectiveFunction<Image3DF>,
//	xSTIR_PoissonLogLikelihoodWithLinearModelForMeanAndProjData3DF>
//	(hp);

//if (!handle)
//	handle = (void*)new DataHandle;
//return handle;
////return status;

//CAST_PTR(DataHandle, ho, ptr_obj);
//DataHandle* handle = new DataHandle;
//try {
//	bool status = xSTIR_setupObjectiveFunction(ho->data());
//	if (status) {
//		ExecutionStatus status("cSTIR_setupObjective failed", __FILE__, __LINE__);
//		handle->set(0, &status);
//		return (void*)handle;
//	}
//}
//catch (StirException& se) {
//	ExecutionStatus status(se);
//	handle->set(0, &status);
//}
//return (void*)handle;

//extern "C"
//void* cSTIR_setupPrior(void* ptr_obj)
//{
//	CAST_PTR(DataHandle, ho, ptr_obj);
//	DataHandle* handle = new DataHandle;
//	try {
//		bool status = xSTIR_setupPrior(ho->data());
//		if (status) {
//			ExecutionStatus status("cSTIR_setupPrior failed", __FILE__, __LINE__);
//			handle->set(0, &status);
//			return (void*)handle;
//		}
//	}
//	catch (StirException& se) {
//		ExecutionStatus status(se);
//		handle->set(0, &status);
//	}
//	return (void*)handle;
//}
//
//extern "C"
//void* cSTIR_setupObjective(void* ptr_obj)
//{
//	return cSTIR_setupObject("objective_function", ptr_obj);
//}

//cSTIR_rayTracingMatrixParameter(const DataHandle* handle, string name)

//return parameterNotFound(name.c_str(), __FILE__, __LINE__);

//cSTIR_rayTracingMatrixParameter(const DataHandle* handle, string name);

void* cSTIR_setupObjective(void* ptr_obj);

//cSTIR_setSubiteration(hp, intDataFromHandle((void*)hv));

//return intDataHandle(cSTIR_getSubiteration(handle));

extern "C"
void cSTIR_setSubiteration(void* ptr_r, int iter)
{
	CAST_PTR(DataHandle, hr, ptr_r);
	xSTIR_subiteration(hr->data()) = iter;
}

extern "C"
int cSTIR_getSubiteration(const void* ptr_r)
{
	CAST_PTR(DataHandle, hr, ptr_r);
	return xSTIR_subiteration(hr->data());
}

//void cSTIR_setSubiteration(void* ptr_r, int i);
//int cSTIR_getSubiteration(const void* ptr_r);

//char buff[32];

//value = intDataFromHandle(hv);
//if (value < 1) {
//	sprintf(buff, "%d", value);
//	return wrongParameterValue
//		("num_tangential_LORs", buff, __FILE__, __LINE__);
//}

//NEW_SPTR(VoxelsOnCartesianGrid<float>, sptr, (IndexRange3D(0, nz - 1,
//	-(ny / 2), -(ny / 2) + ny - 1,
//	-(nx / 2), -(nx / 2) + nx - 1),
//	CartesianCoordinate3D<float>((float)z, (float)y, (float)x),
//	CartesianCoordinate3D<float>((float)sz, (float)sy, (float)sx)));
//boost::shared_ptr<VoxelsOnCartesianGrid<float> >* sptr =
//	new boost::shared_ptr<VoxelsOnCartesianGrid<float> >(
//	new VoxelsOnCartesianGrid<float>(IndexRange3D(0, nz - 1,
//		-(ny / 2), -(ny / 2) + ny - 1,
//		-(nx / 2), -(nx / 2) + nx - 1),
//		CartesianCoordinate3D<float>((float)z, (float)y, (float)x),
//		CartesianCoordinate3D<float>((float)sz, (float)sy, (float)sx)));

//void*
//cSTIR_getProjectorsUsingMatrixParameter
//(const DataHandle* hp, const char* name, DataHandle* hv)
//{
//	NEW(DataHandle, h);
//	ProjectorByBinPairUsingProjMatrixByBin* proj =
//		objectFromHandle<ProjectorByBinPair,
//		ProjectorByBinPairUsingProjMatrixByBin>(hp);
//	if (boost::iequals(name, "matrix_type")) {
//		NEW(boost::shared_ptr<ProjMatrixByBin>, sptr);
//		*sptr = proj->get_proj_matrix_sptr();
//		hv->set((void*)sptr, 0);
//		return (void*)h;
//	}
//	return parameterNotFound(name, __FILE__, __LINE__);
//}

//int
//cSTIR_setTruncateToCylindricalFOVImageProcessorParameter
//(DataHandle* hp, const char* name, const DataHandle* hv)
//{
//	CAST_PTR(boost::shared_ptr<TruncateToCylindricalFOVImageProcessor<float> >, 
//		sptr_p, hp->data());
//	int not_found = 0;
//	not_found = 1;
//	return not_found;
//}

//int

//Method* recon; 

//Method* recon = new Method(parFile);

//Method* recon = new Method();

//boost::shared_ptr<Base>* pp = new boost::shared_ptr<Base>(new Object());

//else if (boost::iequals(name,
//	"PoissonLogLikelihoodWithLinearModelForMeanAndProjData"))
//	deleteObjectHandle<GeneralisedObjectiveFunction<Image3DF> >
//	(handle);

//else if (boost::iequals(name, "ProjectorsUsingMatrix"))
//	deleteObjectHandle<ProjectorByBinPair>(handle);

//else if (boost::iequals(name, "RayTracingMatrix"))
//	deleteObjectHandle<ProjMatrixByBin>(handle);

//else if (boost::iequals(name, "QuadraticPrior"))
//	deleteObjectHandle< GeneralisedPrior<Image3DF> >(handle);

//else if (boost::iequals(name, "TruncateToCylindricalFOVImageProcessor"))
//	deleteObjectHandle< DataProcessor<Image3DF> >(handle);

//else if (boost::iequals(name, "EllipsoidalCylinder"))
//	deleteObjectHandle<Shape3D>(handle);

//CAST_PTR(Reconstruction<Image3DF>, recon, hr->data());

//extern "C"
//void cSTIR_deleteReconstruction(void* ptr_r)
//{
//	if (!ptr_r)
//		return;
//	CAST_PTR(DataHandle, hr, ptr_r);
//	//CAST_PTR(Reconstruction<Image3DF>, recon, hr->data());
//	CAST_PTR(boost::shared_ptr<Reconstruction<Image3DF> >, recon, hr->data());
//	delete recon;
//	delete hr;
//}

//CAST_PTR(IterativeReconstruction<Image3DF>, recon, hr->data());

//template<class Base, class Object = Base>
////static 
//Object*
//objectFromHandle(const DataHandle* handle) {
//	CAST_PTR(boost::shared_ptr<Base>, sptr, handle->data());
//	CAST_PTR(Object, object, sptr->get());
//	return object;
//}

//CAST_PTR(Reconstruction<Image3DF>, recon, hp->data());

//CAST_PTR(IterativeReconstruction<Image3DF>, recon, hp->data());

//xSTIR_set_initial_estimate_file(hp->data(), (char*)hv->data());

//CAST_PTR(IterativeReconstruction<Image3DF>, recon, handle->data());

//CAST_PTR(OSMAPOSLReconstruction<Image3DF>, recon, hp->data());

//CAST_PTR(OSMAPOSLReconstruction<Image3DF>, recon, handle->data());

//CAST_PTR(xSTIR_OSSPSReconstruction3DF, recon, hp->data());

//CAST_PTR(OSSPSReconstruction<Image3DF>, recon, handle->data());

//CAST_PTR(boost::shared_ptr< Reconstruction<Image3DF> >, sptr, hr->data());
//CAST_PTR(Reconstruction<Image3DF>, recon, sptr->get());

//class xSTIR_ProjectorByBinPairUsingProjMatrixByBin : 
//	public ProjectorByBinPairUsingProjMatrixByBin {
//public:
//	bool post_process() {
//// inaccessible
//		return post_processing();
//	}
//};

//CAST_PTR(xSTIR_IterativeReconstruction3DF, recon, ptr);

//CAST_PTR(xSTIR_IterativeReconstruction3DF, recon, ptr);

//CAST_PTR(boost::shared_ptr<xSTIR_IterativeReconstruction3DF>, sptr, ptr);
//xSTIR_IterativeReconstruction3DF* recon = sptr->get();

//CAST_PTR(boost::shared_ptr<xSTIR_IterativeReconstruction3DF>, sptr, ptr);
//xSTIR_IterativeReconstruction3DF* recon = sptr->get();

//template<class Base, class Object = Base>
//Object*
//objectFromHandle(const DataHandle* handle) {
//	CAST_PTR(boost::shared_ptr<Base>, sptr, handle->data());
//	CAST_PTR(Object, object, sptr->get());
//	return object;
//}

//DataProcessor<Image3DF>* proc =
//	objectFromHandle<DataProcessor<Image3DF> >(hp);

//proc->apply(**image);

//Reconstruction<Image3DF>* recon =
//	objectFromHandle< Reconstruction<Image3DF> >(hr);

//if (recon->reconstruct(*image) != Succeeded::yes) {

//CAST_PTR(sptrImage3DF, image, hi->data());

//xSTIR_updateReconstruction(hr->data(), **image);

//CAST_PTR(sptrVoxels3DF, sptr_v, hv->data());

//sptrImage3DF* sptr = new sptrImage3DF((*sptr_v)->clone());

//CAST_PTR(sptrImage3DF, sptr_i, hv->data());

//sptrImage3DF* sptr = new sptrImage3DF((*sptr_i)->clone());

//CAST_PTR(sptrImage3DF, sptr_i, hi->data());
//CAST_PTR(sptrVoxels3DF, sptr_v, hv->data());
CAST_PTR(sptrShape3D, sptr_s, hs->data());

//Image3DF& image = *sptr_i->get();
//Voxels3DF& voxels = *sptr_v->get();

//(*sptr_s)->construct_volume(voxels, num_samples);

//CAST_PTR(sptrImage3DF, sptr_i, hi->data());
//Image3DF& image = *sptr_i->get();

//CAST_PTR(sptrImage3DF, sptrImage, handle->data());
//sptrImage3DF& image = *sptrImage;

//if (!image->get_regular_range(min_indices, max_indices)) {

//CAST_PTR(sptrImage3DF, sptrImage, handle->data());
//sptrImage3DF& image = *sptrImage;

//if (!image->get_regular_range(min_indices, max_indices)) {

//CAST_PTR(sptrImage3DF, first_ptr, first_h->data());
//CAST_PTR(sptrImage3DF, second_ptr, second_h->data());
//sptrImage3DF& first_operand = *first_ptr;
//sptrImage3DF& second_operand = *second_ptr;

//if (!first_operand->has_same_characteristics(*second_operand,

//truncate_rim(*first_operand, rimsize);
//truncate_rim(*second_operand, rimsize);

//float reference_max = first_operand->find_max();
//float reference_min = first_operand->find_min();

//sptrImage3DF* sptr = new sptrImage3DF(first_image.clone());
//Image3DF& image = *(sptr->get());

//*first_operand -= *second_operand;
//const float max_error = first_operand->find_max();
//const float min_error = first_operand->find_min();
//first_image -= second_image;
//const float max_error = first_image.find_max();
//const float min_error = first_image.find_min();

//first_image += second_image;
//*first_operand += *second_operand;

//Shape3D* s = objectFromHandle<Shape3D>(hp);

//Coord3DF origin = s->get_origin();

//s->set_origin(origin);

//Shape3D* s = objectFromHandle<Shape3D>(handle);
//Coord3DF origin = s->get_origin();

//EllipsoidalCylinder* c =
//	objectFromHandle<Shape3D, EllipsoidalCylinder>(hp);

//c->set_length(value);

//c->set_radius_x(value);

//c->set_radius_y(value);

//EllipsoidalCylinder* c =
//	objectFromHandle<Shape3D, EllipsoidalCylinder>(handle);

//return floatDataHandle(c->get_length());

//return floatDataHandle(c->get_radius_x());

//return floatDataHandle(c->get_radius_y());

//CAST_PTR(boost::shared_ptr<TruncateToCylindricalFOVImageProcessor<float> >,
//	sptr_p, hp->data());

//(*sptr_p)->set_strictly_less_than_radius

//CAST_PTR(boost::shared_ptr<TruncateToCylindricalFOVImageProcessor<float> >,
//	sptr_p, handle->data());

//return intDataHandle((*sptr_p)->get_strictly_less_than_radius());

//ProjMatrixByBinUsingRayTracing* matrix =
//	objectFromHandle<ProjMatrixByBin, ProjMatrixByBinUsingRayTracing>(hp);

//matrix->set_num_tangential_LORs(value);

//ProjMatrixByBinUsingRayTracing* matrix =
//	objectFromHandle<ProjMatrixByBin,
//	ProjMatrixByBinUsingRayTracing>(handle);

//return intDataHandle(matrix->get_num_tangential_LORs());

//ProjectorByBinPairUsingProjMatrixByBin* proj =
//	objectFromHandle<ProjectorByBinPair,
//	ProjectorByBinPairUsingProjMatrixByBin>(hp);

//proj->set_proj_matrix_sptr

//ProjectorByBinPairUsingProjMatrixByBin* proj =
//	objectFromHandle<ProjectorByBinPair,
//	ProjectorByBinPairUsingProjMatrixByBin>(handle);

//*sptr = proj->get_proj_matrix_sptr();

//CAST_PTR(boost::shared_ptr<GeneralisedPrior<Image3DF> >, prior, hp->data());

//(*prior)->set_penalisation_factor(floatDataFromHandle((void*)hv));

//CAST_PTR(boost::shared_ptr<GeneralisedPrior<Image3DF> >, 
//	prior, handle->data());

//return floatDataHandle((*prior)->get_penalisation_factor());

//xSTIR_QuadraticPrior3DF* prior =
//	objectFromHandle<GeneralisedPrior<Image3DF>, xSTIR_QuadraticPrior3DF>(hp);

//prior->only2D(intDataFromHandle((void*)hv));

//CAST_PTR(boost::shared_ptr<GeneralisedObjectiveFunction<Image3DF> >, 
//	sptr, hp->data());

//(*sptr)->set_prior_sptr

//PoissonLogLikelihoodWithLinearModelForMean<Image3DF>*
//	obj_fun = objectFromHandle<GeneralisedObjectiveFunction<Image3DF>,
//	PoissonLogLikelihoodWithLinearModelForMean<Image3DF> >
//	(hp);

//obj_fun->set_sensitivity_filename((char*)hv->data());

//obj_fun->set_use_subset_sensitivities

//obj_fun->set_recompute_sensitivity

//xSTIR_PoissonLogLikelihoodWithLinearModelForMeanAndProjData3DF*
//	obj_fun =
//	objectFromHandle<GeneralisedObjectiveFunction<Image3DF>,
//	xSTIR_PoissonLogLikelihoodWithLinearModelForMeanAndProjData3DF>
//	(hp);

//obj_fun->set_input_file((char*)hv->data());

//obj_fun->set_zero_seg0_end_planes

//obj_fun->set_max_segment_num_to_process(intDataFromHandle((void*)hv));

//obj_fun->set_projector_pair_sptr

//Reconstruction<Image3DF>* recon =
//	objectFromHandle< Reconstruction<Image3DF> >(hp);

//recon->set_output_filename_prefix((char*)hv->data());

//IterativeReconstruction<Image3DF>* recon =
//	objectFromHandle< Reconstruction<Image3DF>,
//	IterativeReconstruction<Image3DF> >(hp);

//recon->set_inter_iteration_filter_ptr(*sptr);

//recon->set_objective_function_sptr(*sptr);

//xSTIR_set_initial_estimate_file(recon, (char*)hv->data());

//recon->set_num_subsets(value);

//if (value < 0 || value >= recon->get_num_subsets())

//recon->set_start_subset_num(value);

//recon->set_num_subiterations(value);

//recon->set_start_subiteration_num(value);

//xSTIR_subiteration(recon) = value;

//if (value < 1 || value > recon->get_num_subiterations())

//recon->set_save_interval(value);

//recon->set_inter_iteration_filter_interval(value);

//IterativeReconstruction<Image3DF>* recon =
//	objectFromHandle< Reconstruction<Image3DF>,
//	IterativeReconstruction<Image3DF> >(handle);

//return intDataHandle(recon->get_num_subsets());

//return intDataHandle(recon->get_start_subset_num());

//return intDataHandle(recon->get_num_subiterations());

//return intDataHandle(recon->get_start_subiteration_num());

//return intDataHandle(xSTIR_subiteration(recon));

//*sptr = recon->get_objective_function_sptr();

//OSMAPOSLReconstruction<Image3DF>* recon =
//	objectFromHandle< Reconstruction<Image3DF>,
//	OSMAPOSLReconstruction<Image3DF> >(hp);

//recon->set_MAP_model((char*)hv->data());

//OSMAPOSLReconstruction<Image3DF>* recon =
//	objectFromHandle< Reconstruction<Image3DF>,
//	OSMAPOSLReconstruction<Image3DF> >(handle);

//xSTIR_OSSPSReconstruction3DF* recon =
//	objectFromHandle< Reconstruction<Image3DF>,
//	xSTIR_OSSPSReconstruction3DF >(hp);

//recon->relaxation_parameter_value() = floatDataFromHandle(hv);

//xSTIR_OSSPSReconstruction3DF* recon =
//	objectFromHandle< Reconstruction<Image3DF>,
//	xSTIR_OSSPSReconstruction3DF >(handle);

if (boost::iequals(name, "name"))
return charDataHandle
("PoissonLogLikelihoodWithLinearModelForMeanAndProjData");

if (boost::iequals(name, "name"))
return charDataHandle("OSMAPOSL");

if (boost::iequals(name, "name"))
return charDataHandle("OSSPS");

//NEW(DataHandle, h);
//NEW(boost::shared_ptr<ProjMatrixByBin>, sptr);
//*sptr = proj.get_proj_matrix_sptr();
//h->set((void*)sptr, 0);
//return (void*)h;

//NEW(boost::shared_ptr<GeneralisedPrior<Image3DF> >, sptr);
//*sptr = obj_fun.get_prior_sptr();
//NEW(DataHandle, h);
//h->set((void*)sptr, 0);
//return (void*)h;

//(*(boost::shared_ptr<ProjectorByBinPairUsingProjMatrixByBin>*)hv->data());

//NEW(DataHandle, h);
//NEW(boost::shared_ptr<ProjectorByBinPair>, sptr);
//*sptr = obj_fun.get_projector_pair_sptr();
//h->set((void*)sptr, 0);
//return (void*)h;

//NEW(DataHandle, h);
//NEW(boost::shared_ptr<GeneralisedObjectiveFunction<Image3DF> >, sptr);
//*sptr = recon.get_objective_function_sptr();
//h->set((void*)sptr, 0);
//return (void*)h;

//NEW(DataHandle, h);
//NEW(boost::shared_ptr<DataProcessor<Image3DF> >, sptr);
//*sptr = recon.get_inter_iteration_filter_sptr();
//h->set((void*)sptr, 0);
//return (void*)h;

//NEW(DataHandle, h);
//NEW(boost::shared_ptr<GeneralisedObjectiveFunction<Image3DF> >, sptr);
//*sptr = recon.get_objective_function_sptr();
//h->set((void*)sptr, 0);
//return (void*)h;

//(*(boost::shared_ptr<ProjMatrixByBin>*)hv->data());

//(*(boost::shared_ptr<GeneralisedPrior<Image3DF> >*)hv->data());

//(*(boost::shared_ptr<ProjectorByBinPair>*)hv->data());

//CAST_PTR(boost::shared_ptr<DataProcessor<Image3DF> >, sptr, hv->data());
//recon.set_inter_iteration_filter_ptr(*sptr);

//CAST_PTR(boost::shared_ptr<GeneralisedObjectiveFunction<Image3DF> >,
//	sptr, hv->data());

//(boost::iequals((char*)hv->data(), "true"));

//((char*)hv->data());

//(boost::iequals((char*)hv->data(), "true"));

//(boost::iequals((char*)hv->data(), "true"));

//obj_fun.set_input_file((char*)hv->data());

//(boost::iequals((char*)hv->data(), "true"));

//recon.set_output_filename_prefix((char*)hv->data());

//xSTIR_set_initial_estimate_file(&recon, (char*)hv->data());

//recon.set_MAP_model((char*)hv->data());

//template<class Base, class Object = Base>
//static void*
//newObjHandle()
//{
//	typedef ObjectHandle<Object, Base> Obj;
//	NEW_SPTR1(Obj, ptr_handle);
//	return (void*)ptr_handle;
//}
//
//template<class Base, class Object = Base>
//static void
//deleteObjHandle(void* ptr)
//{
//	typedef ObjectHandle<Object, Base> Obj;
//	if (!ptr)
//		return;
//	CAST_PTR(boost::shared_ptr<Obj>, ptr_sptr, ptr);
//	delete ptr_sptr;
//}

//const DataHandle* ptr_h = (const DataHandle*)ptr;
//void* ptr_d = ptr_h->data();
//if (!ptr_d)
//	return 0;
//else
//	return (char*)ptr_d;

template<class Object, class Base = Object>
class ObjectHandle : public DataHandle {
public:
	ObjectHandle(const ExecutionStatus* status = 0) : DataHandle() {
		ptr_sptr = new boost::shared_ptr<Base>(new Object);
		_data = (void*)ptr_sptr;
		_status = status;
	}
	ObjectHandle(const boost::shared_ptr<Base>& obj,
		const ExecutionStatus* status = 0) : DataHandle() {
		ptr_sptr = new boost::shared_ptr<Base>;
		*ptr_sptr = obj;
		_data = (void*)ptr_sptr;
		_status = new ExecutionStatus(*status);
	}
};

//T* ptr = (T*)malloc(sizeof(T));
//*ptr = x;
//h->set((void*)ptr, 0, GRAB);

//void* refDataHandle(void* ptr);

//class anObjectHandle {
//public:
//	virtual ~anObjectHandle() {}
//	virtual anObjectHandle* copy() = 0;
//};
//
//template<class Object, class Base = Object>
//class ObjectHandle : public DataHandle, public anObjectHandle {
//public:
//	ObjectHandle(const ObjectHandle& obj) {
//		NEW(boost::shared_ptr<Base>, ptr_sptr);
//		*ptr_sptr = *(boost::shared_ptr<Base>*)obj.data();
//		_data = (void*)ptr_sptr;
//		if (obj._status)
//			_status = new ExecutionStatus(*obj._status);
//		_name = obj._name;
//	}
//	ObjectHandle(const char* name, const ExecutionStatus* status = 0) : 
//		DataHandle(name) {
//		//printf("ObjectHandle constructor called...");
//		//boost::shared_ptr<Base>* ptr_sptr = new boost::shared_ptr<Base>(new Object);
//		NEW_SPTR2(Base, ptr_sptr, Object);
//		//printf("ok\n");
//		_data = (void*)ptr_sptr;
//		if (status)
//			_status = new ExecutionStatus(*status);
//	}
//	ObjectHandle(const char* name, const boost::shared_ptr<Base>& obj,
//		const ExecutionStatus* status = 0) : DataHandle(name) {
//		//boost::shared_ptr<Base>* ptr_sptr = new boost::shared_ptr<Base>;
//		NEW(boost::shared_ptr<Base>, ptr_sptr);
//		*ptr_sptr = obj;
//		_data = (void*)ptr_sptr;
//		if (status)
//			_status = new ExecutionStatus(*status);
//	}
//	virtual ~ObjectHandle() {
//		//printf("ObjectHandle destructor called...");
//		CAST_PTR(boost::shared_ptr<Base>, ptr_sptr, _data);
//		delete _status;
//		//printf("ok\n");
//		delete ptr_sptr;
//		//printf("ok\n");
//	}
//	virtual anObjectHandle* copy() {
//		CAST_PTR(boost::shared_ptr<Base>, ptr_sptr, _data);
//		return new ObjectHandle<Object, Base>
//			(_name.c_str(), *ptr_sptr, _status);
//	}
//	virtual const char* name() const {
//		return _name.c_str();
//	}
////private:
////	std::string _name;
//};

#define NEW_OBJ(T, X, Y) T* X = new T(Y)
#define NEW_SPTR1(T, X) \
	boost::shared_ptr< T >* X = new boost::shared_ptr< T >(new T)

//template<class Base, class Object = Base>
//static void*
//newObjHandle()
//{
//	typedef ObjectHandle<Object, Base> Obj;
//	//NEW(Obj, ptr_handle);
//	Obj* ptr_handle = new Obj(Object::registered_name);
//	return (void*)ptr_handle;
//}
//
//static void*
//copyOfObjHandle(void* ptr)
//{
//	CAST_PTR(anObjectHandle, ptr_obj, ptr);
//	return (void*)ptr_obj->copy();
//}
//
//static void
//deleteObjHandle(void* ptr)
//{
//	if (!ptr)
//		return;
//	CAST_PTR(anObjectHandle, ptr_obj, ptr);
//	delete ptr_obj;
//}

//if (boost::iequals(name, "_name")) {
//	CAST_PTR(aDataHandle, h, ptr);
//	return charDataHandle(h->name());
//}

class aDataHandle {
public:
	virtual ~aDataHandle() {}
	virtual const char* name() const = 0;
};

//: public aDataHandle {

//const char* name = "") : 

//, _name(name) {}

//virtual const char* name() const {
//	return _name.c_str();
//}

//std::string _name;

//ObjectHandle<Reconstruction<Image3DF> >* ptr_handle =
//	new ObjectHandle<Reconstruction<Image3DF> >(*ptr_sptr);
//return (void*)ptr_handle;
//NEW_SPTR(Reconstruction<Image3DF>, recon, Method(parFile));
//handle->set((void*)recon, 0);

//ObjectHandle<Reconstruction<Image3DF> >* ptr_handle =
//	new ObjectHandle<Reconstruction<Image3DF> >(*ptr_sptr);
//return (void*)ptr_handle;
//NEW_SPTR(Reconstruction<Image3DF>, recon, Method);
//handle->set((void*)recon, 0);

//if (boost::iequals(name,
//	"PoissonLogLikelihoodWithLinearModelForMeanAndProjData"))
//	return newObjectHandle<GeneralisedObjectiveFunction<Image3DF>,
//	PoissonLogLikelihoodWithLinearModelForMeanAndProjData<Image3DF> >
//	();
//else if (boost::iequals(name, "ProjectorsUsingMatrix"))
//	return newObjectHandle
//	<ProjectorByBinPair, ProjectorByBinPairUsingProjMatrixByBin>();
//if (boost::iequals(name, "RayTracingMatrix"))
//	return newObjectHandle
//	<ProjMatrixByBin, ProjMatrixByBinUsingRayTracing>();
//else if (boost::iequals(name, "QuadraticPrior"))
//	return newObjectHandle
//	<GeneralisedPrior<Image3DF>, QuadraticPrior<float> >();
//else if (boost::iequals(name, "TruncateToCylindricalFOVImageProcessor"))
//	return newObjectHandle<DataProcessor<Image3DF>,
//	TruncateToCylindricalFOVImageProcessor<float> >();
//if (boost::iequals(name, "EllipsoidalCylinder"))
//	return newObjectHandle<Shape3D, EllipsoidalCylinder>();

//printf("%s\n", name);
//printf("%p\n", ptr);
CAST_PTR(DataHandle, handle, ptr);
//if (boost::iequals(name, "ObjectiveFunction"))
//	return copyOfObjectHandle<GeneralisedObjectiveFunction<Image3DF> >
//	(handle);
if (boost::iequals(name, "ObjectiveFunction"))
return copyOfObjHandle(ptr);
//else if (boost::iequals(name, "Projectors"))
//	return copyOfObjectHandle<ProjectorByBinPair>(handle);
//if (boost::iequals(name, "ProjMatrix"))
//	return copyOfObjectHandle<ProjMatrixByBin>(handle);
//else if (boost::iequals(name, "Prior"))
//	return copyOfObjectHandle<GeneralisedPrior<Image3DF> >(handle);
//else if (boost::iequals(name, "DataProcessor"))
//	return copyOfObjectHandle<DataProcessor<Image3DF> >(handle);
else if (boost::iequals(name, "DataProcessor"))
return copyOfObjHandle(ptr);
//if (boost::iequals(name, "Shape"))
//	return copyOfObjectHandle<Shape3D>(handle);
else
return unknownObject("object", name, __FILE__, __LINE__);

CAST_PTR(DataHandle, handle, ptr);
if (boost::iequals(name, "ObjectiveFunction"))
deleteObjHandle(ptr);
//if (boost::iequals(name, "ObjectiveFunction"))
//	deleteObjectHandle<GeneralisedObjectiveFunction<Image3DF> >
//	(handle);
//else if (boost::iequals(name, "Projectors"))
//	deleteObjectHandle<ProjectorByBinPair>(handle);
//else if (boost::iequals(name, "ProjMatrix"))
//	deleteObjectHandle<ProjMatrixByBin>(handle);
//else if (boost::iequals(name, "Prior"))
//	deleteObjectHandle< GeneralisedPrior<Image3DF> >(handle);
//else if (boost::iequals(name, "DataProcessor"))
//	deleteObjectHandle< DataProcessor<Image3DF> >(handle);
//else if (boost::iequals(name, "Shape"))
//	deleteObjectHandle<Shape3D>(handle);
//else if (boost::iequals(name, "Image"))
//	deleteObjectHandle<Image3DF>(handle);
//else if (boost::iequals(name, "Voxels"))
//	deleteObjectHandle<Voxels3DF>(handle);
//else if (boost::iequals(name, "Reconstruction"))
//	deleteObjectHandle< Reconstruction<Image3DF> >(handle);
else
deleteObjHandle(ptr);

//return copyOfObjHandle(ptr);

//deleteObjHandle(ptr);

//DataHandle* handle = new DataHandle;
//handle->set((void*)sptr, 0);

//DataHandle* handle = new DataHandle;
//handle->set((void*)sptr, 0);

//DataHandle* handle = new DataHandle;
//handle->set((void*)sptr, 0);

//handle->set((void*)sptr, 0);

//return newObjHandle<Reconstruction<Image3DF> >(ptr_sptr);

//return newObjHandle<Reconstruction<Image3DF> >(ptr_sptr);

//ObjectHandle<Image3DF>* handle = new ObjectHandle<Image3DF>(*sptr);
//return (void*)handle;

//ObjectHandle<Image3DF>* handle = new ObjectHandle<Image3DF>(*sptr);
//return (void*)handle;

//ObjectHandle<Image3DF>* handle = new ObjectHandle<Image3DF>(*sptr);
//return (void*)handle;

//ObjectHandle<Image3DF>* handle = new ObjectHandle<Image3DF>(*sptr);
//return (void*)handle;

{
	printf("anObjectHandle.copy\n");
	return 0;
}

virtual void print(const char* text) {
	printf("anObjectHandle.prints: %s\n");
}

//class ObjectHandle : public DataHandle, public anObjectHandle {

virtual void print(const char* text) {
	printf("ObjectHandle.prints: %s\n");
}
//protected:
//	void* _data;
//	ExecutionStatus* _status;

//anObjectHandle* ptr = (anObjectHandle*)ptr_handle;
//printf("%p\n", ptr_handle);
//printf("%p\n", ptr);
//printf("%p\n", (void*)ptr_handle);
//printf("%p\n", (void*)ptr);
//ptr_handle->print("in constructor");
//ptr->print();

//static void*
//copyOfObjHandle(void* ptr)
//{
//	//printf("%p\n", ptr);
//	CAST_PTR(anObjectHandle, ptr_obj, ptr);
//	//printf("copying...\n");
//	//return (void*)ptr_obj;
//	//ptr_obj->print("copying");
//	//anObjectHandle* copy = ptr_obj->copy();
//	//return (void*)copy;
//	return (void*)ptr_obj->copy();
//}
//
//static void
//deleteObjHandle(void* ptr)
//{
//	//printf("%p\n", ptr);
//	if (!ptr)
//		return;
//	CAST_PTR(anObjectHandle, ptr_obj, ptr);
//	//ptr_obj->print("in destructor");
//	delete ptr_obj;
//}

//template<class Base, class Object = Base>
//void*
//newObjectHandle()
//{
//	DataHandle* handle = new DataHandle;
//	NEW_SPTR(Base, ptr_sptr, Object);
//	handle->set((void*)ptr_sptr, 0);
//	return (void*)handle;
//}
//
//template<class Base>
//void*
//copyOfObjectHandle(const DataHandle* handle)
//{
//	DataHandle* copy = new DataHandle;
//	boost::shared_ptr<Base>* ptr_sptr =
//		new boost::shared_ptr<Base>(sptrDataFromHandle<Base>(handle));
//	copy->set((void*)ptr_sptr, 0);
//	return (void*)copy;
//}
//
//template<class Object>
//void
//deleteObjectHandle(DataHandle* handle)
//{
//	if (!handle)
//		return;
//	CAST_PTR(boost::shared_ptr<Object>, sptr, handle->data());
//	delete sptr;
//	delete handle;
//}

//template<class T>
//void*
//sptrDataHandle(boost::shared_ptr<T> sptr) {
//	NEW(boost::shared_ptr<T>, ptr_sptr);
//	*ptr_sptr = sptr;
//	NEW(DataHandle, handle);
//	handle->set((void*)ptr_sptr, 0);
//	return (void*)handle;
//}

//if (boost::iequals(name, "matrix_type"))
//	return sptrDataHandle(proj.get_proj_matrix_sptr());

//if (boost::iequals(name, "prior"))
//	return sptrDataHandle(obj_fun.get_prior_sptr());

//if (boost::iequals(name, "projector_pair_type"))
//	return sptrDataHandle(obj_fun.get_projector_pair_sptr());

//if (boost::iequals(name, "objective_function"))
//	return sptrDataHandle(recon.get_objective_function_sptr());

//if (boost::iequals(name, "inter_iteration_filter_type"))
//	return sptrDataHandle(recon.get_inter_iteration_filter_sptr());

//if (boost::iequals(name, "objective_function"))
//	return sptrDataHandle(recon.get_objective_function_sptr());

void* refDataHandle(void* ptr) {
	DataHandle* h = (DataHandle*)ptr;
	DataHandle* hr = new DataHandle;
	hr->set(h->data(), h->status());
	return (void*)hr;
}

//CAST_PTR(sptrImage3DF, ptr_sptr_image, hd->data());
//processor.apply(**ptr_sptr_image);

//CAST_PTR(sptrImage3DF, ptr_sptr_image, hi->data());

//Succeeded s = xSTIR_setupReconstruction(hr->data(), *ptr_sptr_image);

//CAST_PTR(sptrImage3DF, ptr_sptr_image, hi->data());

//if (recon.reconstruct(*ptr_sptr_image) != Succeeded::yes) {

//extern "C"
//void* cSTIR_newReconstruction(const char* method, const char* filename)
//{
//	void* recon;
//	if (boost::iequals(method, "OSMAPOSL"))
//		recon = cSTIR_newReconstructionMethod<OSMAPOSLReconstruction<Image3DF> >
//		(filename);
//	else if (boost::iequals(method, "OSSPS"))
//		recon = cSTIR_newReconstructionMethod<OSSPSReconstruction<Image3DF> >
//		(filename);
//	else
//		recon = unknownObject("method", method, __FILE__, __LINE__);
//	return recon;
//}

//extern "C"
//void* cSTIR_imageFromFile(const char* filename) 
//{
//	try {
//		sptrImage3DF* sptr = new sptrImage3DF
//			(read_from_file<Image3DF>(filename));
//		return newObjectHandle(sptr);
//	}
//	catch (StirException& se) {
//		ExecutionStatus status(se);
//		DataHandle* handle = new DataHandle;
//		handle->set(0, &status);
//		return (void*)handle;
//	}
//	catch (...) {
//		ExecutionStatus status("unhandled exception", __FILE__, __LINE__);
//		DataHandle* handle = new DataHandle;
//		handle->set(0, &status);
//		return (void*)handle;
//	}
//}

//if ((value = floatDataFromHandle(hv)) < 1)
//	return wrongFloatParameterValue("EllipsoidalCylinder::length",
//	value, __FILE__, __LINE__);

//if ((value = floatDataFromHandle(hv)) < 1)
//	return wrongFloatParameterValue("EllipsoidalCylinder::radius_x",
//	value, __FILE__, __LINE__);

//if ((value = floatDataFromHandle(hv)) < 1)
//	return wrongFloatParameterValue("EllipsoidalCylinder::radius_y",
//	value, __FILE__, __LINE__);

//if ((value = intDataFromHandle(hv)) < 1)
//	return wrongIntParameterValue
//	("ProjMatrixByBinUsingRayTracing::num_tangential_LORs",
//	value, __FILE__, __LINE__);

//if (value < 1)
//	return wrongIntParameterValue
//	("IterativeReconstruction::num_subsets", 
//	value, __FILE__, __LINE__);

//if (value < 0 || value >= recon.get_num_subsets())
//	return wrongIntParameterValue
//	("IterativeReconstruction::start_subset_num",
//	value, __FILE__, __LINE__);

//if (value < 1)
//	return wrongIntParameterValue
//	("IterativeReconstruction::num_subiterations",
//	value, __FILE__, __LINE__);

//if (value < 1)
//	return wrongIntParameterValue
//	("IterativeReconstruction::start_subiteration_num",
//	value, __FILE__, __LINE__);

//if (value < 1)
//	return wrongIntParameterValue
//	("IterativeReconstruction::subiteration_num",
//	value, __FILE__, __LINE__);

//if (value < 1 || value > recon.get_num_subiterations())
//	return wrongIntParameterValue
//	("IterativeReconstruction::save_interval",
//	value, __FILE__, __LINE__);

//if (value < 0)
//	return wrongIntParameterValue
//	("IterativeReconstruction::inter_iteration_filter_interval",
//	value, __FILE__, __LINE__);

//return newObjectHandle<GeneralisedObjectiveFunction<Image3DF>,
//PoissonLogLikelihoodWithLinearModelForMeanAndProjData<Image3DF> >

//ProjectorByBinPairUsingProjMatrixByBin>();

//ProjMatrixByBinUsingRayTracing>();

//		<GeneralisedPrior<Image3DF>, QuadraticPrior<float> >();

//<DataProcessor<Image3DF>,
//TruncateToCylindricalFOVImageProcessor<float> >();

//TruncateToCylindricalFOVImageProcessor<float>& proc =
//	objectFromHandle<DataProcessor<Image3DF>,
//	TruncateToCylindricalFOVImageProcessor<float> >(hp);

//TruncateToCylindricalFOVImageProcessor<float>& proc =
//	objectFromHandle<DataProcessor<Image3DF>,
//	TruncateToCylindricalFOVImageProcessor<float> >(handle);

//ProjMatrixByBinUsingRayTracing& matrix =
//	objectFromHandle<ProjMatrixByBin, ProjMatrixByBinUsingRayTracing>
//	(hp);

//ProjMatrixByBinUsingRayTracing& matrix =
//	objectFromHandle<ProjMatrixByBin,
//	ProjMatrixByBinUsingRayTracing>(handle);

//ProjectorByBinPairUsingProjMatrixByBin& proj =
//	objectFromHandle<ProjectorByBinPair,
//	ProjectorByBinPairUsingProjMatrixByBin>(hp);

//ProjectorByBinPairUsingProjMatrixByBin& proj =
//	objectFromHandle<ProjectorByBinPair,
//	ProjectorByBinPairUsingProjMatrixByBin>(handle);

//GeneralisedPrior<Image3DF>& prior =
//	objectFromHandle<GeneralisedPrior<Image3DF> >(hp);

//GeneralisedPrior<Image3DF>& prior =
//	objectFromHandle<GeneralisedPrior<Image3DF> >(handle);

//objectFromHandle<GeneralisedPrior<Image3DF>, xSTIR_QuadraticPrior3DF>

//GeneralisedObjectiveFunction<Image3DF>& obj_fun =
//	objectFromHandle<GeneralisedObjectiveFunction<Image3DF> >(hp);

//GeneralisedObjectiveFunction<Image3DF>& obj_fun =
//	objectFromHandle<GeneralisedObjectiveFunction<Image3DF> >(handle);

//PoissonLogLikelihoodWithLinearModelForMean<Image3DF>&
//	obj_fun = objectFromHandle<GeneralisedObjectiveFunction<Image3DF>,
//	PoissonLogLikelihoodWithLinearModelForMean<Image3DF> >
//	(hp);

//IterativeReconstruction<Image3DF>& recon =
//	objectFromHandle< Reconstruction<Image3DF>,
//	IterativeReconstruction<Image3DF> >(hp);

//(sptrDataFromHandle<DataProcessor<Image3DF> >(hv));

//(sptrDataFromHandle<GeneralisedObjectiveFunction<Image3DF> >(hv));

//IterativeReconstruction<Image3DF>& recon =
//	objectFromHandle< Reconstruction<Image3DF>,
//	IterativeReconstruction<Image3DF> >(handle);

//catch (StirException& se) {
//	ExecutionStatus status(se);
//	DataHandle* handle = new DataHandle;
//	handle->set(0, &status);
//	return (void*)handle;
//}
//catch (...) {
//	ExecutionStatus status("unhandled exception", __FILE__, __LINE__);
//	DataHandle* handle = new DataHandle;
//	handle->set(0, &status);
//	return (void*)handle;
//}

//catch (StirException& se) {
//	ExecutionStatus status(se);
//	handle->set(0, &status);
//}
//catch (...) {
//	ExecutionStatus status("unhandled exception", __FILE__, __LINE__);
//	handle->set(0, &status);
//}
//return (void*)handle;

//catch (StirException& se) {
//	ExecutionStatus status(se);
//	handle->set(0, &status);
//}
////catch (exception& e) {
////	cout << "exception " << e.what() << endl;
////}
//catch (...) {
//	ExecutionStatus status("unhandled exception", __FILE__, __LINE__);
//	handle->set(0, &status);
//}

//return (void*)handle;

//CAST_PTR(DataHandle, hi, ptr_i);
//Image3DF& image = objectFromHandle<Image3DF>(hi);

//CAST_PTR(DataHandle, handle, ptr);
//Image3DF& image = objectFromHandle<Image3DF>(handle);

//data[i] = (*image)[z][y][x];

void* cSTIR_acquisitionModelDataFromTemplate(void* ptr_t);

// obsolete model methods
void* cSTIR_acquisitionDataFromTemplate(const char* file, const char* tmpl);
void* cSTIR_setupAcquisitionModel(void* ptr_am, void* ptr_ad, void* ptr_im);
void* cSTIR_acquisitionModelFwd(void* ptr_am, void* ptr_im, void* ptr_ad);
void* cSTIR_acquisitionModelBwd(void* ptr_am, void* ptr_ad, void* ptr_im);

extern "C"
void* cSTIR_acquisitionDataFromTemplate(const char* datafile, const char* templ)
{
	try {
		sptrProjData sptr_t = ProjData::read_from_file(templ);
		NEW_SPTR(ProjData, ptr_sptr,
			ProjDataInterfile(sptr_t->get_exam_info_sptr(),
			sptr_t->get_proj_data_info_sptr(), datafile));
		return newObjectHandle(ptr_sptr);
	}
	CATCH
}

extern "C"
void* cSTIR_setupAcquisitionModel(void* ptr_am, void* ptr_ad, void* ptr_im)
{
	try {
		CAST_PTR(DataHandle, ha, ptr_am);
		CAST_PTR(DataHandle, hd, ptr_ad);
		CAST_PTR(DataHandle, hi, ptr_im);
		sptrProjPair& sptr_am = objectSptrFromHandle<ProjectorByBinPair>(ha);
		sptrProjData& sptr_ad = objectSptrFromHandle<ProjData>(hd);
		sptrImage3DF& sptr_im = objectSptrFromHandle<Image3DF>(hi);
		Succeeded s =
			sptr_am->set_up(sptr_ad->get_proj_data_info_sptr(), sptr_im);
		DataHandle* handle = new DataHandle;
		if (s != Succeeded::yes) {
			ExecutionStatus status("cSTIR_setupAcquisitionModel failed",
				__FILE__, __LINE__);
			handle->set(0, &status);
		}
		return (void*)handle;
	}
	CATCH
}

extern "C"
void* cSTIR_acquisitionModelFwd(void* ptr_am, void* ptr_im, void* ptr_ad)
{
	try {
		CAST_PTR(DataHandle, ha, ptr_am);
		CAST_PTR(DataHandle, hd, ptr_ad);
		CAST_PTR(DataHandle, hi, ptr_im);
		sptrProjPair& sptr_am = objectSptrFromHandle<ProjectorByBinPair>(ha);
		sptrProjData& sptr_ad = objectSptrFromHandle<ProjData>(hd);
		sptrImage3DF& sptr_im = objectSptrFromHandle<Image3DF>(hi);
		sptr_am->get_forward_projector_sptr()->forward_project
			(*sptr_ad, *sptr_im);
		return new DataHandle;
	}
	CATCH
}

extern "C"
void* cSTIR_acquisitionModelBwd(void* ptr_am, void* ptr_ad, void* ptr_im)
{
	try {
		CAST_PTR(DataHandle, ha, ptr_am);
		CAST_PTR(DataHandle, hd, ptr_ad);
		CAST_PTR(DataHandle, hi, ptr_im);
		sptrProjPair& sptr_am = objectSptrFromHandle<ProjectorByBinPair>(ha);
		sptrProjData& sptr_ad = objectSptrFromHandle<ProjData>(hd);
		sptrImage3DF& sptr_im = objectSptrFromHandle<Image3DF>(hi);
		sptr_am->get_back_projector_sptr()->back_project(*sptr_im, *sptr_ad);
		return new DataHandle;
	}
	CATCH
}

//std::cerr << "\nINFO: "
// << string
// << std::endl;

//sptrImage3DF* ptr_sptr = (sptrImage3DF*)handle->data();
if (ptr_sptr == 0)
return 0;

//NEW(boost::shared_ptr<T>, ptr_sptr);
//*ptr_sptr = sptr;
//ObjectHandle<T>* ptr_handle = new ObjectHandle<T>(*ptr_sptr);

//// DataHandle methods
//void* newDataHandle();
//void* charDataHandle(const char* s);
//void* intDataHandle(int i);
//void* floatDataHandle(float i);
//void* doubleDataHandle(double i);
//char* charDataFromHandle(const void* ptr);
//int intDataFromHandle(const void* ptr);
//float floatDataFromHandle(const void* ptr);
//double doubleDataFromHandle(const void* ptr);
//void deleteDataHandle(void* ptr);

//// ExecutionStatus methods
//int executionStatus(const void* ptr);
//const char* executionError(const void* ptr);
//const char* executionErrorFile(const void* ptr);
//int executionErrorLine(const void* ptr);

//extern "C" {

//#ifndef CSTIR_FOR_MATLAB
//	void cSTIR_getImageDimensions(const void* ptr, size_t pd);
//	void cSTIR_getImageData(const void* ptr, size_t pd);
//#else
//	void cSTIR_getImageDimensions(const void* ptr, int* pd);
//	void cSTIR_getImageData(const void* ptr, double* pd);
//#endif

//void cSTIR_deleteObject(void* ptr);

//void* cSTIR_copyOfObject(void* ptr);

// create empty-copy in memory

//extern "C"
//void* cSTIR_copyOfObject(void* ptr)
//{
//	try {
//		CAST_PTR(anObjectHandle, ptr_obj, ptr);
//		return (void*)ptr_obj->copy();
//	}
//	CATCH
//}
//
//extern "C"
//void cSTIR_deleteObject(void* ptr)
//{
//	if (!ptr)
//		return;
//	CAST_PTR(anObjectHandle, ptr_obj, ptr);
//	delete ptr_obj;
//}

//filename = path + "Utahscat600k_ca_seg4.hs";
//boost::shared_ptr<ProjData> sptr_t = ProjData::read_from_file(filename);
//Succeeded sppm =
//	sptr_ppm->set_up(sptr_t->get_proj_data_info_sptr(), sptr_image);
//boost::shared_ptr<ProjData> sptr_fd(
//	new ProjDataInMemory(sptr_t->get_exam_info_sptr(),
//	sptr_t->get_proj_data_info_sptr()));
//sptr_ppm->get_forward_projector_sptr()->forward_project
//	(*sptr_fd, *sptr_image);

//size_t size = sptr_background_->size_all();
//size_t foreground_size = sptr_fd->size_all();
//if (size != foreground_size) {
//	std::cout << "wrong background term size " << size
//		<< ", must be " << foreground_size << ", skipping\n";
//	return sptr_fd;
//}
//std::cout << "adding background term...\n";
//double* fdata = new double[size];
//double* bdata = new double[size];
//sptr_fd->copy_to(fdata);
//sptr_background_->copy_to(bdata);
//add_(size, fdata, bdata);
//sptr_fd->fill_from(fdata);
//delete[] fdata;
//delete[] bdata;

//#define CAST_PTR(T, X, Y) T* X = (T*)Y
//#define NEW_SPTR(Base, X, Object) \
//	boost::shared_ptr< Base >* X = new boost::shared_ptr< Base >(new Object)
//#define CATCH \
//	catch (LocalisedException& se) {\
//		ExecutionStatus status(se);\
//		DataHandle* handle = new DataHandle;\
//		handle->set(0, &status);\
//		return (void*)handle;\
//	}\
//	catch (...) {\
//		ExecutionStatus status("unhandled exception", __FILE__, __LINE__);\
//		DataHandle* handle = new DataHandle;\
//		handle->set(0, &status);\
//		return (void*)handle;\
//	}\

//((ProjectorPairUsingMatrix*)sptr_ppm_.get())->
//	set_proj_matrix_sptr(sptr_matrix);

//boost::shared_ptr<ProjectorByBinPair> sptr_ppm_;

//PETAcquisitionModel(

//boost::shared_ptr<ProjectorByBinPair> sptr_projectors,

//sptr_projectors_ = sptr_projectors;

//void set_multiplicative_term(boost::shared_ptr<ProjData> sptr)
//{
//	sptr_mult_ = sptr;
//}

//boost::shared_ptr<ProjData> sptr_mult_;

//sptr_ppm_->set_proj_matrix_sptr(sptr_matrix);

//void set_up(
//	boost::shared_ptr<ProjData> sptr_acq,
//	boost::shared_ptr<Image> sptr_image)
//{
//	PETAcquisitionModel::set_up(sptr_ppm_, sptr_acq, sptr_image);
//}

//return sptr_ppm_->get_proj_matrix_sptr();

//boost::shared_ptr<ProjectorByBinPair> projectors_sptr()
//{
//	return sptr_ppm_;
//}

//sptr_ppm_.reset(new ProjectorPairUsingMatrix);
//set_projectors(sptr_ppm_);

//boost::shared_ptr<ProjectorPairUsingMatrix> sptr_ppm_;

size_t segments, sinos, views, tangs;
segments = sptr_acq->get_num_segments();
sinos = sptr_acq->get_num_sinograms();
views = sptr_acq->get_num_views();
tangs = sptr_acq->get_num_tangential_poss();
std::cout << "segments: " << segments << '\n';
std::cout << "sinograms: " << sinos << '\n';
std::cout << "views: " << views << '\n';
std::cout << "tangential positions: " << tangs << '\n';
Coordinate3D<int> min_indices;
Coordinate3D<int> max_indices;
sptr_image->get_regular_range(min_indices, max_indices);
for (int i = 0; i < 3; i++)
	std::cout << max_indices[i + 1] - min_indices[i + 1] + 1 << '\n';

std::cout << "ok\n";
int dim[3];
cSTIR_getAcquisitionsDimensions(ptr_dt, (size_t)&dim[0]);
std::cout << dim[0] << ' ' << dim[1] << ' ' << dim[2] << '\n';

//std::cout << "ok\n";

//sptrImage3DF sptr_im = objectSptrFromHandle<Image3DF>(hi);

//*ptr_sptr = sptr_am->forward(*sptr_im, datafile);

//sptrAcqMod3DF sptr_am = objectSptrFromHandle<AcqMod3DF>(ha);

//Succeeded s = sptr_am->set_up(sptr_dt, sptr_im);

//sptrAcqMod3DF sptr_am = objectSptrFromHandle<AcqMod3DF>(ha);

//*ptr_sptr = sptr_am->forward(im, datafile);

//sptrAcqMod3DF sptr_am = objectSptrFromHandle<AcqMod3DF>(ha);
//sptrProjData sptr_ad = objectSptrFromHandle<ProjData>(hd);
//sptrImage3DF* ptr_sptr = new sptrImage3DF(sptr_am->backward(*sptr_ad));

//ProjDataInMemory adc(ad.get_exam_info_sptr(), ad.get_proj_data_info_sptr());
//adc.fill(ad);

//double* data_a = new double[size_a];
//double* data_b = new double[size_b];

//void add_(size_t n, double* u, double* v)

//void* cSTIR_acquisitionModelSetup
//	(void* ptr_am, const char* templ, void* ptr_im);
//void* cSTIR_acquisitionModelForward
//	(void* ptr_am, const char* datafile, void* ptr_dt, void* ptr_im);
//void* cSTIR_acquisitionModelBackward
//	(void* ptr_am, void* ptr_ad, void* ptr_im);

//extern "C"
//void* cSTIR_acquisitionModelSetup(void* ptr_am, const char* templ, void* ptr_im)
//{
//	try {
//		sptrProjData* ptr_sptr =
//			new sptrProjData(ProjData::read_from_file(templ));
//		sptrProjData& sptr_t = *ptr_sptr;
//		CAST_PTR(DataHandle, ha, ptr_am);
//		CAST_PTR(DataHandle, hi, ptr_im);
//		sptrProjPair& sptr_am = objectSptrFromHandle<ProjectorByBinPair>(ha);
//		sptrImage3DF& sptr_im = objectSptrFromHandle<Image3DF>(hi);
//		Succeeded s =
//			sptr_am->set_up(sptr_t->get_proj_data_info_sptr(), sptr_im);
//		DataHandle* handle = new DataHandle;
//		if (s != Succeeded::yes) {
//			ExecutionStatus status("cSTIR_acquisitionModelSetup failed",
//				__FILE__, __LINE__);
//			handle->set(0, &status);
//		}
//		else
//			handle->set((void*)ptr_sptr);
//		return (void*)handle;
//	}
//	CATCH
//}
//
//extern "C"
//void* cSTIR_acquisitionModelForward
//(void* ptr_am, const char* datafile, void* ptr_dt, void* ptr_im)
//{
//	try {
//		CAST_PTR(DataHandle, ha, ptr_am);
//		CAST_PTR(DataHandle, ht, ptr_dt);
//		CAST_PTR(DataHandle, hi, ptr_im);
//		sptrProjPair& sptr_am = objectSptrFromHandle<ProjectorByBinPair>(ha);
//		sptrProjData& sptr_dt = objectSptrFromHandle<ProjData>(ht);
//		sptrImage3DF& sptr_im = objectSptrFromHandle<Image3DF>(hi);
//		DataHandle* handle = new DataHandle;
//		if (strlen(datafile) < 1) {
//			NEW_SPTR(ProjData, ptr_sptr,
//				ProjDataInMemory(sptr_dt->get_exam_info_sptr(),
//				sptr_dt->get_proj_data_info_sptr()));
//			sptrProjData& sptr_t = *ptr_sptr;
//			sptr_am->get_forward_projector_sptr()->forward_project
//				(*sptr_t, *sptr_im);
//			handle->set((void*)ptr_sptr);
//		}
//		else {
//			NEW_SPTR(ProjData, ptr_sptr,
//				ProjDataInterfile(sptr_dt->get_exam_info_sptr(),
//				sptr_dt->get_proj_data_info_sptr(), datafile));
//			sptrProjData& sptr_t = *ptr_sptr;
//			sptr_am->get_forward_projector_sptr()->forward_project
//				(*sptr_t, *sptr_im);
//			handle->set((void*)ptr_sptr);
//		}
//		return (void*)handle;
//	}
//	CATCH
//}
//
//extern "C"
//void* cSTIR_acquisitionModelBackward(void* ptr_am, void* ptr_ad, void* ptr_im)
//{
//	try {
//		CAST_PTR(DataHandle, ha, ptr_am);
//		CAST_PTR(DataHandle, hd, ptr_ad);
//		CAST_PTR(DataHandle, hi, ptr_im);
//		sptrProjPair& sptr_am = objectSptrFromHandle<ProjectorByBinPair>(ha);
//		sptrProjData& sptr_ad = objectSptrFromHandle<ProjData>(hd);
//		Image3DF& image = objectFromHandle<Image3DF>(hi);
//		sptrImage3DF* ptr_sptr = new sptrImage3DF(image.clone());
//		sptrImage3DF& sptr_im = *ptr_sptr;
//		sptr_im->fill(0.0F);
//		sptr_am->get_back_projector_sptr()->back_project(*sptr_im, *sptr_ad);
//		return newObjectHandle(ptr_sptr);
//	}
//	CATCH
//}

//inline bool xSTIR_setupPrior(void* ptr)
//{
//	CAST_PTR(boost::shared_ptr< GeneralisedPrior<Image3DF> >, sptr_obj, ptr);
//	CAST_PTR(xSTIR_GeneralisedPrior3DF, obj, sptr_obj->get());
//	bool status = obj->post_process();
//	return status;
//}
//
//inline bool xSTIR_setupObjectiveFunction(void* ptr)
//{
//	CAST_PTR(boost::shared_ptr< GeneralisedObjectiveFunction<Image3DF> >,
//		sptr_obj, ptr);
//	CAST_PTR(xSTIR_GeneralisedObjectiveFunction3DF, obj, sptr_obj->get());
//	bool status = obj->post_process();
//	return status;
//}

//inline Succeeded xSTIR_setupReconstruction(void* ptr, sptrImage3DF const& image)
//{
//	CAST_PTR(boost::shared_ptr<xSTIR_IterativeReconstruction3DF>, sptr, ptr);
//	xSTIR_IterativeReconstruction3DF* recon = sptr->get();
//	// not needed - default is non-zero string ("1") anyway
//	//recon->set_initial_estimate_file("dummy.hv");
//	Succeeded s = Succeeded::no;
//	if (recon->post_process())
//		return s;
//	s = recon->setup(image);
//	recon->subiteration() = recon->get_start_subiteration_num();
//	return s;
//}
//
//inline void xSTIR_updateReconstruction(void* ptr, Image3DF& image) 
//{
//	CAST_PTR(boost::shared_ptr<xSTIR_IterativeReconstruction3DF>, sptr, ptr);
//	xSTIR_IterativeReconstruction3DF* recon = sptr->get();
//	recon->update(image);
//}

//inline int& xSTIR_subiteration(void* ptr) 
//{
//	CAST_PTR(xSTIR_IterativeReconstruction3DF, recon, ptr);
//	return recon->subiteration();
//}
//
//inline void xSTIR_set_initial_estimate_file(void* ptr, const char* filename) 
//{
//	CAST_PTR(xSTIR_IterativeReconstruction3DF, recon, ptr);
//	recon->set_initial_estimate_file(filename);
//}

//			< ObjectiveFunction3DF, PoissonLogLhLinModMeanProjData3DF >();

//CAST_PTR(DataHandle, ho, ptr_obj);

//status = xSTIR_setupObjectiveFunction(ho->data());

//CAST_PTR(DataHandle, hp, ptr_p);
//CAST_PTR(DataHandle, hi, ptr_i);

//CAST_PTR(DataHandle, ha, ptr_am);
//CAST_PTR(DataHandle, ht, ptr_dt);
//CAST_PTR(DataHandle, hi, ptr_im);

//CAST_PTR(DataHandle, ha, ptr_am);
//CAST_PTR(DataHandle, hi, ptr_im);

//CAST_PTR(DataHandle, ha, ptr_am);
//CAST_PTR(DataHandle, hd, ptr_ad);

//sptrProjData& sptr_t =
//	objectSptrFromHandle<ProjData>((DataHandle*)ptr_t);

//sptrProjData& sptr_ad = objectSptrFromHandle<ProjData>((DataHandle*)ptr_acq);

//sptrProjData& sptr_ad = objectSptrFromHandle<ProjData>((DataHandle*)ptr_acq);

//sptrProjData& sptr_ad = objectSptrFromHandle<ProjData>((DataHandle*)ptr_acq);

//sptrProjData& sptr_ad = objectSptrFromHandle<ProjData>((DataHandle*)ptr_acq);

//sptrProjData& sptr_from = objectSptrFromHandle<ProjData>((DataHandle*)ptr_from);

//sptrProjData& sptr_ad = objectSptrFromHandle<ProjData>((DataHandle*)ptr_acq);

//CAST_PTR(DataHandle, hr, ptr_r);
//CAST_PTR(DataHandle, hi, ptr_i);

//Succeeded s = xSTIR_setupReconstruction(hr->data(), sptr_image);

//CAST_PTR(DataHandle, hr, ptr_r);
//CAST_PTR(DataHandle, hi, ptr_i);

//CAST_PTR(DataHandle, hr, ptr_r);
//CAST_PTR(DataHandle, hi, ptr_i);

//xSTIR_updateReconstruction(hr->data(), image);

//CAST_PTR(DataHandle, hf, ptr_f);
//CAST_PTR(DataHandle, hi, ptr_i);

//CAST_PTR(DataHandle, hf, ptr_f);
//CAST_PTR(DataHandle, hi, ptr_i);

//CAST_PTR(DataHandle, hv, ptr_v);

//CAST_PTR(DataHandle, hi, ptr_i);

//sptrProjData& sptr_ad = objectSptrFromHandle<ProjData>((DataHandle*)ptr_ad);

//CAST_PTR(DataHandle, hi, ptr_i);
//CAST_PTR(DataHandle, hv, ptr_v);
//CAST_PTR(DataHandle, hs, ptr_s);

//Image3DF* ptr_image = objectPtrFromHandle<Image3DF>((DataHandle*)ptr_i);

//Image3DF* ptr_image = objectPtrFromHandle<Image3DF>((DataHandle*)ptr_im);

//Image3DF* ptr_image = objectPtrFromHandle<Image3DF>((DataHandle*)ptr_im);

//Image3DF* ptr_image = objectPtrFromHandle<Image3DF>((DataHandle*)ptr_im);

//CAST_PTR(DataHandle, first_h, first);
//CAST_PTR(DataHandle, second_h, second);

//PoissonLogLhLinModMeanProjData3DF& obj_fun = objectFromHandle
//	< ObjectiveFunction3DF, PoissonLogLhLinModMeanProjData3DF >(hp);

//{	
//AcqMod3DF& am = objectFromHandle<AcqMod3DF>(hv);
//obj_fun.set_projector_pair_sptr(am.projectors_sptr());
//if (am.additive_term_sptr().get())
//	obj_fun.set_additive_proj_data_sptr(am.additive_term_sptr());
//if (am.normalisation_sptr().get())
//	obj_fun.set_normalisation_sptr(am.normalisation_sptr());
//}

//xSTIR_set_initial_estimate_file(&recon, charDataFromHandle(hv));

//xSTIR_subiteration(&recon) = value;

//return intDataHandle(xSTIR_subiteration(&recon));

