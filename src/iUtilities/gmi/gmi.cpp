#include <stdio.h>

#include <iostream>

using namespace std;

#include "SIRF/common/envar.h"

int generate_matlab_interface(
	const char* library,
	const char* prefix,
	const string& path_in,
	const string& chfile,
	const string& path_out,
	const string& mhfile,
	const string& mcfile,
	int wp = 0);

int main(int argc, char **argv)
{
	int status;
	string path_in;
	string path_out;
	string SIRF_path;
	string name("SIRF_PATH");
	SIRF_path = EnvironmentVariable(name.c_str());
	if (SIRF_path.length() < 1) {
		cout << name.c_str() << " not defined" << endl;
		return 1;
	}
	path_in = SIRF_path + "/src/iUtilities/";
	path_out = path_in;
	status = generate_matlab_interface\
		("IUTILITIES", "---", path_in, "iutilities.h", path_out, "miutilities.h", "miutilities.c");
	if (status)
		cout << "wrong input file format" << endl;
}