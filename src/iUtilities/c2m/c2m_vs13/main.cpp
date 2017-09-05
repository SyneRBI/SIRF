#include <iostream>
#include <string>
using namespace std;

#ifdef _MSC_VER
#include <conio.h>
#endif

#include <tchar.h>
#include <Windows.h>

#define BUFSIZE 4096

string EnvironmentVariable(const char* name)
{
	DWORD dwRet, dwErr;
	LPTSTR ptr_value;
	string value;
	ptr_value = (LPTSTR)malloc(BUFSIZE*sizeof(TCHAR));
	dwRet = GetEnvironmentVariable(TEXT(name), ptr_value, BUFSIZE);
	if (dwRet)
		value = string(ptr_value);
	else
		value = string("");
	free(ptr_value);
	return value;
}

int c2m(
	const char* library,
	const char* prefix,
	const string& path,
	const string& chfile, 
	const string& mhfile, 
	const string& mcfile,
	int wp = 0);

int main(int argc, char **argv)
{
	int status;
	string path;

	string SIRF_path;
	SIRF_path = EnvironmentVariable("SIRF_PATH");
	if (SIRF_path.length() < 1) {
		cout << "SIRF_PATH not defined" << endl;
		return 1;
	}

	if (argc < 2 || argv[1][0] == 'u') {
		//path = "../";
		path = SIRF_path + "/src/iUtilities/";
		status = c2m\
			("IUTILITIES", "---", path, "iutilities.h", "miutilities.h", "miutilities.c");
		if (status)
			cout << "wrong input file format" << endl;
	}

	if (argc < 2 || argv[1][0] == 'g') {
		path = SIRF_path + "/src/xGadgetron/cGadgetron/";
		//path = "../../xGadgetron/cGadgetron/";
		status = c2m\
			("CGADGETRON", "cGT", path, "cgadgetron.h", "mgadgetron.h", "mgadgetron.c");
		if (status)
			cout << "wrong input file format" << endl;
	}

	if (argc < 2 || argv[1][0] == 's') {
		path = SIRF_path + "/src/xSTIR/cSTIR/";
		//path = "../../xSTIR/cSTIR/";
		status = c2m("CSTIR", "cSTIR", path, "cstir.h", "mstir.h", "mstir.c", 1);
		if (status)
			cout << "wrong input file format" << endl;
	}

#ifdef _MSC_VER
	cout << "Press any key to continue ";
	_getch();
#endif
	return 0;
}

