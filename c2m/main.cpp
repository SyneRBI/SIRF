#include <iostream>
#include <string>
using namespace std;

#include <conio.h>

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
	path = "../";
	status = c2m\
		("IUTILITIES", "---", path, "iutilities.h", "mutilities.h", "mutilities.c");
	if (status)
		cout << "wrong input file format" << endl;

	path = "../../xGadgetron/cGadgetron/";
	status = c2m\
		("CGADGETRON", "cGT", path, "cgadgetron.h", "mgadgetron.h", "mgadgetron.c");
	if (status)
		cout << "wrong input file format" << endl;

	path = "../../xSTIR/cSTIR/";
	status = c2m("CSTIR", "cSTIR", path, "cstir.h", "mstir.h", "mstir.c", 1);
	if (status)
		cout << "wrong input file format" << endl;

	cout << "Press any key to continue ";
	_getch();
	return 0;
}

