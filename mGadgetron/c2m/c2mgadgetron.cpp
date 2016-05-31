#include <iostream>
#include <string>
using namespace std;

#include <conio.h>

int c2m(
	const string& path, 
	const string& chfile, 
	const string& mhfile, 
	const string& mcfile,
	int wp = 0);

int main(int argc, char **argv)
{
	string path = "../../cGadgetron/";
	int status = c2m(path, "cgadgetron.h", "mgadgetron.h", "mgadgetron.c");
	if (status)
		cout << "wrong input file format" << endl;
	path = "../../../iUtilities/";
	status = c2m(path, "iutilities.h", "mutilities.h", "mutilities.c");
	if (status)
		cout << "wrong input file format" << endl;
	cout << "Press any key to continue ";
	_getch();
	return 0;
}

