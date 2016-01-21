#include <iostream>
#include <string>
using namespace std;

#include <conio.h>

int c2m(
	const string& path, 
	const string& chfile, 
	const string& mhfile, 
	const string& mcfile);

int main(int argc, char **argv)
{
	string path = "../../src/";
	int status = c2m(path, "cgadgetron.h", "mgadgetron.h", "mgadgetron.c");
	if (status)
		cout << "wrong input file format" << endl;
	cout << "Press any key to continue ";
	_getch();
	return 0;
}

