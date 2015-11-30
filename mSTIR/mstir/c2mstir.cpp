#include <conio.h>
#include <iostream>
#include <string>
using namespace std;

int c2m(
	const string& path, 
	const string& chfile, 
	const string& mhfile, 
	const string& mcfile);

int main(int argc, char **argv)
{
	string path = "C:/Users/wps46139/Codes/PETMR/xSTIR/cSTIR/";
	int status = c2m(path, "cstir.h", "mstir.h", "mstir.c");
	if (status)
		cout << "wrong input file format" << endl;
	//cout << "Press any key to continue";
	//_getch();
	return 0;
}

