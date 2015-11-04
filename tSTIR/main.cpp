#include <conio.h>
#include <iostream>
#include <string>
using namespace std;

void stir_test0();
void stir_test1();
void stir_test2();
void cstir_test1();
void cstir_test2();
void cstir_test3();
int c2m(const string& path, const string& chfile, const string& mhfile, const string& mcfile);

#include "try_static.h"

void set_i(int i) {
	A a;
	a.set_i(i);
}

int get_i() {
	A a;
	return a.get_i();
}

void print_i() {
	A a;
	a.print_i();
}

void test_static() {
	cout << get_i() << endl;
	print_i();
	set_i(100);
	print_i();
}

int main(int argc, char **argv)
{
	//test_static();
	//stir_test1();
	cstir_test1();
	//string path = "C:/Users/wps46139/Codes/PETMR/xSTIR/cSTIR/";
	//int status = c2m(path, "cstir.h", "mstir.h", "mstir.c");
	//if (status)
	//	cout << "wrong input file format" << endl;
	cout << "Press any key to continue";
	_getch();
	return 0;
}

