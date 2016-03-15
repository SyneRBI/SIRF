#include <conio.h>
#include <fstream>
#include <iostream>
#include <string>
using namespace std;

int convert_h(const string& in, string& out) {
	out = "EXPORTED_FUNCTION ";
	size_t i = in.find("cGT");
	if (i != string::npos) {
		out += in;
		out[i + 18] = 'm';
	}
	else {
		size_t j = in.find("(", 0);
		if (j == string::npos)
			return 1;
		j--;
		while (isspace(in[j]))
			j--;
		i = j - 1;
		while (!isspace(in[i]) && in[i] != '*')
			i--;
		i++;
		out += in.substr(0, i);
		out += "m";
		out += toupper(in[i]);
		out += in.substr(i + 1);
	}
	return 0;
}

int convert_c(const string& in, string& out) {

	int p, q;
	size_t i, j, k, l;

	out = "EXPORTED_FUNCTION ";

	k = in.find(";");
	if (k == string::npos)
		return 3;

	i = in.find("cGT");
	if (i == string::npos) {
		j = in.find("(", 0);
		if (j == string::npos)
			return 2;
		l = j - 1;
		while (isspace(in[l]))
			l--;
		i = l - 1;
		while (!isspace(in[i]) && in[i] != '*')
			i--;
		i++;
		out += in.substr(0, i);
		out += "m";
		out += toupper(in[i]);
		out += in.substr(i + 1, k - i - 1);
	}
	else {
		j = in.find("(", i);
		if (j == string::npos)
			return 2;
		out += in.substr(0, k);
		out[i + 18] = 'm';
	}
	out += " {\n\t";

	q = (int)i - 1;
	while (isspace(in[q]))
		q--;
	p = q - 1;
	while (p >= 0 && !isspace(in[p]))
		p--;
	p++;

	if (in.substr(p, q - p + 1).compare("void") != 0)
		out += "return ";

	out += in.substr(i, j - i + 1);

	for (;;) {
		k = in.find(',', j);
		if (k == string::npos)
			k = in.find(')', j);
		if (k == string::npos)
			return 4;
		j = k - 1;
		while (isspace(in[j]))
			j--;
		i = j;
		while (!isspace(in[i]) && in[i] != '*' && in[i] != '(')
			i--;
		i++;
		l = j - i + 1;
		out += in.substr(i, l);
		out += in[k];
		if (in[k] == ')') {
			out += ";";
			break;
		}
		else
			out += " ";
		j = k + 1;
	}

	out += "\n}";
	return 0;

}

int print_hash(int ifcount, int hashcount)
{
	if (hashcount < 2 || ifcount < 1)
		return 0;
	return 1;
}

int c2m(
	const string& path, 
	const string& chfile, 
	const string& mhfile, 
	const string& mcfile,
	int with_print = 0)
{

	ifstream fin;
	//ofstream fout;
	ofstream fh;
	ofstream fc;
	string in;
	string out;
	string line;
	int hashcount;
	int ifcount;
	int status;
	bool quit;

	quit = false;

	fin.open(path + chfile, ios::in);
	if (!fin) {
		cout << "file " << path + chfile << " not found" << endl;
		return -1;
	}

	size_t i = string::npos;
	size_t m;

	fh.open(mhfile, ios::out);
	fc.open(mcfile, ios::out);

	fh << "#ifndef CGADGETRON_TO_MATLAB_INTERFACE" << endl;
	fh << "#define CGADGETRON_TO_MATLAB_INTERFACE" << endl << endl;
	fh << "#define CGADGETRON_FOR_MATLAB" << endl;
	fh << "#include \"shrhelp.h\"" << endl << endl;

	fc << "#include <mex.h>" << endl;
	fc << "#define EXPORT_FCNS" << endl;
	fc << "#define CGADGETRON_FOR_MATLAB" << endl;
	fc << "#include \"matrix.h\"" << endl;
	fc << "#include \"shrhelp.h\"" << endl;
	fc << "#include \"" << chfile << '"' << endl << endl;

	hashcount = 0;
	ifcount = 0;

	for (;;) {
		if (fin.eof())
			break;
		getline(fin, line);
		cout << line << endl;
		i = line.find_first_not_of(" \t\n\v\f\r");
		if (i == string::npos || line[i] == '/' && line[i + 1] == '/')
			continue;
		if (line[i] == '#') {
			if (line[i + 1] == 'e' && line[i + 2] == 'n' && line[i + 3] == 'd')
				ifcount--;
			if (print_hash(ifcount, hashcount)) {
				fh << line << endl;
				fc << line << endl;
			}
			hashcount++;
			if (line[i + 1] == 'i' && line[i + 2] == 'f')
				ifcount++;
			continue;
		}
		in = line;
		i = in.find(';');
		while (i == string::npos) {
			if (fin.eof())
				break;
			getline(fin, line);
			cout << line << endl;
			m = line.find("//");
			if (m != string::npos)
				line.erase(m);
			i = line.find_first_not_of(" \t\n\v\f\r");
			if (i != string::npos && line[i] == '#') {
				fh << in << endl;
				fc << in << endl;
				if (line[i + 1] == 'e' && line[i + 2] == 'n' && line[i + 3] == 'd')
					ifcount--;
				if (print_hash(ifcount, hashcount)) {
					fh << line << endl;
					fc << line << endl;
				}
				hashcount++;
				if (line[i + 1] == 'i' && line[i + 2] == 'f')
					ifcount++;
				i = string::npos;
				in.clear();
				continue;
			}
			if (i != string::npos)
				in += line.substr(i);
			i = line.find(';');
		}
		if (i == string::npos)
			break;
		//cout << in << endl;
		status = convert_h(in, out);
		if (status) {
			cout << in << endl;
			fin.close();
			fh.close();
			fc.close();
			return status;
		}
		fh << out << endl;
		status = convert_c(in, out);
		if (status) {
			cout << in << endl;
			fin.close();
			fh.close();
			fc.close();
			return status;
		}
		fc << out << endl;
		if (fin.eof())
			break;
		in = line.substr(i + 1);
		i = in.find_first_not_of(" \t\n\v\f\r");
		if (i == string::npos)
			in.clear();
	}

	fh << endl;
	if (with_print) {
		fh << "EXPORTED_FUNCTION void* mNewMexPrinter();" << endl;
		fh << "EXPORTED_FUNCTION void mDeleteMexPrinter(void* ptr);" << endl;
	}
	fh << "#endif" << endl;

	fc << endl;
	if (with_print) {
		fc << "EXPORTED_FUNCTION void* mNewMexPrinter() {" << endl;
		fc << "  return newMexPrinter();" << endl;
		fc << "}" << endl;
		fc << "EXPORTED_FUNCTION void mDeleteMexPrinter(void* ptr) {" << endl;
		fc << "  deleteMexPrinter(ptr);" << endl;
		fc << "}" << endl;
	}
	fc <<
		"void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])";
	fc << " {}" << endl;

	fin.close();
	fh.close();
	fc.close();
	return 0;
}