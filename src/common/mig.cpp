/*
CCP PETMR Synergistic Image Reconstruction Framework (SIRF)
Copyright 2015 - 2017 Rutherford Appleton Laboratory STFC

This is software developed for the Collaborative Computational
Project in Positron Emission Tomography and Magnetic Resonance imaging
(http://www.ccppetmr.ac.uk/).

Licensed under the Apache License, Version 2.0 (the "License");
you may not use this file except in compliance with the License.
You may obtain a copy of the License at
http://www.apache.org/licenses/LICENSE-2.0
Unless required by applicable law or agreed to in writing, software
distributed under the License is distributed on an "AS IS" BASIS,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
See the License for the specific language governing permissions and
limitations under the License.

*/

#include <fstream>
#include <iostream>
#include <string>
using namespace std;

int convert_h(const char* prefix, const string& in, string& out) {
	out = "EXPORTED_FUNCTION ";
	size_t i = in.find(prefix);
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

int convert_c(const char* prefix, const string& in, string& out) {

	int p, q;
	size_t i, j, k, l;

	out = "EXPORTED_FUNCTION ";

	k = in.find(";");
	if (k == string::npos)
		return 3;

	i = in.find(prefix);
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

int generate_matlab_interface(
	const char* library,
	const char* prefix,
	const string& path_in, 
	const string& chfile, 
	const string& path_out,
	const string& mhfile,
	const string& mcfile,
	int with_print = 0,
	bool to_cout = false)
{

	ifstream fin;
	ofstream fh;
	ofstream fc;
	string in;
	string out;
	string line;
	int comment;
	int hashcount;
	int ifcount;
	int head_printed;
	int status;
	bool quit;

	quit = false;

	fin.open((path_in + chfile).c_str(), ios::in);
	if (!fin) {
		cout << "file " << path_in + chfile << " not found" << endl;
		return -1;
	}

	size_t i = string::npos;
	size_t m;

	fh.open((path_out + mhfile).c_str(), ios::out);
	fc.open((path_out + mcfile).c_str(), ios::out);

	head_printed = 0;
	comment = 0;
	hashcount = 0;
	ifcount = 0;

	for (;;) {

		if (fin.eof())
			break;

		getline(fin, line);
		if (to_cout)
			cout << line << endl;
		i = line.find_first_not_of(" \t\n\v\f\r");

		if (i == string::npos || line[i] == '/' && line[i + 1] == '/')
			continue;

		if (line[i] == '/' && line[i + 1] == '*')
			comment = 1;
		if (comment) {
			fh << line << endl;
			fc << line << endl;
		}
		if (line[i] == '*' && line[i + 1] == '/') {
			line.erase(0);
			comment = 0;
		}
		if (comment)
			continue;

		if (!comment && !head_printed) {
			fh << "#ifndef " << library << "_TO_MATLAB_INTERFACE" << endl;
			fh << "#define " << library << "_TO_MATLAB_INTERFACE" << endl << endl;
			fh << "#define " << library << "_FOR_MATLAB" << endl;
			fh << "#ifdef _WIN32" << endl;
			fh << "#define EXPORTED_FUNCTION __declspec(dllexport)" << endl;
			fh << "#else" << endl;
			fh << "#define EXPORTED_FUNCTION" << endl;
			fh << "#endif" << endl << endl;
			//fh << "#include \"shrhelp.h\"" << endl << endl;
			//fc << "#define EXPORT_FCNS" << endl;
			fc << "#define " << library << "_FOR_MATLAB" << endl;
			fc << "#ifdef _WIN32" << endl;
			fc << "#define EXPORTED_FUNCTION __declspec(dllexport)" << endl;
			fc << "#else" << endl;
			fc << "#define EXPORTED_FUNCTION" << endl;
			fc << "#endif" << endl << endl;
			fc << "#include <mex.h>" << endl;
			fc << "#include \"matrix.h\"" << endl;
			//fc << "#include \"shrhelp.h\"" << endl;
			fc << "#include \"" << chfile << '"' << endl << endl;
			head_printed = 1;
		}

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
			if (to_cout)
				cout << line << endl;
			m = line.find("//");
			if (m != string::npos)
				line.erase(m);
			i = line.find_first_not_of(" \t\n\v\f\r");
			if (i != string::npos && line[i] == '#') {
				if (in.length()) {
					fh << in << endl;
					fc << in << endl;
				}
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
			if (i != string::npos) {
				in += " ";
				in += line.substr(i);
			}
			i = line.find(';');
		}
		if (i == string::npos)
			break;
		//cout << in << endl;
		status = convert_h(prefix, in, out);
		if (status) {
			if (to_cout)
				cout << in << endl;
			fin.close();
			fh.close();
			fc.close();
			return status;
		}
		fh << out << endl;
		status = convert_c(prefix, in, out);
		if (status) {
			if (to_cout)
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

	if (with_print) {
		fh << "EXPORTED_FUNCTION void* mNewMexPrinter();" << endl;
		fh << "EXPORTED_FUNCTION void* mDeleteMexPrinter(void* ptr);" << endl;
	}
	fh << endl;
	fh << "#endif" << endl;

	if (with_print) {
		fc << "void* newMexPrinter();" << endl;
		fc << "void* deleteMexPrinter(void* ptr);" << endl;
		fc << "EXPORTED_FUNCTION void* mNewMexPrinter() {" << endl;
		fc << "  return newMexPrinter();" << endl;
		fc << "}" << endl;
		fc << "EXPORTED_FUNCTION void* mDeleteMexPrinter(void* ptr) {" << endl;
		fc << "  return deleteMexPrinter(ptr);" << endl;
		fc << "}" << endl;
	}
	fc << endl;
	fc <<
		"void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])";
	fc << " {}" << endl;

	fin.close();
	fh.close();
	fc.close();
	return 0;
}
