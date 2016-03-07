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

int c2m(
	const string& path, 
	const string& chfile, 
	const string& mhfile, 
	const string& mcfile,
	int with_print = 0)
{

	ifstream fin;
	ofstream fout;
	ofstream fh;
	ofstream fc;
	string in;
	string out;
	string line;
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
	while (fin) {
		getline(fin, line);
		//cout << line << endl;
		m = line.find("//");
		if (m != string::npos)
			line.erase(m);
		i = line.find('{');
		if (i != string::npos)
			break;
	}

	if (i == string::npos) {
		fin.close();
		return -2;
	}

	fout.open(chfile, ios::out);
	fh.open(mhfile, ios::out);
	fc.open(mcfile, ios::out);

	fout << "#ifndef CGADGETRON_INTERFACE" << endl;
	fout << "#define CGADGETRON_INTERFACE" << endl << endl;

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
	//fc << "#include \"cgadgetron.h\"" << endl << endl;

	in = line.substr(i + 1);
	i = in.find_first_not_of(" \t\n\v\f\r");
	if (i == string::npos)
		in.clear();
	else
		fout << in << endl;

	for (;;) {
		i = in.find('}');
		if (i != string::npos)
			break;
		i = in.find(';');
		while (i == string::npos) {
			getline(fin, line);
			//cout << line << endl;
			if (fin.eof())
				break;
			i = in.find('}');
			if (i != string::npos) {
				quit = true;
				break;
			}
			m = line.find("//");
			if (m != string::npos)
				line.erase(m);
			i = line.find_first_not_of(" \t\n\v\f\r");
			if (i != string::npos && line[i] == '#') {
				fout << line << endl;
				fh << line << endl;
				fc << line << endl;
				i = string::npos;
				break;
			}
			if (i != string::npos)
				in += line.substr(i);
			i = line.find(';');
		}
		if (quit)
			break;
		if (i == string::npos) {
			if (fin.eof())
				break;
			in.clear();
			continue;
		}
		//cout << in << endl;
		fout << in << endl;
		if (i != string::npos) {
			status = convert_h(in, out);
			if (status) {
				cout << in << endl;
				fin.close();
				fout.close();
				fh.close();
				fc.close();
				return status;
			}
			fh << out << endl;
			status = convert_c(in, out);
			if (status) {
				cout << in << endl;
				fin.close();
				fout.close();
				fh.close();
				fc.close();
				return status;
			}
			fc << out << endl;
		}
		if (fin.eof())
			break;
		in = line.substr(i + 1);
		i = in.find_first_not_of(" \t\n\v\f\r");
		if (i == string::npos)
			in.clear();
	}

	if (with_print) {
		fout << "void* newMexPrinter();" << endl;
		fout << "void deleteMexPrinter(void* ptr);" << endl;
	}
	fout << endl << "#endif" << endl;

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
	fout.close();
	fh.close();
	fc.close();
	return 0;
}