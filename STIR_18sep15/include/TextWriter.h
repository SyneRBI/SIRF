#ifndef TEXT_WRITER_TYPES
#define TEXT_WRITER_TYPES

#include <fstream>
#include <iostream>
#include <string>
using namespace std;

#define DEFAULT_STREAM std::cerr

class aTextWriter {
public:
	virtual ~aTextWriter() {}
	virtual void write(const char* text) const = 0;
};

class TextPrinter : public aTextWriter {
public:
	TextPrinter(const char* s = 0) : _stream(0) {
		if (s) {
			if (strcmp(s, "stdout") == 0 || strcmp(s, "cout") == 0)
				_stream = 1;
			else if (strcmp(s, "stderr") == 0 || strcmp(s, "cerr") == 0)
				_stream = 2;
		}
	}
	virtual void write(const char* text) const {
		switch (_stream) {
		case 1:
			cout << text;
			break;
		case 2:
			cerr << text;
			break;
		default:
			DEFAULT_STREAM << text;
		}
	}
private:
	int _stream;
};

class TextWriter : public aTextWriter {
public:
	ostream* out;
	TextWriter(ostream* os = 0) : out(os) {}
	virtual void write(const char* text) const {
		if (out) {
			(*out) << text;
			(*out).flush();
		}
		else
			DEFAULT_STREAM << text;
	}
};

class TextWriterHandle {
public:
	//void set_writer(aTextWriter * my_writer = 0);
	void write(const char* text) const {
		//static bool first_call = true;
		//if (first_call)
		//{
		//	if (!writer)
		//	    this->set_writer(0);
		//	first_call = false;
		//}
		if (writer)
			writer->write(text);
		else
			DEFAULT_STREAM << text;
	}
//private:
	static aTextWriter* writer;
};

class TextOutputHandler {
public:
	void set_output
		(aTextWriter* info = 0, aTextWriter* warn = 0, aTextWriter* errr = 0) {
		init_();
		info_ = info;
		warn_ = warn;
		errr_ = errr;
	}
	void print_info(const char* text) {
		init_();
		if (info_)
			info_->write(text);
	}
private:
	static aTextWriter* info_;
	static aTextWriter* warn_;
	static aTextWriter* errr_;
	void init_() {
		static bool initialized = false;
		if (!initialized) {
			info_ = 0;
			warn_ = 0;
			errr_ = 0;
			initialized = true;
		}
	}
};

void writeText(const char* text);

class AltCout {
public:
	AltCout& operator<<(int i) {
		char buff[32];
		sprintf_s(buff, 32, "%d", i);
		//printf("%s", buff);
		writeText(buff);
		return *this;
	}
	AltCout& operator<<(double x) {
		char buff[32];
		sprintf_s(buff, 32, "%f", x);
		//printf("%s", buff);
		writeText(buff);
		return *this;
	}
	AltCout& operator<<(char a) {
		//printf("%c", a);
		char buff[2];
		buff[0] = a;
		buff[1] = '\0';
		writeText(buff);
		return *this;
	}
	AltCout& operator<<(const char* s) {
		//printf("%s", s);
		writeText(s);
		return *this;
	}
	AltCout& operator<<(string s) {
		//printf("%s", s);
		writeText(s.c_str());
		return *this;
	}
};

class Stir {
public:
	static char endl;
	static AltCout cout;
	static AltCout cerr;
};

#endif