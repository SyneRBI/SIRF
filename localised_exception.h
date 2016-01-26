#ifndef LOCALISED_EXCEPTION
#define LOCALISED_EXCEPTION

#include <string.h>

#include <exception>
#include <iostream>

class LocalisedException : public std::exception {
public:
	LocalisedException(const char* reason, const char* file, int line) {
		size_t len = strlen(reason) + 1;
		reason_ = new char[len];
		memcpy(reason_, reason, len);
		len = strlen(file) + 1;
		file_ = new char[len];
		memcpy(file_, file, len);
		line_ = line;
	}
	virtual ~LocalisedException() throw() {
		delete[] reason_;
		delete[] file_;
	}
	virtual const char* what() const throw()
	{
		return reason_;
	}
	const char* file() const throw()
	{
		return file_;
	}
	int line() const throw() {
		return line_;
	}
private:
	char* reason_;
	char* file_;
	int line_;
};

#endif