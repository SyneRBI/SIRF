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

#ifndef LOCALISED_EXCEPTION
#define LOCALISED_EXCEPTION

#include <string.h>

#include <exception>
#include <iostream>

#define THROW(msg) throw LocalisedException(msg, __FILE__, __LINE__)
#define ASSERT(condition, msg) if (!(condition)) THROW(msg)

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