/*
SyneRBI Synergistic Image Reconstruction Framework (SIRF)
Copyright 2015 - 2020 Rutherford Appleton Laboratory STFC

This is software developed for the Collaborative Computational
Project in Synergistic Reconstruction for Biomedical Imaging (formerly CCP PETMR)
(http://www.ccpsynerbi.ac.uk/).

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
	LocalisedException(const std::string& reason, const std::string& file, int line)
        : reason_(reason), file_(file), line_(line)
    { }
	virtual ~LocalisedException() throw()
    { }
	virtual const char* what() const throw()
	{
		return reason_.c_str();
	}
	const std::string& file() const throw()
	{
		return file_;
	}
	int line() const throw() {
		return line_;
	}
private:
	const std::string reason_;
	const std::string file_;
	int line_;
};

#endif
