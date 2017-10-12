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

#ifndef GET_ENVIRONMENT_VARIABLES
#define GET_ENVIRONMENT_VARIABLES

#ifdef _MSC_VER
#if _MSC_VER >= 1900
using namespace System;
static string toStandardString(System::String^ var)
{
	using System::Runtime::InteropServices::Marshal;
	System::IntPtr pointer = Marshal::StringToHGlobalAnsi(var);
	char* charPointer = reinterpret_cast<char*>(pointer.ToPointer());
	string returnString(charPointer, var->Length);
	Marshal::FreeHGlobal(pointer);
	return returnString;
}
string EnvironmentVariable(const char* name)
{
	try {
		System::String^ var = gcnew String(name);
		return toStandardString(Environment::GetEnvironmentVariable(var));
	}
	catch (...) {
		return "";
	}
}
#else
#include <tchar.h>
#include <Windows.h>
#define BUFSIZE 4096
string EnvironmentVariable(const char* name)
{
	DWORD dwRet;
	LPTSTR ptr_value;
	string value;
	ptr_value = (LPTSTR)malloc(BUFSIZE*sizeof(TCHAR));
	dwRet = GetEnvironmentVariable(TEXT(name), ptr_value, BUFSIZE);
	if (dwRet)
		value = string(ptr_value);
	else
		value = string("");
	free(ptr_value);
	return value;
}
#endif
#else
#include <cstdlib>
string EnvironmentVariable(const char* name)
{
	const char* val = ::getenv(name);
	if (val)
		return val;
	else
		return "";
}
#endif

#endif
