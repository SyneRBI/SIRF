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

/*!
\file
\ingroup Registration
\brief NiftyReg's aladin class for rigid and affine registrations.

\author Richard Brown
\author CCP PETMR
*/

#ifndef _SIRFRegParserKEY_H_
#define _SIRFRegParserKEY_H_

#include <string>
#include <iostream>
#include <algorithm>

namespace sirf {
/// Get template type
template<class A>
std::string get_typeid(A)
{
    if      (typeid(A) == typeid(const char *))  return "const char *";
    else if (typeid(A) == typeid(int))           return "int";
    else if (typeid(A) == typeid(float))         return "float";
    else if (typeid(A) == typeid(double))        return "double";
    else if (typeid(A) == typeid(unsigned int))  return "unsigned int";
    else if (typeid(A) == typeid(bool))          return "bool";
    else                                         return "unknown";
}

/// Abstract base class for a parser key
template<class Z>
class SIRFRegParserKeyBase
{
public:

    /// Set object
    void set_object(std::shared_ptr<Z> object) { _object = object; }

    /// Get arguments
    virtual void get_arguments(std::string &line) = 0;

    /// Call function
    virtual void call_function() = 0;

    /// Print number of arguments expected
    virtual void print_num_arguments_expected() = 0;

protected:

    /// Check that the number of arguments matches the expected amount
    void check_num_arguments(std::string &line, int num_args_expected)
    {
        int num_args;

        // If last character is "=", then num args == 0
        if (line.back() == '=') num_args = 0;

        // If not... Parameters are separated by commas. so num parameters is num commas + 1
        else {
            num_args = int(std::count(line.begin(), line.end(), ',')) + 1;
        }

        if (num_args != num_args_expected) {
            throw std::runtime_error("Error. Num arguments expected: " + std::to_string(num_args_expected) + ", num arguments found: " + std::to_string(num_args) + ".");
        }

        this->print_num_arguments_expected();
        std::cout << "\tNumber arguments found: "    << num_args          << "\n";
    }


    /// Get the nth argument as a string
    std::string get_arg_as_string(std::string line, int arg_num)
    {
        // Start index is always to right of the ":="
        // End index is up to ","
        int index_start = line.find(":=") + 2;
        int index_end   = line.find(",");

        // For subsequent arguments
        for (int i=0; i<arg_num; i++) {
            index_start = line.find(",",index_start) + 1; // add 1 to go to right of ","
            index_end   = line.find(",",index_start);
        }

        // Error check
        if (index_end == index_start) {
            throw std::runtime_error("Error. Argument number " + std::to_string(arg_num) + " was not found.");
        }

        std::string arg_as_string = line.substr(index_start, index_end-index_start);
        std::cout << "\t\tArg as string: " << arg_as_string << "\n";

        // Return the substring
        return arg_as_string;
    }

    /// Get argument - const char *
    void get_argument(std::string line, int arg_num, const char *&arg  ) { arg = get_arg_as_string(line, arg_num).c_str();     }

    /// Get argument - bool
    void get_argument(std::string line, int arg_num, bool &arg         ) { arg = std::stoi(get_arg_as_string(line, arg_num));  }

    /// Get argument - int
    void get_argument(std::string line, int arg_num, int &arg          ) { arg = std::stoi(get_arg_as_string(line, arg_num));  }

    /// Get argument - unsigned int
    void get_argument(std::string line, int arg_num, unsigned int &arg ) { arg = std::stoll(get_arg_as_string(line, arg_num)); }

    /// Get argument - float
    void get_argument(std::string line, int arg_num, float &arg        ) { arg = std::stof(get_arg_as_string(line, arg_num));  }

    /// Get argument - double
    void get_argument(std::string line, int arg_num, double &arg       ) { arg = std::stod(get_arg_as_string(line, arg_num));  }

    /// Get argument - unsigned long
    void get_argument(std::string line, int arg_num, unsigned long &arg) { arg = std::stoul(get_arg_as_string(line, arg_num)); }

    /// Object to call the function on
    std::shared_ptr<Z> _object;
};

/// Abstract base class for a parser key with 0 arguments
template<class Z>
class SIRFRegParserKey0Arg : public SIRFRegParserKeyBase<Z>
{
public:

    /// Set function
    void set_function(void (Z::*function)()) { _function = function; }

    /// Get arguments - don't need to do anything
    virtual void get_arguments(std::string &line)
    {
        this->check_num_arguments(line,0);
    }

    /// Call the function
    virtual void call_function()
    {
        std::cout << "\tCalling function with no arguments..." << std::flush;
        (this->_object.get()->*_function)();
        std::cout << "done." << std::endl;
    }

    /// Print number of arguments expected
    virtual void print_num_arguments_expected() { std::cout << "\tNumber arguments expected: 0\n"; }


protected:

    /// Function
    void (Z::*_function)();
};

/// Abstract base class for a parser key with 1 argument
template<class Z, class A>
class SIRFRegParserKey1Arg : public SIRFRegParserKeyBase<Z>
{
public:

    /// Set function
    void set_function(void (Z::*function)(A)) { _function = function; }

    /// Get arguments
    virtual void get_arguments(std::string &line)
    {
        this->check_num_arguments(line, 1);

        // Set the argument
        std::cout << "\tSetting argument...\n";
        this->get_argument(line, 0, _arg1);
        std::cout << "\t\tExpected type: " << get_typeid(_arg1) << "\n";
        std::cout << "\t\tValue:         " << _arg1 << "\n";
    }

    /// Call the function
    virtual void call_function()
    {
        std::cout << "\tCalling function with one argument..." << std::flush;
        (this->_object.get()->*_function)(_arg1);
        std::cout << "done." << std::endl;
    }

    /// Print number of arguments expected
    virtual void print_num_arguments_expected() { std::cout << "\tNumber arguments expected: 1 (" << get_typeid(_arg1) << ")\n"; }

protected:

    /// Function
    void (Z::*_function)(A);
    /// Argument
    A _arg1;
};

/// Abstract base class for a parser key with 2 arguments
template<class Z, class A, class B>
class SIRFRegParserKey2Arg : public SIRFRegParserKeyBase<Z>
{
public:

    /// Set function
    void set_function(void (Z::*function)(A,B)) { _function = function; }

    /// Get arguments
    virtual void get_arguments(std::string &line)
    {
        this->check_num_arguments(line, 2);

        // Set the first argument
        std::cout << "\tSetting first argument...\n";
        this->get_argument(line, 0, _arg1);
        std::cout << "\t\tExpected type: " << get_typeid(_arg1) << "\n";
        std::cout << "\t\tValue:         " << _arg1 << "\n";

        // Set the second argument
        std::cout << "\tSetting second argument...\n";
        this->get_argument(line, 1, _arg2);
        std::cout << "\t\tExpected type: " << get_typeid(_arg2) << "\n";
        std::cout << "\t\tValue:         " << _arg2 << "\n";
    }

    /// Call the function
    virtual void call_function()
    {
        std::cout << "\tCalling function with two arguments..." << std::flush;
        (this->_object.get()->*_function)(_arg1,_arg2);
        std::cout << "done." << std::endl;
    }

    /// Print number of arguments expected
    virtual void print_num_arguments_expected() { std::cout << "\tNumber arguments expected: 2 (" << get_typeid(_arg1) << " and " << get_typeid(_arg2) << ")\n"; }

protected:

    /// Function
    void (Z::*_function)(A,B);
    /// First argument
    A _arg1;
    /// Second argument
    B _arg2;
};
}

#endif
