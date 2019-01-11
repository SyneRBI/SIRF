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

#pragma once

#include "sirf/cReg/ParserKey.h"
#include <memory>
#include <map>
#include <fstream>
#include <boost/filesystem.hpp>
#include <type_traits>

namespace sirf {
/// Reg parser
template<class Z>
class Parser
{
public:

    /// Set object
    void set_object(std::shared_ptr<Z> object) { _object = object; }

    /// Set filename
    void set_filename(boost::filesystem::path filename) { _filename = filename; }

    /// Add key with 0 arguments
    void add_key(const std::string &keyword, void (Z::*function)())
    {
        std::cout << "Adding a keyword to the list of possible parameters: " << keyword << "..." << std::flush;

        // Create a new key, then add keyword and function
        std::shared_ptr<ParserKey0Arg<Z> > key = std::make_shared<ParserKey0Arg<Z> >();
        key->set_object(_object);
        key->set_function(function);
        _keys[keyword] = key;

        std::cout << " done.\n" << std::flush;
    }

    /// Add key with 1 argument - we need to template Y in case Z is the derived class and Y is a parent class
    template<class Y, class A>
    void add_key(const std::string &keyword, void (Y::*function)(A))
    {
        std::cout << "Adding a keyword to the list of possible parameters: " << keyword << "..." << std::flush;

        // Create a new key, then add keyword and function
        std::shared_ptr<ParserKey1Arg<Z,A> > key = std::make_shared<ParserKey1Arg<Z,A> >();
        key->set_object(_object);
        key->set_function(function);
        _keys[keyword] = key;

        std::cout << " done.\n" << std::flush;
    }

    /// Add key with 2 arguments - we need to template Y in case Z is the derived class and Y is a parent class
    template<class Y, class A, class B>
    void add_key(const std::string &keyword, void (Y::*function)(A, B))
    {
        std::cout << "Adding a keyword to the list of possible parameters: " << keyword << "..." << std::flush;

        // Create a new key, then add keyword and function
        std::shared_ptr<ParserKey2Arg<Z,A,B> > key = std::make_shared<ParserKey2Arg<Z,A,B> >();
        key->set_object(_object);
        key->set_function(function);
        _keys[keyword] = key;

        std::cout << " done.\n" << std::flush;
    }

    /// Parse
    void parse()
    {
        // Open file
        std::ifstream file;
        this->open_file(file);

        std::string line;

        // Do nothing with the first line (as this is the title)
        this->get_line(file,line);

        // Keep in the while loop until break
        while(true) {

            // Get next line
            this->get_line(file,line);

            // If end was found, break
            if (this->is_end(line)) break;

            // If no parameter was found (blank line or comment), continue
            if (!this->contains_parameter(line)) continue;

            // Give info to user
            std::cout << "\nReading the line: \"" << line << "\"...\n";

            // Get the keyword
            std::string keyword = this->get_keyword(line);

            // Get the key that corresponds to the keyword
            std::shared_ptr<ParserKeyBase<Z> > key;
            key = this->get_key_sptr(keyword);

            // Figure out the variables from the line of text
            key->get_arguments(line);

            // Call the function
            key->call_function();
        }

        // Close the file
        file.close();
    }

    /// Open file
    void open_file(std::ifstream &file)
    {
        // Open file
        file.open(_filename.string());

        // Check file is open
        if (!file) {
            file.close();
            throw std::runtime_error("Unable to open registration parameter file:\n\t" + _filename.string());
        }

        std::cout << "\nSuccessfully opened the parameter file: " << _filename << "\n";
    }

    /// Get line
    void get_line(std::ifstream &file, std::string &line)
    {
        // If end of file was reached (this shouldn't happen, as the end statement in the
        // parameter file should cause a break in the while loop before this gets called
        if (!getline(file,line)) {
            file.close();
            throw std::runtime_error("Error reading parameter file. Maybe the end of the file was reached but no \"End\" statement was found.");
        }

        // Take out any spaces and tabs
        line.erase(remove(line.begin(),line.end(),' '),line.end());
        line.erase(remove(line.begin(),line.end(),'\t'),line.end());

        // Remove everything to the right of the comment (if there is one)
        line = line.substr(0, line.find(';'));

    }

    /// Is the end of the file?
    bool is_end(std::string line)
    {
        // Convert to lower case
        std::transform(line.begin(), line.end(), line.begin(), ::tolower);

        // If first three characters are "end", return true
        if (line.substr(0,3) == "end") {
            std::cout << "\nThe end line was reached. Parsing will stop.\n\n";
            return true;
        }

        // Else, return false
        return false;
    }

    /// Does the line contain a parameter?
    bool contains_parameter(std::string &line)
    {
        // If no ":=" is found, return false
        if (line.find(":=") == -1) return false;

        return true;
    }

    /// Get the keyword from a line
    std::string get_keyword(std::string &line)
    {
        // get index of the ":="
        int index_end = line.find(":=");

        std::string keyword = line.substr(0, index_end);

        std::cout << "\tKeyword: " << keyword << "\n";

        // return the substring
        return line.substr(0, index_end);
    }

    /// Get a key from its keyword
    std::shared_ptr<ParserKeyBase<Z> > get_key_sptr(std::string &keyword)
    {
        // Compare the keyword to the list of keys
        typename std::map<std::string, std::shared_ptr<ParserKeyBase<Z> > >::iterator it;
        it = _keys.find(keyword);

        // If the keyword is not found in the list of keys, display error & go to next line.
        if (it == _keys.end()) {
            throw std::runtime_error("Error. Key not found in list of possible keywords.");
        }

        // If the keyword was found, print it
        std::cout << "\tMatching key was found.\n";

        // Get the key
        return it->second;
    }

protected:
    /// Object
    std::shared_ptr<Z>                                                _object;
    /// Filename
    boost::filesystem::path                                           _filename;
    /// Map of keys
    std::map<std::string, std::shared_ptr<ParserKeyBase<Z> > > _keys;
};
}
