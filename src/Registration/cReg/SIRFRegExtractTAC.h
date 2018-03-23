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
\brief Base class for all SIRF registration.

The parameter file should have interfile-like syntax.
The variables will be stored as a vector of floats and converted into the required type (int, unsigned int, etc) if necessary.
Multiple variables for a given parameter should be comma separated.
Spaces and tabs will be ignored.
For the title, it doesn't matter what is written as it will be ignored, but something has to be there (otherwise the first parameter will be ignored).
Possible parameters are all the Set<something> methods for each class (e.g., nifty_aladin::SetPerformRigid) and should be written in the parameter file without the "Set" (e.g., PerformRigid).

An example is given below:
    SomeTitle :=
        ReferenceTimePoint := 1
        FloatingTimePoint := 2
        LinearEnergyWeights := 1.5,1
        AdditiveMC :=
    end :=

More examples can be found in // Need to give path

\author Richard Brown
\author CCP PETMR
*/

#ifndef _SIRFREGEXTRACTTAC_H_
#define _SIRFREGEXTRACTTAC_H_

#include <stdexcept>
#include <nifti1_io.h>
#include <string>
#include <vector>
#include <iostream>
#include <boost/filesystem.hpp>

/// Extract the TAC from a set of images
class SIRFRegExtractTAC
{
public:

    /// Constructor
    SIRFRegExtractTAC() {}

    /// Destructor
    virtual ~SIRFRegExtractTAC() {}

    /// Add input filename
    void add_input_filename(const std::string filename) { _filenames.push_back(filename); }

    /// Clear input filenames
    void clear_input_filenames() { _filenames.clear(); }

    /// Set segmentation filename
    void set_segmentation_filename(const std::string filename) { _segmentation_filename = filename; }

    /// Add VOI (all labels with the same string will be joined together as one VOI)
    void add_VOI(const int index, const std::string name) { _VOIs.push_back(std::pair<int, std::string>(index,name)); }

    /// Clear VOIs
    void clear_VOIs() { _VOIs.clear(); }

    /// Update
    virtual void update();

    /// Print TAC

    /// Save TAC to file

protected:

    /// Check parameters
    virtual void check_parameters();

    /// Vector of filenames
    std::vector<boost::filesystem::path>     _filenames;
    /// Segmentation filename
    boost::filesystem::path                  _segmentation_filename;
    /// Vector of VOIs
    std::vector<std::pair<int,std::string> > _VOIs;
};

#endif
