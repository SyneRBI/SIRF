# Table of contents

1. [Overview](#Overview)
2. [SIRF structure](#SIRF_structure)
    1. [Rationale](#Rationale)
    2. [Software layers](#Software_layers)
    3. [Data handling](#Data_handling)
    4. [Illustration](#Illustration)
        1. [Python](#Illustration_Python)
        2. [Matlab](#Illustration_Matlab)

# Overview <a name="Overview"></a>

The SIRF (Synergistic Image Reconstruction Framework) software is an Open Source toolkit for the reconstruction of PET and MRI raw data. The aim is to provide code simple enough to easily perform a reconstruction, yet powerful enough to be able to handle real, full-size datasets. Our strategy in achieving this aim is to employ available Open Source reconstruction software written in advanced programming languages such as C++ and provide basic-user-friendly interfaces to it written in script languages, primarily Matlab and Python. This document describes in detail the interfacing principles and structure in order to facilitate the contributions to SIRF from any interested developer.
This Developer's Guide describes version 0.9 of SIRF. The software can be found on [https://github.com/CCPPETMR](https://github.com/CCPPETMR).

# SIRF structure <a name="SIRF_structure"></a>

## Rationale <a name="Rationale"></a>

SIRF builds upon existing Open Source software packages for medical image reconstruction, referred to as *reconstruction engines*, or simply *engines*, in this document. At present, we use STIR (Software for Tomographic Image Reconstruction) as a PET engine and Gadgetron as MR engine.

The way we use reconstruction engines differs from the one intended by their developers, which necessitated various additions to the engines' functionality. Since the engines we use at present are coded in C\+\+, these additions, which essentially form the first layer of SIRF interface software into which we wrap the engines, are in C\+\+ as well.

To avoid difficulties with the interfacing of a C\+\+ code into the script languages we are targeting, we add another layer by wrapping C\+\+ code into C. This extra step (pretty routine as C\+\+ is built upon C) dramatically simplifies interfacing and widens the scope of our targets to practically any programming language in use.

Interfacing of C into Matlab is pretty straightforward, as Matlab has facilities (the function calllib etc.) for calling functions from C libraries, with just a little 'syntactic sugar' needed under Windows. For Python we use [SWIG](http://www.swig.org) (Simple Wrapper and Interface Generator), which can also be used for a number of other prospective languages (follow SWIG link to see the list of currently covered languages).

Our topmost interface layer are Object-Oriented Matlab and Python modules that cover up in a nice user-friendly way all the technicalities involved in the inner interface layers.

## Software layers <a name="Software_layers"></a>

To summarise, SIRF software is structured as the following set of layers (from top to bottom):

SIRF software layers |
-------------------------|
Matlab/Python modules |
Matlab/Python interfaces to C |
C interface to C\+\+ code |
Extended engine functionality |
Reconstruction engines |

## Data handling <a name="Data_handling"></a>

Data processed by reconstruction engines are often of complicated structure and generally cannot be efficiently handled by script languages. For this reason, in SIRF all data processing is performed by the engines and their extensions, and scripts do not have direct access to data. In order to be passed around by scripts, an engine data item is encapsulated in a special interface object ObjectHandle, and scripts can only pass a void pointer to it between engine-calling methods. In addition to encapsulating the engine data, ObjectHandle also stores information on the current execution status, which can be inspected (but not changed) by scripts.

The only way for user's scripts to work with the engine data is to get a copy of it in a script array format via as_array() methods of SIRF data container classes and pass an array to the engine via their fill() methods.

## Illustration <a name="Illustration"></a>

### Python <a name="Illustration_Python"></a>

A reconstruction script normally would contain the following line indicating the source of raw acquisition data:

    acq_data = AcquisitionData(filename)

(mind that the actual data reading from the indicated file happens later on). The above line creates a Python object of AcquisitionData class. Focusing, for the sake of illustration, on MR reconstruction by Gadgetron, the definition of this class is provided in the Python module pGadgetron.py, and the constructor method that is called to create this object calls, in turn, a C function via SWIG-generated C-to-Python interface module pygadgetron.py:

    self.handle = pygadgetron.cGT_ISMRMRDAcquisitionsFromFile(filename)

The C function that is called is as follows:

    extern "C"
    void*
    cGT_ISMRMRDAcquisitionsFromFile(const char* file)
    {
        if (!boost::filesystem::exists(file))
            return fileNotFound(file, __FILE__, __LINE__);
        try {
            boost::shared_ptr<AcquisitionsContainer> 
                acquisitions(new AcquisitionsFile(file));
            return sptrObjectHandle<AcquisitionsContainer>(acquisitions);
        }
        CATCH;
    }

As can be seen from the above C source, the function checks if the acquisition data file exists, and if it does, creates an AcquisitionsContainer object of derived type AcquisitionsFile. The return value of this function, which ends up in `acq_data.handle`, is a C void pointer to a C++ object ObjectHandle that encapsulates a shared pointer to an AcquisitionsContainer object (cf. [Data handling](#Data_handling)). We note that AcquisitionsData and AcquisitionsFile types are part of extended Gadgetron functionality, the first layer above the Gadgetron engine, and the ObjectHandle type is part of the second layer. Finally, the constructor of AcquisitionsFile has the following line

    dataset_ = boost::shared_ptr<ISMRMRD::Dataset>
        (new ISMRMRD::Dataset(filename.c_str(), "/dataset"));

showing that the acquisition data file is handled by the Dataset object of ISMRMRD library employed by Gadgetron engine.

### Matlab <a name="Illustration_Matlab"></a>

In the case of a Matlab script, the line

    acq_data = AcquisitionData(filename);

creates a Matlab object of AcquisitionData class, the constructor of which contains a C function call (via Matlab's calllib function)

    self.handle_ = calllib('mgadgetron', ...
        'mGT_ISMRMRDAcquisitionsFromFile', filename);

to the following C function

    #ifdef _WIN32
    #define EXPORTED_FUNCTION __declspec(dllexport)
    #else
    #define EXPORTED_FUNCTION
    #endif
    EXPORTED_FUNCTION void* mGT_ISMRMRDAcquisitionsFile(const char* file) {
        return cGT_ISMRMRDAcquisitionsFile(file);
    }


