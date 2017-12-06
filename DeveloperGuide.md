# Table of contents

1. [Overview](#Overview)
2. [SIRF structure](#SIRF_structure)
    1. [Rationale](#Rationale)
    2. [Software layers](#Software_layers)
    	1. [Reconstruction engines](#Reconstruction_engines)
    	2. [Extended engine functionality](#Extended_engine_functionality)
    	3. [C interface](#C_interface)
    	4. [Matlab and Pyton interfaces to C](#Matlab_Python_interfaces)
    	5. [Matlab and Pyton OO interfaces](#Matlab_Python_OO_interfaces)
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
Matlab/Python OO interface modules | +mSTIR pSTIR.py | +mGadgetron pGadgetron.py
Matlab/Python interfaces to C | mstir.\* pystir.* | mgadgetron.\* pygadgetron.*
C interface to C\++ code | cstir.* | cgadgetron.*
Extended engine functionality | xSTIR/cSTIR/* | xGadgetron/cGadgetron/*
Reconstruction engines | STIR | Gadgetron

As you start to explore SIRF code, you may notice files and folders preceded by the characters "x", "c", "p" and "m". These correspond to the various layers of SIRF: eXtended engine functionality, C-interface, Python and Matlab.

### Reconstruction engines <a name="Reconstruction_engines"></a>

AAt present, SIRF uses software package STIR for PET reconstruction and Gadgetron for MR reconstruction.

STIR implements a library of C++ classes for performing PET reconstruction and related tasks such as data input/output. Parameters of a STIR reconstructor object are normally defined by the user in an Interfile rather than set by directly calling its mutator (set) methods. A set of executables is provided that read the parameter Interfile and raw data file specified in the command line and perform the required tasks, so that the user does not need to know C\+\+ or any other programming languages at all.

With Gadgetron, reconstruction is performed by a chain of gadgets, pieces of code implementing specific tasks. The chain of gadgets runs on the server, which can be just a command line window, or it can be another computer or a VM. In order to set up the chain, the server needs to receive an xml text describing it from the client, which again can be another command line window on the same or another computer. The first gadget in the chain then starts waiting for acquisition data to arrive from the client in chunks of a certain size. Having processed a chunk of data, the first gadget passes the result to the second and starts processing the next chunk and so on. The last gadget sends the reconstructed images back to the client. As with STIR, the user is not required to have any knowledge of C\++. Instead, the reconstruction tasks are essentially programmed in xml language, which is very similar to STIR, where the reconstruction tasks are programmed in Interfile language.

### Extended engine functionality <a name="Extended_engine_functionality"></a>

The intended usage of STIR and Gadgetron requires minimal participation from the user - normally, the composition of an xml or Interfile document at most. SIRF allows the users to actually code the reconstruction tasks by working with reconstruction objects and the data objects that they operate on. Since neither STIR nor Gadgetron developers were concerned with this kind of usage, their code was not well-suited for this kind of interaction and therefore required some extensions. The extended engine functionality layer of SIRF is a C++ implementation of various data types that provide the necessary extensions.

[//]: <> (extend the functionality of STIR and Gadgetron.)

[//]: <> (The following data types extend the functionality of STIR and Gadgetron.)

#### Extended STIR functionality <a name="Extended STIR functionality"></a>

Extended STIR functionality types are essentially wrappers around relevant STIR classes that hide STIR implementation specifics, so that a different PET reconstruction engine may be used instead without changing anything in the interface layers further above.

###### PETAcquisitionData <a name="PETAcquisitionData"></a>

A class for handling PET acquisition data. Based on STIR `ProjData` class and retains most of its functionality. Has additional algebraic operations functionality. *Files:* `SIRF/src/xSTIR/cSTIR/stir_data_containers.*`.

###### PETImageData <a name="PETImageData"></a>

A class for handling PET image data. Based on STIR `DiscretisedDensity` class and retains most of its functionality. Has additional algebraic operations functionality. *Files:* `SIRF/src/xSTIR/cSTIR/stir_data_containers.*`.

###### PETAcquisitionModel <a name="PETAcquisitionModel"></a>

A class for PET acquisition process simulation. Has method `forward` for simulating acquisition process in a PET scanner and method `backward` for the adjoint (transposed) operation. *Files:* `SIRF/src/xSTIR/cSTIR/stir_x.*`.

###### PETAcquisitionModelUsingMatrix <a name="PETAcquisitionModelUsingMatrix"></a>

A class derived from `PETAcquisitionModel`. Employs STIR object `ProjMatrixByBin` in forward- and backprojection. *Files:* `SIRF/src/xSTIR/cSTIR/stir_x.*`.

#### Extended Gadgetron functionality <a name="Extended Gadgetron functionality"></a>

Extended Gadgetron functionality in SIRF is built entirely upon the Gadgetron client `gadgetron_ismrmrd_client.cpp`: SIRF essentially provides an alternative client. The server-side parts of the Gadgetron code are not used at present. This, for example, allows the use of Gadgetron without its installation under Windows - instead, Gadgetron can be run on a Linux Virtual Machine.

###### MRAcquisitionData <a name="MRAcquisitionData"></a>

A container class for storing and handling ISMRMRD acquisitions. *Files:* `SIRF/src/xGadgetron/cGadgetron/gadgetron_data_containers.*`.

###### MRImageData <a name="MRImageData"></a>

A container class for storing and handling ISMRMRD images. *Files:* `SIRF/src/xGadgetron/cGadgetron/gadgetron_data_containers.*`.

###### MRAcquisitionModel <a name="MRAcquisitionModel"></a>

A class for MR acquisition process simulation. Has method `forward` for simulating acquisition process in an MR scanner and method `backward` for the adjoint (transposed) operation. *Files:* `SIRF/src/xGadgetron/cGadgetron/gadgetron_x.*`.

###### Gadget <a name="Gadget"></a>

A class for storing properties and generating xml description of a Gadgetron gadget. *Files:* `SIRF/src/xGadgetron/cGadgetron/gadget_lib.h`.

###### SIRF gadgets library <a name="Gadgets_library"></a>

At present, the set of classes that have been derived from Gadgets (representing a subset of Gadgetron gadgets) contains:

	NoiseAdjustGadget
    AsymmetricEchoAdjustROGadget
    RemoveROOversamplingGadget
    AcquisitionAccumulateTriggerGadget
    BucketToBufferGadget
    GenericReconCartesianReferencePrepGadget
    SimpleReconGadget
    GenericReconCartesianGrappaGadget
    GenericReconFieldOfViewAdjustmentGadget
    GenericReconImageArrayScalingGadget
    ImageArraySplitGadget
    ExtractGadget
    ComplexToFloatGadget
    FloatToShortGadget

and some other (reader/writer/finish) gadgets not accessible from SIRF scripts. The library also has a gadget `SimpleReconGadgetSet` that actually represents a small gadget chain performing fully sampled reconstruction. *Files:* `SIRF/src/xGadgetron/cGadgetron/gadget_lib.h`.

### C interface <a name="C_interface"></a>

The C interface layer is a set of C functions that wrap SIRF C++ code.

Each interface function has arguments of basic C types: `void*`, `int`, `float` and `char*`, and some functions have either arguments of types `int*` and `float*` (for interfacing to Matlab) or of type `size_t` (for interfacing to Python). The return value of each interface function is either `void` (to be deprecated) or `void*`. *Files:* `SIRF/src/xSTIR/cSTIR/cstir.*`, `SIRF/src/xGadgetron/cGadgetron/cgadgetron.*`.

SIRF Python and Matlab objects do not have direct access to C\++ objects and their data. Their role is to obtain a reference to a C\++ object or its data by calling a C interface function and pass it to another C\++ object by calling another C interface function. A C\++ object or data reference is wrapped into an object of the class DataHandle, which additionally has a property that records any exceptions thrown by the C++ code, and the pointer to this DataHandle object is passed as a `void*` argument or return value. To simplify/shorten the coding of wrapping/unwrapping, a class template ObjectHandle is derived from DataHandle. *Files:* `SIRF/src/iUtilities/data_handle.h`, `SIRF/scr/common/include/SIRF/common/object_handle.inl`.

### Matlab and Python interfaces to C <a name="Matlab_Python_interfaces"></a>

Wrapping C\++ into C dramatically simplifies the interfacing into any programming language. In the case of Matlab, no interfacing is required under Linux, where Matlab can call C library functions directly via `calllib()`, whereas under Windows one just need to add  `__declspec(dllexport)` in front of every C function, which is done by executables `gmi_xstir.exe` and `gmi_xgadgetron.exe`.

For Python, we use SWIG, which requires just these 5 lines to generate the interface for STIR (and similar 5 lines for Gadgetron):

	%module pystir
	%{
	#include "cstir.h"
	%}
	%include "cstir.h"


### Matlab and Python OO interfaces <a name="Matlab_Python_OO_interfaces"></a>

Matlab and Python interfaces of the previous section are not user-friendly and not Object-Oriented, which is why on top of them we have Object-Oriented modules +mSTIR and +mGadgetron in Matlab and pSTIR.py and pGadgetron.py in Python. These modules are described in User Guide.

## Data handling principles <a name="Data_handling"></a>

Data processed by reconstruction engines are often of complicated structure and generally cannot be efficiently handled by script languages. For this reason, in SIRF all data processing is performed by the engines and their extensions, and scripts do not have direct access to data.

[//]: <> (In order to be passed around by scripts, an engine data item is encapsulated in a special interface object ObjectHandle, and scripts can only pass a void pointer to it between engine-calling methods. In addition to encapsulating the engine data, ObjectHandle also stores information on the current execution status, which can be inspected, but not changed, by scripts.)

The only way for users' scripts to work with the engine data is to get a copy of it in a script array format via as_array() methods of SIRF data container classes and pass an array to the engine via their fill() methods.

# Illustration <a name="Illustration"></a>

### Python <a name="Illustration_Python"></a>

A reconstruction script would normally contain the following line indicating the source of raw acquisition data:

    acq_data = AcquisitionData(filename)

(mind that the actual data reading from the indicated file happens later on). The above line creates a Python object of AcquisitionData class. Focusing, for the sake of illustration, on MR reconstruction by Gadgetron, the definition of this class is provided in the Python module pGadgetron.py, and the constructor method that is called to create this object calls, in turn, a C function via the SWIG-generated C-to-Python interface module pygadgetron.py:

    self.handle = pygadgetron.cGT_ISMRMRDAcquisitionsFromFile(filename)

The C function that is called is as follows:

    extern "C"
    void*
    cGT_ISMRMRDAcquisitionsFromFile(const char* file)
    {
        if (!boost::filesystem::exists(file))
            return fileNotFound(file, __FILE__, __LINE__);
        try {
            boost::shared_ptr<MRAcquisitionData> 
                acquisitions(new AcquisitionsFile(file));
            return sptrObjectHandle<MRAcquisitionData>(acquisitions);
        }
        CATCH;
    }

As can be seen from the above C source, the function checks if the acquisition data file exists, and if it does, creates an MRAcquisitionData object of derived type AcquisitionsFile. The return value of this function, which ends up in `acq_data.handle`, is a C void pointer to a C++ object ObjectHandle that encapsulates a shared pointer to an MRAcquisitionData object (cf. [C interface](#C_interface)). We note that MRAcquisitionData and AcquisitionsFile types are part of extended Gadgetron functionality, the first layer above the Gadgetron engine, and the ObjectHandle type is part of the second layer. Finally, the constructor of AcquisitionsFile has the following line

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


