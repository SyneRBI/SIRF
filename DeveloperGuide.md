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
3. [Adding new functionality to SIRF](#Adding_functionality)
	1. [Step 1: create C wrappers](#Create_wrappers)
	2. [Step 2: create Matlab and Python Object-Oriented interface for your additions](#Create_OO_interface)
4. [Adding a gadget to SIRF gadget library](#Adding_gadget)

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

### Reconstruction engines <a name="Reconstruction_engines"></a>

At present, SIRF uses software package STIR for PET reconstruction and Gadgetron for MR reconstruction.

STIR implements a library of C++ classes for performing PET reconstruction and related tasks such as data input/output. Parameters of a STIR reconstructor object normally are defined by the user in an Interfile rather than set by calling directly its mutator (set) methods. A set of executables is provided that would read the parameter Interfile and raw data file specified in the command line and perform the required tasks, so that the user do not have to know C\+\+ or any other programming language at all.

With Gadgetron, reconstruction is performed by a chain of gadgets, pieces of code implementing specific tasks. The chain of gadgets runs on the server, which can be just a command line window, or it can be another computer or a VM. In order to set up the chain, the server needs to receive an xml text describing it from the client, which again can be another command line window on the same or another computer. The first gadget in the chain then starts waiting for acquisition data to arrive from the client in chunks of certain size. Having processed a chunk of data, the first gadget passes the result to the second and starts processing the next chunk and so on. The last gadget sends the reconstructed images back to the client. Just like with STIR, the user is not required to have any knowledge of C\++. Instead, the reconstruction tasks are essentially programmed in xml language, very much like for STIR, reconstruction tasks are programmed in Interfile language.

### Extended engine functionality <a name="Extended_engine_functionality"></a>

The intended usage of STIR and Gadgetron requires minimal participation from the user - normally, the composition of an xml or Interfile document at most. SIRF allows the users to actually code the reconstruction tasks by working with reconstruction objects and data objects they operate on. Since neither STIR nor Gadgetron developers were concerned with this kind of usage, their code was not well-suited for it and required some extensions. The extended engine functionality layer of SIRF is a C++ implementation of various data types that provide the necessary extensions.

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

A set of classes derived from Gadget representing a subset of Gadgetron gadgets, at present contains

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

C interface layer is a set of C functions that wrap SIRF C++ code.

Each interface function has arguments of basic C types: `void*`, `int`, `float` and `char*`, and some functions have either arguments of types `int*` and `float*` (for interfacing to Matlab) or of type `size_t` (for interfacing to Python). The return value of each interface function is either `void` (to be deprecated) or `void*`. *Files:* `SIRF/src/xSTI/cSTIR/cstir.*`, `SIRF/src/xGadgetron/cGadgetron/cgadgetron.*`.

SIRF Python and Matlab objects do not have direct access to C\++ objects and their data. Their role is to obtain a reference to a C\++ object or its data by calling a C interface function and pass it to another C\++ object by calling another C interface function. A C\++ object or data reference is wrapped into an object of the class DataHandle, which additionally has a property that records any exceptions thrown by the C++ code, and the pointer to this DataHandle object is passed as a `void*` argument or return value. To simplify/shorten the coding of wrapping/unwrapping, a class template ObjectHandle is derived from DataHandle. *Files:* `SIRF/src/iUtilities/data_handle.h`, `SIRF/scr/common/include/SIRF/common/object_handle.inl`.

### Matlab and Python interfaces to C <a name="Matlab_Python_interfaces"></a>

Wrapping C\++ into C dramatically simplifies the interfacing into any programming language. In the case of Matlab, no interfacing is meeded under Linux, where Matlab can call C library functions directly via `calllib()`, whereas under Windows one just need to add  `__declspec(dllexport)` in front of every C function, which is done by executables `gmi_xstir.exe` and `gmi_xgadgetron.exe`.

For Python, we use SWIG, which requires just these 5 lines to generate the interface for STIR (and similar 5 lines for Gadgetron):

	%module pystir
	%{
	#include "cstir.h"
	%}
	%include "cstir.h"


### Matlab and Python OO interfaces <a name="Matlab_Python_OO_interfaces"></a>

Matlab and Python interfaces of the previous section are not user-friendly and not Object-Oriented, which is why on top of them we have Object-Oriented modules +mSTIR and +mGadgetron in Matlab and pSTIR.py and pGadgetron.py in Python. These modules are described in User Guide.

## Source files folder structure

    SIRF/src
        common         : code common to all engines
        iUtilities     : interface utilities
        xGadgetron     : Gadgetron extesions and interfaces
        	cGadgetron : Gadgetron extensions and C wrappers
            mGadgetron : Object-Oriented Matlab interface
            pGadgetron : Object-Oriented Python interface
        xSTIR          : STIR extensions and interfaces
        	cSTIR      : STIR extensions and C wrappers
            mSTIR      : Object-Oriented Matlab interface
            pSTIR      : Object-Oriented Python interface

## Data handling principles <a name="Data_handling"></a>

Data processed by reconstruction engines are often of complicated structure and generally cannot be efficiently handled by script languages. For this reason, in SIRF all data processing is performed by the engines and their extensions, and scripts do not have direct access to data.

[//]: <> (In order to be passed around by scripts, an engine data item is encapsulated in a special interface object ObjectHandle, and scripts can only pass a void pointer to it between engine-calling methods. In addition to encapsulating the engine data, ObjectHandle also stores information on the current execution status, which can be inspected, but not changed, by scripts.)

The only way for user's scripts to work with the engine data is to get a copy of it in a script array format via as_array() methods of SIRF data container classes and pass an array to the engine via their fill() methods.

# Illustration <a name="Illustration"></a>

In this section we illustrate how the SIRF software layers interact. To avoid duplication, we focus on the case of MR reconstruction by Gadgetron. The PET case is totally parallel.

### Python <a name="Illustration_Python"></a>

A reconstruction script normally would contain the following line indicating the source of raw acquisition data:

    acq_data = AcquisitionData(filename)

(mind that the actual data reading from the indicated file happens later on). The above line creates a Python object of AcquisitionData class. The definition of this class is provided in the Python module pGadgetron.py, and the constructor method that is called to create this object calls, in turn, a C function via SWIG-generated C-to-Python interface module pygadgetron.py:

    self.handle = pygadgetron.cGT_ISMRMRDAcquisitionsFromFile(filename)

The C function that is called is as follows:

    extern "C"
    void* cGT_ISMRMRDAcquisitionsFromFile(const char* file)
    {
        if (!boost::filesystem::exists(file))
            return fileNotFound(file, __FILE__, __LINE__);
        try {
            shared_ptr<MRAcquisitionData> 
                acquisitions(new AcquisitionsFile(file));
            return newObjectHandle<MRAcquisitionData>(acquisitions);
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

# Adding new functionality to SIRF <a name="Adding_new_functionality"></a>

The preferred way for adding new functionality to STIR is by implementing it in C\++. Having implemented and tested your C++ addition, put your class definitions to `gadgetron_x.h` and your implementation code to `gadgetron_x.cpp` in folder `SIRF/src/xGadgetron/cGadgetron`. Alternatively, add new `*.h` and `*.cpp` files, in which case you will need to list them in `add_library` statement in `CMakeLists.txt` too. Then follow the two steps described below. Again, to reduce duplication, instructions are for Gadgetron-related additions; those for STIR-related are symmetric.

We stress that the instructions below are for classes that will be exposed to the user only. We note that currently, only a tiny fraction of STIR classes and no Gadgetron classes are exposed.

### Step 1: create C wrappers. <a name="Create_wrappers"></a>

If a new class has default constructor, you may like to add the following two lines to the function `cGT_newObject()` in `cgadgetron.cpp` before the last `return` statement there:

	if (boost::iequals(name, "YourClassName"))
		return newObjectHandle<YourClassName>();

For all other constructors and other methods of your class you will need to write a C wrapper, the prototype of which **must** be placed in `cgadgetron.h`. The definition of the wrapper may be added to `cgadgetron.cpp` or a new file in `SIRF/src/xGadgetron/cGadgetron`, which will need to be mentioned in `add_library` statement in `CMakeLists.txt`.

It is advised to use existing C wrappers as templates. An example of a C wrapper for a constructor of an object of class MRAcquisitionData was given above in section Illustration. It exposes some general principles in creating C wrappers:

- wrapper should return a pointer to `DataHandle` object cast into `void*`;
- wrapper should use `try{} CATCH;` bracket for passing C++ exceptions to Python and Matlab;
- a pointer to newly constructed object should be stored in a shared pointer (file `cgadgetron_shared_ptr.h` specifies which shared pointer to use, `std::shared_ptr` or `boost::shared_ptr`);
- a `DataHandle` wrapper for a shared pointer can be created by the function template `newObjectHandle`, which returns a pointer to the `DataHandle` object it creates.

Some further principles are illustrated by the following example of a wrapper:

	extern "C"
	void* cGT_norm(const void* ptr_x)
	{
		try {
			aDataContainer<complex_float_t>& x =
				objectFromHandle<aDataContainer<complex_float_t> >(ptr_x);
			return dataHandle(x.norm());
		}
		CATCH;
	}

This wrapper returns the return value of the method `norm()` of an object of the class template `aDataContainer<complex_float_t>`. The argument of the wrapper is actually a pointer to a `DataHandle` object that stores a shared pointer to an `aDataContainer<complex_float_t>` object. The function template `objectFromHandle` obtains a reference to this object, thus enabling the call to its method `norm()`. The function template `dataHandle`, which has one argument of arbitrary scalar type, wraps the return value of `norm()` into `DataHandle` object, and returns the pointer to this object as `void*`.

The next two wrappers demonstrate how data is exchanged between C++ and Matlab/Python arrays.

	extern "C"
	void* cGT_getAcquisitionsData
	(void* ptr_acqs, unsigned int slice, size_t ptr_re, size_t ptr_im)
	{
		try {
			float* re = (float*)ptr_re;
			float* im = (float*)ptr_im;
			MRAcquisitionData& acqs =
				objectFromHandle<MRAcquisitionData>(ptr_acqs);
			acqs.get_acquisitions_data(slice, re, im);
			return (void*)new DataHandle;
		}
		CATCH;
	}

	extern "C"
	void* cGT_setAcquisitionsData
	(void* ptr_acqs, unsigned int na, unsigned int nc, unsigned int ns,
	size_t ptr_re, size_t ptr_im)
	{
		try {
			float* re = (float*)ptr_re;
			float* im = (float*)ptr_im;
			MRAcquisitionData& acqs =
				objectFromHandle<MRAcquisitionData>(ptr_acqs);
			int err = acqs.set_acquisition_data(na, nc, ns, re, im);
			DataHandle* handle = new DataHandle;
            if (err)
                handle->set_status("Mismatching acquisition dimensions",
                    __FILE__, __LINE__);
			return (void*)handle;
		}
		CATCH;
	}

Both wrappers obtain a pointer to `DataHandle` containing the reference `acqs` to an `MRAcquisitionData` object and pointers to the data of two Matlab/Python arrays storing real and imaginary part of MR acquisitions. The reference `acqs` is used to call methods `get_acquisitions_data` and `set_acquisitions_data` which perform the data exchange.

The second wrapper also shows how to deal with errors. If `set_acquisitions_data` returns non-zero error flag, method `set_status` of a `DataHandle` object is called that records the error message and location. This error message and location are used by SIRF Python and Matlab OO interfaces to throw an exception.

### Step 2: create Matlab and Python Object-Oriented interface for your additions <a name="Create_OO_interface"></a>

This section demonstrates how C wrappers of the previous section are used by the objects of SIRF Object-Oriented interface modules.

#### Python

The SIRF build creates a Python module `pygadgetron.py` with SWIG-generated interface to Python for C wrappers declared in `cgadgetron.h` and a Python module `pyutilities.py` with various interface utilities. An example of a Python interface function from `pygadgetron.py` is the function `pygadgetron.cGT_ISMRMRDAcquisitionsFromFile` featuring in Illustration section. This function is called by the constructor of the Python class `AcquisitionData` defined in SIRF Object-Oriented interface module `pGadgetron.py`:

    self.handle = pygadgetron.cGT_ISMRMRDAcquisitionsFromFile(filename)

and the destructor of this class calls a Python interface function `deleteDataHandle` from `pyutilities.py` to destruct the `DataHandle` object created by the above line:

    pyiutilities.deleteDataHandle(self.handle)

Class `AcquisitionData` has method `as_array()` that creates and returns a Python array containing acquisition data copied from `MRAcquisitionData` as follows

    re = numpy.ndarray((ny, nc, ns), dtype = numpy.float32)
    im = numpy.ndarray((ny, nc, ns), dtype = numpy.float32)
    try_calling(pygadgetron.cGT_getAcquisitionsData\
        (self.handle, n, re.ctypes.data, im.ctypes.data))
    return re + 1j*im

where `try_calling` is a function that checks the execution status stored by `DataHandle` object, the pointer to which is returned by `cGT_getAcquisitionsData`, and throws Python exception if the execution status is abnormal. In a similar way, the data stored in a complex Python array `data` is copied to `MRAcquisitionData` as follows

    na, nc, ns = data.shape
    re = numpy.real(data).astype(numpy.float32)
    im = numpy.imag(data).astype(numpy.float32)
    try_calling(pygadgetron.cGT_setAcquisitionsData\
        (self.handle, na, nc, ns, re.ctypes.data, im.ctypes.data))

The hierarchy of SIRF interface classes mirrors that of its C\++ classes: the C++ class `MRAcquisitionData` inherits from `aDataContainer`, and, accordingly, `AcquisitionData` inherits from `DataContainer` class of the interface module `pGadgetron.py`. The latter class has method `norm()`, which has the following implementation:

    handle = pygadgetron.cGT_norm(self.handle)
    check_status(handle)
    r = pyiutility.floatDataFromHandle(handle)
    pyiutility.deleteDataHandle(handle)
    return r

where `check_status` is the status-checking function called by `try_calling`, and the function `floatDataFromHandle` retrieves a single precision `float` value stored in `DataHandle` object of `handle`.

#### Matlab

SIRF build at present does not generate Matlab interface. to generate them, the user must execute `gmi_xstir` and `gmi_xgadgetron` in subfolder `bin` of the installation folder. This will create files `mstir.c` and `mgadgetron.c` containing Matlab callable C wrappers for SIRF C interface functions (cf. `mGT_ISMRMRDAcquisitionsFile` in section Illustration), and corresponding specification files `mstir.h` and `mgadgetron.h`.

The Matlab counterparts of the code pieces of the previous section are:

    self.handle_ = calllib('mgadgetron', ...
        'mGT_ISMRMRDAcquisitionsFromFile', filename);

    if ~isempty(self.handle_)
        mUtilities.delete(self.handle_)
        self.handle_ = [];
    end

(note that `self.handle_` needs to be set to empty to avoid its deletion by the base class)

    ptr_re = libpointer('singlePtr', zeros(na*nc*ns, 1));
    ptr_im = libpointer('singlePtr', zeros(na*nc*ns, 1));
    calllib('mgadgetron', 'mGT_getAcquisitionsData', ...
        self.handle_, n, ptr_re, ptr_im);
    re = reshape(ptr_re.Value, ns, nc, na);
    im = reshape(ptr_im.Value, ns, nc, na);
    data = re + 1i*im;

    [ns, nc, na] = size(data);
    re = real(data);
    im = imag(data);
    if isa(re, 'single')
        ptr_re = libpointer('singlePtr', re);
        ptr_im = libpointer('singlePtr', im);
    else
        ptr_re = libpointer('singlePtr', single(re));
        ptr_im = libpointer('singlePtr', single(im));
    end
    h = calllib('mgadgetron', 'mGT_setAcquisitionsData', ...
        self.handle_, na, nc, ns, ptr_re, ptr_im);
    mUtilities.check_status('AcquisitionData', h);
    mUtilities.delete(h)

    handle = calllib('mgadgetron', 'mGT_norm', self.handle_);
    mUtilities.check_status('DataContainer', handle);
    r = calllib('miutilities', 'mFloatDataFromHandle', handle);
    mUtilities.delete(handle)

# Adding a gadget to SIRF gadget library <a name="Adding_gadget"></a>

To add a Gadgetron gadget to SIRF gadgets library, follow the steps below.

1. Locate the subfolder of `Gadgetron/gadgets` contaning the gadget definition (below referred to as `GADGET_FOLDER`).

2. Add a new class to `gadget_lib.h` in folder `Gadgetron` using the following template (see also other gadgetrs declarations in `gadget_lib.h` for guidance):

    class NEW_GADGET_CLASS: public Gadget {
    public:
        NEW_GADGET_CLASS() :
            Gadget(GADGET_NAME, GADGET_FOLDER, GADGET_CLASS)
        {
            add_property(PROPERTY1, VALUE1);
            add_property(PROPERTY2, VALUE2);
            ...
        }
        static const char* class_name()
        {
            return NEW_GADGET_CLASS_STRING;
        }
    };
Here `GADGET_CLASS` is the name of the gadget class, `GADGET_NAME` and `NEW_GADGET_CLASS` are arbitrary names for the class and gadget (it is recommended that `NEW_GADGET_CLASS = GADGET_CLASS`), `PROPERTY1` and `VALUE1` are C strings containing the first property name and value etc. (these lines are not needed if the gadget does not have properties), and `NEW_GADGET_CLASS_STRING` is a C string containing the name `NEW_GADGET_CLASS`.

3. Add a line

		NEW_GADGET(NEW_GADGET_CLASS);
to the function `cGT_newObject()` in the file `cgadgetron.cpp`.

4. Build SIRF.
