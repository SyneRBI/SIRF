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
5. [SIRF programming conventions](#SIRF_programming_conventions)
	1. [Naming](#conventions_naming)
	2. [File conventions](#conventions_file)
	3. [Others](#conventions_others)
6. [Cloning data container](#clone_datacontainer)

# Overview <a name="Overview"></a>

The SIRF (Synergistic Image Reconstruction Framework) software is an Open Source toolkit for the reconstruction of PET and MRI raw data. The aim is to provide code simple enough to easily perform a reconstruction, yet powerful enough to be able to handle real, full-size datasets. Our strategy in achieving this aim is to employ available Open Source reconstruction software written in advanced programming languages such as C++ and provide basic-user-friendly interfaces to it written in script languages, primarily Matlab and Python. This document describes in detail the interfacing principles and structure in order to facilitate the contributions to SIRF from any interested developer.
This Developer's Guide is for version 2.1 of SIRF. The software can be found on [https://github.com/CCPPETMR](https://github.com/CCPPETMR).

# SIRF structure <a name="SIRF_structure"></a>

## Rationale <a name="Rationale"></a>

SIRF builds upon existing Open Source software packages for medical image reconstruction, referred to as *reconstruction engines*, or simply *engines*, in this document. At present, we use STIR (Software for Tomographic Image Reconstruction) as a PET engine and Gadgetron as MR engine.

The way we use reconstruction engines differs from the one intended by their developers, which necessitated various additions to the engines' functionality. Since the engines we use at present are coded in C\+\+, these additions, which essentially form the first layer of SIRF interface software into which we wrap the engines, are in C\+\+ as well.

To avoid difficulties with the interfacing of a C\+\+ code into the script languages we are targeting, we add another layer by wrapping C\+\+ code into C. This extra step (pretty routine as C\+\+ is built upon C) dramatically simplifies interfacing and widens the scope of our targets to practically any programming language in use.

Interfacing of C into Matlab is pretty straightforward, as Matlab has facilities (the function calllib etc.) for calling functions from C libraries, with just a little 'syntactic sugar' needed under Windows. For Python we use [SWIG](http://www.swig.org) (Simple Wrapper and Interface Generator), which can also be used for a number of other prospective languages (follow SWIG link to see the list of currently covered languages).

Our topmost interface layer are Object-Oriented Matlab and Python modules that cover up in a nice user-friendly way all the technicalities involved in the inner interface layers.

## Software layers <a name="Software_layers"></a>

To summarise, SIRF software is structured as the following set of layers (from top to bottom):

| SIRF software layers | PET | MR | Registration | Common
| --- | --- | --- |
| Matlab/Python OO interfaces | `+STIR STIR.py` | `+Gadgetron Gadgetron.py`| `+Reg Reg.py` | `+SIRF SIRF.py`
| Matlab/Python interfaces to C | `mstir.* pystir.*` | `mgadgetron.* pygadgetron.*` | `mreg.* pyreg.py`| `msirf.m pysirf.py`
| C interface to C\++ code | `cstir.*` | `cgadgetron.*` | `cReg.*`| `csirf.*`
| Extended engine functionality | `xSTIR/cSTIR/*` | `xGadgetron/cGadgetron/*` |
| Reconstruction engines | `STIR/*` | `Gadgetron/*` |

As you start to explore SIRF code, you may notice files and folders preceded by the characters "x", "c", "p" and "m". These correspond to the various layers of SIRF: eXtended engine functionality, C-interface, Python and Matlab.

### Reconstruction engines <a name="Reconstruction_engines"></a>

At present, SIRF uses software package STIR for PET reconstruction and Gadgetron for MR reconstruction.

STIR implements a library of C++ classes for performing PET reconstruction and related tasks such as data input/output. Parameters of a STIR reconstructor object are normally defined by the user in an Interfile rather than set by directly calling its mutator (set) methods. A set of executables is provided that read the parameter Interfile and raw data file specified in the command line and perform the required tasks, so that the user does not need to know C\+\+ or any other programming languages at all.

With Gadgetron, reconstruction is performed by a chain of gadgets, pieces of code implementing specific tasks. The chain of gadgets runs on the server, which can be just a command line window, or it can be another computer or a VM. In order to set up the chain, the server needs to receive an xml text describing it from the client, which again can be another command line window on the same or another computer. The first gadget in the chain then starts waiting for acquisition data to arrive from the client in chunks of a certain size. Having processed a chunk of data, the first gadget passes the result to the second and starts processing the next chunk and so on. The last gadget sends the reconstructed images back to the client. As with STIR, the user is not required to have any knowledge of C\++. Instead, the reconstruction tasks are essentially programmed in xml language, which is very similar to STIR, where the reconstruction tasks are programmed in Interfile language.

### Extended engine functionality <a name="Extended_engine_functionality"></a>

The intended usage of STIR and Gadgetron requires minimal participation from the user - normally, the composition of an xml or Interfile document at most. SIRF allows the users to actually code the reconstruction tasks by working with reconstruction objects and the data objects that they operate on. Since neither STIR nor Gadgetron developers were concerned with this kind of usage, their code was not well-suited for this kind of interaction and therefore required some extensions. The extended engine functionality layer of SIRF is a C++ implementation of various data types that provide the necessary extensions.

#### Extended STIR functionality <a name="Extended STIR functionality"></a>

Extended STIR functionality types are essentially wrappers around relevant STIR classes that hide STIR implementation specifics, so that a different PET reconstruction engine may be used instead without changing anything in the interface layers further above.

###### PETAcquisitionData <a name="PETAcquisitionData"></a>

A class for handling PET acquisition data. Based on STIR `ProjData` class and retains most of its functionality. Has additional algebraic operations functionality. *Files:* `SIRF/src/xSTIR/cSTIR/stir_data_containers.*`.

###### PETImageData <a name="PETImageData"></a>

A class for handling PET image data. Based on STIR `DiscretisedDensity` class and retains most of its functionality. Has additional algebraic operations functionality. *Files:* `SIRF/src/xSTIR/cSTIR/stir_data_containers.*`.

###### ListmodeToSinograms <a name="ListmodeToSinograms"></a>

A class for listmode-to-sinogram conversion and randoms estimation.

###### PETAcquisitionModel <a name="PETAcquisitionModel"></a>

A class for PET acquisition process simulation. Has method `forward` for simulating acquisition process in a PET scanner and method `backward` for the adjoint (transposed) operation. *Files:* `SIRF/src/xSTIR/cSTIR/stir_x.*`.

###### PETAcquisitionModelUsingMatrix <a name="PETAcquisitionModelUsingMatrix"></a>

A class derived from `PETAcquisitionModel`. Employs STIR object `ProjMatrixByBin` in forward- and backprojection. *Files:* `SIRF/src/xSTIR/cSTIR/stir_x.*`.

###### PETAcquisitionSensitivityModel <a name="PETAcquisitionSensitivityModel"></a>

A class for PET acquisition sensitivity simulation, responsible for modelling detector efficiencies and attenuation.

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
    FatWaterGadget
    PhysioInterpolationGadget

    ImageArraySplitGadget
    ExtractGadget
    ComplexToFloatGadget
    FloatToShortGadget

and some other (reader/writer/finish) gadgets not accessible from SIRF scripts. The library also has a gadget `SimpleReconGadgetSet` that actually represents a small gadget chain performing fully sampled reconstruction. *Files:* `SIRF/src/xGadgetron/cGadgetron/gadget_lib.h`.

### C interface <a name="C_interface"></a>

The C interface layer is a set of C functions that wrap SIRF C++ code.

Each interface function has arguments of basic C types: `void*`, `int`, `float` and `char*`, and some functions have either arguments of types `int*` and `float*` (for interfacing to Matlab) or of type `size_t` (for interfacing to Python). The return value of each interface function is either `void` (to be deprecated) or `void*`. *Files:* `SIRF/src/xSTIR/cSTIR/cstir.*`, `SIRF/src/xGadgetron/cGadgetron/cgadgetron.*`.

SIRF Python and Matlab objects do not have direct access to C\++ objects and their data. Their role is to obtain a reference to a C\++ object or its data by calling a C interface function and pass it to another C\++ object by calling another C interface function. A C\++ object or data reference is wrapped into an object of the class `DataHandle`, which additionally has a property that records any exceptions thrown by the C++ code, and the pointer to this `DataHandle` object is passed as a `void*` argument or return value. To simplify/shorten the coding of wrapping/unwrapping, a class template `ObjectHandle` is derived from `DataHandle`. *Files:* `SIRF/src/iUtilities/data_handle.h`, `SIRF/scr/common/include/SIRF/common/object_handle.inl`.

### Matlab and Python interfaces to C <a name="Matlab_Python_interfaces"></a>

Wrapping C\++ into C dramatically simplifies the interfacing into any programming language. In the case of Matlab, no interfacing is required under Linux, where Matlab can call C library functions directly via `calllib()`, whereas under Windows one just need to add  `__declspec(dllexport)` in front of every C function, which is done by executables `gmi_xstir.exe`, `gmi_xgadgetron.exe` etc.

For Python, we use SWIG, which requires just these 5 lines to generate the interface for STIR (and similar lines for Gadgetron and Registration):

    %module pystir
    %{
    #include "cstir.h"
    %}
    %include "cstir.h"


### Matlab and Python OO interfaces <a name="Matlab_Python_OO_interfaces"></a>

Matlab and Python interfaces of the previous section are not user-friendly and not Object-Oriented, which is why on top of them we have Object-Oriented modules `+STIR`, `+Gadgetron` and `+Reg` in Matlab and `STIR.py`, `Gadgetron.py` and `Reg.py` in Python. These modules are described in User Guide.

## Source files folder structure

    SIRF/src
        common         : code common to all engines
        iUtilities     : interface utilities
        Registration   
            cReg       : C wrappers for NiftyReg
            mReg       : Object-Oriented Matlab interface
            pReg       : Object-Oriented Python interface
        xGadgetron     
            cGadgetron : Gadgetron extensions and C wrappers
            mGadgetron : Object-Oriented Matlab interface
            pGadgetron : Object-Oriented Python interface
        xSTIR          
            cSTIR      : STIR extensions and C wrappers
            mSTIR      : Object-Oriented Matlab interface
            pSTIR      : Object-Oriented Python interface

## Data handling principles <a name="Data_handling"></a>

Data processed by reconstruction engines are often of complicated structure and generally cannot be efficiently handled by script languages. For this reason, in SIRF all data processing is performed by the engines and their extensions, and scripts do not have direct access to data.

The only way for user's scripts to work with the engine data is to get a copy of it in a script array format via `as_array()` methods of SIRF data container classes and pass an array to the engine via their `fill()` methods.

# Illustration <a name="Illustration"></a>

In this section we illustrate how the SIRF software layers interact. To avoid duplication, we focus on the case of MR reconstruction by Gadgetron. The PET and Registration cases are totally parallel.

### Python <a name="Illustration_Python"></a>

A reconstruction script would normally contain the following line indicating the source of raw acquisition data:

    acq_data = AcquisitionData(filename)

(mind that the actual data reading from the indicated file happens later on). The above line creates a Python object of `AcquisitionData` class. The definition of this class is provided in the Python module `Gadgetron.py`, and the constructor method that is called to create this object calls, in turn, a C function via SWIG-generated C-to-Python interface module pygadgetron.py:

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

As can be seen from the above C source, the function checks if the acquisition data file exists, and if it does, creates an `MRAcquisitionData` object of derived type `AcquisitionsFile`. The return value of this function, which ends up in `acq_data.handle`, is a C void pointer to a C++ object `ObjectHandle` that encapsulates a shared pointer to an `MRAcquisitionData` object (cf. [C interface](#C_interface)). We note that `MRAcquisitionData` and `AcquisitionsFile` types are part of extended Gadgetron functionality, the first layer above the Gadgetron engine, and the `ObjectHandle` type is part of the second layer. Finally, the constructor of `AcquisitionsFile` has the following line:

    dataset_ = boost::shared_ptr<ISMRMRD::Dataset>
        (new ISMRMRD::Dataset(filename.c_str(), "/dataset"));

showing that the acquisition data file is handled by the Dataset object of ISMRMRD library employed by Gadgetron engine.

If exception is thrown, it is caught by the `CATCH` macro, which records it in the `DataHandle` object that this function will return, so that it can be reported by Python and Matlab interface module `AcquistionData` objects.

### Matlab <a name="Illustration_Matlab"></a>

In the case of a Matlab script, the line

    acq_data = AcquisitionData(filename);

creates a Matlab object of `AcquisitionData` class, the constructor of which contains a C function call (via Matlab's `calllib` function)

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

The preferred way for adding new functionality to SIRF is by implementing it in C\++. Having implemented and tested your C++ addition, put your class definitions to `gadgetron_x.h` and your implementation code to `gadgetron_x.cpp` in folder `SIRF/src/xGadgetron/cGadgetron`. Alternatively, add new `*.h` and `*.cpp` files, in which case you will need to list them in `add_library` statement in `CMakeLists.txt` too. Then follow the two steps described below. Again, to reduce duplication, instructions are for Gadgetron-related additions; those for STIR-related are symmetric.

We stress that the instructions below are for classes that will be exposed to the user only. We note that currently, only a tiny fraction of STIR classes and no Gadgetron classes are exposed.

### Step 1: create C wrappers. <a name="Create_wrappers"></a>

If a new class has default constructor, you may like to add the following two lines to the function `cGT_newObject()` in `cgadgetron.cpp` before the last `return` statement there:

	if (boost::iequals(name, "YourClassName"))
		return newObjectHandle<YourClassName>();

For all other constructors and other methods of your class you will need to write a C wrapper, the prototype of which **must** be placed in `cgadgetron.h`. The definition of the wrapper may be added to `cgadgetron.cpp` or a new file in `SIRF/src/xGadgetron/cGadgetron`, which will need to be mentioned in `add_library` statement in `CMakeLists.txt`.

It is advised to use existing C wrappers as templates. An example of a C wrapper for a constructor of an object of class `MRAcquisitionData` was given above in section Illustration. It exposes some general principles in creating C wrappers:

- wrapper should return a pointer to `DataHandle` object cast into `void*`;
- wrapper should use `try{} CATCH;` bracket for passing C++ exceptions to Python and Matlab;
- a pointer to newly constructed object should be stored in a shared pointer (file `cgadgetron_shared_ptr.h` specifies which shared pointer to use, `std::shared_ptr` or `boost::shared_ptr`);
- a `DataHandle` wrapper for a shared pointer can be created by the function template `newObjectHandle`, which returns a pointer to the `DataHandle` object it creates.

Some further principles are illustrated by the following example of a wrapper:

	extern "C"
	void*
	cSIRF_norm(const void* ptr_x)
	{
		try {
			DataContainer& x =
				objectFromHandle<DataContainer >(ptr_x);
			return dataHandle(x.norm());
		}
		CATCH;
	}

This wrapper returns the return value of the method `norm()` of an object of the class `DataContainer`, from which all SIRF data container classes (`MRImageData`, `PETAcquisitionData` etc.) are derived. The argument of the wrapper is actually a pointer to a `DataHandle` object that stores a shared pointer to a `DataContainer` object. The function template `objectFromHandle` obtains a reference to this object, thus enabling the call to its method `norm()`. The function template `dataHandle`, which has one argument of arbitrary scalar type, wraps the return value of `norm()` into a `DataHandle` object, and returns the pointer to this object as `void*`. 

The next two wrappers demonstrate how data is exchanged between C++ and Matlab/Python arrays.

	extern "C"
	void*
	cGT_acquisitionsDataAsArray(void* ptr_acqs, size_t ptr_z, int all)
	{
		try {
			complex_float_t* z = (complex_float_t*)ptr_z;
			CAST_PTR(DataHandle, h_acqs, ptr_acqs);
			MRAcquisitionData& acqs =
				objectFromHandle<MRAcquisitionData>(h_acqs);
			acqs.get_data(z, all);
			return new DataHandle;
		}
		CATCH;
	}
	
	extern "C"
	void*
	cGT_fillAcquisitionsData(void* ptr_acqs, size_t ptr_z, int all)
	{
		try {
			complex_float_t* z = (complex_float_t*)ptr_z;
			CAST_PTR(DataHandle, h_acqs, ptr_acqs);
			MRAcquisitionData& acqs =
				objectFromHandle<MRAcquisitionData>(h_acqs);
			acqs.set_data(z, all);
			return new DataHandle;
		}
		CATCH;
	}

Both wrappers obtain a pointer `ptr_acqs` to a `DataHandle` object containing the reference `acqs` to an `MRAcquisitionData` object and pointer `ptr_z` to the data of a Matlab/Python array storing MR acquisitions data. The reference `acqs` is used to call methods `get_data` and `set_data` which perform the data exchange. The argument `all` specifies whether all acquisition data must be put in the receiving `ptr_z`-referenced array or whether data not directly used in the reconstruction (e.g. noise calibration data) must be ignored. Exceptions are handled by the `CATCH` macro, successful execution returns a default empty `DataHandle` object, the only role of which in this case is to report successful return.

### Step 2: create Matlab and Python Object-Oriented interface for your additions <a name="Create_OO_interface"></a>

This section demonstrates how C wrappers of the previous section are used by the objects of SIRF Object-Oriented interface modules.

#### Python

The SIRF build creates a Python module `pygadgetron.py` with SWIG-generated interface to Python for C wrappers declared in `cgadgetron.h` and a Python module `pyutilities.py` with various interface utilities. An example of a Python interface function from `pygadgetron.py` is the function `pygadgetron.cGT_ISMRMRDAcquisitionsFromFile` featuring in Illustration section. This function is called by the constructor of the Python class `AcquisitionData` defined in SIRF Object-Oriented interface module `Gadgetron.py`:

    self.handle = pygadgetron.cGT_ISMRMRDAcquisitionsFromFile(filename)

and the destructor of this class calls a Python interface function `deleteDataHandle` from `pyutilities.py` to destruct the `DataHandle` object created by the above line:

    pyiutilities.deleteDataHandle(self.handle)

Class `AcquisitionData` has method `as_array()` that creates and returns a Python array containing acquisition data copied from `MRAcquisitionData` as follows

    z = numpy.ndarray((ny, nc, ns), dtype = numpy.complex64)
    try_calling(pygadgetron.cGT_acquisitionsDataAsArray\
        (self.handle, z.ctypes.data, return_all))
    return z

where `try_calling` is a function that checks the execution status stored by `DataHandle` object, the pointer to which is returned by `cGT_getAcquisitionsData`, and throws Python exception if the execution status is abnormal. In a similar way, the data stored in a complex Python array `data` is copied to `MRAcquisitionData` by the method `fill()` as follows

    try_calling(pygadgetron.cGT_fillAcquisitionsData\
        (self.handle, data.ctypes.data, fill_all))

The hierarchy of SIRF interface classes mirrors that of its C\++ classes: the C++ class `MRAcquisitionData` inherits from C++ `DataContainer`, and, accordingly, `AcquisitionData` inherits from `DataContainer` class of the interface module `SIRF.py`. The latter class has method `norm()`, which has the following implementation:

    handle = pysirf.cSIRF_norm(self.handle)
    check_status(handle)
    r = pyiutility.floatDataFromHandle(handle)
    pyiutility.deleteDataHandle(handle)
    return r

where `check_status` is the status-checking function (called by `try_calling` in the previous examples), and the function `floatDataFromHandle` retrieves a single precision `float` value stored in `DataHandle` object of `handle`.

#### Matlab

SIRF build at present does not generate Matlab interface. to generate them, the user must execute `gmi_*` executables in subfolder `bin` of the installation folder. This will create files `m*.c` containing Matlab callable C wrappers for SIRF C interface functions (cf. `mGT_ISMRMRDAcquisitionsFile` in section Illustration), and corresponding specification files `m*.h`.

The Matlab counterparts of the code pieces of the previous section are:

    % AcquisitionData constructor contains
    self.handle_ = calllib('mgadgetron', ...
        'mGT_ISMRMRDAcquisitionsFromFile', filename);

    % AcquisitionData destructor contains
    if ~isempty(self.handle_)
        sirf.Utilities.delete(self.handle_)
        self.handle_ = []; % to avoid deletion by the base class
    end

    % AcquisitionData.as_array() contains
    n = ns*nc*na;
    ptr_z = libpointer('singlePtr', zeros(2, n));
    calllib('mgadgetron', 'mGT_acquisitionsDataAsArray', ...
        self.handle_, ptr_z, all);
    data = reshape(ptr_z.Value(1:2:end) + 1i*ptr_z.Value(2:2:end), ...
        ns, nc, na);

    % AcquisitionData.fill() contains
    z = [real(data(:))'; imag(data(:))'];
    h = calllib('mgadgetron', 'mGT_fillAcquisitionsData', ...
        self.handle_, ptr_z, all);
    sirf.Utilities.check_status('AcquisitionData', h);
    sirf.Utilities.delete(h)

    % DataContainer.norm() contains
    handle = calllib('msirf', 'mSIRF_norm', self.handle_);
    self.Utilities.check_status('DataContainer', handle);
    r = calllib('miutilities', 'mFloatDataFromHandle', handle);
    sirf.Utilities.delete(handle)

# Adding a gadget to SIRF gadget library <a name="Adding_gadget"></a>

To add a Gadgetron gadget to SIRF gadgets library, follow the steps below.

* Locate the subfolder of `Gadgetron/gadgets` contaning the gadget definition (below referred to as `GADGET_FOLDER`).

* Add a new class to `gadget_lib.h` in folder `Gadgetron` using the following template (see also other gadgetrs declarations in `gadget_lib.h` for guidance):

~~~
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
~~~ 
Here `GADGET_CLASS` is the name of the gadget class, `GADGET_NAME` and `NEW_GADGET_CLASS` are arbitrary names for the gadget and new class and (it is recommended that `NEW_GADGET_CLASS = GADGET_CLASS`), `PROPERTY1` and `VALUE1` are C strings containing the first property name and value etc. (these lines are not needed if the gadget does not have properties), and `NEW_GADGET_CLASS_STRING` is a C string containing the name `NEW_GADGET_CLASS`.

* Add a line

		NEW_GADGET(NEW_GADGET_CLASS);
to the function `cGT_newObject()` in the file `cgadgetron.cpp`.

* Build SIRF.

# SIRF programming conventions <a name="SIRF_programming_conventions"></a>
The programming style used in SIRF resembles closely that used in STIR. When implementing new code in C++, try to follow the following conventions.

## Naming <a name="conventions_naming"></a>
* Classes/types: Use camel caps (also known as medial capitals), for example `ProjMatrixByBin`.
* Variables, methods and members: lower case, with underscores between each word, for example `set_bin_efficiency`.
	* Variable for a number of something: `num_var` 
	* Number of an item in a sequence: `var_num`
	* Pointers, shared pointers and auto pointers should have the respective suffixes: `ptr`, `sptr` and `aptr`.

## File conventions <a name="conventions_file"></a>
* `.inl` files should contain all inline code for functions. This keeps the `h` files clean and short.
* For "supporting" classes in `.h` files, forward declarations should be preferred over `#include`, as this reduces build times when `.h` files are modified.

## Others <a name="conventions_others"></a>
* Access to all class data should be achieved via `get` and `set` methods.
* If outputs are passed as an argument (via reference or pointer), put them first in the list of arguments.
* Error handling should be performed via try/catch statements.

# Cloning data container <a name="clone_datacontainer"></a>
This section should hopefully clear up a bit of confusion around the way cloning an image/acquisition data in SIRF.

We often store an image/acquisition data as a more general version of itself. For example, we might store a `STIRImageData` as an `ImageData`. However, we need to be able to create a clone of the original object. This requires covariant return types, which is tricky with shared/unique pointers. We followed the advice in [this Fluent C++ blog](https://www.fluentcpp.com/2017/09/12/how-to-return-a-smart-pointer-and-use-covariance/). Therefore any new classes should contain a `clone` method (which returns a unique pointer), and a `clone_impl` (which returns a bare pointer) that should only be used by `clone`.