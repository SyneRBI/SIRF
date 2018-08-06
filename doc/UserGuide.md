# Table of Contents

1. [Overview](#Overview)
    1. [General architecture](#General_architecture)
	2. [Supported scanners and file formats](#Supported_scanners_and_file_formats)
	    1. [MRI](#MRI)
		2. [PET](#PET)
2. [Where to find further information](#Further_information)
3. [General notes of usage](#General_notes)
4. [Framework basic functionality](#Basic_functionality)
    1. [General conventions](#General_conventions)
	    1. [Object-oriented paradigm](#Object-oriented_paradigm)
	    2. [Error handling](#Error_handling)
	    3. [Naming conventions](#Naming_conventions)
	    4. [Units and index ordering](#Units_and_index_ordering)
	    5. [Handles](#Handles)
    2. [Library components](#Library_components)
        1. [Getting help on SIRF library modules](#Getting_help_on_SIRF_library_modules)
        2. [General structure of the classes](#General_structure_of_the_classes)
        3. [Basic classes](#Basic_classes)
        4. [Other classes](#Other_classes)
        5. [Functions](#Functions)
5. [Appendix](#Appendix)
    1. [Acquisition data storage scheme management](#storage_management)
    2. [Programming chains of Gadgetron gadgets](#programming_Gadgetron_chains)
        1. [Creating and running gadget chains by SIRF script](#creating_and_running_gadget_chains)
        2. [SIRF gadget library](#SIRF_gadget_library)

# Overview <a name="Overview"></a>

The SIRF (Synergistic Image Reconstruction Framework) software is an Open Source toolkit for the reconstruction of PET and MRI raw data. The aim is to provide code simple enough to easily perform a reconstruction, yet powerful enough to be able to handle real, full-size datasets. Our strategy in achieving this aim is to employ available Open Source reconstruction software written in advanced programming languages such as C++ and provide basic-user-friendly interfaces to it written in script languages, primarily Matlab and Python. The interface style permits a reconstruction to be performed in stages, allowing the user to inspect or modify data, or insert their own code. 

This User’s Guide describes version 1.1 of SIRF. The software can be found on [https://github.com/CCPPETMR](https://github.com/CCPPETMR).

## General architecture <a name="General_architecture"></a>

The code builds upon existing Open Source software packages for medical image reconstruction. At the outset, these packages are STIR for PET reconstruction, and Gadgetron for MRI. SIRF provides MATLAB and Python interfaces to these underlying reconstruction engines. Underlying SIRF are C interfaces to these reconstruction engines. These interfaces are called from higher level MATLAB or Python interfaces. 

At present, you should only use the MATLAB and Python interfaces. The underlying C library is internal and likely to change over the next few releases. 

## Supported scanners and file formats <a name="Supported_scanners_and_file_formats"></a>

### MRI <a name="MRI"></a>

SIRF expects raw MR data in the ISMRMRD format. We currently provide a clone of the siemens_to_ismrmrd Git repository. This enables raw data from Siemens mMR Biograph PET-MR scanners to be converted to ISMRMRD format. For more details of how to export the raw MR data from Siemens PET-MR scanners and how to convert the data to ISMRMRD please see the wiki: [https://github.com/CCPPETMR/SIRF/wiki/MR-raw-data](https://github.com/CCPPETMR/SIRF/wiki/MR-raw-data).  

Converters for data from other scanners are available from [https://github.com/ismrmrd](https://github.com/ismrmrd) but we have not tried these yet. 

### PET <a name="PET"></a>

STIR can handle data from the Siemens mMR Biograph with progress being made for the GE Signa PET/MR. However, STIR currently still relies on some bash scripts for file format conversion, estimation of randoms and scatter etc. Therefore, in the current SIRF release, we do not yet support measured data from any scanner. This will be fixed for version 1.0. 

# Where to find further information <a name="Further_information"></a>

- CCPPETMR website [http://www.ccppetmr.ac.uk](http://www.ccppetmr.ac.uk) for links to project overview, meeting notes, design documents etc 

- CCPPETMR SIRF Wiki [https://github.com/CCPPETMR/SIRF/wiki](https://github.com/CCPPETMR/SIRF/wiki)  for detailed instructions and user-contributed content 

- CCPPETMR Virtual Machine Wiki [https://github.com/CCPPETMR/CCPPETMR_VM/wiki](https://github.com/CCPPETMR/CCPPETMR_VM/wiki) with information on how to use the Virtual Machine that we supply with pre-installed software. 

- Inline documentation within MATLAB and Python functions, see below for examples. 

- Demo functions to demonstrate SIRF features. After installing SIRF, these will be available in SIRF/examples. 

- Our plan for future releases and additional features is available from the Software Documents tab [http://www.ccppetmr.ac.uk/softwareframework.html](http://www.ccppetmr.ac.uk/softwareframework.html)  on our website.

- Installation instructions can be found on our Wiki at [https://github.com/CCPPETMR/SIRF/wiki/Installation-instructions](https://github.com/CCPPETMR/SIRF/wiki/Installation-instructions). Note that on the [Virtual machine](https://github.com/CCPPETMR/CCPPETMR_VM/wiki), this has all has been done for you and you can just use SIRF. 

# General notes of usage <a name="General_notes"></a>

Please note that with the installation set-up, you will normally have two copies of the Matlab/Python module files: the original ones in the SIRF clone and the installed ones. This only matters if you want to debug or modify the files. The installation instructions point Python and Matlab to the “installed” files. 

The MR module and the demos create temporary files during operation. They are normally created in the same folder as the input data, but are cleaned up afterwards. Therefore,  the data cannot reside in a read-only folder. 
	
# Framework basic functionality <a name="Basic_functionality"></a>

## General conventions <a name="General_conventions"></a> 

### Object-oriented paradigm <a name="Object-oriented_paradigm"></a>

SIRF library modules are interfaces to object-oriented C++, which makes it reasonable for them to follow the object-oriented programming paradigm as well. This means that instead of having data containers (arrays, files etc.) and functions that operate on them, we employ objects, which contain data and come with sets of functions, called their _methods_, that operate on data. Each object contains a special method called constructor, which has the same name as the object class name and must be called to create that object. For example, to create an object of class ImageData that handles MR image data and fill it with data stored in the HDF5 file 'my_image.h5' one needs to do assignment 

    image = ImageData('my_image.h5'); 

We note that an MR ImageData object contains not only the voxel values, but also a number of parameters specified by the ISMRMRD format of MR image data. The object data is encapsulated, i.e. is not directly accessible from the user's code (being handled mostly by the underpinning C++ code) and is processed by the object methods. For example, to display the data encapsulated by image, one needs to call its method show(): 

    image.show(); 

and to copy the data into a Matlab array one uses method as_array(): 

    image_data_array = image.as_array(); 

Parameters of objects are modified/accessed via set/get methods (mutators and accessors). For example, the value of an objective function handled by object named obj_fun on an image data object image is computed by its method  get_value() as  

    obj_fun_value = obj_fun.get_value(image); 

The mutators are also responsible for basic error checking. 

Some classes are derived from other classes, which means that they have all the methods of the classes they are derived from: we say that these methods are _inherited_. For example, class `AcquisitionModelUsingRayTracingMatrix` is derived from `AcquisitionModelUsingMatrix`, which in turn is derived from `AcquisitionModel`, and so it inherits all the methods of the latter two.

### Error handling <a name="Error_handling"></a>

Error handling is via exceptions, i.e. functions do not return an error status, but throw an error if something did not work. The user can catch these exceptions if required as illustrated in the demos. 

### Naming conventions <a name="Naming_conventions"></a>

- Types/classes start with capitals, every word is capitalised, no underscores, e.g. AcquisitionModel. 

- Class methods are lower case, underscores between different words, e.g. get_voxel_size(). 

- Methods indicating  

    - a number of things start with num, e.g. num_gates. 

    - the number of an item in a sequence end with num, e.g. gate_num. 

### Units and index ordering <a name="Units_and_index_ordering"></a>

Distances are expressed in mm. 

For arrays in the target language, we use “native” ordering of indices in Python and Matlab. These are unfortunately opposite, so we would write  

    image_array[z,y,x] # Python 

    image_array(x,y,z) % Matlab 

For images, the meaning of z,y,x is currently acquisition dependent. Geometric information will be added later. 

### Handles <a name="Handles"></a>

In both Matlab and Python, SIRF operates with handles to objects, which affects the meaning of the assignment x = y: instead of creating a separate copy of y stored in x, x simply points to the same underlying data. As the result, any changes in x simultaneously change y. 

In order to have a true (i.e. independent) copy of a SIRF object, the user must call the object methods that create copies of them (see below). 
	
## Library components <a name="Library_components"></a>

At present, the SIRF library provides two Python interface modules pSTIR and pGadgetron for STIR and Gadgetron respectively, and two respective Matlab modules mSTIR and mGadgetron. 

### Getting help on SIRF library modules <a name="Getting_help_on_SIRF_library_modules"></a>

We remind that to see the contents of a Python module, the user needs to import it and use Python's help, and in Matlab one needs to use doc. For example,  

    # Python  
    import pSTIR 
    help(pSTIR) 

will show the components of pSTIR, and similarly 

    % Matlab 
    doc mSTIR 

will show the components of mGadgetron. In the same way,   

    # Python  
    help(pGadgetron.ImageData) 

will provide information on pGadgetron ImageData class, and  

    % Matlab 
    doc mSTIR.AcquisitionData  

on the mSTIR.AcquisitionData class. Regrettably, help and doc show all methods, including some common built-in methods such as __weakref__ method in Python or addlistener method in Matlab. Methods that are not related to SIRF is relatively easy to identify in Python (built-in methods have underscores in names). In Matlab they are difficult to identify, which is why we mark relevant Matlab methods other than constructors with \*\*\*SIRF\*\*\*. Methods not marked this way should be ignored. 

In order to understand the functionality of a derived class (see [Object-oriented paradigm](#Object-oriented_paradigm)), you are advised to first get help on the classes it is derived from. In Python, you can see that a class is derived by the presence of "Method resolution order" section in Python help output, which lists all classes it is derived from. You are advised to get help on all these classes except Python's class `builtins.object`. In Matlab, look at "Superclasses" item in "Class Details", and get help on the classes listed there except Matlab's class `handle`.

### General structure of the classes <a name="General_structure_of_the_classes"></a>

Most classes have a constructor to create an object from a file 

    image_data = ImageData(filename) 

and a method to create a copy of the object 

    a_copy = image_data.clone() 

“Processing” classes normally use the following pattern 

    recon.set_input(acquisition_data); 
    recon.setup(image_data); 
    recon.process(); 
    output_image_data=recon.get_output(); 

Classes follow a simple hierarchy, where top-level describes the generic functionality, and derived classes add/specify functionality. To see an example, look up Reconstructor and IterativeReconstructor classes in pSTIR or mSTIR using help or doc. We note that help(pSTIR.IterativeReconstructor) and doc mSTIR.IterativeReconstructor will show all the functionality of this class, i.e. including that of Reconstructor (and also some built-in functionality common to Python/Matlab classes). 

In what follows we use PET instead of pSTIR/mSTIR and MR instead of pGadgetron/mGadgetron to cover both Python and Matlab and also prospective alternative reconstruction engines. 

In the rest of the document we give basic information on the SIRF classes, including brief descriptions of the methods that are of interest to the user. Please use the inline help facility discussed above for more information. 

### Basic classes <a name="Basic_classes"></a>

#### Data Containers 

Reconstructed data are represented by ImageData objects. Currently they represent 3D volumes discretised using voxels.  

Measured data (either raw or after some pre-processing) are represented by AcquisitionData objects. These contain everything what is needed to be able to reconstruct the data (including scanner information and geometry). 

##### AcquisitionData

Class for acquisition data.

###### Methods:

    AcquisitionData  
               (PET/MR) Constructor. If no arguments are present, creates an
                        empty object, otherwise:
                        PET: Specifies the file containing raw data or
                        creates new AcquisitionData based on scanner information 
                        that comes either from a template AcquisitionData object
                        or discerned from the scanner name and parameters given
                        in the arguments; 
                        MR: Specifies the file containing raw data. 
    set_storage_scheme
               (PET/MR) Specifies whether the intermediate data should be kept in 
                        files or in RAM (see Acquisition data storage scheme 
                        management in Appendix).
    get_storage_scheme
               (PET/MR) Returns currently used storage scheme. 
    create_uniform_image  
                  (PET) Returns new compatible ImageData object. 
    as_array   (PET/MR) Returns the object data as an array. 
    fill       (PET/MR) Replaces the object data with user-supplied data. 
    clone      (PET/MR) Returns a copy of this object. 
    write      (PET/MR) Writes the object data to a file. 
    sort           (MR) Sort the object data. 
    is_sorted      (MR) Returns true if and only if the object data is sorted. 
    get_info       (MR) Returns information on the object data. 
    process        (MR) Processes the object data by a chain of gadgets. 
    dimensions     (MR) Returns the object data dimensions

##### ImageData

Class for data representing 3D objects.

###### Methods:

    ImageData   (PET)  Constructor. Reads data from a file or creates empty object. 
    initialise  (PET)  Sets the image size in voxels, voxel sizes and the origin. 
    fill     (PET/MR)  Replaces the object data with user-supplied data. 
    as_array (PET/MR)  Returns the object data as an array. 
    clone    (PET/MR)  Returns a copy of this object. 
    read_from_file
             (PET/MR)  Reads the image data from file.
    get_uniform_copy   
                (PET)  Returns a copy of this image filled with a constant value. 
    add_shape   (PET)  Adds a uniform shape to the image. 
    show     (PET/MR)  Interactively displays the image. 
    write    (PET/MR)  Writes the object data to a file. 
    dimensions  (PET)  Returns the object data dimensions
    voxel_sizes (PET)  Returns the voxel sizes
	
##### CoilSensitivityData (MR)

Class for storing coil sensitivity maps.

###### Methods:

    CoilSensitivityData  Constructor. Creates empty object.
    calculate            Calculates coil sensitivities from the acquisition data 
                         Specified by the argument. 
    csm_as_array         Returns the coil sensitivity map for the slice/repetition 
                         specified by the argument as an array. 	
					 
##### Examples:

    PET_image = ImageData('image.hv'); % read image data from a file 
    PET_image0 = ImageData(); % create empty image object 
    PET_image0.initialise([128,128,31], [3,3,3.375]); % in Python: (128,128,31) etc. 
    PET_image0.fill(1.0); % assign value 1.0 at each voxel 
    PET_image_array = PET_image.as_array(); % copy image data to a Matlab array 
 
    MR_image_array = MR_image.as_array(); % copy image data to a Matlab array 
    MR_acquisition_data = AcquisitionData('mr_raw_data.h5'); 
    cs_data = CoilSensitivityData(); % create empty object 
    cs_data.calculate(MR_acquisition_data); % calculate coil sensitivities 
    csm0 = cs_data.csm_as_array(0); % obtain coil sensitivities for slice 0 as array 

#### Data Processors

##### ImageDataProcessor

Class for objects that process ImageData objects. 

###### Methods:

    ImageDataProcessor  
               (PET/MR) Constructor. Creates new ImageDataProcessor object 
                        (PET: empty, MR: defined by the argument). 
    set_input  (PET/MR) Sets the processor input. 
    process    (PET/MR) Processes image data on input. 
    get_output (PET/MR) Retrieves the processed image data. 
    apply         (PET) Processes the ImageData argument. 

##### TruncateToCylinderProcessor (PET)

Class for the image processor that zeroes the image outside a cylinder. Inherits the methods of ImageDataProcessor. 

###### Methods (in addition to those of ImageDataProcessor): 

    set_strictly_less_than_radius  Defines the behaviour on the cylinder boundary. 
    get_strictly_less_than_radius  Exposes the behaviour on the cylinder boundary. 

##### AcquisitionDataProcessor

Class for objects that process AcquisitionData objects. 

###### Methods:

    AcquisitionDataProcessor  
               (MR) Constructor. Creates new processor object defined by the argument. 
    set_input  (MR) Sets the processor input. 
    process    (MR) Processes the image data on input. 
    get_output (MR) Retrieves the processed image data. 

###### Examples: 

    filter = TruncateToCylinderProcessor(); 
    filter.apply(PET_image); 
    img_proc = ImageDataProcessor({'ExtractGadget'}); % Python: ['ExtractGadget'] 
    img_proc.set_input(MR_image);  
    img_proc.process();  
    MR_image_magnitude = img_proc.get_output(); 
    acq_proc.set_input(MR_acquired_data);  
    acq_proc.process();  
    preprocessed_data = acq_proc.get_output(); 

#### Reconstructors 

##### Reconstructor 

Class for a generic image reconstructor. 

###### Methods: 

    set_input  (PET/MR) Sets the input (AquisitionData object). 
    process    (PET/MR) Runs the reconstruction. 
    get_output (PET/MR) Returns the output (ImageData object). 
    set_output_filename_prefix  
               (PET/MR) Specifies the naming for the output files.  

##### IterativeReconstructor (PET) 

Class for PET reconstruction algorithms that use Ordered Subsets technique whereby the acquisition data is split into subsets, and the objective function and its gradient are represented as the sums of components corresponding to subsets. Typically, one iteration of such algorithm would deal with one subset, and is therefore referred to as sub-iteration. Inherits the methods of Reconstructor. 

###### Methods (in addition to those of Reconstructor): 

    set_num_subsets            Sets the number of subsets, 
    get_num_subsets            Returns the number of subsets. 
    set_num_subiterations      Sets the number of subiterations. 
    get_num_subiterations      Returns the number of subiterations. 
    get_subiterations_num      Returns the current subiteration number. 
    set_save_interval          Specifies how often to save image estimates. 
    set_objective_function     Specifies the objective function. 
    set_up                     Prepares the reconstructor for use. 
    set_current_estimate       Sets the current image estimate. 
    get_current_estimate       Returns the current image estimate.  
    update_current_estimate    Updates the current image estimate. 
    set_current_subset_num     Specifies the current subset number. 
    get_subset_sensitivity     Returns sensitivity image for the current subset. 
    reconstruct                Reconstructs using the argument as initial image. 
    process                    Reconstructs using current image estimate as initial. 
    update                     Updates using the argument as current image estimate. 

##### OSMAPOSLReconstructor (PET) 

Class for reconstructor objects using Ordered Subsets Maximum A Posteriori One Step Late reconstruction algorithm, see [http://stir.sourceforge.net/documentation/doxy/html/classstir_1_1OSMAPOSLReconstruction.html](http://stir.sourceforge.net/documentation/doxy/html/classstir_1_1OSMAPOSLReconstruction.html). Inherits the methods of IterativeReconstructor. 

###### Methods (in addition to those of IterativeReconstructor): 

    OSMAPOSLReconstructor  Constructor. Creates new reconstructor object.  

##### OSSPSReconstructor (PET) 

Class for reconstructor objects using Ordered Subsets Separable Paraboloidal Surrogate reconstruction algorithm, see [http://stir.sourceforge.net/documentation/doxy/html/classstir_1_1OSSPSReconstruction.html](http://stir.sourceforge.net/documentation/doxy/html/classstir_1_1OSSPSReconstruction.html). Inherits the methods of IterativeReconstructor.  

###### Methods (in addition to those of IterativeReconstructor): 

    OSSPSReconstructor   Constructor. Creates new reconstructor object. 

##### FullySampledReconstructor (MR) 

Class for a reconstructor from fully sampled Cartesian raw data. Inherits the methods of Reconstructor. 

###### Methods (in addition to those of Reconstructor): 

    FullySampledReconstructor  Constructor. Creates new reconstructor object.  

##### CartesianGRAPPAReconstructor (MR) 

Class for a reconstructor from undersampled Cartesian raw data. Inherits the methods of Reconstructor. 

###### Methods (in addition to those of Reconstructor): 

    CartesianGRAPPAReconstructor  Constructor. Creates new reconstructor object. 

### Other classes <a name="Other_classes"></a>

##### ListmodeToSinograms (PET)

Class for converting raw data from listmode format into *sinograms*,
i.e. histogrammed data in the format of PETAcquisitionData.

It has 2 main functions:
  - `process()` can be used to read prompts and/or delayed coincidences to produce a single
    PETAcquisitionData. 2 variables decide what done with 3 possible cases:
       - `store_prompts`=`true`, `store_delayeds`=`false`: only prompts are stored
       - `store_prompts`=`false`, `store_delayeds`=`true`: only delayeds are stored
       - `store_prompts`=`true`, `store_delayeds`=`true`: prompts-delayeds are stored
    Clearly, enabling the `store_delayeds` option only makes sense if the data was 
    acquired accordingly.
  - estimate_randoms() can be used to get a relatively noiseless estimate of the 
    random coincidences.

###### Methods:

    ListmodeToSinograms  Constructor. Takes an optional text string argument with
                         the name of a STIR parameter file defining the conversion options.
                         If no argument is given, default settings apply except
                         for the names of input raw data file, template file and
                         output filename prefix, which must be set by the user by
                         calling respective methods.
    set_input            Specifies the input raw data file.
    set_output_prefix    Specifies the prefix for the output file(s), which will
                         be appended by `_g1f1d0b0.hs`.
    set_template         Specifies the file containing acquisition data to be
                         used as a source of information about the scanner.
    set_time_interval    Specifies the scanning time sub-interval to be converted
                         (an empty interval indicates that all raw data must be converted)
    flag_on              Turns on (i.e. assigns value true to) a conversion flag.
    flag_off             Turns off (i.e. assigns value false to) a conversion flag.
    set_up               Sets up the converter.
    process              Performs the conversion.
    get_output           Returns AcquisitionData object containing converted data.
    estimate_randoms     Estimates randoms. (Currently via a Maximum Likelihood estimate
                         of the singles, based on the delayed coincidences).

###### Examples: 

    lm2sino = ListmodeToSinograms()
    lm2sino.set_input(list_file)
    lm2sino.set_output_prefix(sino_file)
    lm2sino.set_template(tmpl_file)
    lm2sino.set_time_interval(0, 10)
    lm2sino.flag_on('store_prompts')

    lm2sino.set_up()
    lm2sino.process()
    acq_data = lm2sino.get_output()
    randoms_acq_data = lm2sino.estimate_randoms()

##### AcquisitionModel 

Class for the acquisition process modelling. Main component is the forward projection operation F that for a given image data x estimates the data y = F(x) to be acquired by the scanner (simulated acquisition data). The transpose B of the Frechet derivative of F is referred to as backprojection (if F is linear, e.g. a matrix, then B is the transpose of F).

For PET, F(x) is the right-hand side of the following equation:

    (F)    y = (1/n)(G x + a) + b 

where  

G, ray tracing matrix, is a matrix whose columns correspond to the image voxels and rows to pairs of scanner's detectors (bins), each column simulating the impact of this voxel's radiation on the data acquired by the bins; 

a and b, additive and background terms, represent the effects of accidental coincidences and scattering; 

n, bin normalization, is the inverse of bin efficiencies. 

Accordingly, the backprojection B is given by

    (B)    x = G' (1/n) y 

where G' is the transpose of G. 

###### Methods: 

    AcquisitionModel (PET/MR) Constructor. Creates an acquisition model 
                              (PET: empty, MR: based on the image and 
                              acquisition data templates specified by the 
                              arguments). 
    forward          (PET/MR) Returns F(x) for the image data x specified 
                              by the argument. 
    backward         (PET/MR) Returns B(y) for the acquisition data y specified 
                              by the argument. 
    set_up              (PET) Sets up the model based on acquisition and image data  
                              templates provided by the arguments. 
    set_additive_term   (PET) Sets term a in (F). 
    set_acquisition_sensitivity   
                        (PET) Defines n in (F) via AcquisitionSensitivityModel
                              (see below). 
    set_coil_sensitivity_maps  
                         (MR) Sets coil sensitivity maps to be used.  

###### Examples: 

    MR_model = AcquisitionModel(acq_template, image_template); 
    MR_model.set_coil_sensitivity_maps(cs_data); 
    sim_data = MR_model.forward(MR_image);

##### AcquisitionModelUsingRayTracingMatrix (PET) 

Class for the PET acquisition process model that uses (implicitly) a sparse matrix for G in (F). This class inherits the methods of PET AcquisitionModel class, with forward projection defined by (F) and backprojection by (B).

###### Methods (in addition to those of AcquisitionModel): 

    AcquisitionModelUsingRayTracingMatrix  
                          Constructor. Creates an acquisition model. 

###### Examples: 

    acq_model = AcquisitionModelUsingRayTracingMatrix(); 
    acq_model.set_up(acq_template, image_template) 
    sim_data = acq_model.forward(image); 

##### AcquisitionSensitivityModel (PET)

Class for a part of AcquisitionModel that accounts for bin efficiencies and attenuation.
Provides methods for for applying (1/n) factor in (F) and (B) or its inverse. 

###### Methods: 

    AcquisitionSensitivityModel 
                     Constructor. Creates a new object of this class
                     - from an ECAT8 file or
                     - from an attenuation image (ImageData object) or
                     - from bin efficiencies (AcquisitionData object) or
                     - by chaining two objects of this class.
                     In the last case, the normalisation n is the product 
                     of the two objects' normalisations.

    set_up           Sets up the object.
    normalise        Multiplies the argument (AcquisitionData) by n.
    unnormalise      Multiplies the argument (AcquisitionData) by (1/n).
    forward          Returns the argument multiplied by (1/n). The argument
                     is not changed.
    invert           Returns the argument multiplied by n. The argument
                     is not changed.

###### Examples: 

    # obtain an acquisition data template
    template = AcquisitionData(temp_file)

    # create acquisition sensitivity model from ECAT8 normalization data
    asm = AcquisitionSensitivityModel(norm_file)

    # create acquisition sensitivity model from attenuation image
    attn_image = ImageData(attn_file)
    am = AcquisitionModelUsingRayTracingMatrix()
    am.set_up(template, attn_image)
    asm = AcquisitionSensitivityModel(attn_image, am)

    # create acquisition sensitivity model from bin efficiencies
    asm = AcquisitionSensitivityModel(bin_eff)

##### ObjectiveFunction (PET) 

Class for objective functions maximized by iterative Ordered Subsets reconstruction algorithms. At present we use Poisson logarithmic likelihood function with linear model for mean and a specific arrangement of the acquisition data. To make our interface more user-friendly, we provide a convenience function make_PoissonLogLikelihood that creates objects of this class (instead of the usual constructor) based on the acquisition data to be used. 

The user have an option of adding a penalty term (referred to as prior) to the objective function. At present, we have just one particular kind of prior implemented, the quadratic prior described in the next section. 

###### Methods: 

    set_prior        Specifies the prior. 
    set_num_subsets  Specifies the number of subsets. 
    set_up           Prepares this object for use. 
    get_value        Returns the value of the objective function. 
    get_gradient     Returns the gradient of the objective function. 
    get_subset_gradient 
                     Returns the component of the gradient for the specified subset. 
    get_backprojection_of_acquisition_ratio 
                     Returns the backprojection of the ratio of measured to estimated 
                     acquisition data. 
    set_acquisition_model 
                     Specifies the acquisition model to be used. 
    set_acquisition_data 
                     Specifies the acquisition data to be used.  

##### QuadraticPrior (PET)

Class for a penalty term to be added to the objective function. 

###### Methods: 

    set_penalisation_factor  Specifies the prior's scaling factor. 
    get_gradient             Returns the prior gradient.  

### Functions <a name="Functions"></a>

    preprocess_acquisition_data (MR)  Preprocesses the MR acquisition data.  

    make_Poisson_loglikelihood (PET)  Returns Poisson objective function.

# Appendix <a name="Appendix"></a>

## Acquisition data storage scheme management <a name="storage_management"></a>

SIRF offers the users two options in handling acquisition data generated by SIRF scripts. The default one keeps all acquisition data generated by the script in scratch files deleted after the script terminates. An alternative one stores acquisition data in memory. The user can specify the storage scheme to be employed by calling a static method `set_storage_scheme` of `AcquisitionData` class:

    AcquisitionData.set_storage_scheme(scheme)
    
where `scheme` is either `"memory"` or `"file"`. To see which scheme is currently used, use

    scheme = AcquisitionData.get_storage_scheme()
    
A particular setting of storage scheme by a Matlab script or a Python script run from Spyder is persistent: any script run afterwards will use the same storage scheme unless a different storage scheme is explicitly set by `set_storage_scheme` or Matlab/Spyder is re-started.

## Programming chains of Gadgetron gadgets <a name="programming_Gadgetron_chains"></a>

Gadgetron is a MR reconstruction framework which was designed to process a datastream, i.e. rather than waiting for a complete 3D k-space to be acquired, each readout (frequency encoding line) is processed immidiately (if possible, e.g. Fourier transform along phase encoding can only be applied once all phase encoding lines have been acquired) (https://github.com/gadgetron/gadgetron/wiki/What-Is-The-Gadgetron). The reconstruction is performed by a chain of gadgets, i.e. pieces of code implementing specific tasks. The chain of gadgets runs on the server, which can be just a command line window, or it can be another computer or a VM. In order to set up the chain, the server needs to receive an xml text describing it from the client, which again can be another command line window on the same or another computer. The first gadget in the chain then starts waiting for acquisition data to arrive from the client in chunks of certain size. Having processed a chunk of data, the first gadget passes the result to the second and starts processing the next chunk and so on. The last gadget sends the reconstructed images back to the client.

### Creating and running gadget chains by SIRF script  <a name="creating_and_running_gadget_chains"></a>

The standard way of using Gadgetron is to run `gadgetron_ismrmrd_client` from a command line (with Gadgetron running in another terminal window), providing the name of the raw data file (in HDF5 format) and the name of the xml file containing the description of the gadget chain via command-line options. SIRF offers an equivalent alternative whereby the data and the gadget chain are defined in a Python or Matlab script. The gadget chain is defined by creating a Reconstructor object and providing the list of gadgets descriptions as an argument:

    my_recon = Reconstructor(my_gadget_list);

Here `my_gadget_list` is a list of strings in Python or a cell array of strings in Matlab, each string describing a gadget in the following format:

    [label:]gadget_name[(property1=value1[,property2=value2,...])]

where the names of the gadget and its properties are same as in Gadgetron xml files, and an optional label can be used to change the labelled gadget properties at any time by using `set_gadget_property` method:

    my_recon.set_gadget_property(label, property, value);

The following example of a gadget chain definition is taken from demo script `fully_sampled_recon_single_chain.py`:

    recon = Reconstructor(['RemoveROOversamplingGadget', \
        'AcquisitionAccumulateTriggerGadget(trigger_dimension=repetition)', \
        'BucketToBufferGadget(split_slices=true, verbose=false)', \
        'SimpleReconGadget', 'ImageArraySplitGadget', 'ex:ExtractGadget'])
    # ExtractGadget defines which type of image should be returned:
    # none      0
    # magnitude 1
    # real      2
    # imag      4
    # phase     8
    # max       16
    # in this example '5' returns both magnitude and imag
    recon.set_gadget_property('ex', 'extract_mask', 5)

The input data is defined by creating an AcquisitionData object and passing it to the reconstruction object via its method `set_input`:

    acq_data = AcquisitionData(input_file_name);
    my_recon.set_input(acq_data);

and the reconstruction is performed by calling the method `process`, and the reconstructed images are returned as an ImageData object by the method `get_output`:

    my_recon.process();
    image_data = my_recon.get_output();

While the way to use Gadgetron just described is the most efficient performance-wise, some users may like to get more involved in the reconstruction process. SIRF offers such users an opportunity to split a standard reconstruction chain into sub-chains and process intermediate data. For example, the chain defined above can be split into an acquisition processing chain that removes oversampling, a shortened reconstruction chain and an image processing chain:

    acq_proc = AcquisitionDataProcessor(['RemoveROOversamplingGadget'])
    acq_proc.set_input(acq_data)
    acq_proc.process()
    preprocessed_data = acq_proc.get_output()

    # do something with preprocessed data here

    recon = Reconstructor\
        (['AcquisitionAccumulateTriggerGadget(trigger_dimension=repetition)', \
        'BucketToBufferGadget(split_slices=true, verbose=false)',
        'SimpleReconGadget', 'ImageArraySplitGadget'])
    recon.set_input(preprocessed_data)
    recon.process()
    complex_image_data = recon.get_output()

    # do something with the complex image data here

    img_proc = ImageDataProcessor(['ExtractGadget(extract_mask=1)'])
    img_proc.set_input(complex_image_data)
    img_proc.process()
    real_image_data = img_proc.get_output()


### SIRF gadget library <a name="SIRF_gadget_library"></a>

This section provides a concise description of Gadgetron gadgets that can be used by current SIRF release scripts in a way described in the previous section. For further information consult Gadgetron documentation.

Below `internal<N>` refers to Gadgetron data objects to which SIRF does not provide interface at present. We emphasize that splitting Gadgetron chains into sub-chains in a way described in the previous section only makes sense if the input of the first gadget and the output of the last gadget of each sub-chain are either AcquisitionData or ImageData.

#### RemoveROOversamplingGadget

| input | output | parameters |
| - | - | - |
| AcquisitionData | AcquisitionData | none |

Removes the oversampling along the readout direction.

#### NoiseAdjustGadget

| input | output | parameters |
| - | - | - |
| AcquisitionData | AcquisitionData | none |

Ensures that the noise between different receiver coils is not correlated and that each receiver coils has a similar noise level.

#### AsymmetricEchoAdjustGadget

| input | output | parameters |
| - | - | - |
| AcquisitionData | AcquisitionData | none |

Pads each readout with zeros to compensate for partial echo acquisitions.

#### AcquisitionAccumulateTriggerGadget

| input | output | parameters | default values |
| - | - | - | - |
| AcquisitionData | internal1 | trigger_dimension | "repetition" |
| | | sorting_dimension | "slice" |

Collects lines of k-space until a certain trigger condition is encountered, i.e., when there is enough data to reconstruct an image.

#### BucketToBufferGadget

| input | output | parameters | default values |
| - | - | - | - |
| internal1 | internal2 | N_dimension | "" |
| | | S_dimension | "" |
| | | split_slices | "true" |
| | | ignore_segment | "true" |
| | | verbose | "true" |

Inserts the collected data into a buffer more suitable for recon processing.

#### SimpleReconGadget

| input | output | parameters |
| - | - | - |
| internal2 | internal3 | none |

Performs simple fast Fouriertransforms to transform acquired k-space data to image space.

#### GenericReconCartesianReferencePrepGadget

| input | output | parameters | default values |
| - | - | - | - |
| internal2 | internal4 | debug_folder | "" |
| | | perform_timing | "true" |
| | | verbose | "true" |
| | | average_all_ref_N | "true" |
| | | average_all_ref_S | "true" |
| | | prepare_ref_always | "true" |

Selects the reference data used to calculate the GRAPPA kernel

#### GenericReconCartesianGrappaGadget

| input | output | parameters | default values |
| - | - | - | - |
| internal4 | internal5 | debug_folder | ""
| | | perform_timing | "true" |
| | | verbose | "true" |
| | | image_series | "0" |
| | | coil_map_algorithm | "Inati" |
| | | send_out_gfactor | "true" |
| | | downstream_coil_compression | "true" |
| | | downstream_coil_compression_thres | "0.01" |
| | | downstream_coil_compression_num_modesKept | "0" |

Performs GRAPPA kernel calibration, calculate coil sensitivity maps and carry out unfolding.

#### GenericReconFieldOfViewAdjustmentGadget

| input | output | parameters | default values |
| - | - | - | - |
| internal5 | internal6 | debug_folder | "" |
| | | perform_timing | "false" |
| | | verbose | "false" |

Adjusts FOV and image resolution according to the parameters given for the reconstructed image in the file header.

#### GenericReconImageArrayScalingGadget

| input | output | parameters | default values |
| - | - | - | - |
| internal6 | internal3 | perform_timing | "false" |
| | | verbose | "false" |
| | | min_intensity_value | "64" |
| | | max_intensity_value | "4095" |
| | | scalingFactor | "10.0" |
| | | scalingFactor_dedicated | "100.0" |
| | | use_constant_scalingFactor | "true" |
| | | auto_scaling_only_once | "true" |

Applie scaling to image, g-factor map, SNR map and/or SNR standard deviation map.

#### ImageArraySplitGadget

| input | output | parameters |
| - | - | - |
| internal3 | ImageData | none |

Splits array of images (7D [X, Y, Z, CHA, N, S, LOC]) into individual images (4D [X, Y, Z, CHA])

#### ExtractGadget

| input | output | parameters | default values |
| - | - | - | - |
| ImageData | ImageData | extract_mask | "1" |

Extracts a certain type of image data from the reconstructed image stream, i.e. extract_mask=1 yields magnitude images, extract_mask=2 yields readl images.

#### ComplexToFloatGadget

| input | output | parameters |
| - | - | - |
| ImageData | ImageData | none |

Depending on the image type, the magnitude, real, imaginary or phase of the complex image is returned.

#### FloatToShortGadget

| input | output | parameters | default values |
| - | - | - | - |
| ImageData | ImageData | min_intensity | "0" |
| | | max_intensity | "32767" |
| | | intensity_offset | "0" |

Scales and transforms images from float to short.

#### SimpleReconGadgetSet

A chain of Gadgetron gadgets
~~~
AcquisitionAccumulateTriggerGadget
BucketToBufferGadget
SimpleReconGadget
ImageArraySplitGadget
~~~
for fully sampled reconstruction.

| input | output | parameters | default values |
| - | - | - | - |
| AcquisitionData | ImageData | N_dimension | "" |
| | | S_dimension | "" |
| | | sorting_dimension | "slice" |
| | | trigger_dimension | "repetition" |
| | | split_slices | "true" |

