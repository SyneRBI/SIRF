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

# Overview <a name="Overview"></a>

The SIRF (Synergistic Image Reconstruction Framework) software is an Open Source toolkit for the reconstruction of PET and MRI raw data. The aim is to provide code simple enough to easily perform a reconstruction, yet powerful enough to be able to handle real, full-size datasets. Our strategy in achieving this aim is to employ available Open Source reconstruction software written in advanced programming languages such as C++ and provide basic-user-friendly interfaces to it written in script languages, primarily Matlab and Python. The interface style permits a reconstruction to be performed in stages, allowing the user to inspect or modify data, or insert their own code. 

This User’s Guide describes version 0.9 of SIRF. The software can be found on [https://github.com/CCPPETMR](https://github.com/CCPPETMR).

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

# Installation instructions <a name="Installation_instructions"></a>

Brief installation instructions can be found on [https://raw.githubusercontent.com/CCPPETMR/SIRF/master/INSTALL.txt](https://raw.githubusercontent.com/CCPPETMR/SIRF/master/INSTALL.txt). More detailed instructions for specific Operating Systems etc are on our Wiki at [https://github.com/CCPPETMR/SIRF/wiki/Installation-instructions](https://github.com/CCPPETMR/SIRF/wiki/Installation-instructions). 

Note that on the [Virtual machine](https://github.com/CCPPETMR/CCPPETMR_VM/wiki), this has all be done for you and you can just use SIRF. 

# General notes of usage <a name="General_notes"></a>

Please note that with the installation set-up, you will normally have two copies of the Matlab/Python module files: the original ones in the SIRF clone and the installed ones. This only matters if you want to debug or modify the files. The installation instructions point Python and Matlab to the “installed” files. 

The MR module and the demos create temporary files during operation. They are normally created in the same folder as the input data, but are cleaned up afterwards. Therefore,  the data cannot reside in a read-only folder. 
	
# Framework basic functionality <a name="Basic_functionality"></a>

## General conventions <a name="General_conventions"></a> 

### Object-oriented paradigm <a name="Object-oriented_paradigm"></a>

SIRF library modules are interfaces to object-oriented C++, which makes it reasonable for them to follow the object-oriented programming paradigm as well. This means that instead of having data containers (arrays, files etc.) and functions that operate on them, we employ objects, which contain data and come with sets of functions, called their methods, that operate on data. Each object contains a special method called constructor, which has the same name as the object class name and must be called to create that object. For example, to create an object of class ImageData that handles MR image data and fill it with data stored in the HDF5 file 'my_image.h5' one needs to do assignment 

    image = ImageData('my_image.h5'); 

We note that an MR ImageData object contains not only the voxel values, but also a number of parameters specified by the ISMRMRD format of MR image data. The object data is encapsulated, i.e. is not directly accessible from the user's code (being handled mostly by the underpinning C++ code) and is processed by the object methods. For example, to display the data encapsulated by image, one needs to call its method show(): 

    image.show(); 

and to copy the data into a Matlab array one uses method as_array(): 

    image_data_array = image.as_array(); 

Parameters of objects are modified/accessed via set/get methods (mutators and accessors). For example, the value of an objective function handled by object named obj_fun on an image data object image is computed by its method  get_value() as  

    obj_fun_value = obj_fun.get_value(image); 

The mutators are also responsible for basic error checking. 

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

At present, the SIRF library provides two Python interface modules pStir and pGadgetron for STIR and Gadgetron respectively, and two respective Matlab modules mStir and mGadgetron. 

### Getting help on SIRF library modules <a name="Getting_help_on_SIRF_library_modules"></a>

We remind that to see the contents of a Python module, the user needs to import it and use Python's help, and in Matlab one needs to use doc. For example,  

    # Python  
    import pStir 
    help(pStir) 

will show the components of pStir, and similarly 

    % Matlab 
    doc mStir 

will show the components of mGadgetron. In the same way,   

    # Python  
    help(pGadgetron.ImageData) 

will provide information on pGadgetron ImageData class, and  

    % Matlab 
    doc mStir.AcquisitionData  

on the mStir.AcquisitionData class. Regrettably, help and doc show all methods, including some common built-in methods such as __weakref__ method in Python or addlistener method in Matlab. Methods that are not related to SIRF is relatively easy to identify in Python (built-in methods have underscores in names). In Matlab they are difficult to identify, which is why we mark relevant Matlab methods other than constructors with \*\*\*SIRF\*\*\*. Methods not marked this way should be ignored. 

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

Classes follow a simple hierarchy, where top-level describes the generic functionality, and derived classes add/specify functionality. To see an example, look up Reconstructor and IterativeReconstructor classes in pStir or mStir using help or doc. We note that help(pStir.IterativeReconstructor) and doc mStir.IterativeReconstructor will show all the functionality of this class, i.e. including that of Reconstructor (and also some built-in functionality common to Python/Matlab classes). 

In what follows we use PET instead of pStir/mStir and MR instead of pGadgetron/mGadgetron to cover both Python and Matlab and also prospective alternative reconstruction engines. 

In the rest of the document we give basic information on the SIRF classes, including brief descriptions of the methods that are of interest to the user. Please use the inline help facility discussed above for more information. 

### Basic classes <a name="Basic_classes"></a>

#### Data Containers 

Reconstructed data are represented by ImageData objects. Currently they represent 3D volumes discretised using voxels.  

Measured data (either raw or after some pre-processing) are represented by AcquisitionData objects. These contain everything what is needed to be able to reconstruct the data (including scanner information and geometry). 

##### AcquisitionData

Class for acquisition data.

###### Methods:

    AcquisitionData  
               (PET/MR) Constructor. Reads data from a file or creates empty object. 
    create_uniform_image  
                  (PET) Returns new compatible ImageData object. 
    as_array   (PET/MR) Returns the object data as an array. 
    fill       (PET/MR) Replaces the object data with user-supplied data. 
    clone      (PET/MR) Returns a copy of this object. 
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
    add_shape   (PET)  Adds a uniform shape to the image. 
    show     (PET/MR)  Interactively displays the image. 
    write    (PET/MR)  Writes the object data to a file. 
	
##### CoilSensitivityData (MR)

Class for storing coil sensitivity maps.

###### Methods:

    CoilSensitivityData  Constructor. Creates empty object.
    calculate            Calculates coil sensitivities from the acquisition data 
                         Specified by the argument. 
    csm_as_array         Returns the coil sensitivity map for the slice/repetition 
                         specified by the argument as an array. 	
					 
###### Examples:

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

##### AcquisitionModel 

Class for the acquisition process modelling. Main component is the forward projection operation F that for a given image data x estimates the data y = F(x) to be acquired by the scanner (simulated acquisition data). The transpose B of the Frechet derivative of F is referred to as backprojection (if F is linear, e.g. a matrix, then B is the transpose of F). 

###### Methods: 

    AcquisitionModel  
             (PET/MR) Constructor. Creates an acquisition model (PET: empty, MR: based 
                      on the image and acquisition data templates specified by the 
                      arguments). 
    forward  (PET/MR) Returns F(x) for the image data x specified by the argument. 
    backward (PET/MR) Returns B(y) for the acquisition data y specified by the  
                      argument. 
    set_coil_sensitivity_maps  
                 (MR) Sets coil sensitivity maps to be used.  

###### Examples: 

    MR_model = AcquisitionModel(acq_template, image_template); 
    MR_model.set_coil_sensitivity_maps(cs_data); 
    sim_data = MR_model.forward(MR_image);

##### AcquisitionModelUsingRayTracingMatrix (PET) 

Class for the PET acquisition process model whereby the image data x is related to the acquisition data y as 

    (F)    y = (1/n)(G x + a) + b 

where  

G, ray tracing matrix, is a matrix whose columns correspond to the image voxels and rows to pairs of scanner's detectors (bins), each column simulating the impact of this voxel's radiation on the data acquired by the bins; 

a and b, additive and background terms, represent the effects of accidental coincidences and scattering (b is not implemented yet); 

n, bin normalization, is the inverse of bin efficiencies. 

This class inherits the methods of PET AcquisitionModel class, with forward projection defined by (F) and backprojection by 

    (B)    x = G' (1/n) y 

where G' is the transpose of G. 

###### Methods (in addition to those of AcquisitionModel): 

    AcquisitionModelUsingRayTracingMatrix  
                          Constructor. Creates an acquisition model. 
    set_up                Sets up the model based on acquisition and image data  
                          templates provided by the arguments. 
    set_additive_term     Sets term a in (F). 
    set_normalisation     Sets n in (F).  
    set_bin_efficiency    Sets 1/n in (F). This is useful if some bin efficiencies 
                          are zero. 

###### Examples: 

    acq_model = AcquisitionModelUsingRayTracingMatrix(); 
    acq_model.set_up(acq_template, image_template) 
    sim_data = acq_model.forward(image); 

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
