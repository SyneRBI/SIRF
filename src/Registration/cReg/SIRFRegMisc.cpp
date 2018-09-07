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
\brief Generic tools (e.g., opening files). Currently dependent on NiftyReg, could cut that dependence?

\author Richard Brown
\author CCP PETMR
*/

#include "SIRFRegMisc.h"
#include "SIRFImageData.h"
#include "SIRFImageDataDeformation.h"
#include "SIRFRegTransformation.h"
#include <_reg_tools.h>
#if NIFTYREG_VER_1_5
#include <_reg_globalTrans.h>
#include <_reg_localTrans.h>
#include <_reg_ReadWriteMatrix.h>
#elif NIFTYREG_VER_1_3
#include <_reg_localTransformation.h>
#endif
#include <iostream>
#include <iomanip>
#include <sstream>

using namespace std;

namespace SIRFRegMisc {

/// Open nifti image
void open_nifti_image(std::shared_ptr<nifti_image> &image, const boost::filesystem::path filename)
{
    // If no filename has been set, return
    if (filename == "") {
        throw runtime_error("Empty filename has been supplied, cannot open nifti image.");
    }

    // Check that the file exists
    if (!boost::filesystem::exists(filename)) {
        throw runtime_error("Cannot find the file: " + filename.string() + ".");
    }

    // Check that file is nifti
    if (is_nifti_file(filename.c_str()) == -1) {
        throw runtime_error("Attempting to open a file that is not a NIFTI image.\n\tFilename: " + filename.string());
    }

    // Open file
    image = make_shared<nifti_image>(*nifti_image_read(filename.c_str(), 1));

    // Ensure the image has all the values correctly set
    reg_checkAndCorrectDimension(image.get());
}

/// Save nifti image
void save_nifti_image(nifti_image *image, const string filename)
{
    if (!image)
        throw runtime_error("Cannot save image to file.");

    cout << "\nSaving image to file (" << filename << ")..." << flush;

    boost::filesystem::path filename_boost(filename);

    // If the folder doesn't exist, create it
    if (!boost::filesystem::exists(filename_boost.parent_path())) {
        if (filename_boost.parent_path().string() != "") {
            cout << "\n\tCreating folder: \"" << filename_boost.parent_path().string() << "\"\n" << flush;
            boost::filesystem::create_directory(filename_boost.parent_path());
        }
    }

    nifti_set_filenames(image, filename.c_str(), 0, 0);
    nifti_image_write(image);
    cout << "done.\n\n";
}

/// Save nifti image
void save_nifti_image(std::shared_ptr<nifti_image> image, const string filename)
{
    save_nifti_image(image.get(), filename);
}

/// Split multi-component image
vector<std::shared_ptr<nifti_image> >
    split_multicomponent_nifti_image(std::shared_ptr<nifti_image> input_sptr)
{
    // Only works for ndim==5
    if (input_sptr->ndim != 5)
        throw runtime_error("Splitting only currently works for ndim==5.");

    // Create the vector to store the single component images
    vector<std::shared_ptr<nifti_image> > output;

    // Loop over all of the components
    int num_components = input_sptr->dim[5];
    for (int component=0; component<num_components; component++) {

        // Create new image
        nifti_image *image_ptr;

        // Copy the input image
        image_ptr = nifti_copy_nim_info(input_sptr.get());

        // Alter the info to change the number of dimensions
        image_ptr->dim[0] = image_ptr->ndim = 3;
        image_ptr->dim[5] = image_ptr->nu   = 1;
        image_ptr->nvox   = image_ptr->nvox / input_sptr->nu;

        // How much memory do we need?
        size_t mem = image_ptr->nvox*image_ptr->nbyper;

        // Allocate the data
        image_ptr->data=static_cast<void *>(malloc(mem));

        // Start index
        size_t index = mem*component;

        // Copy the data, assuming that the highest dimension are stored first.
        char *src  = static_cast<char*>(input_sptr->data);
        char *dest = static_cast<char*>(image_ptr->data);
        memcpy(dest, src+index, mem);

        // The code was vector. now set to none
        image_ptr->intent_code = NIFTI_INTENT_NONE;

        // Add to vector of single-component images
        output.push_back(make_shared<nifti_image>(*image_ptr));
    }

    return output;
}

/// Save a multicomponent nifti image
void save_multicomponent_nifti_image(std::shared_ptr<nifti_image> input_sptr, const string &filename, const bool &split_xyz)
{
    // If the user wants it saved as multicomponent image
    if (!split_xyz) {
        SIRFRegMisc::save_nifti_image(input_sptr,filename);
        return;
    }

    // Else, the user wants the multicomponent image split into 3 separate images.
    vector<std::shared_ptr<nifti_image> > components = split_multicomponent_nifti_image(input_sptr);

    // Loop over each component
    for (int i=0; i<components.size(); i++) {

        // Edit the filename
        string appended_filename;

        if (components.size() == 3) {
            if      (i == 0) appended_filename = filename + "_x";
            else if (i == 1) appended_filename = filename + "_y";
            else if (i == 2) appended_filename = filename + "_z";
        }
        else appended_filename = filename + "_" + to_string(i+1);

        // And save it
        save_nifti_image(components[i],appended_filename);
    }
}

/// Copy nifti image
void copy_nifti_image(std::shared_ptr<nifti_image> &output_image_sptr, const std::shared_ptr<nifti_image> &image_to_copy_sptr)
{
    cout << "\nPerforming hard copy of nifti image..." << flush;

    // Copy the info
    nifti_image *output_ptr;

    output_ptr = nifti_copy_nim_info(image_to_copy_sptr.get());
    output_image_sptr = make_shared<nifti_image>(*output_ptr);

    // How much memory do we need to copy?
    size_t mem = output_image_sptr->nvox * output_image_sptr->nbyper;

    // Allocate the memory
    output_image_sptr->data=static_cast<void *>(malloc(mem));

    // Copy!
    memcpy(output_image_sptr->data, image_to_copy_sptr->data, mem);

    cout << "done.\n\n";
}

/// Flip multicomponent image along a given axis
void flip_multicomponent_image(SIRFImageDataDeformation &im, int dim)
{
    cout << "\nFlipping multicomponent image in dim number: " << dim << "..." << flush;

    std::shared_ptr<nifti_image> im_sptr = im.get_image_as_nifti();

    // Check the dimension to flip, that dims==5 and nu==3
    if (dim < 0 || dim > 2)
        throw runtime_error("\n\tDimension to flip should be between 0 and 2.");
    if (im_sptr->dim[0] != 5)
        throw runtime_error("\n\tNifti image is not a multicomponent image.");
    if (im_sptr->nu != 3)
        throw runtime_error("\n\tMulticomponent aspect should contain three values (x,y,z).");

    // Data is ordered such that the multicomponent info is last.
    // So, the first third of the data is the x-values, second third is y and last third is z.
    // Data is therefore = dim_number * num_voxels/3
    int start_index =   dim   * int(im_sptr->nvox/3);
    // End index is one before the start of the next dimension (thus the minus 1)
    int end_index   = (dim+1) * int(im_sptr->nvox/3 - 1);

    // Check whether single or double precision
    if (im_sptr->datatype == DT_FLOAT32) {
        // Get data
        float *data_ptr = static_cast<float*>(im_sptr->data);
        for (int i=start_index; i<=end_index; i++) {
            data_ptr[i] = -data_ptr[i];
        }
    }
    else if (im_sptr->datatype == DT_FLOAT64) {
        // Get data
        double *data_ptr = static_cast<double*>(im_sptr->data);
        for (int i=start_index; i<=end_index; i++)
            data_ptr[i] = -data_ptr[i];
    }
    else
        throw runtime_error("\n\tOnly double and float images are supported at the moment. (This would be easy to change.)");

    cout << "done.\n\n";
}

#if NIFTYREG_VER_1_3
/// Get cpp from transformation matrix
void get_cpp_from_transformation_matrix(std::shared_ptr<nifti_image> &cpp_sptr, const mat44 &TM_sptr, const std::shared_ptr<nifti_image> &warped_sptr)
{
    // Copy info from the reference image
    nifti_image *cpp_ptr = cpp_sptr.get();
    cpp_ptr = nifti_copy_nim_info(warped_sptr.get());

    // Edit some of the information to make it a cpp image
    cpp_ptr->dim[0]      = cpp_ptr->ndim = 5;
    cpp_ptr->dim[5]      = cpp_ptr->nu   = 3;
    cpp_ptr->datatype    = 16;
    cpp_ptr->nbyper      = 4;
    cpp_ptr->nvox       *= 3;
    cpp_ptr->intent_code = NIFTI_INTENT_VECTOR;

    // Allocate memory
    cpp_ptr->data=static_cast<void *>(malloc(cpp_ptr->nvox*cpp_ptr->nbyper));

    // Convert affine transformation to cpp
    reg_bspline_initialiseControlPointGridWithAffine(TM_sptr.get(), cpp_ptr);

    // Need to correct the control point position image (otherwise nv=0 and you can't read with matlab)
    reg_checkAndCorrectDimension(cpp_ptr);

    // Copy output
    cpp_sptr = make_shared<nifti_image>(*cpp_ptr);
}
#endif
/// Get def from cpp
void get_def_from_cpp(SIRFImageDataDeformation &def, const SIRFImageDataDeformation &cpp, const SIRFImageData &ref)
{
    def.create_from_3D_image(ref);
    
    reg_spline_getDeformationField(cpp.get_image_as_nifti().get(),
#if NIFTYREG_VER_1_3
                                   ref.get_image_as_nifti().get(),
#endif
                                   def.get_image_as_nifti().get(),
                                   NULL,
                                   false, //composition
                                   true // bspline
                                   );
}

/// Convert from deformation to displacement field image
void convert_from_def_to_disp(SIRFImageDataDeformation &im)
{
    // Get the disp field from the def field
    reg_getDisplacementFromDeformation(im.get_image_as_nifti().get());
    im.get_image_as_nifti()->intent_p1 = DISP_FIELD;
}

/// Convert from displacement to deformation field image
void convert_from_disp_to_def(SIRFImageDataDeformation &im)
{
    // Get the def field from the disp field
    reg_getDeformationFromDisplacement(im.get_image_as_nifti().get());
    im.get_image_as_nifti()->intent_p1 = DEF_FIELD;
}

/// Multiply image
void multiply_image(SIRFImageData &output, const SIRFImageData &input, const float &value)
{
#if NIFTYREG_VER_1_5
    reg_tools_multiplyValueToImage(input.get_image_as_nifti().get(), output.get_image_as_nifti().get(), value);
#elif NIFTYREG_VER_1_3
    // for last argument: 0=add, 1=subtract, 2=multiply, 3=divide
    reg_tools_addSubMulDivValue(input.get_image_as_nifti().get(), output.get_image_as_nifti().get(), value, 2);
#endif
}

/// Do nifti image metadatas match?
bool do_nifti_image_metadata_match(const SIRFImageData &im1, const SIRFImageData &im2)
{
    const std::shared_ptr<nifti_image> im1_sptr = im1.get_image_as_nifti();
    const std::shared_ptr<nifti_image> im2_sptr = im2.get_image_as_nifti();

    bool images_match = true;
    if (!do_nifti_image_metadata_elements_match("analyze75_orient",im1_sptr->analyze75_orient,im2_sptr->analyze75_orient)) images_match = false;
    if (!do_nifti_image_metadata_elements_match("byteorder",       im1_sptr->byteorder,       im2_sptr->byteorder       )) images_match = false;
    if (!do_nifti_image_metadata_elements_match("cal_max",         im1_sptr->cal_max,         im2_sptr->cal_max         )) images_match = false;
    if (!do_nifti_image_metadata_elements_match("cal_min",         im1_sptr->cal_min,         im2_sptr->cal_min         )) images_match = false;
    if (!do_nifti_image_metadata_elements_match("datatype",        im1_sptr->datatype,        im2_sptr->datatype        )) images_match = false;
    if (!do_nifti_image_metadata_elements_match("du",              im1_sptr->du,              im2_sptr->du              )) images_match = false;
    if (!do_nifti_image_metadata_elements_match("dv",              im1_sptr->dv,              im2_sptr->dv              )) images_match = false;
    if (!do_nifti_image_metadata_elements_match("dw",              im1_sptr->dw,              im2_sptr->dw              )) images_match = false;
    if (!do_nifti_image_metadata_elements_match("dx",              im1_sptr->dx,              im2_sptr->dx              )) images_match = false;
    if (!do_nifti_image_metadata_elements_match("dy",              im1_sptr->dy,              im2_sptr->dy              )) images_match = false;
    if (!do_nifti_image_metadata_elements_match("dz",              im1_sptr->dz,              im2_sptr->dz              )) images_match = false;
    if (!do_nifti_image_metadata_elements_match("ext_list",        im1_sptr->ext_list,        im2_sptr->ext_list        )) images_match = false;
    if (!do_nifti_image_metadata_elements_match("freq_dim",        im1_sptr->freq_dim,        im2_sptr->freq_dim        )) images_match = false;
    if (!do_nifti_image_metadata_elements_match("iname_offset",    im1_sptr->iname_offset,    im2_sptr->iname_offset    )) images_match = false;
    if (!do_nifti_image_metadata_elements_match("intent_code",     im1_sptr->intent_code,     im2_sptr->intent_code     )) images_match = false;
    if (!do_nifti_image_metadata_elements_match("intent_p1",       im1_sptr->intent_p1,       im2_sptr->intent_p1       )) images_match = false;
    if (!do_nifti_image_metadata_elements_match("intent_p2",       im1_sptr->intent_p2,       im2_sptr->intent_p2       )) images_match = false;
    if (!do_nifti_image_metadata_elements_match("intent_p3",       im1_sptr->intent_p3,       im2_sptr->intent_p3       )) images_match = false;
    if (!do_nifti_image_metadata_elements_match("nbyper",          im1_sptr->nbyper,          im2_sptr->nbyper          )) images_match = false;
    if (!do_nifti_image_metadata_elements_match("ndim",            im1_sptr->ndim,            im2_sptr->ndim            )) images_match = false;
    if (!do_nifti_image_metadata_elements_match("nifti_type",      im1_sptr->nifti_type,      im2_sptr->nifti_type      )) images_match = false;
    if (!do_nifti_image_metadata_elements_match("nt",              im1_sptr->nt,              im2_sptr->nt              )) images_match = false;
    if (!do_nifti_image_metadata_elements_match("nu",              im1_sptr->nu,              im2_sptr->nu              )) images_match = false;
    if (!do_nifti_image_metadata_elements_match("num_ext",         im1_sptr->num_ext,         im2_sptr->num_ext         )) images_match = false;
    if (!do_nifti_image_metadata_elements_match("nv",              im1_sptr->nv,              im2_sptr->nv              )) images_match = false;
    if (!do_nifti_image_metadata_elements_match("nvox",            im1_sptr->nvox,            im2_sptr->nvox            )) images_match = false;
    if (!do_nifti_image_metadata_elements_match("nw",              im1_sptr->nw,              im2_sptr->nw              )) images_match = false;
    if (!do_nifti_image_metadata_elements_match("nx",              im1_sptr->nx,              im2_sptr->nx              )) images_match = false;
    if (!do_nifti_image_metadata_elements_match("ny",              im1_sptr->ny,              im2_sptr->ny              )) images_match = false;
    if (!do_nifti_image_metadata_elements_match("nz",              im1_sptr->nz,              im2_sptr->nz              )) images_match = false;
    if (!do_nifti_image_metadata_elements_match("phase_dim",       im1_sptr->phase_dim,       im2_sptr->phase_dim       )) images_match = false;
    if (!do_nifti_image_metadata_elements_match("qfac",            im1_sptr->qfac,            im2_sptr->qfac            )) images_match = false;
    if (!do_nifti_image_metadata_elements_match("qform_code",      im1_sptr->qform_code,      im2_sptr->qform_code      )) images_match = false;
    if (!do_nifti_image_metadata_elements_match("qoffset_x",       im1_sptr->qoffset_x,       im2_sptr->qoffset_x       )) images_match = false;
    if (!do_nifti_image_metadata_elements_match("qoffset_y",       im1_sptr->qoffset_y,       im2_sptr->qoffset_y       )) images_match = false;
    if (!do_nifti_image_metadata_elements_match("qoffset_z",       im1_sptr->qoffset_z,       im2_sptr->qoffset_z       )) images_match = false;
    if (!do_nifti_image_metadata_elements_match("quatern_b",       im1_sptr->quatern_b,       im2_sptr->quatern_b       )) images_match = false;
    if (!do_nifti_image_metadata_elements_match("quatern_c",       im1_sptr->quatern_c,       im2_sptr->quatern_c       )) images_match = false;
    if (!do_nifti_image_metadata_elements_match("quatern_d",       im1_sptr->quatern_d,       im2_sptr->quatern_d       )) images_match = false;
    if (!do_nifti_image_metadata_elements_match("scl_inter",       im1_sptr->scl_inter,       im2_sptr->scl_inter       )) images_match = false;
    if (!do_nifti_image_metadata_elements_match("scl_slope",       im1_sptr->scl_slope,       im2_sptr->scl_slope       )) images_match = false;
    if (!do_nifti_image_metadata_elements_match("sform_code",      im1_sptr->sform_code,      im2_sptr->sform_code      )) images_match = false;
    if (!do_nifti_image_metadata_elements_match("slice_code",      im1_sptr->slice_code,      im2_sptr->slice_code      )) images_match = false;
    if (!do_nifti_image_metadata_elements_match("slice_dim",       im1_sptr->slice_dim,       im2_sptr->slice_dim       )) images_match = false;
    if (!do_nifti_image_metadata_elements_match("slice_duration",  im1_sptr->slice_duration,  im2_sptr->slice_duration  )) images_match = false;
    if (!do_nifti_image_metadata_elements_match("slice_end",       im1_sptr->slice_end,       im2_sptr->slice_end       )) images_match = false;
    if (!do_nifti_image_metadata_elements_match("slice_start",     im1_sptr->slice_start,     im2_sptr->slice_start     )) images_match = false;
    if (!do_nifti_image_metadata_elements_match("swapsize",        im1_sptr->swapsize,        im2_sptr->swapsize        )) images_match = false;
    if (!do_nifti_image_metadata_elements_match("time_units",      im1_sptr->time_units,      im2_sptr->time_units      )) images_match = false;
    if (!do_nifti_image_metadata_elements_match("toffset",         im1_sptr->toffset,         im2_sptr->toffset         )) images_match = false;
    if (!do_nifti_image_metadata_elements_match("xyz_units",       im1_sptr->xyz_units,       im2_sptr->xyz_units       )) images_match = false;
    if (!do_nifti_image_metadata_elements_match("qto_ijk",         im1_sptr->qto_ijk,         im2_sptr->qto_ijk         )) images_match = false;
    if (!do_nifti_image_metadata_elements_match("qto_xyz",         im1_sptr->qto_xyz,         im2_sptr->qto_xyz         )) images_match = false;
    if (!do_nifti_image_metadata_elements_match("sto_ijk",         im1_sptr->sto_ijk,         im2_sptr->sto_ijk         )) images_match = false;
    if (!do_nifti_image_metadata_elements_match("sto_xyz",         im1_sptr->sto_xyz,         im2_sptr->sto_xyz         )) images_match = false;

    for (int i=0; i<8; i++) {
        if (!do_nifti_image_metadata_elements_match("dim["+std::to_string(i)+"]",    im1_sptr->dim[i],    im2_sptr->dim[i] ))   images_match = false;
        if (!do_nifti_image_metadata_elements_match("pixdim["+std::to_string(i)+"]", im1_sptr->pixdim[i], im2_sptr->pixdim[i])) images_match = false;
    }

    if (images_match) cout << "\tOK!\n";

    return images_match;
}

bool do_nifti_image_metadata_elements_match(const std::string &name, const mat44 &elem1, const mat44 &elem2)
{
    if(do_mat44_match(elem1, elem1))
        return true;
    std::cout << "mismatch in " << name << "\n";
    std::vector<mat44>vec;
    vec.push_back(elem1);
    vec.push_back(elem2);
    print_mat44(vec);
    std::cout << "\n";
    return false;
}

bool do_nifti_images_match(const SIRFImageData &im1, const SIRFImageData &im2, const float accuracy_percentage_of_max)
{
    if(!im1.is_initialised())
        throw runtime_error("do_nifti_images_match: Image 1 not initialised.");
    if(!im2.is_initialised())
        throw runtime_error("do_nifti_images_match: Image 2 not initialised.");

    if (!do_nifti_image_metadata_match(im1,im2)) {
        cout << "\nImage metadata does not match.\n";
        return false;
    }

    if (im1.get_image_as_nifti()->datatype == DT_BINARY)   return do_arrays_match<bool>              (im1, im2, accuracy_percentage_of_max);
    if (im1.get_image_as_nifti()->datatype == DT_INT8)     return do_arrays_match<signed char>       (im1, im2, accuracy_percentage_of_max);
    if (im1.get_image_as_nifti()->datatype == DT_INT16)    return do_arrays_match<signed short>      (im1, im2, accuracy_percentage_of_max);
    if (im1.get_image_as_nifti()->datatype == DT_INT32)    return do_arrays_match<signed int>        (im1, im2, accuracy_percentage_of_max);
    if (im1.get_image_as_nifti()->datatype == DT_FLOAT32)  return do_arrays_match<float>             (im1, im2, accuracy_percentage_of_max);
    if (im1.get_image_as_nifti()->datatype == DT_FLOAT64)  return do_arrays_match<double>            (im1, im2, accuracy_percentage_of_max);
    if (im1.get_image_as_nifti()->datatype == DT_UINT8)    return do_arrays_match<unsigned char>     (im1, im2, accuracy_percentage_of_max);
    if (im1.get_image_as_nifti()->datatype == DT_UINT16)   return do_arrays_match<unsigned short>    (im1, im2, accuracy_percentage_of_max);
    if (im1.get_image_as_nifti()->datatype == DT_UINT32)   return do_arrays_match<unsigned int>      (im1, im2, accuracy_percentage_of_max);
    if (im1.get_image_as_nifti()->datatype == DT_INT64)    return do_arrays_match<signed long long>  (im1, im2, accuracy_percentage_of_max);
    if (im1.get_image_as_nifti()->datatype == DT_UINT64)   return do_arrays_match<unsigned long long>(im1, im2, accuracy_percentage_of_max);
    if (im1.get_image_as_nifti()->datatype == DT_FLOAT128) return do_arrays_match<long double>       (im1, im2, accuracy_percentage_of_max);

    stringstream ss;
    ss << "do_nifti_images_match not implemented for your data type: ";
    ss << nifti_datatype_string(im1.get_image_as_nifti()->datatype);
    ss << " (bytes per voxel: ";
    ss << im1.get_image_as_nifti()->nbyper << ").";
    throw std::runtime_error(ss.str());
}

/// Dump info of nifti image
void dump_nifti_info(const string &im_filename)
{
    SIRFImageData image(im_filename);
    SIRFRegMisc::dump_nifti_info(image);
}

/// Dump info of nifti image
void dump_nifti_info(const SIRFImageData &im)
{
    vector<SIRFImageData> image;
    image.push_back(im);
    dump_nifti_info(image);
}

/// Dump info of multiple nifti images
void dump_nifti_info(const vector<SIRFImageData> &ims)
{
    cout << "\nPrinting info for " << ims.size() << " nifti image(s):\n";
    dump_nifti_element(ims, "analyze_75_orient", &nifti_image::analyze75_orient);
    dump_nifti_element(ims, "analyze75_orient",  &nifti_image::analyze75_orient);
    dump_nifti_element(ims, "byteorder",         &nifti_image::byteorder);
    dump_nifti_element(ims, "cal_max",           &nifti_image::cal_max);
    dump_nifti_element(ims, "cal_min",           &nifti_image::cal_min);
    dump_nifti_element(ims, "datatype",          &nifti_image::datatype);
    dump_nifti_element(ims, "dt",                &nifti_image::dt);
    dump_nifti_element(ims, "du",                &nifti_image::du);
    dump_nifti_element(ims, "dv",                &nifti_image::dv);
    dump_nifti_element(ims, "dw",                &nifti_image::dw);
    dump_nifti_element(ims, "dx",                &nifti_image::dx);
    dump_nifti_element(ims, "dy",                &nifti_image::dy);
    dump_nifti_element(ims, "dz",                &nifti_image::dz);
    dump_nifti_element(ims, "ext_list",          &nifti_image::ext_list);
    dump_nifti_element(ims, "freq_dim",          &nifti_image::freq_dim);
    dump_nifti_element(ims, "iname_offset",      &nifti_image::iname_offset);
    dump_nifti_element(ims, "intent_code",       &nifti_image::intent_code);
    dump_nifti_element(ims, "intent_p1",         &nifti_image::intent_p1);
    dump_nifti_element(ims, "intent_p2",         &nifti_image::intent_p2);
    dump_nifti_element(ims, "intent_p3",         &nifti_image::intent_p3);
    dump_nifti_element(ims, "nbyper",            &nifti_image::nbyper);
    dump_nifti_element(ims, "ndim",              &nifti_image::ndim);
    dump_nifti_element(ims, "nifti_type",        &nifti_image::nifti_type);
    dump_nifti_element(ims, "nt",                &nifti_image::nt);
    dump_nifti_element(ims, "nu",                &nifti_image::nu);
    dump_nifti_element(ims, "num_ext",           &nifti_image::num_ext);
    dump_nifti_element(ims, "nv",                &nifti_image::nv);
    dump_nifti_element(ims, "nvox",              &nifti_image::nvox);
    dump_nifti_element(ims, "nw",                &nifti_image::nw);
    dump_nifti_element(ims, "nx",                &nifti_image::nx);
    dump_nifti_element(ims, "ny",                &nifti_image::ny);
    dump_nifti_element(ims, "nz",                &nifti_image::nz);
    dump_nifti_element(ims, "phase_dim",         &nifti_image::phase_dim);
    dump_nifti_element(ims, "qfac",              &nifti_image::qfac);
    dump_nifti_element(ims, "qform_code",        &nifti_image::qform_code);
    dump_nifti_element(ims, "qoffset_x",         &nifti_image::qoffset_x);
    dump_nifti_element(ims, "qoffset_y",         &nifti_image::qoffset_y);
    dump_nifti_element(ims, "qoffset_z",         &nifti_image::qoffset_z);
    dump_nifti_element(ims, "quatern_b",         &nifti_image::quatern_b);
    dump_nifti_element(ims, "quatern_c",         &nifti_image::quatern_c);
    dump_nifti_element(ims, "quatern_d",         &nifti_image::quatern_d);
    dump_nifti_element(ims, "scl_inter",         &nifti_image::scl_inter);
    dump_nifti_element(ims, "scl_slope",         &nifti_image::scl_slope);
    dump_nifti_element(ims, "sform_code",        &nifti_image::sform_code);
    dump_nifti_element(ims, "slice_code",        &nifti_image::slice_code);
    dump_nifti_element(ims, "slice_dim",         &nifti_image::slice_dim);
    dump_nifti_element(ims, "slice_duration",    &nifti_image::slice_duration);
    dump_nifti_element(ims, "slice_end",         &nifti_image::slice_end);
    dump_nifti_element(ims, "slice_start",       &nifti_image::slice_start);
    dump_nifti_element(ims, "swapsize",          &nifti_image::swapsize);
    dump_nifti_element(ims, "time_units",        &nifti_image::time_units);
    dump_nifti_element(ims, "toffset",           &nifti_image::toffset);
    dump_nifti_element(ims, "xyz_units",         &nifti_image::xyz_units);
    dump_nifti_element(ims, "dim",               &nifti_image::dim,    8);
    dump_nifti_element(ims, "pixdim",            &nifti_image::pixdim, 8);

    vector<std::shared_ptr<nifti_image> > images;
    for(int i=0;i<ims.size();i++)
        images.push_back(ims[i].get_image_as_nifti());

    // Print transformation matrices
    vector<mat44> qto_ijk_vec, qto_xyz_vec, sto_ijk_vec, sto_xyz_vec;
    for(int j=0; j<int(images.size()); j++) {
        qto_ijk_vec.push_back(images[j]->qto_ijk);
        qto_xyz_vec.push_back(images[j]->qto_xyz);
        sto_ijk_vec.push_back(images[j]->sto_ijk);
        sto_xyz_vec.push_back(images[j]->sto_xyz);
    }
    cout << "\t" << left << setw(19) << "qto_ijk:" << "\n";
    SIRFRegMisc::print_mat44(qto_ijk_vec);
    cout << "\t" << left << setw(19) << "qto_xyz:" << "\n";
    SIRFRegMisc::print_mat44(qto_xyz_vec);
    cout << "\t" << left << setw(19) << "sto_ijk:" << "\n";
    SIRFRegMisc::print_mat44(sto_ijk_vec);
    cout << "\t" << left << setw(19) << "sto_xyz:" << "\n";
    SIRFRegMisc::print_mat44(sto_xyz_vec);

    cout << "\n";
}

/// Save transformation matrix to file
void save_transformation_matrix(const mat44 &transformation_matrix, const string &filename)
{
    // Check that input isn't blank
    if (filename == "")
        throw runtime_error("Error, cannot write transformation matrix to file because filename is blank");

    // Have to create copy as the following function isn't const
    mat44 temp = transformation_matrix;
    reg_tool_WriteAffineFile(&temp, filename.c_str());
}

/// Read transformation matrix from file
void open_transformation_matrix(mat44 &transformation_matrix, const string& filename)
{
    // Check that the file exists
    if (!boost::filesystem::exists(filename))
        throw runtime_error("Cannot find the file: " + filename + ".");

    cout << "\n\nReading transformation matrix from file...\n\n";
    reg_tool_ReadAffineFile(&transformation_matrix, (char*)filename.c_str());
    cout << "\n\nSuccessfully read transformation matrix from file:\n";

    print_mat44(transformation_matrix);
}

/// Print mat44
void print_mat44(const mat44 &mat)
{
    std::vector<mat44> mat_vec;
    mat_vec.push_back(mat);
    print_mat44(mat_vec);
}

/// Print mat44
void print_mat44(const vector<mat44> &mats)
{
    for(int i=0;i<4;i++) {
        cout << "\t\t\t   ";
        for(int j=0;j<mats.size();j++) {
            ostringstream ss;
            ss << "[" <<
                  setprecision(3) << mats[j].m[i][0] << "," <<
                  setprecision(3) << mats[j].m[i][1] << "," <<
                  setprecision(3) << mats[j].m[i][2] << "," <<
                  setprecision(3) << mats[j].m[i][3] << "] ";
            cout << setw(19) << ss.str();
        }
        cout << "\n";
    }
}

bool do_mat44_match( const mat44& mat1, const mat44& mat2 )
{
    for (int i=0; i<4; i++)
        for (int j=0; j<4; j++)
            if( fabs(mat1.m[j][i] - mat2.m[j][i]) > 1.e-7 )
                return false;
    return true;
}

/// Mat44 multiplier
mat44 multiply_mat44(const mat44 &x, const mat44 &y)
{
    // Print info
    cout << "\nMultiplying two matrices...\n";
    cout << "Matrix 1:\n";
    print_mat44(x);
    cout << "Matrix 2:\n";
    print_mat44(y);

    // Create result, set to zero
    mat44 res;
    for (int i=0;i<4;i++) {
        for (int j=0;j<4;j++) {
            res.m[i][j] = 0.;
        }
    }

    for (int i=0;i<4;i++) {
        for (int j=0;j<4;j++) {
            for (int k=0;k<4;k++) {
                res.m[i][j] += x.m[i][k] * y.m[k][j];
            }
        }
    }

    cout << "Result:\n";
    print_mat44(res);

    return res;
}

/// Compose multiple transformations into single deformation field
void compose_transformations_into_single_deformation(SIRFRegTransformationDeformation &def, const std::vector<SIRFRegTransformation*> &transformations, const SIRFImageData &ref)
{
    if (transformations.size() == 0)
        throw std::runtime_error("SIRFRegMisc::compose_transformations_into_single_deformation no transformations given.");

    SIRFImageDataDeformation deformation = transformations.at(0)->get_as_deformation_field(ref).deep_copy();

    for (unsigned i=1; i<transformations.size(); ++i) {

        SIRFImageDataDeformation temp = transformations.at(i)->get_as_deformation_field(ref);
        reg_defField_compose(temp.get_image_as_nifti().get(),deformation.get_image_as_nifti().get(),nullptr);
    }

    def = SIRFRegTransformationDeformation(deformation);
}

/// Compose multiple transformations into single deformation field
void compose_transformations_into_single_deformation(SIRFRegTransformationDeformation &def, const std::vector<std::shared_ptr<SIRFRegTransformation> > &transformations, const SIRFImageData &ref)
{
    std::vector<SIRFRegTransformation*> vec;
    for (unsigned i=0; i<transformations.size(); ++i)
        vec.push_back(transformations.at(i).get());
    compose_transformations_into_single_deformation(def, vec, ref);
}

/// Get identity matrix
mat44 get_identity_matrix()
{
    mat44 mat;
    for (int i=0; i<4; ++i) {
        for (int j=0; j<4; ++j) {
            if (i==j)
                mat.m[i][j] = 1.F;
            else
                mat.m[i][j] = 0.F;
        }
    }
    return mat;
}

}
