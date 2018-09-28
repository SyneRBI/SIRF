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
#include "NiftiImage3D.h"
#include "NiftiImage3DTensor.h"
#include "SIRFRegTransformation.h"
#include "NiftiImage3DDeformation.h"
#include "NiftiImage3DDisplacement.h"
#include "SIRFRegMat44.h"
#include <_reg_tools.h>
#if NIFTYREG_VER_1_5
#include <_reg_globalTrans.h>
#include <_reg_localTrans.h>
#elif NIFTYREG_VER_1_3
#include <_reg_localTransformation.h>
#endif
#include <iostream>
#include <iomanip>
#include <sstream>

using namespace std;
using namespace sirf;

namespace SIRFRegMisc {

/// Open nifti image
void open_nifti_image(shared_ptr<nifti_image> &image, const boost::filesystem::path &filename)
{
    // If no filename has been set, return
    if (filename.size() == 0) {
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
    nifti_image *im = nifti_image_read(filename.c_str(), 1);
    image = shared_ptr<nifti_image>(im, nifti_image_free);

    // Ensure the image has all the values correctly set
    reg_checkAndCorrectDimension(image.get());
}

/// Save nifti image
void save_nifti_image(const NiftiImage &image, const string &filename)
{
    if (!image.is_initialised())
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

    nifti_set_filenames(image.get_raw_nifti_sptr().get(), filename.c_str(), 0, 0);
    nifti_image_write(image.get_raw_nifti_sptr().get());
    cout << "done.\n\n";
}

/// Split multi-component image
vector<NiftiImage3D>
    split_multicomponent_nifti_image(const NiftiImage3DTensor &input)
{
    // Only works for ndim==5
    if (input.get_raw_nifti_sptr()->ndim != 5 || input.get_raw_nifti_sptr()->nu != 3)
        throw runtime_error("Splitting only currently works for ndim==5, nu==3.");

    // Create the vector to store the single component images
    vector<NiftiImage3D> output;

    // Loop over all of the components
    for (int component=0; component<3; component++) {

        // Create new image
        nifti_image *image_ptr;

        // Copy the input image
        image_ptr = nifti_copy_nim_info(input.get_raw_nifti_sptr().get());

        // Alter the info to change the number of dimensions
        image_ptr->dim[0] = image_ptr->ndim = 3;
        image_ptr->dim[5] = image_ptr->nu   = 1;
        image_ptr->nvox   = image_ptr->nvox / 3;

        // How much memory do we need?
        size_t mem = image_ptr->nvox*image_ptr->nbyper;

        // Allocate the data
        image_ptr->data=static_cast<void *>(malloc(mem));

        // Start index
        size_t index = mem*component;

        // Copy the data, assuming that the highest dimension are stored first.
        char *src  = static_cast<char*>(input.get_raw_nifti_sptr()->data);
        char *dest = static_cast<char*>(image_ptr->data);
        memcpy(dest, src+index, mem);

        // The code was vector. now set to none
        image_ptr->intent_code = NIFTI_INTENT_NONE;

        // Create NiftiImage3D from this
        NiftiImage3D im(*image_ptr);

        // Add to vector of single-component images
        output.push_back(im);
    }

    return output;
}

/// Save a multicomponent nifti image
void save_multicomponent_nifti_image_split_xyz(const NiftiImage3DTensor &input, const string &filename)
{
    vector<string> filenames;

    // Loop over each component
    for (unsigned i=0; i<3; ++i) {
        if      (i == 0) filenames.push_back(filename + "_x");
        else if (i == 1) filenames.push_back(filename + "_y");
        else if (i == 2) filenames.push_back(filename + "_z");
    }
    save_multicomponent_nifti_image_split_xyz(input,filenames[0],filenames[1],filenames[2]);
}

/// Save a multicomponent nifti image
void save_multicomponent_nifti_image_split_xyz(const NiftiImage3DTensor &input, const string &filename_x, const string &filename_y, const string &filename_z)
{
    // Split into 3 separate images.
    vector<NiftiImage3D> components = split_multicomponent_nifti_image(input);
    save_nifti_image(components.at(0),filename_x);
    save_nifti_image(components.at(1),filename_y);
    save_nifti_image(components.at(2),filename_z);
}

/// Copy nifti image
void copy_nifti_image(shared_ptr<nifti_image> &output_image_sptr, const shared_ptr<nifti_image> &image_to_copy_sptr)
{
    cout << "\nPerforming hard copy of nifti image..." << flush;

    // Copy the info
    nifti_image *output_ptr;

    output_ptr = nifti_copy_nim_info(image_to_copy_sptr.get());
    output_image_sptr = shared_ptr<nifti_image>(output_ptr, nifti_image_free);

    // How much memory do we need to copy?
    size_t mem = output_image_sptr->nvox * unsigned(output_image_sptr->nbyper);

    // Allocate the memory
    output_image_sptr->data=static_cast<void *>(malloc(mem));

    // Copy!
    memcpy(output_image_sptr->data, image_to_copy_sptr->data, mem);

    // Check everything is ok
    reg_checkAndCorrectDimension(output_image_sptr.get());

    cout << "done.\n\n";
}

/// Flip multicomponent image along a given axis
void flip_multicomponent_image(NiftiImage3DTensor &im, int dim)
{
    cout << "\nFlipping multicomponent image in dim number: " << dim << "..." << flush;

    shared_ptr<nifti_image> im_sptr = im.get_raw_nifti_sptr();

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
void get_cpp_from_transformation_matrix(shared_ptr<nifti_image> &cpp_sptr, const mat44 &TM_sptr, const shared_ptr<nifti_image> &warped_sptr)
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
    cpp_sptr = shared_ptr<nifti_image>(cpp_ptr, nifti_image_free);
}
#endif
/// Get def from cpp
void get_def_from_cpp(NiftiImage3DDeformation &def, const NiftiImage3DTensor &cpp, const NiftiImage3D &ref)
{
    def.create_from_3D_image(ref);
    
    reg_spline_getDeformationField(cpp.get_raw_nifti_sptr().get(),
#if NIFTYREG_VER_1_3
                                   ref.get_raw_nifti_sptr().get(),
#endif
                                   def.get_raw_nifti_sptr().get(),
                                   NULL,
                                   false, //composition
                                   true // bspline
                                   );
}

/// Convert from deformation to displacement field image
void convert_from_def_to_disp(NiftiImage3DDisplacement &disp, const NiftiImage3DDeformation &def)
{
    // Get the disp field from the def field
    NiftiImage3DTensor temp = def.deep_copy();
    reg_getDisplacementFromDeformation(temp.get_raw_nifti_sptr().get());
    temp.get_raw_nifti_sptr()->intent_p1 = DISP_FIELD;
    disp = temp;
}

/// Convert from displacement to deformation field image
void convert_from_disp_to_def(NiftiImage3DDeformation &def, const NiftiImage3DDisplacement &disp)
{
    // Get the def field from the disp field
    NiftiImage3DTensor temp = disp.deep_copy();
    reg_getDeformationFromDisplacement(temp.get_raw_nifti_sptr().get());
    temp.get_raw_nifti_sptr()->intent_p1 = DEF_FIELD;
    def = temp.deep_copy();
}

/// Do nifti image metadatas match?
bool do_nifti_image_metadata_match(const NiftiImage &im1, const NiftiImage &im2)
{
    const shared_ptr<nifti_image> im1_sptr = im1.get_raw_nifti_sptr();
    const shared_ptr<nifti_image> im2_sptr = im2.get_raw_nifti_sptr();

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
        if (!do_nifti_image_metadata_elements_match("dim["+to_string(i)+"]",    im1_sptr->dim[i],    im2_sptr->dim[i] ))   images_match = false;
        if (!do_nifti_image_metadata_elements_match("pixdim["+to_string(i)+"]", im1_sptr->pixdim[i], im2_sptr->pixdim[i])) images_match = false;
    }

    if (images_match) cout << "\tOK!\n";

    return images_match;
}

template<typename T>
bool do_nifti_image_metadata_elements_match(const string &name, const T &elem1, const T &elem2)
{
    if(float(fabs(elem1-elem2)) < 1.e-7F)
        return true;
    cout << "mismatch in " << name << " , (values: " <<  elem1 << " and " << elem2 << ")\n";
    return false;
}
template bool do_nifti_image_metadata_elements_match<float> (const string &name, const float &elem1, const float &elem2);

bool do_nifti_image_metadata_elements_match(const string &name, const mat44 &elem1, const mat44 &elem2)
{
    if(SIRFRegMat44(elem1)==SIRFRegMat44(elem2))
        return true;
    cout << "mismatch in " << name << "\n";
    SIRFRegMat44::print({elem1, elem2});
    cout << "\n";
    return false;
}

bool do_nifti_images_match(const NiftiImage &im1, const NiftiImage &im2, const float accuracy_percentage_of_max)
{
    if(!im1.is_initialised())
        throw runtime_error("do_nifti_images_match: Image 1 not initialised.");
    if(!im2.is_initialised())
        throw runtime_error("do_nifti_images_match: Image 2 not initialised.");

    if (!do_nifti_image_metadata_match(im1,im2)) {
        cout << "\nImage metadata does not match.\n";
        return false;
    }

    if (im1.get_raw_nifti_sptr()->datatype == DT_BINARY)   return do_arrays_match<bool>              (im1, im2, accuracy_percentage_of_max);
    if (im1.get_raw_nifti_sptr()->datatype == DT_INT8)     return do_arrays_match<signed char>       (im1, im2, accuracy_percentage_of_max);
    if (im1.get_raw_nifti_sptr()->datatype == DT_INT16)    return do_arrays_match<signed short>      (im1, im2, accuracy_percentage_of_max);
    if (im1.get_raw_nifti_sptr()->datatype == DT_INT32)    return do_arrays_match<signed int>        (im1, im2, accuracy_percentage_of_max);
    if (im1.get_raw_nifti_sptr()->datatype == DT_FLOAT32)  return do_arrays_match<float>             (im1, im2, accuracy_percentage_of_max);
    if (im1.get_raw_nifti_sptr()->datatype == DT_FLOAT64)  return do_arrays_match<double>            (im1, im2, accuracy_percentage_of_max);
    if (im1.get_raw_nifti_sptr()->datatype == DT_UINT8)    return do_arrays_match<unsigned char>     (im1, im2, accuracy_percentage_of_max);
    if (im1.get_raw_nifti_sptr()->datatype == DT_UINT16)   return do_arrays_match<unsigned short>    (im1, im2, accuracy_percentage_of_max);
    if (im1.get_raw_nifti_sptr()->datatype == DT_UINT32)   return do_arrays_match<unsigned int>      (im1, im2, accuracy_percentage_of_max);
    if (im1.get_raw_nifti_sptr()->datatype == DT_INT64)    return do_arrays_match<signed long long>  (im1, im2, accuracy_percentage_of_max);
    if (im1.get_raw_nifti_sptr()->datatype == DT_UINT64)   return do_arrays_match<unsigned long long>(im1, im2, accuracy_percentage_of_max);
    if (im1.get_raw_nifti_sptr()->datatype == DT_FLOAT128) return do_arrays_match<long double>       (im1, im2, accuracy_percentage_of_max);

    stringstream ss;
    ss << "do_nifti_images_match not implemented for your data type: ";
    ss << nifti_datatype_string(im1.get_raw_nifti_sptr()->datatype);
    ss << " (bytes per voxel: ";
    ss << im1.get_raw_nifti_sptr()->nbyper << ").";
    throw runtime_error(ss.str());
}

template<typename T>
bool do_arrays_match(const NiftiImage &im1, const NiftiImage &im2, const float &required_accuracy_compared_to_max)
{
    if(!im1.is_initialised())
        throw runtime_error("do_arrays_match: Image 1 not initialised.");
    if(!im2.is_initialised())
        throw runtime_error("do_arrays_match: Image 2 not initialised.");

    // Get norm between two images
    float norm = im1.get_norm(im2);
    float epsilon = (im1.get_max()+im2.get_max())/2.F;
    epsilon *= required_accuracy_compared_to_max;

    if (norm > epsilon) {
        cout << "\nImages are not equal (norm > epsilon).\n";
        cout << "\tmax1                              = " << im1.get_max() << "\n";
        cout << "\tmax2                              = " << im1.get_max() << "\n";
        cout << "\tmin1                              = " << im1.get_min() << "\n";
        cout << "\tmin2                              = " << im2.get_min() << "\n";
        cout << "\trequired accuracy compared to max = " << required_accuracy_compared_to_max << "\n";
        cout << "\tepsilon                           = " << epsilon << "\n";
        cout << "\tnorm                              = " << norm << "\n";
        return false;
    }
    return true;
}
template bool do_arrays_match<bool>              (const NiftiImage &im1, const NiftiImage &im2, const float &required_accuracy_compared_to_max);
template bool do_arrays_match<signed char>       (const NiftiImage &im1, const NiftiImage &im2, const float &required_accuracy_compared_to_max);
template bool do_arrays_match<signed short>      (const NiftiImage &im1, const NiftiImage &im2, const float &required_accuracy_compared_to_max);
template bool do_arrays_match<signed int>        (const NiftiImage &im1, const NiftiImage &im2, const float &required_accuracy_compared_to_max);
template bool do_arrays_match<float>             (const NiftiImage &im1, const NiftiImage &im2, const float &required_accuracy_compared_to_max);
template bool do_arrays_match<double>            (const NiftiImage &im1, const NiftiImage &im2, const float &required_accuracy_compared_to_max);
template bool do_arrays_match<unsigned char>     (const NiftiImage &im1, const NiftiImage &im2, const float &required_accuracy_compared_to_max);
template bool do_arrays_match<unsigned short>    (const NiftiImage &im1, const NiftiImage &im2, const float &required_accuracy_compared_to_max);
template bool do_arrays_match<unsigned int>      (const NiftiImage &im1, const NiftiImage &im2, const float &required_accuracy_compared_to_max);
template bool do_arrays_match<signed long long>  (const NiftiImage &im1, const NiftiImage &im2, const float &required_accuracy_compared_to_max);
template bool do_arrays_match<unsigned long long>(const NiftiImage &im1, const NiftiImage &im2, const float &required_accuracy_compared_to_max);
template bool do_arrays_match<long double>       (const NiftiImage &im1, const NiftiImage &im2, const float &required_accuracy_compared_to_max);

template<typename T>
float get_array_max(const NiftiImage &im)
{
    if(!im.is_initialised())
        throw runtime_error("get_array_max: Image not initialised.");

    // Check sizes
    if (im.get_raw_nifti_sptr()->nbyper != sizeof(T))
        throw runtime_error("get_array_max: Datatype does not match desired cast type (" + to_string(im.get_raw_nifti_sptr()->nbyper) + " versus " + to_string(sizeof(T)) + ").");

    // Get data
    T *data = static_cast<T*>(im.get_raw_nifti_sptr()->data);
    return *max_element(data, data + im.get_raw_nifti_sptr()->nvox);
}
template float get_array_max<bool>              (const NiftiImage &im);
template float get_array_max<signed char>       (const NiftiImage &im);
template float get_array_max<signed short>      (const NiftiImage &im);
template float get_array_max<signed int>        (const NiftiImage &im);
template float get_array_max<float>             (const NiftiImage &im);
template float get_array_max<double>            (const NiftiImage &im);
template float get_array_max<unsigned char>     (const NiftiImage &im);
template float get_array_max<unsigned short>    (const NiftiImage &im);
template float get_array_max<unsigned int>      (const NiftiImage &im);
template float get_array_max<signed long long>  (const NiftiImage &im);
template float get_array_max<unsigned long long>(const NiftiImage &im);
template float get_array_max<long double>       (const NiftiImage &im);

template<typename T>
float get_array_min(const NiftiImage &im)
{
    if(!im.is_initialised())
        throw runtime_error("get_array_min: Image not initialised.");

    // Check sizes
    if (im.get_raw_nifti_sptr()->nbyper != sizeof(T))
        throw runtime_error("get_array_min: Datatype does not match desired cast type (" + to_string(im.get_raw_nifti_sptr()->nbyper) + " versus " + to_string(sizeof(T)) + ").");

    // Get data
    T *data = static_cast<T*>(im.get_raw_nifti_sptr()->data);
    return *min_element(data, data + im.get_raw_nifti_sptr()->nvox);
}
template float get_array_min<bool>              (const NiftiImage &im);
template float get_array_min<signed char>       (const NiftiImage &im);
template float get_array_min<signed short>      (const NiftiImage &im);
template float get_array_min<signed int>        (const NiftiImage &im);
template float get_array_min<float>             (const NiftiImage &im);
template float get_array_min<double>            (const NiftiImage &im);
template float get_array_min<unsigned char>     (const NiftiImage &im);
template float get_array_min<unsigned short>    (const NiftiImage &im);
template float get_array_min<unsigned int>      (const NiftiImage &im);
template float get_array_min<signed long long>  (const NiftiImage &im);
template float get_array_min<unsigned long long>(const NiftiImage &im);
template float get_array_min<long double>       (const NiftiImage &im);

template<typename T>
float get_array_sum(const NiftiImage &im)
{
    if(!im.is_initialised())
        throw runtime_error("get_array_min: Image not initialised.");

    // Check sizes
    if (im.get_raw_nifti_sptr()->nbyper != sizeof(T))
        throw runtime_error("get_array_min: Datatype does not match desired cast type (" + to_string(im.get_raw_nifti_sptr()->nbyper) + " versus " + to_string(sizeof(T)) + ").");

    // Get data
    T *data = static_cast<T*>(im.get_raw_nifti_sptr()->data);
    float sum = 0.F;
    for (unsigned i=0; i<im.get_raw_nifti_sptr()->nvox; ++i)
        sum += float(data[i]);
    return sum;
}
template float get_array_sum<bool>              (const NiftiImage &im);
template float get_array_sum<signed char>       (const NiftiImage &im);
template float get_array_sum<signed short>      (const NiftiImage &im);
template float get_array_sum<signed int>        (const NiftiImage &im);
template float get_array_sum<float>             (const NiftiImage &im);
template float get_array_sum<double>            (const NiftiImage &im);
template float get_array_sum<unsigned char>     (const NiftiImage &im);
template float get_array_sum<unsigned short>    (const NiftiImage &im);
template float get_array_sum<unsigned int>      (const NiftiImage &im);
template float get_array_sum<signed long long>  (const NiftiImage &im);
template float get_array_sum<unsigned long long>(const NiftiImage &im);
template float get_array_sum<long double>       (const NiftiImage &im);

template<typename T>
float get_array_element(const NiftiImage &im, const int idx[])
{
    if(!im.is_initialised())
        throw runtime_error("get_3D_array_element: Image not initialised.");

    // Check sizes
    if (im.get_raw_nifti_sptr()->nbyper != sizeof(T))
        throw runtime_error("get_3D_array_element: Datatype does not match desired cast type (" + to_string(im.get_raw_nifti_sptr()->nbyper) + " versus " + to_string(sizeof(T)) + ").");

    // Get dims and spacing
    int   *dim    = im.get_raw_nifti_sptr()->dim;

    // Check it's in bounds
    for (int i=0; i<dim[0]; ++i) {
        if (idx[i]<0 || idx[i]>=dim[i+1]) {
            stringstream ss;
            ss << "get_array_element: Element out of bounds.\n";
            ss << "\tRequested = ( "; for (int i=0;i<dim[0];++i) ss << idx[i] << " ";
            ss << ")\n\tBounds = ( "; for (int i=0;i<dim[0];++i) ss << dim[i+1] << " "; ss << " ).";
            throw runtime_error(ss.str());
        }
    }

    int idx_1d = 0;
    for (int i=0; i<dim[0]; ++i) {
        int component = idx[i];
        for (int j=i+2; j<dim[0]+1; ++j)
            component *= dim[j];
        idx_1d += component;
    }

    // Get data
    T *data = static_cast<T*>(im.get_raw_nifti_sptr()->data);

    cout << "\nget_array_element: Be careful, I made this quickly for debugging and haven't thought about data order.\n";

    return float(data[idx_1d]);
}
template float get_array_element<bool>              (const NiftiImage &im, const int idx[]);
template float get_array_element<signed char>       (const NiftiImage &im, const int idx[]);
template float get_array_element<signed short>      (const NiftiImage &im, const int idx[]);
template float get_array_element<signed int>        (const NiftiImage &im, const int idx[]);
template float get_array_element<float>             (const NiftiImage &im, const int idx[]);
template float get_array_element<double>            (const NiftiImage &im, const int idx[]);
template float get_array_element<unsigned char>     (const NiftiImage &im, const int idx[]);
template float get_array_element<unsigned short>    (const NiftiImage &im, const int idx[]);
template float get_array_element<unsigned int>      (const NiftiImage &im, const int idx[]);
template float get_array_element<signed long long>  (const NiftiImage &im, const int idx[]);
template float get_array_element<unsigned long long>(const NiftiImage &im, const int idx[]);
template float get_array_element<long double>       (const NiftiImage &im, const int idx[]);

template<typename T>
NiftiImage sum_arrays(const NiftiImage &im1, const NiftiImage &im2)
{
    if(!im1.is_initialised())
        throw runtime_error("sum_arrays: Image 1 not initialised.");
    if(!im2.is_initialised())
        throw runtime_error("sum_arrays: Image 2 not initialised.");
    if(!do_nifti_image_metadata_match(im1,im2))
        throw runtime_error("sum_arrays: Cannot add images as metadata does not match.");
    // Check sizes
    if (im1.get_raw_nifti_sptr()->nbyper != sizeof(T))
        throw runtime_error("sum_arrays: Datatype of image 1 does not match desired cast type (" + to_string(im1.get_raw_nifti_sptr()->nbyper) + " versus " + to_string(sizeof(T)) + ").");

    NiftiImage result = im1.deep_copy();

    // Get data
    T *im2_data = static_cast<T*>(im2.get_raw_nifti_sptr()->data);
    T *res_data = static_cast<T*>(result.get_raw_nifti_sptr()->data);

    for (unsigned i=0; i<im1.get_raw_nifti_sptr()->nvox; ++i) {
        res_data[i] += im2_data[i];
    }
    return result;
}
template NiftiImage sum_arrays<bool>              (const NiftiImage &im1, const NiftiImage &im2);
template NiftiImage sum_arrays<signed char>       (const NiftiImage &im1, const NiftiImage &im2);
template NiftiImage sum_arrays<signed short>      (const NiftiImage &im1, const NiftiImage &im2);
template NiftiImage sum_arrays<signed int>        (const NiftiImage &im1, const NiftiImage &im2);
template NiftiImage sum_arrays<float>             (const NiftiImage &im1, const NiftiImage &im2);
template NiftiImage sum_arrays<double>            (const NiftiImage &im1, const NiftiImage &im2);
template NiftiImage sum_arrays<unsigned char>     (const NiftiImage &im1, const NiftiImage &im2);
template NiftiImage sum_arrays<unsigned short>    (const NiftiImage &im1, const NiftiImage &im2);
template NiftiImage sum_arrays<unsigned int>      (const NiftiImage &im1, const NiftiImage &im2);
template NiftiImage sum_arrays<signed long long>  (const NiftiImage &im1, const NiftiImage &im2);
template NiftiImage sum_arrays<unsigned long long>(const NiftiImage &im1, const NiftiImage &im2);
template NiftiImage sum_arrays<long double>       (const NiftiImage &im1, const NiftiImage &im2);

template<typename T>
NiftiImage sub_arrays(const NiftiImage &im1, const NiftiImage &im2)
{
    if(!im1.is_initialised())
        throw runtime_error("sub_arrays: Image 1 not initialised.");
    if(!im2.is_initialised())
        throw runtime_error("sub_arrays: Image 2 not initialised.");
    if(!do_nifti_image_metadata_match(im1,im2))
        throw runtime_error("sub_arrays: Cannot subtract images as metadata does not match.");
    // Check sizes
    if (im1.get_raw_nifti_sptr()->nbyper != sizeof(T))
        throw runtime_error("sub_arrays: Datatype of image 1 does not match desired cast type (" + to_string(im1.get_raw_nifti_sptr()->nbyper) + " versus " + to_string(sizeof(T)) + ").");

    NiftiImage result = im1.deep_copy();

    // Get data
    T *im2_data = static_cast<T*>(im2.get_raw_nifti_sptr()->data);
    T *res_data = static_cast<T*>(result.get_raw_nifti_sptr()->data);

    for (unsigned i=0; i<im1.get_raw_nifti_sptr()->nvox; ++i) {
        res_data[i] -= im2_data[i];
    }
    return result;
}
template NiftiImage sub_arrays<bool>              (const NiftiImage &im1, const NiftiImage &im2);
template NiftiImage sub_arrays<signed char>       (const NiftiImage &im1, const NiftiImage &im2);
template NiftiImage sub_arrays<signed short>      (const NiftiImage &im1, const NiftiImage &im2);
template NiftiImage sub_arrays<signed int>        (const NiftiImage &im1, const NiftiImage &im2);
template NiftiImage sub_arrays<float>             (const NiftiImage &im1, const NiftiImage &im2);
template NiftiImage sub_arrays<double>            (const NiftiImage &im1, const NiftiImage &im2);
template NiftiImage sub_arrays<unsigned char>     (const NiftiImage &im1, const NiftiImage &im2);
template NiftiImage sub_arrays<unsigned short>    (const NiftiImage &im1, const NiftiImage &im2);
template NiftiImage sub_arrays<unsigned int>      (const NiftiImage &im1, const NiftiImage &im2);
template NiftiImage sub_arrays<signed long long>  (const NiftiImage &im1, const NiftiImage &im2);
template NiftiImage sub_arrays<unsigned long long>(const NiftiImage &im1, const NiftiImage &im2);
template NiftiImage sub_arrays<long double>       (const NiftiImage &im1, const NiftiImage &im2);

template<typename T>
float arrays_norm(const NiftiImage &im1, const NiftiImage &im2)
{
    if(!im1.is_initialised())
        throw runtime_error("sum_arrays: Image 1 not initialised.");
    if(!im2.is_initialised())
        throw runtime_error("sum_arrays: Image 2 not initialised.");
    for (int i=0; i<8; ++i)
        if (im1.get_raw_nifti_sptr()->dim[i] != im2.get_raw_nifti_sptr()->dim[i])
            throw runtime_error("SIRFRegMisc::arrays_norm: dimensions do not match.");
    // Check sizes
    if (im1.get_raw_nifti_sptr()->nbyper != sizeof(T))
        throw runtime_error("arrays_norm: Datatype of image 1 does not match desired cast type (" + to_string(im1.get_raw_nifti_sptr()->nbyper) + " versus " + to_string(sizeof(T)) + ").");
    if (im2.get_raw_nifti_sptr()->nbyper != sizeof(T))
        throw runtime_error("arrays_norm: Datatype of image 1 does not match desired cast type (" + to_string(im2.get_raw_nifti_sptr()->nbyper) + " versus " + to_string(sizeof(T)) + ").");

    // Get data
    T *im1_data = static_cast<T*>(im1.get_raw_nifti_sptr()->data);
    T *im2_data = static_cast<T*>(im2.get_raw_nifti_sptr()->data);

    // Use double precision to minimise rounding errors
    double result(0);
    for (unsigned i=0; i<im1.get_raw_nifti_sptr()->nvox; ++i)
        result += double(pow(im1_data[i] - im2_data[i], 2));
    return float(sqrt(result));
}
template float arrays_norm<bool>              (const NiftiImage &im1, const NiftiImage &im2);
template float arrays_norm<signed char>       (const NiftiImage &im1, const NiftiImage &im2);
template float arrays_norm<signed short>      (const NiftiImage &im1, const NiftiImage &im2);
template float arrays_norm<signed int>        (const NiftiImage &im1, const NiftiImage &im2);
template float arrays_norm<float>             (const NiftiImage &im1, const NiftiImage &im2);
template float arrays_norm<double>            (const NiftiImage &im1, const NiftiImage &im2);
template float arrays_norm<unsigned char>     (const NiftiImage &im1, const NiftiImage &im2);
template float arrays_norm<unsigned short>    (const NiftiImage &im1, const NiftiImage &im2);
template float arrays_norm<unsigned int>      (const NiftiImage &im1, const NiftiImage &im2);
template float arrays_norm<signed long long>  (const NiftiImage &im1, const NiftiImage &im2);
template float arrays_norm<unsigned long long>(const NiftiImage &im1, const NiftiImage &im2);
template float arrays_norm<long double>       (const NiftiImage &im1, const NiftiImage &im2);

/// Dump info of multiple nifti images
void dump_headers_actual(const vector<NiftiImage> &ims)
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
    dump_nifti_element(ims, "num_ext",           &nifti_image::num_ext);
    dump_nifti_element(ims, "nvox",              &nifti_image::nvox);
    dump_nifti_element(ims, "nx",                &nifti_image::nx);
    dump_nifti_element(ims, "ny",                &nifti_image::ny);
    dump_nifti_element(ims, "nz",                &nifti_image::nz);
    dump_nifti_element(ims, "nt",                &nifti_image::nt);
    dump_nifti_element(ims, "nu",                &nifti_image::nu);
    dump_nifti_element(ims, "nv",                &nifti_image::nv);
    dump_nifti_element(ims, "nw",                &nifti_image::nw);
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

    vector<shared_ptr<nifti_image> > images;
    for(unsigned i=0;i<ims.size();i++)
        images.push_back(ims[i].get_raw_nifti_sptr());

    // Print transformation matrices
    vector<SIRFRegMat44> qto_ijk_vec, qto_xyz_vec, sto_ijk_vec, sto_xyz_vec;
    for(unsigned j=0; j<images.size(); j++) {
        qto_ijk_vec.push_back(images[j]->qto_ijk);
        qto_xyz_vec.push_back(images[j]->qto_xyz);
        sto_ijk_vec.push_back(images[j]->sto_ijk);
        sto_xyz_vec.push_back(images[j]->sto_xyz);
    }
    cout << "\t" << left << setw(19) << "qto_ijk:" << "\n";
    SIRFRegMat44::print(qto_ijk_vec);
    cout << "\t" << left << setw(19) << "qto_xyz:" << "\n";
    SIRFRegMat44::print(qto_xyz_vec);
    cout << "\t" << left << setw(19) << "sto_ijk:" << "\n";
    SIRFRegMat44::print(sto_ijk_vec);
    cout << "\t" << left << setw(19) << "sto_xyz:" << "\n";
    SIRFRegMat44::print(sto_xyz_vec);

    cout << "\n";
}

template<typename T>
void dump_nifti_element(const vector<NiftiImage> &ims, const string &name, const T &call_back)
{
    string header = name + ": ";
    cout << "\t" << left << setw(19) << header;
    for(unsigned i=0; i<ims.size(); i++)
        cout << setw(19) << ims[i].get_raw_nifti_sptr().get()->*call_back;
    cout << "\n";
}

template<typename T>
void dump_nifti_element(const vector<NiftiImage> &ims, const string &name, const T &call_back, const unsigned num_elems)
{
    for(unsigned i=0; i<num_elems; i++) {
        string header = name + "[" + to_string(i) + "]: ";
        cout << "\t" << left << setw(19) << header;
        for(unsigned j=0; j<ims.size(); j++)
            cout << setw(19) << (ims[j].get_raw_nifti_sptr().get()->*call_back)[i];
        cout << "\n";
    }
}

template<typename T>
void fill_array(const NiftiImage &im, const float &v)
{
    if(!im.is_initialised())
        throw runtime_error("fill_array: Image not initialised.");

    // Check sizes
    if (im.get_raw_nifti_sptr()->nbyper != sizeof(T))
        throw runtime_error("fill_array: Datatype does not match desired cast type (" + to_string(im.get_raw_nifti_sptr()->nbyper) + " versus " + to_string(sizeof(T)) + ").");

    // Get data
    T *data = static_cast<T*>(im.get_raw_nifti_sptr()->data);
    for (unsigned i=0; i<im.get_raw_nifti_sptr()->nvox; i++) data[i] = T(v);
}
template void fill_array<bool>              (const NiftiImage &im, const float &v);
template void fill_array<signed char>       (const NiftiImage &im, const float &v);
template void fill_array<signed short>      (const NiftiImage &im, const float &v);
template void fill_array<signed int>        (const NiftiImage &im, const float &v);
template void fill_array<float>             (const NiftiImage &im, const float &v);
template void fill_array<double>            (const NiftiImage &im, const float &v);
template void fill_array<unsigned char>     (const NiftiImage &im, const float &v);
template void fill_array<unsigned short>    (const NiftiImage &im, const float &v);
template void fill_array<unsigned int>      (const NiftiImage &im, const float &v);
template void fill_array<signed long long>  (const NiftiImage &im, const float &v);
template void fill_array<unsigned long long>(const NiftiImage &im, const float &v);
template void fill_array<long double>       (const NiftiImage &im, const float &v);
}
