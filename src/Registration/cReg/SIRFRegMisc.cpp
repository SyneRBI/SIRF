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

    std:shared_ptr<nifti_image> im_sptr = im.get_image_as_nifti();

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
    double start_index =   dim   * im_sptr->nvox/3;
    // End index is one before the start of the next dimension (thus the minus 1)
    double end_index   = (dim+1) * im_sptr->nvox/3 - 1;

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
void get_cpp_from_transformation_matrix(std::shared_ptr<nifti_image> &cpp_sptr, const std::shared_ptr<mat44> &TM_sptr, const std::shared_ptr<nifti_image> &warped_sptr)
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
}

/// Convert from displacement to deformation field image
void convert_from_disp_to_def(SIRFImageDataDeformation &im)
{
    // Get the def field from the disp field
    reg_getDeformationFromDisplacement(im.get_image_as_nifti().get());
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

/// Do nifti images match?
bool do_nift_image_match(const SIRFImageData &im1, const SIRFImageData &im2)
{
    const std::shared_ptr<nifti_image> im1_sptr = im1.get_image_as_nifti();
    const std::shared_ptr<nifti_image> im2_sptr = im2.get_image_as_nifti();

    bool images_match = true;
    if( im1_sptr->analyze75_orient  != im2_sptr->analyze75_orient   ) { images_match = false; cout << "mismatch in analyze75_orient , (values: " <<  im1_sptr->analyze75_orient << " and " << im2_sptr->analyze75_orient << ")\n"; }
    if( im1_sptr->byteorder         != im2_sptr->byteorder          ) { images_match = false; cout << "mismatch in byteorder , (values: " <<  im1_sptr->byteorder << " and " << im2_sptr->byteorder << ")\n"; }
    if( im1_sptr->cal_max           != im2_sptr->cal_max            ) { images_match = false; cout << "mismatch in cal_max , (values: " <<  im1_sptr->cal_max << " and " << im2_sptr->cal_max << ")\n"; }
    if( im1_sptr->cal_min           != im2_sptr->cal_min            ) { images_match = false; cout << "mismatch in cal_min , (values: " <<  im1_sptr->cal_min << " and " << im2_sptr->cal_min << ")\n"; }
    if( im1_sptr->datatype          != im2_sptr->datatype           ) { images_match = false; cout << "mismatch in datatype , (values: " <<  im1_sptr->datatype << " and " << im2_sptr->datatype << ")\n"; }
    if( im1_sptr->du                != im2_sptr->du                 ) { images_match = false; cout << "mismatch in du , (values: " <<  im1_sptr->du << " and " << im2_sptr->du << ")\n"; }
    if( im1_sptr->dv                != im2_sptr->dv                 ) { images_match = false; cout << "mismatch in dv , (values: " <<  im1_sptr->dv << " and " << im2_sptr->dv << ")\n"; }
    if( im1_sptr->dw                != im2_sptr->dw                 ) { images_match = false; cout << "mismatch in dw , (values: " <<  im1_sptr->dw << " and " << im2_sptr->dw << ")\n"; }
    if( im1_sptr->dx                != im2_sptr->dx                 ) { images_match = false; cout << "mismatch in dx , (values: " <<  im1_sptr->dx << " and " << im2_sptr->dx << ")\n"; }
    if( im1_sptr->dy                != im2_sptr->dy                 ) { images_match = false; cout << "mismatch in dy , (values: " <<  im1_sptr->dy << " and " << im2_sptr->dy << ")\n"; }
    if( im1_sptr->dz                != im2_sptr->dz                 ) { images_match = false; cout << "mismatch in dz , (values: " <<  im1_sptr->dz << " and " << im2_sptr->dz << ")\n"; }
    if( im1_sptr->ext_list          != im2_sptr->ext_list           ) { images_match = false; cout << "mismatch in ext_list , (values: " <<  im1_sptr->ext_list << " and " << im2_sptr->ext_list << ")\n"; }
    if( im1_sptr->freq_dim          != im2_sptr->freq_dim           ) { images_match = false; cout << "mismatch in freq_dim , (values: " <<  im1_sptr->freq_dim << " and " << im2_sptr->freq_dim << ")\n"; }
    if( im1_sptr->iname_offset      != im2_sptr->iname_offset       ) { images_match = false; cout << "mismatch in iname_offset , (values: " <<  im1_sptr->iname_offset << " and " << im2_sptr->iname_offset << ")\n"; }
    if( im1_sptr->intent_code       != im2_sptr->intent_code        ) { images_match = false; cout << "mismatch in intent_code , (values: " <<  im1_sptr->intent_code << " and " << im2_sptr->intent_code << ")\n"; }
    if( im1_sptr->intent_p1         != im2_sptr->intent_p1          ) { images_match = false; cout << "mismatch in intent_p1 , (values: " <<  im1_sptr->intent_p1 << " and " << im2_sptr->intent_p1 << ")\n"; }
    if( im1_sptr->intent_p2         != im2_sptr->intent_p2          ) { images_match = false; cout << "mismatch in intent_p2 , (values: " <<  im1_sptr->intent_p2 << " and " << im2_sptr->intent_p2 << ")\n"; }
    if( im1_sptr->intent_p3         != im2_sptr->intent_p3          ) { images_match = false; cout << "mismatch in intent_p3 , (values: " <<  im1_sptr->intent_p3 << " and " << im2_sptr->intent_p3 << ")\n"; }
    if( im1_sptr->nbyper            != im2_sptr->nbyper             ) { images_match = false; cout << "mismatch in nbyper , (values: " <<  im1_sptr->nbyper << " and " << im2_sptr->nbyper << ")\n"; }
    if( im1_sptr->ndim              != im2_sptr->ndim               ) { images_match = false; cout << "mismatch in ndim , (values: " <<  im1_sptr->ndim << " and " << im2_sptr->ndim << ")\n"; }
    if( im1_sptr->nifti_type        != im2_sptr->nifti_type         ) { images_match = false; cout << "mismatch in nifti_type , (values: " <<  im1_sptr->nifti_type << " and " << im2_sptr->nifti_type << ")\n"; }
    if( im1_sptr->nt                != im2_sptr->nt                 ) { images_match = false; cout << "mismatch in nt , (values: " <<  im1_sptr->nt << " and " << im2_sptr->nt << ")\n"; }
    if( im1_sptr->nu                != im2_sptr->nu                 ) { images_match = false; cout << "mismatch in nu , (values: " <<  im1_sptr->nu << " and " << im2_sptr->nu << ")\n"; }
    if( im1_sptr->num_ext           != im2_sptr->num_ext            ) { images_match = false; cout << "mismatch in num_ext , (values: " <<  im1_sptr->num_ext << " and " << im2_sptr->num_ext << ")\n"; }
    if( im1_sptr->nv                != im2_sptr->nv                 ) { images_match = false; cout << "mismatch in nv , (values: " <<  im1_sptr->nv << " and " << im2_sptr->nv << ")\n"; }
    if( im1_sptr->nvox              != im2_sptr->nvox               ) { images_match = false; cout << "mismatch in nvox , (values: " <<  im1_sptr->nvox << " and " << im2_sptr->nvox << ")\n"; }
    if( im1_sptr->nw                != im2_sptr->nw                 ) { images_match = false; cout << "mismatch in nw , (values: " <<  im1_sptr->nw << " and " << im2_sptr->nw << ")\n"; }
    if( im1_sptr->nx                != im2_sptr->nx                 ) { images_match = false; cout << "mismatch in nx , (values: " <<  im1_sptr->nx << " and " << im2_sptr->nx << ")\n"; }
    if( im1_sptr->ny                != im2_sptr->ny                 ) { images_match = false; cout << "mismatch in ny , (values: " <<  im1_sptr->ny << " and " << im2_sptr->ny << ")\n"; }
    if( im1_sptr->nz                != im2_sptr->nz                 ) { images_match = false; cout << "mismatch in nz , (values: " <<  im1_sptr->nz << " and " << im2_sptr->nz << ")\n"; }
    if( im1_sptr->phase_dim         != im2_sptr->phase_dim          ) { images_match = false; cout << "mismatch in phase_dim , (values: " <<  im1_sptr->phase_dim << " and " << im2_sptr->phase_dim << ")\n"; }
    if( im1_sptr->qfac              != im2_sptr->qfac               ) { images_match = false; cout << "mismatch in qfac , (values: " <<  im1_sptr->qfac << " and " << im2_sptr->qfac << ")\n"; }
    if( im1_sptr->qform_code        != im2_sptr->qform_code         ) { images_match = false; cout << "mismatch in qform_code , (values: " <<  im1_sptr->qform_code << " and " << im2_sptr->qform_code << ")\n"; }
    if( im1_sptr->qoffset_x         != im2_sptr->qoffset_x          ) { images_match = false; cout << "mismatch in qoffset_x , (values: " <<  im1_sptr->qoffset_x << " and " << im2_sptr->qoffset_x << ")\n"; }
    if( im1_sptr->qoffset_y         != im2_sptr->qoffset_y          ) { images_match = false; cout << "mismatch in qoffset_y , (values: " <<  im1_sptr->qoffset_y << " and " << im2_sptr->qoffset_y << ")\n"; }
    if( im1_sptr->qoffset_z         != im2_sptr->qoffset_z          ) { images_match = false; cout << "mismatch in qoffset_z , (values: " <<  im1_sptr->qoffset_z << " and " << im2_sptr->qoffset_z << ")\n"; }
    if( im1_sptr->quatern_b         != im2_sptr->quatern_b          ) { images_match = false; cout << "mismatch in quatern_b , (values: " <<  im1_sptr->quatern_b << " and " << im2_sptr->quatern_b << ")\n"; }
    if( im1_sptr->quatern_c         != im2_sptr->quatern_c          ) { images_match = false; cout << "mismatch in quatern_c , (values: " <<  im1_sptr->quatern_c << " and " << im2_sptr->quatern_c << ")\n"; }
    if( im1_sptr->quatern_d         != im2_sptr->quatern_d          ) { images_match = false; cout << "mismatch in quatern_d , (values: " <<  im1_sptr->quatern_d << " and " << im2_sptr->quatern_d << ")\n"; }
    if( im1_sptr->scl_inter         != im2_sptr->scl_inter          ) { images_match = false; cout << "mismatch in scl_inter , (values: " <<  im1_sptr->scl_inter << " and " << im2_sptr->scl_inter << ")\n"; }
    if( im1_sptr->scl_slope         != im2_sptr->scl_slope          ) { images_match = false; cout << "mismatch in scl_slope , (values: " <<  im1_sptr->scl_slope << " and " << im2_sptr->scl_slope << ")\n"; }
    if( im1_sptr->sform_code        != im2_sptr->sform_code         ) { images_match = false; cout << "mismatch in sform_code , (values: " <<  im1_sptr->sform_code << " and " << im2_sptr->sform_code << ")\n"; }
    if( im1_sptr->slice_code        != im2_sptr->slice_code         ) { images_match = false; cout << "mismatch in slice_code , (values: " <<  im1_sptr->slice_code << " and " << im2_sptr->slice_code << ")\n"; }
    if( im1_sptr->slice_dim         != im2_sptr->slice_dim          ) { images_match = false; cout << "mismatch in slice_dim , (values: " <<  im1_sptr->slice_dim << " and " << im2_sptr->slice_dim << ")\n"; }
    if( im1_sptr->slice_duration    != im2_sptr->slice_duration     ) { images_match = false; cout << "mismatch in slice_duration , (values: " <<  im1_sptr->slice_duration << " and " << im2_sptr->slice_duration << ")\n"; }
    if( im1_sptr->slice_end         != im2_sptr->slice_end          ) { images_match = false; cout << "mismatch in slice_end , (values: " <<  im1_sptr->slice_end << " and " << im2_sptr->slice_end << ")\n"; }
    if( im1_sptr->slice_start       != im2_sptr->slice_start        ) { images_match = false; cout << "mismatch in slice_start , (values: " <<  im1_sptr->slice_start << " and " << im2_sptr->slice_start << ")\n"; }
    if( im1_sptr->swapsize          != im2_sptr->swapsize           ) { images_match = false; cout << "mismatch in swapsize , (values: " <<  im1_sptr->swapsize << " and " << im2_sptr->swapsize << ")\n"; }
    if( im1_sptr->time_units        != im2_sptr->time_units         ) { images_match = false; cout << "mismatch in time_units , (values: " <<  im1_sptr->time_units << " and " << im2_sptr->time_units << ")\n"; }
    if( im1_sptr->toffset           != im2_sptr->toffset            ) { images_match = false; cout << "mismatch in toffset , (values: " <<  im1_sptr->toffset << " and " << im2_sptr->toffset << ")\n"; }
    if( im1_sptr->xyz_units         != im2_sptr->xyz_units          ) { images_match = false; cout << "mismatch in xyz_units , (values: " <<  im1_sptr->xyz_units << " and " << im2_sptr->xyz_units << ")\n"; }
    if( !do_mat44_match(im1_sptr->qto_ijk, im2_sptr->qto_ijk)       ) { images_match = false; cout << "mismatch in qto_ijk\n"; vector<mat44>vec; vec.push_back(im1_sptr->qto_ijk); vec.push_back(im2_sptr->qto_ijk); print_mat44(vec); cout << "\n"; }
    if( !do_mat44_match(im1_sptr->qto_xyz, im2_sptr->qto_xyz)       ) { images_match = false; cout << "mismatch in qto_xyz\n"; vector<mat44>vec; vec.push_back(im1_sptr->qto_xyz); vec.push_back(im2_sptr->qto_xyz); print_mat44(vec); cout << "\n"; }
    if( !do_mat44_match(im1_sptr->sto_ijk, im2_sptr->sto_ijk)       ) { images_match = false; cout << "mismatch in sto_ijk\n"; vector<mat44>vec; vec.push_back(im1_sptr->sto_ijk); vec.push_back(im2_sptr->sto_ijk); print_mat44(vec); cout << "\n"; }
    if( !do_mat44_match(im1_sptr->sto_xyz, im2_sptr->sto_xyz)       ) { images_match = false; cout << "mismatch in sto_xyz\n"; vector<mat44>vec; vec.push_back(im1_sptr->sto_xyz); vec.push_back(im2_sptr->sto_xyz); print_mat44(vec); cout << "\n"; }

    for (int i=0; i<8; i++) {
        if (   im1_sptr->dim[i] !=   im2_sptr->dim[i]  ) { images_match = false; cout << "mismatch in dim[" <<i<<  "] , (values: " <<   im1_sptr->dim[i]   << " and " <<  im2_sptr->dim[i]   << ")\n"; }
        if ( im1_sptr->pixdim[i] != im2_sptr->pixdim[i]) { images_match = false; cout << "mismatch in pixdim["<<i<<"] , (values: " <<  im1_sptr->pixdim[i] << " and " << im2_sptr->pixdim[i] << ")\n"; }
    }

    if (images_match) cout << "\tOK!\n";

    return images_match;
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
    vector<std::shared_ptr<nifti_image> > images;
    for(int i=0;i<ims.size();i++)
        images.push_back(ims[i].get_image_as_nifti());

    cout << "\nPrinting info for " << images.size() <<" nifti image:\n";
    cout << "\t" << left << setw(19) << "analyze_75_orient:"; for(int i=0;i<images.size();i++) { cout << setw(19) << images[i]->analyze75_orient; } cout << "\n";
    cout << "\t" << left << setw(19) << "analyze75_orient:"; for(int i=0;i<images.size();i++) { cout << setw(19) << images[i]->analyze75_orient; } cout << "\n";
    cout << "\t" << left << setw(19) << "byteorder:"; for(int i=0;i<images.size();i++) { cout << setw(19) << images[i]->byteorder; } cout << "\n";
    cout << "\t" << left << setw(19) << "cal_max:"; for(int i=0;i<images.size();i++) { cout << setw(19) << images[i]->cal_max; } cout << "\n";
    cout << "\t" << left << setw(19) << "cal_min:"; for(int i=0;i<images.size();i++) { cout << setw(19) << images[i]->cal_min; } cout << "\n";
    cout << "\t" << left << setw(19) << "datatype:"; for(int i=0;i<images.size();i++) { cout << setw(19) << images[i]->datatype; } cout << "\n";
    cout << "\t" << left << setw(19) << "du:"; for(int i=0;i<images.size();i++) { cout << setw(19) << images[i]->du; } cout << "\n";
    cout << "\t" << left << setw(19) << "dv:"; for(int i=0;i<images.size();i++) { cout << setw(19) << images[i]->dv; } cout << "\n";
    cout << "\t" << left << setw(19) << "dw:"; for(int i=0;i<images.size();i++) { cout << setw(19) << images[i]->dw; } cout << "\n";
    cout << "\t" << left << setw(19) << "dx:"; for(int i=0;i<images.size();i++) { cout << setw(19) << images[i]->dx; } cout << "\n";
    cout << "\t" << left << setw(19) << "dy:"; for(int i=0;i<images.size();i++) { cout << setw(19) << images[i]->dy; } cout << "\n";
    cout << "\t" << left << setw(19) << "dz:"; for(int i=0;i<images.size();i++) { cout << setw(19) << images[i]->dz; } cout << "\n";
    cout << "\t" << left << setw(19) << "ext_list:"; for(int i=0;i<images.size();i++) { cout << setw(19) << images[i]->ext_list; } cout << "\n";
    cout << "\t" << left << setw(19) << "fname:"; for(int i=0;i<images.size();i++) { if(images[i]->fname) cout << setw(19) << images[i]->fname << std::flush; } cout << "\n";
    cout << "\t" << left << setw(19) << "freq_dim:"; for(int i=0;i<images.size();i++) { cout << setw(19) << images[i]->freq_dim; } cout << "\n";
    cout << "\t" << left << setw(19) << "iname_offset:"; for(int i=0;i<images.size();i++) { cout << setw(19) << images[i]->iname_offset; } cout << "\n";
    cout << "\t" << left << setw(19) << "intent_code:"; for(int i=0;i<images.size();i++) { cout << setw(19) << images[i]->intent_code; } cout << "\n";
    cout << "\t" << left << setw(19) << "intent_p1:"; for(int i=0;i<images.size();i++) { cout << setw(19) << images[i]->intent_p1; } cout << "\n";
    cout << "\t" << left << setw(19) << "intent_p2:"; for(int i=0;i<images.size();i++) { cout << setw(19) << images[i]->intent_p2; } cout << "\n";
    cout << "\t" << left << setw(19) << "intent_p3:"; for(int i=0;i<images.size();i++) { cout << setw(19) << images[i]->intent_p3; } cout << "\n";
    cout << "\t" << left << setw(19) << "nbyper:"; for(int i=0;i<images.size();i++) { cout << setw(19) << images[i]->nbyper; } cout << "\n";
    cout << "\t" << left << setw(19) << "ndim:"; for(int i=0;i<images.size();i++) { cout << setw(19) << images[i]->ndim; } cout << "\n";
    cout << "\t" << left << setw(19) << "nifti_type:"; for(int i=0;i<images.size();i++) { cout << setw(19) << images[i]->nifti_type; } cout << "\n";
    cout << "\t" << left << setw(19) << "nt:"; for(int i=0;i<images.size();i++) { cout << setw(19) << images[i]->nt; } cout << "\n";
    cout << "\t" << left << setw(19) << "nu:"; for(int i=0;i<images.size();i++) { cout << setw(19) << images[i]->nu; } cout << "\n";
    cout << "\t" << left << setw(19) << "num_ext:"; for(int i=0;i<images.size();i++) { cout << setw(19) << images[i]->num_ext; } cout << "\n";
    cout << "\t" << left << setw(19) << "nv:"; for(int i=0;i<images.size();i++) { cout << setw(19) << images[i]->nv; } cout << "\n";
    cout << "\t" << left << setw(19) << "nvox:"; for(int i=0;i<images.size();i++) { cout << setw(19) << images[i]->nvox; } cout << "\n";
    cout << "\t" << left << setw(19) << "nw:"; for(int i=0;i<images.size();i++) { cout << setw(19) << images[i]->nw; } cout << "\n";
    cout << "\t" << left << setw(19) << "nx:"; for(int i=0;i<images.size();i++) { cout << setw(19) << images[i]->nx; } cout << "\n";
    cout << "\t" << left << setw(19) << "ny:"; for(int i=0;i<images.size();i++) { cout << setw(19) << images[i]->ny; } cout << "\n";
    cout << "\t" << left << setw(19) << "nz:"; for(int i=0;i<images.size();i++) { cout << setw(19) << images[i]->nz; } cout << "\n";
    cout << "\t" << left << setw(19) << "phase_dim:"; for(int i=0;i<images.size();i++) { cout << setw(19) << images[i]->phase_dim; } cout << "\n";
    cout << "\t" << left << setw(19) << "qfac:"; for(int i=0;i<images.size();i++) { cout << setw(19) << images[i]->qfac; } cout << "\n";
    cout << "\t" << left << setw(19) << "qform_code:"; for(int i=0;i<images.size();i++) { cout << setw(19) << images[i]->qform_code; } cout << "\n";
    cout << "\t" << left << setw(19) << "qoffset_x:"; for(int i=0;i<images.size();i++) { cout << setw(19) << images[i]->qoffset_x; } cout << "\n";
    cout << "\t" << left << setw(19) << "qoffset_y:"; for(int i=0;i<images.size();i++) { cout << setw(19) << images[i]->qoffset_y; } cout << "\n";
    cout << "\t" << left << setw(19) << "qoffset_z:"; for(int i=0;i<images.size();i++) { cout << setw(19) << images[i]->qoffset_z; } cout << "\n";
    cout << "\t" << left << setw(19) << "quatern_b:"; for(int i=0;i<images.size();i++) { cout << setw(19) << images[i]->quatern_b; } cout << "\n";
    cout << "\t" << left << setw(19) << "quatern_c:"; for(int i=0;i<images.size();i++) { cout << setw(19) << images[i]->quatern_c; } cout << "\n";
    cout << "\t" << left << setw(19) << "quatern_d:"; for(int i=0;i<images.size();i++) { cout << setw(19) << images[i]->quatern_d; } cout << "\n";
    cout << "\t" << left << setw(19) << "scl_inter:"; for(int i=0;i<images.size();i++) { cout << setw(19) << images[i]->scl_inter; } cout << "\n";
    cout << "\t" << left << setw(19) << "scl_slope:"; for(int i=0;i<images.size();i++) { cout << setw(19) << images[i]->scl_slope; } cout << "\n";
    cout << "\t" << left << setw(19) << "sform_code:"; for(int i=0;i<images.size();i++) { cout << setw(19) << images[i]->sform_code; } cout << "\n";
    cout << "\t" << left << setw(19) << "slice_code:"; for(int i=0;i<images.size();i++) { cout << setw(19) << images[i]->slice_code; } cout << "\n";
    cout << "\t" << left << setw(19) << "slice_dim:"; for(int i=0;i<images.size();i++) { cout << setw(19) << images[i]->slice_dim; } cout << "\n";
    cout << "\t" << left << setw(19) << "slice_duration:"; for(int i=0;i<images.size();i++) { cout << setw(19) << images[i]->slice_duration; } cout << "\n";
    cout << "\t" << left << setw(19) << "slice_end:"; for(int i=0;i<images.size();i++) { cout << setw(19) << images[i]->slice_end; } cout << "\n";
    cout << "\t" << left << setw(19) << "slice_start:"; for(int i=0;i<images.size();i++) { cout << setw(19) << images[i]->slice_start; } cout << "\n";
    cout << "\t" << left << setw(19) << "swapsize:"; for(int i=0;i<images.size();i++) { cout << setw(19) << images[i]->swapsize; } cout << "\n";
    cout << "\t" << left << setw(19) << "time_units:"; for(int i=0;i<images.size();i++) { cout << setw(19) << images[i]->time_units; } cout << "\n";
    cout << "\t" << left << setw(19) << "toffset:"; for(int i=0;i<images.size();i++) { cout << setw(19) << images[i]->toffset; } cout << "\n";
    cout << "\t" << left << setw(19) << "xyz_units:"; for(int i=0;i<images.size();i++) { cout << setw(19) << images[i]->xyz_units; } cout << "\n";
    for(int i=0;i<8;i++) { cout << "\tdim[" << i << "]:\t\t   "; for(int j=0;j<images.size();j++) { cout << setw(19) << images[j]->dim[i]; } cout << "\n"; }
    for(int i=0;i<8;i++) { cout << "\tpixdim[" << i << "]:\t   "; for(int j=0;j<images.size();j++) { cout << setw(19) << images[j]->pixdim[i]; } cout << "\n"; }

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
void save_transformation_matrix(const std::shared_ptr<mat44> &transformation_matrix_sptr, const string &filename)
{
    // Check that the matrix exists
    if (!transformation_matrix_sptr)
        throw runtime_error("Transformation matrix is null pointer, can't save.");

    // Check that input isn't blank
    if (filename == "")
        throw runtime_error("Error, cannot write transformation matrix to file because filename is blank");

    reg_tool_WriteAffineFile(transformation_matrix_sptr.get(), filename.c_str());
}

/// Read transformation matrix from file
void open_transformation_matrix(std::shared_ptr<mat44> &transformation_matrix_sptr, const string& filename)
{
    // Check that the file exists
    if (!boost::filesystem::exists(filename))
        throw runtime_error("Cannot find the file: " + filename + ".");

    cout << "\n\nReading transformation matrix from file...\n\n";
    transformation_matrix_sptr = make_shared<mat44>();
    reg_tool_ReadAffineFile(transformation_matrix_sptr.get(), (char*)filename.c_str());
    cout << "\n\nSuccessfully read transformation matrix from file:\n";

    print_mat44(*transformation_matrix_sptr);
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

}
