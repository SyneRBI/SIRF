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
#include "_reg_globalTransformation.h"
#include <_reg_tools.h>
#include <iostream>

using namespace std;

namespace SIRFRegMisc {

/// Open nifti image
void open_nifti_image(shared_ptr<nifti_image> &image, const boost::filesystem::path filename)
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
void save_nifti_image(shared_ptr<nifti_image> image, const string filename)
{
    if (!image) {
        throw runtime_error("Cannot save image to file.");
    }

    cout << "\nSaving image to file (" << filename << ")..." << flush;

    boost::filesystem::path filename_boost(filename);

    // If the folder doesn't exist, create it
    if (!boost::filesystem::exists(filename_boost.parent_path())) {
        cout << "\n\tCreating folder: \"" << filename_boost.parent_path().string() << "\"\n" << flush;
        boost::filesystem::create_directory(filename_boost.parent_path());
    }

    nifti_set_filenames(image.get(), filename.c_str(), 0, 0);
    nifti_image_write(image.get());
    cout << "done.\n\n";
}

/// Copy nifti image
void copy_nifti_image(const string input_filename, const string output_filename)
{
    shared_ptr<nifti_image> image;
    open_nifti_image(image, input_filename);
    save_nifti_image(image, output_filename);
}

/// Copy nifti image
void copy_nifti_image(shared_ptr<nifti_image> &output_image_sptr, const shared_ptr<nifti_image> &image_to_copy_sptr)
{
    const char *fname = image_to_copy_sptr->fname;
    open_nifti_image(output_image_sptr, fname);
}

/// Do nifti images match?
bool do_nift_image_match(const shared_ptr<nifti_image> &im1_sptr, const shared_ptr<nifti_image> &im2_sptr)
{
    bool images_match = true;
    if( im1_sptr->analyze75_orient  != im2_sptr->analyze75_orient   ) { images_match = false; cout << "mismatch in analyze75_orient , (values: " <<  im1_sptr->analyze75_orient << " and " << im2_sptr->analyze75_orient << ")\n"; }
    if( im1_sptr->byteorder         != im2_sptr->byteorder          ) { images_match = false; cout << "mismatch in byteorder , (values: " <<  im1_sptr->byteorder << " and " << im2_sptr->byteorder << ")\n"; }
    if( im1_sptr->cal_max           != im2_sptr->cal_max            ) { images_match = false; cout << "mismatch in cal_max , (values: " <<  im1_sptr->cal_max << " and " << im2_sptr->cal_max << ")\n"; }
    if( im1_sptr->cal_min          != im2_sptr->cal_min           ) { images_match = false; cout << "mismatch in cal_min , (values: " <<  im1_sptr->cal_min << " and " << im2_sptr->cal_min << ")\n"; }
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

    for (int i=0; i<4; i++) {
        for (int j=0; j<4; j++) {
            if( im1_sptr->qto_ijk.m[i][j] != im2_sptr->qto_ijk.m[i][j] ) { images_match = false; cout << "mismatch in qto_ijk["<<i<<"]["<<j<<"] , (values: " <<  im1_sptr->qto_ijk.m[i][j] << " and " << im2_sptr->qto_ijk.m[i][j] << ")\n"; }
            if( im1_sptr->qto_xyz.m[i][j] != im2_sptr->qto_xyz.m[i][j] ) { images_match = false; cout << "mismatch in qto_xyz["<<i<<"]["<<j<<"] , (values: " <<  im1_sptr->qto_xyz.m[i][j] << " and " << im2_sptr->qto_xyz.m[i][j] << ")\n"; }
            if( im1_sptr->sto_ijk.m[i][j] != im2_sptr->sto_ijk.m[i][j] ) { images_match = false; cout << "mismatch in sto_ijk["<<i<<"]["<<j<<"] , (values: " <<  im1_sptr->sto_ijk.m[i][j] << " and " << im2_sptr->sto_ijk.m[i][j] << ")\n"; }
            if( im1_sptr->sto_xyz.m[i][j] != im2_sptr->sto_xyz.m[i][j] ) { images_match = false; cout << "mismatch in sto_xyz["<<i<<"]["<<j<<"] , (values: " <<  im1_sptr->sto_xyz.m[i][j] << " and " << im2_sptr->sto_xyz.m[i][j] << ")\n"; }
        }
    }

    for (int i=0; i<8; i++) {
        if (   im1_sptr->dim[i] !=   im2_sptr->dim[i]  ) { images_match = false; cout << "mismatch in dim[" <<i<<  "] , (values: " <<   im1_sptr->dim[i]   << " and " <<  im2_sptr->dim[i]   << ")\n"; }
        if ( im1_sptr->pixdim[i] != im2_sptr->pixdim[i]) { images_match = false; cout << "mismatch in pixdim["<<i<<"] , (values: " <<  im1_sptr->pixdim[i] << " and " << im2_sptr->pixdim[i] << ")\n"; }
    }

    if (images_match) cout << "\tOK!\n";

    return images_match;
}

/// Dump info of nifti image
void dump_nifti_info(const std::string &im_filename)
{
    shared_ptr<nifti_image> image;
    SIRFRegMisc::open_nifti_image(image,im_filename);
    SIRFRegMisc::dump_nifti_info(image);
}

/// Dump info of nifti image
void dump_nifti_info(const shared_ptr<nifti_image> &im1_sptr)
{
    vector<shared_ptr<nifti_image> > image;
    image.push_back(im1_sptr);
    dump_nifti_info(image);
}

/// Dump info of multiple nifti images
void dump_nifti_info(const vector<shared_ptr<nifti_image> > &images)
{
    cout << "\nPrinting info for " << images.size() <<" nifti image:\n";
    cout << "\tanalyze_75_orient: \t"; for(int i=0;i<images.size();i++) { cout << images[i]->analyze75_orient << " \t"; } cout << "\n";
    cout << "\tanalyze75_orient: \t"; for(int i=0;i<images.size();i++) { cout << images[i]->analyze75_orient << " \t"; } cout << "\n";
    cout << "\tbyteorder: \t"; for(int i=0;i<images.size();i++) { cout << images[i]->byteorder << " \t"; } cout << "\n";
    cout << "\tcal_max: \t"; for(int i=0;i<images.size();i++) { cout << images[i]->cal_max << " \t"; } cout << "\n";
    cout << "\tcal_min: \t"; for(int i=0;i<images.size();i++) { cout << images[i]->cal_min << " \t"; } cout << "\n";
    cout << "\tdatatype: \t"; for(int i=0;i<images.size();i++) { cout << images[i]->datatype << " \t"; } cout << "\n";
    cout << "\tdu: \t"; for(int i=0;i<images.size();i++) { cout << images[i]->du << " \t"; } cout << "\n";
    cout << "\tdv: \t"; for(int i=0;i<images.size();i++) { cout << images[i]->dv << " \t"; } cout << "\n";
    cout << "\tdw: \t"; for(int i=0;i<images.size();i++) { cout << images[i]->dw << " \t"; } cout << "\n";
    cout << "\tdx: \t"; for(int i=0;i<images.size();i++) { cout << images[i]->dx << " \t"; } cout << "\n";
    cout << "\tdy: \t"; for(int i=0;i<images.size();i++) { cout << images[i]->dy << " \t"; } cout << "\n";
    cout << "\tdz: \t"; for(int i=0;i<images.size();i++) { cout << images[i]->dz << " \t"; } cout << "\n";
    cout << "\text_list: \t"; for(int i=0;i<images.size();i++) { cout << images[i]->ext_list << " \t"; } cout << "\n";
    cout << "\tfreq_dim: \t"; for(int i=0;i<images.size();i++) { cout << images[i]->freq_dim << " \t"; } cout << "\n";
    cout << "\tiname_offset: \t"; for(int i=0;i<images.size();i++) { cout << images[i]->iname_offset << " \t"; } cout << "\n";
    cout << "\tintent_code: \t"; for(int i=0;i<images.size();i++) { cout << images[i]->intent_code << " \t"; } cout << "\n";
    cout << "\tintent_p1: \t"; for(int i=0;i<images.size();i++) { cout << images[i]->intent_p1 << " \t"; } cout << "\n";
    cout << "\tintent_p2: \t"; for(int i=0;i<images.size();i++) { cout << images[i]->intent_p2 << " \t"; } cout << "\n";
    cout << "\tintent_p3: \t"; for(int i=0;i<images.size();i++) { cout << images[i]->intent_p3 << " \t"; } cout << "\n";
    cout << "\tnbyper: \t"; for(int i=0;i<images.size();i++) { cout << images[i]->nbyper << " \t"; } cout << "\n";
    cout << "\tndim: \t"; for(int i=0;i<images.size();i++) { cout << images[i]->ndim << " \t"; } cout << "\n";
    cout << "\tnifti_type: \t"; for(int i=0;i<images.size();i++) { cout << images[i]->nifti_type << " \t"; } cout << "\n";
    cout << "\tnt: \t"; for(int i=0;i<images.size();i++) { cout << images[i]->nt << " \t"; } cout << "\n";
    cout << "\tnu: \t"; for(int i=0;i<images.size();i++) { cout << images[i]->nu << " \t"; } cout << "\n";
    cout << "\tnum_ext: \t"; for(int i=0;i<images.size();i++) { cout << images[i]->num_ext << " \t"; } cout << "\n";
    cout << "\tnv: \t"; for(int i=0;i<images.size();i++) { cout << images[i]->nv << " \t"; } cout << "\n";
    cout << "\tnvox: \t"; for(int i=0;i<images.size();i++) { cout << images[i]->nvox << " \t"; } cout << "\n";
    cout << "\tnw: \t"; for(int i=0;i<images.size();i++) { cout << images[i]->nw << " \t"; } cout << "\n";
    cout << "\tnx: \t"; for(int i=0;i<images.size();i++) { cout << images[i]->nx << " \t"; } cout << "\n";
    cout << "\tny: \t"; for(int i=0;i<images.size();i++) { cout << images[i]->ny << " \t"; } cout << "\n";
    cout << "\tnz: \t"; for(int i=0;i<images.size();i++) { cout << images[i]->nz << " \t"; } cout << "\n";
    cout << "\tphase_dim: \t"; for(int i=0;i<images.size();i++) { cout << images[i]->phase_dim << " \t"; } cout << "\n";
    cout << "\tqfac: \t"; for(int i=0;i<images.size();i++) { cout << images[i]->qfac << " \t"; } cout << "\n";
    cout << "\tqform_code: \t"; for(int i=0;i<images.size();i++) { cout << images[i]->qform_code << " \t"; } cout << "\n";
    cout << "\tqoffset_x: \t"; for(int i=0;i<images.size();i++) { cout << images[i]->qoffset_x << " \t"; } cout << "\n";
    cout << "\tqoffset_y: \t"; for(int i=0;i<images.size();i++) { cout << images[i]->qoffset_y << " \t"; } cout << "\n";
    cout << "\tqoffset_z: \t"; for(int i=0;i<images.size();i++) { cout << images[i]->qoffset_z << " \t"; } cout << "\n";
    cout << "\tquatern_b: \t"; for(int i=0;i<images.size();i++) { cout << images[i]->quatern_b << " \t"; } cout << "\n";
    cout << "\tquatern_c: \t"; for(int i=0;i<images.size();i++) { cout << images[i]->quatern_c << " \t"; } cout << "\n";
    cout << "\tquatern_d: \t"; for(int i=0;i<images.size();i++) { cout << images[i]->quatern_d << " \t"; } cout << "\n";
    cout << "\tscl_inter: \t"; for(int i=0;i<images.size();i++) { cout << images[i]->scl_inter << " \t"; } cout << "\n";
    cout << "\tscl_slope: \t"; for(int i=0;i<images.size();i++) { cout << images[i]->scl_slope << " \t"; } cout << "\n";
    cout << "\tsform_code: \t"; for(int i=0;i<images.size();i++) { cout << images[i]->sform_code << " \t"; } cout << "\n";
    cout << "\tslice_code: \t"; for(int i=0;i<images.size();i++) { cout << images[i]->slice_code << " \t"; } cout << "\n";
    cout << "\tslice_dim: \t"; for(int i=0;i<images.size();i++) { cout << images[i]->slice_dim << " \t"; } cout << "\n";
    cout << "\tslice_duration: \t"; for(int i=0;i<images.size();i++) { cout << images[i]->slice_duration << " \t"; } cout << "\n";
    cout << "\tslice_end: \t"; for(int i=0;i<images.size();i++) { cout << images[i]->slice_end << " \t"; } cout << "\n";
    cout << "\tslice_start: \t"; for(int i=0;i<images.size();i++) { cout << images[i]->slice_start << " \t"; } cout << "\n";
    cout << "\tswapsize: \t"; for(int i=0;i<images.size();i++) { cout << images[i]->swapsize << " \t"; } cout << "\n";
    cout << "\ttime_units: \t"; for(int i=0;i<images.size();i++) { cout << images[i]->time_units << " \t"; } cout << "\n";
    cout << "\ttoffset: \t"; for(int i=0;i<images.size();i++) { cout << images[i]->toffset << " \t"; } cout << "\n";
    cout << "\txyz_units: \t"; for(int i=0;i<images.size();i++) { cout << images[i]->xyz_units << " \t"; } cout << "\n";
    cout << "\tqto_ijk:\n"; for(int i=0;i<4;i++) { cout << "\t\t"; for(int j=0;j<images.size();j++) { cout << "[" << images[j]->qto_ijk.m[0][i] << "," << images[j]->qto_ijk.m[1][i] << "," << images[j]->qto_ijk.m[2][i] << "," << images[j]->qto_ijk.m[3][i] << "]\t"; }; cout <<"\n"; }
    cout << "\tqto_xyz:\n"; for(int i=0;i<4;i++) { cout << "\t\t"; for(int j=0;j<images.size();j++) { cout << "[" << images[j]->qto_xyz.m[0][i] << "," << images[j]->qto_xyz.m[1][i] << "," << images[j]->qto_xyz.m[2][i] << "," << images[j]->qto_xyz.m[3][i] << "]\t"; }; cout <<"\n"; }
    cout << "\tsto_ijk:\n"; for(int i=0;i<4;i++) { cout << "\t\t"; for(int j=0;j<images.size();j++) { cout << "[" << images[j]->sto_ijk.m[0][i] << "," << images[j]->sto_ijk.m[1][i] << "," << images[j]->sto_ijk.m[2][i] << "," << images[j]->sto_ijk.m[3][i] << "]\t"; }; cout <<"\n"; }
    cout << "\tsto_xyz:\n"; for(int i=0;i<4;i++) { cout << "\t\t"; for(int j=0;j<images.size();j++) { cout << "[" << images[j]->sto_xyz.m[0][i] << "," << images[j]->sto_xyz.m[1][i] << "," << images[j]->sto_xyz.m[2][i] << "," << images[j]->sto_xyz.m[3][i] << "]\t"; }; cout <<"\n"; }
    for(int i=0;i<8;i++) { cout << "\tdim[" << i << "]: \t"; for(int j=0;j<images.size();j++) { cout << images[j]->dim[i] << "\t"; } cout << "\n"; }
    for(int i=0;i<8;i++) { cout << "\tpixdim[" << i << "]: \t"; for(int j=0;j<images.size();j++) { cout << images[j]->pixdim[i] << "\t"; } cout << "\n"; }
    cout << "\n";
}

/// Print info of element of nifti image
void print_nifti_info(const string &im_filename, const string keyword)
{
    shared_ptr<nifti_image> image;
    SIRFRegMisc::open_nifti_image(image,im_filename);
    SIRFRegMisc::print_nifti_info(image,keyword);
}

/// Print info of element of nifti image
void print_nifti_info(const shared_ptr<nifti_image> &im1_sptr, const string keyword)
{
    if (keyword == "analyze75_orient") cout << "analyze75_orient: " << im1_sptr->analyze75_orient << "\n";
    if (keyword == "analyze75_orient") cout << "analyze75_orient: " << im1_sptr->analyze75_orient << "\n";
    if (keyword == "byteorder") cout << "byteorder: " << im1_sptr->byteorder << "\n";
    if (keyword == "cal_max") cout << "cal_max: " << im1_sptr->cal_max << "\n";
    if (keyword == "cal_min") cout << "cal_min: " << im1_sptr->cal_min << "\n";
    if (keyword == "datatype") cout << "datatype: " << im1_sptr->datatype << "\n";
    if (keyword == "du") cout << "du: " << im1_sptr->du << "\n";
    if (keyword == "dv") cout << "dv: " << im1_sptr->dv << "\n";
    if (keyword == "dw") cout << "dw: " << im1_sptr->dw << "\n";
    if (keyword == "dx") cout << "dx: " << im1_sptr->dx << "\n";
    if (keyword == "dy") cout << "dy: " << im1_sptr->dy << "\n";
    if (keyword == "dz") cout << "dz: " << im1_sptr->dz << "\n";
    if (keyword == "ext_list") cout << "ext_list: " << im1_sptr->ext_list << "\n";
    if (keyword == "freq_dim") cout << "freq_dim: " << im1_sptr->freq_dim << "\n";
    if (keyword == "iname_offset") cout << "iname_offset: " << im1_sptr->iname_offset << "\n";
    if (keyword == "intent_code") cout << "intent_code: " << im1_sptr->intent_code << "\n";
    if (keyword == "intent_p1") cout << "intent_p1: " << im1_sptr->intent_p1 << "\n";
    if (keyword == "intent_p2") cout << "intent_p2: " << im1_sptr->intent_p2 << "\n";
    if (keyword == "intent_p3") cout << "intent_p3: " << im1_sptr->intent_p3 << "\n";
    if (keyword == "nbyper") cout << "nbyper: " << im1_sptr->nbyper << "\n";
    if (keyword == "ndim") cout << "ndim: " << im1_sptr->ndim << "\n";
    if (keyword == "nifti_type") cout << "nifti_type: " << im1_sptr->nifti_type << "\n";
    if (keyword == "nt") cout << "nt: " << im1_sptr->nt << "\n";
    if (keyword == "nu") cout << "nu: " << im1_sptr->nu << "\n";
    if (keyword == "num_ext") cout << "num_ext: " << im1_sptr->num_ext << "\n";
    if (keyword == "nv") cout << "nv: " << im1_sptr->nv << "\n";
    if (keyword == "nvox") cout << "nvox: " << im1_sptr->nvox << "\n";
    if (keyword == "nw") cout << "nw: " << im1_sptr->nw << "\n";
    if (keyword == "nx") cout << "nx: " << im1_sptr->nx << "\n";
    if (keyword == "ny") cout << "ny: " << im1_sptr->ny << "\n";
    if (keyword == "nz") cout << "nz: " << im1_sptr->nz << "\n";
    if (keyword == "phase_dim") cout << "phase_dim: " << im1_sptr->phase_dim << "\n";
    if (keyword == "qfac") cout << "qfac: " << im1_sptr->qfac << "\n";
    if (keyword == "qform_code") cout << "qform_code: " << im1_sptr->qform_code << "\n";
    if (keyword == "qoffset_x") cout << "qoffset_x: " << im1_sptr->qoffset_x << "\n";
    if (keyword == "qoffset_y") cout << "qoffset_y: " << im1_sptr->qoffset_y << "\n";
    if (keyword == "qoffset_z") cout << "qoffset_z: " << im1_sptr->qoffset_z << "\n";
    if (keyword == "quatern_b") cout << "quatern_b: " << im1_sptr->quatern_b << "\n";
    if (keyword == "quatern_c") cout << "quatern_c: " << im1_sptr->quatern_c << "\n";
    if (keyword == "quatern_d") cout << "quatern_d: " << im1_sptr->quatern_d << "\n";
    if (keyword == "scl_inter") cout << "scl_inter: " << im1_sptr->scl_inter << "\n";
    if (keyword == "scl_slope") cout << "scl_slope: " << im1_sptr->scl_slope << "\n";
    if (keyword == "sform_code") cout << "sform_code: " << im1_sptr->sform_code << "\n";
    if (keyword == "slice_code") cout << "slice_code: " << im1_sptr->slice_code << "\n";
    if (keyword == "slice_dim") cout << "slice_dim: " << im1_sptr->slice_dim << "\n";
    if (keyword == "slice_duration") cout << "slice_duration: " << im1_sptr->slice_duration << "\n";
    if (keyword == "slice_end") cout << "slice_end: " << im1_sptr->slice_end << "\n";
    if (keyword == "slice_start") cout << "slice_start: " << im1_sptr->slice_start << "\n";
    if (keyword == "swapsize") cout << "swapsize: " << im1_sptr->swapsize << "\n";
    if (keyword == "time_units") cout << "time_units: " << im1_sptr->time_units << "\n";
    if (keyword == "toffset") cout << "toffset: " << im1_sptr->toffset << "\n";
    if (keyword == "xyz_units") cout << "xyz_units: " << im1_sptr->xyz_units << "\n";

    if (keyword == "qto_ijk") cout << "qto_ijk:\n"; SIRFRegMisc::print_mat44(&im1_sptr->qto_ijk);
    if (keyword == "qto_xyz") cout << "qto_xyz:\n"; SIRFRegMisc::print_mat44(&im1_sptr->qto_xyz);
    if (keyword == "sto_ijk") cout << "sto_ijk:\n"; SIRFRegMisc::print_mat44(&im1_sptr->sto_ijk);
    if (keyword == "sto_xyz") cout << "sto_xyz:\n"; SIRFRegMisc::print_mat44(&im1_sptr->sto_xyz);
    if (keyword == "dim")     for (int i=0;i<8;i++) cout << "dim["    << i << "]: " << im1_sptr->dim[i]    << "\n";
    if (keyword == "pixdim")  for (int i=0;i<8;i++) cout << "pixdim[" << i << "]: " << im1_sptr->pixdim[i] << "\n";
}

/// Save transformation matrix to file
void save_transformation_matrix(shared_ptr<mat44> &transformation_matrix_sptr, const string filename)
{
    // Check that the matrix exists
    if (!transformation_matrix_sptr) {
        throw std::runtime_error("Transformation matrix is null pointer. Have you run the registration?");
    }

    // Check that input isn't blank
    if (filename == "") {
        throw std::runtime_error("Error, cannot write transformation matrix to file because filename is blank");
    }

    reg_tool_WriteAffineFile(transformation_matrix_sptr.get(), filename.c_str());
}

/// Read transformation matrix from file
void open_transformation_matrix(shared_ptr<mat44> &transformation_matrix_sptr, const string filename)
{
    // Check that the file exists
    if (!boost::filesystem::exists(filename)) {
        throw runtime_error("Cannot find the file: " + filename + ".");
    }

    cout << "\n\nReading transformation matrix from file...\n\n";

    // Initialise the matrix
    transformation_matrix_sptr = make_shared<mat44>();

    reg_tool_ReadAffineFile(transformation_matrix_sptr.get(), (char*)filename.c_str());

    cout << "\n\nSuccessfully read transformation matrix from file:\n";

    for (int i=0; i<4; i++) {
        cout << "\t" << flush;
        cout << transformation_matrix_sptr->m[i][0] << ", " << flush;
        cout << transformation_matrix_sptr->m[i][1] << ", " << flush;
        cout << transformation_matrix_sptr->m[i][2] << ", " << flush;
        cout << transformation_matrix_sptr->m[i][3] << "\n" << flush;
    }
}

/// Print mat44
void print_mat44(const mat44 *mat_ptr)
{
    for (int i=0; i<4; i++) {
        cout << "[";
        cout << mat_ptr->m[i][0] << " "   << flush;
        cout << mat_ptr->m[i][1] << " "   << flush;
        cout << mat_ptr->m[i][2] << " "   << flush;
        cout << mat_ptr->m[i][3] << "]\n" << flush;
    }
}

/// Mat44 multiplier
mat44 multiply_mat44(const mat44 &x, const mat44 &y)
{
    // Print info
    cout << "\nMultiplying two matrices...\n";
    cout << "Matrix 1:\n";
    print_mat44(&x);
    cout << "Matrix 2:\n";
    print_mat44(&y);

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
    print_mat44(&res);

    return res;
}

}
