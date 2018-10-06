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
#include <_reg_globalTrans.h>
#include <_reg_localTrans.h>
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

/// Copy nifti image
void copy_nifti_image(shared_ptr<nifti_image> &output_image_sptr, const shared_ptr<nifti_image> &image_to_copy_sptr)
{
#ifndef NDEBUG
    cout << "\nPerforming hard copy of nifti image..." << flush;
#endif

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

#ifndef NDEBUG
    cout << "done.\n\n";
#endif
}

/// Do nifti image metadatas match?
bool do_nifti_image_metadata_match(const NiftiImage &im1, const NiftiImage &im2)
{
#ifndef NDEBUG
    cout << "\nChecking if metadata of two images match..." << flush;
#endif

    const shared_ptr<nifti_image> im1_sptr = im1.get_raw_nifti_sptr();
    const shared_ptr<nifti_image> im2_sptr = im2.get_raw_nifti_sptr();

    bool images_match =
            do_nifti_image_metadata_elements_match("analyze75_orient",im1_sptr->analyze75_orient,im2_sptr->analyze75_orient) &&
            do_nifti_image_metadata_elements_match("byteorder",       im1_sptr->byteorder,       im2_sptr->byteorder       ) &&
            do_nifti_image_metadata_elements_match("cal_max",         im1_sptr->cal_max,         im2_sptr->cal_max         ) &&
            do_nifti_image_metadata_elements_match("cal_min",         im1_sptr->cal_min,         im2_sptr->cal_min         ) &&
            do_nifti_image_metadata_elements_match("datatype",        im1_sptr->datatype,        im2_sptr->datatype        ) &&
            do_nifti_image_metadata_elements_match("du",              im1_sptr->du,              im2_sptr->du              ) &&
            do_nifti_image_metadata_elements_match("dv",              im1_sptr->dv,              im2_sptr->dv              ) &&
            do_nifti_image_metadata_elements_match("dw",              im1_sptr->dw,              im2_sptr->dw              ) &&
            do_nifti_image_metadata_elements_match("dx",              im1_sptr->dx,              im2_sptr->dx              ) &&
            do_nifti_image_metadata_elements_match("dy",              im1_sptr->dy,              im2_sptr->dy              ) &&
            do_nifti_image_metadata_elements_match("dz",              im1_sptr->dz,              im2_sptr->dz              ) &&
            do_nifti_image_metadata_elements_match("ext_list",        im1_sptr->ext_list,        im2_sptr->ext_list        ) &&
            do_nifti_image_metadata_elements_match("freq_dim",        im1_sptr->freq_dim,        im2_sptr->freq_dim        ) &&
            do_nifti_image_metadata_elements_match("iname_offset",    im1_sptr->iname_offset,    im2_sptr->iname_offset    ) &&
            do_nifti_image_metadata_elements_match("intent_code",     im1_sptr->intent_code,     im2_sptr->intent_code     ) &&
            do_nifti_image_metadata_elements_match("intent_p1",       im1_sptr->intent_p1,       im2_sptr->intent_p1       ) &&
            do_nifti_image_metadata_elements_match("intent_p2",       im1_sptr->intent_p2,       im2_sptr->intent_p2       ) &&
            do_nifti_image_metadata_elements_match("intent_p3",       im1_sptr->intent_p3,       im2_sptr->intent_p3       ) &&
            do_nifti_image_metadata_elements_match("nbyper",          im1_sptr->nbyper,          im2_sptr->nbyper          ) &&
            do_nifti_image_metadata_elements_match("ndim",            im1_sptr->ndim,            im2_sptr->ndim            ) &&
            do_nifti_image_metadata_elements_match("nifti_type",      im1_sptr->nifti_type,      im2_sptr->nifti_type      ) &&
            do_nifti_image_metadata_elements_match("nt",              im1_sptr->nt,              im2_sptr->nt              ) &&
            do_nifti_image_metadata_elements_match("nu",              im1_sptr->nu,              im2_sptr->nu              ) &&
            do_nifti_image_metadata_elements_match("num_ext",         im1_sptr->num_ext,         im2_sptr->num_ext         ) &&
            do_nifti_image_metadata_elements_match("nv",              im1_sptr->nv,              im2_sptr->nv              ) &&
            do_nifti_image_metadata_elements_match("nvox",            im1_sptr->nvox,            im2_sptr->nvox            ) &&
            do_nifti_image_metadata_elements_match("nw",              im1_sptr->nw,              im2_sptr->nw              ) &&
            do_nifti_image_metadata_elements_match("nx",              im1_sptr->nx,              im2_sptr->nx              ) &&
            do_nifti_image_metadata_elements_match("ny",              im1_sptr->ny,              im2_sptr->ny              ) &&
            do_nifti_image_metadata_elements_match("nz",              im1_sptr->nz,              im2_sptr->nz              ) &&
            do_nifti_image_metadata_elements_match("phase_dim",       im1_sptr->phase_dim,       im2_sptr->phase_dim       ) &&
            do_nifti_image_metadata_elements_match("qfac",            im1_sptr->qfac,            im2_sptr->qfac            ) &&
            do_nifti_image_metadata_elements_match("qform_code",      im1_sptr->qform_code,      im2_sptr->qform_code      ) &&
            do_nifti_image_metadata_elements_match("qoffset_x",       im1_sptr->qoffset_x,       im2_sptr->qoffset_x       ) &&
            do_nifti_image_metadata_elements_match("qoffset_y",       im1_sptr->qoffset_y,       im2_sptr->qoffset_y       ) &&
            do_nifti_image_metadata_elements_match("qoffset_z",       im1_sptr->qoffset_z,       im2_sptr->qoffset_z       ) &&
            do_nifti_image_metadata_elements_match("quatern_b",       im1_sptr->quatern_b,       im2_sptr->quatern_b       ) &&
            do_nifti_image_metadata_elements_match("quatern_c",       im1_sptr->quatern_c,       im2_sptr->quatern_c       ) &&
            do_nifti_image_metadata_elements_match("quatern_d",       im1_sptr->quatern_d,       im2_sptr->quatern_d       ) &&
            do_nifti_image_metadata_elements_match("scl_inter",       im1_sptr->scl_inter,       im2_sptr->scl_inter       ) &&
            do_nifti_image_metadata_elements_match("scl_slope",       im1_sptr->scl_slope,       im2_sptr->scl_slope       ) &&
            do_nifti_image_metadata_elements_match("sform_code",      im1_sptr->sform_code,      im2_sptr->sform_code      ) &&
            do_nifti_image_metadata_elements_match("slice_code",      im1_sptr->slice_code,      im2_sptr->slice_code      ) &&
            do_nifti_image_metadata_elements_match("slice_dim",       im1_sptr->slice_dim,       im2_sptr->slice_dim       ) &&
            do_nifti_image_metadata_elements_match("slice_duration",  im1_sptr->slice_duration,  im2_sptr->slice_duration  ) &&
            do_nifti_image_metadata_elements_match("slice_end",       im1_sptr->slice_end,       im2_sptr->slice_end       ) &&
            do_nifti_image_metadata_elements_match("slice_start",     im1_sptr->slice_start,     im2_sptr->slice_start     ) &&
            do_nifti_image_metadata_elements_match("swapsize",        im1_sptr->swapsize,        im2_sptr->swapsize        ) &&
            do_nifti_image_metadata_elements_match("time_units",      im1_sptr->time_units,      im2_sptr->time_units      ) &&
            do_nifti_image_metadata_elements_match("toffset",         im1_sptr->toffset,         im2_sptr->toffset         ) &&
            do_nifti_image_metadata_elements_match("xyz_units",       im1_sptr->xyz_units,       im2_sptr->xyz_units       ) &&
            do_nifti_image_metadata_elements_match("qto_ijk",         im1_sptr->qto_ijk,         im2_sptr->qto_ijk         ) &&
            do_nifti_image_metadata_elements_match("qto_xyz",         im1_sptr->qto_xyz,         im2_sptr->qto_xyz         ) &&
            do_nifti_image_metadata_elements_match("sto_ijk",         im1_sptr->sto_ijk,         im2_sptr->sto_ijk         ) &&
            do_nifti_image_metadata_elements_match("sto_xyz",         im1_sptr->sto_xyz,         im2_sptr->sto_xyz         );

    for (int i=0; i<8; i++) {
        if (!do_nifti_image_metadata_elements_match("dim["+to_string(i)+"]",    im1_sptr->dim[i],    im2_sptr->dim[i] ))   images_match = false;
        if (!do_nifti_image_metadata_elements_match("pixdim["+to_string(i)+"]", im1_sptr->pixdim[i], im2_sptr->pixdim[i])) images_match = false;
    }

#ifndef NDEBUG
    if (images_match) cout << "\tOK!\n";
#endif

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

/// Dump info of multiple nifti images
void dump_headers(const vector<NiftiImage> &ims)
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

    // Print min
    string min_header = "min: ";
    cout << "\t" << left << setw(19) << min_header;
    for(unsigned i=0; i<ims.size(); i++)
        cout << setw(19) << ims[i].get_min();

    // Print max
    cout << "\n\t" << left << setw(19) << "max: ";
    for(unsigned i=0; i<ims.size(); i++)
        cout << setw(19) << ims[i].get_max();

    cout << "\n\n";
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

template<typename newType>
void change_datatype(NiftiImage &im)
{
    if (im.get_raw_nifti_sptr()->datatype == DT_BINARY)   return SIRFRegMisc::change_datatype<newType,bool>              (im);
    if (im.get_raw_nifti_sptr()->datatype == DT_INT8)     return SIRFRegMisc::change_datatype<newType,signed char>       (im);
    if (im.get_raw_nifti_sptr()->datatype == DT_INT16)    return SIRFRegMisc::change_datatype<newType,signed short>      (im);
    if (im.get_raw_nifti_sptr()->datatype == DT_INT32)    return SIRFRegMisc::change_datatype<newType,signed int>        (im);
    if (im.get_raw_nifti_sptr()->datatype == DT_FLOAT32)  return SIRFRegMisc::change_datatype<newType,float>             (im);
    if (im.get_raw_nifti_sptr()->datatype == DT_FLOAT64)  return SIRFRegMisc::change_datatype<newType,double>            (im);
    if (im.get_raw_nifti_sptr()->datatype == DT_UINT8)    return SIRFRegMisc::change_datatype<newType,unsigned char>     (im);
    if (im.get_raw_nifti_sptr()->datatype == DT_UINT16)   return SIRFRegMisc::change_datatype<newType,unsigned short>    (im);
    if (im.get_raw_nifti_sptr()->datatype == DT_UINT32)   return SIRFRegMisc::change_datatype<newType,unsigned int>      (im);
    if (im.get_raw_nifti_sptr()->datatype == DT_INT64)    return SIRFRegMisc::change_datatype<newType,signed long long>  (im);
    if (im.get_raw_nifti_sptr()->datatype == DT_UINT64)   return SIRFRegMisc::change_datatype<newType,unsigned long long>(im);
    if (im.get_raw_nifti_sptr()->datatype == DT_FLOAT128) return SIRFRegMisc::change_datatype<newType,long double>       (im);

    stringstream ss;
    ss << "NiftImage::get_max not implemented for your data type: ";
    ss << nifti_datatype_string(im.get_raw_nifti_sptr()->datatype);
    ss << " (bytes per voxel: ";
    ss << im.get_raw_nifti_sptr()->nbyper << ").";
    throw runtime_error(ss.str());
}
template void change_datatype<bool>              (NiftiImage &im);
template void change_datatype<signed char>       (NiftiImage &im);
template void change_datatype<signed short>      (NiftiImage &im);
template void change_datatype<signed int>        (NiftiImage &im);
template void change_datatype<float>             (NiftiImage &im);
template void change_datatype<double>            (NiftiImage &im);
template void change_datatype<unsigned char>     (NiftiImage &im);
template void change_datatype<unsigned short>    (NiftiImage &im);
template void change_datatype<unsigned int>      (NiftiImage &im);
template void change_datatype<signed long long>  (NiftiImage &im);
template void change_datatype<unsigned long long>(NiftiImage &im);
template void change_datatype<long double>       (NiftiImage &im);

}
