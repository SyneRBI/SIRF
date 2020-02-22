/*
CCP PETMR Synergistic Image Reconstruction Framework (SIRF)
Copyright 2018-2020 University College London

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
\brief Perform registration with any types of SIRF image and any SIRF registration algorithm.

\author Richard Brown
\author CCP PETMR
*/

#include "sirf/Reg/NiftyAladinSym.h"
#include "sirf/Reg/NiftyF3dSym.h"
#include "sirf/Reg/AffineTransformation.h"
#include "sirf/Reg/NiftiImageData3D.h"
#include "sirf/Gadgetron/gadgetron_data_containers.h"
#include "sirf/STIR/stir_data_containers.h"
#include <boost/filesystem.hpp>
#ifdef SIRF_SPM
#include "sirf/Reg/SPMRegistration.h"
#endif

using namespace sirf;

enum Algorithm {
    Aladin,
    F3d,
    SPM
};

static std::shared_ptr<const ImageData> image_as_sptr(const std::string &filename, const std::string &engine)
{
    if      (strcmp(engine.c_str(), "Nifti") == 0)
        return std::make_shared<const NiftiImageData3D<float> >(filename);
    else if (strcmp(engine.c_str(), "STIR") == 0)
        return std::make_shared<const STIRImageData>(filename);
    else if (strcmp(engine.c_str(), "Gadgetron") == 0) {
        std::shared_ptr<GadgetronImageData> sptr_img(new GadgetronImagesVector);
		sptr_img->read(filename);
        return std::move(sptr_img);
    }
    else
        throw std::runtime_error("sirf_registration: unknown image engine - " + engine + ".\n");
}
//reg,algo,is_affine_or_rigid,algo_str
static void algo_as_sptr(std::shared_ptr<Registration<float> > &algo_sptr, Algorithm &algo, const std::string &algorithm)
{
    if (strcmp(algorithm.c_str(), "aladin") == 0) {
        algo_sptr = std::make_shared<NiftyAladinSym<float> >();
        algo = Aladin;
    }
    else if (strcmp(algorithm.c_str(), "f3d") == 0) {
        algo_sptr = std::make_shared<NiftyF3dSym<float> >();
        algo = F3d;
    }
    else if (strcmp(algorithm.c_str(), "spm") == 0) {
#ifdef SIRF_SPM
        algo_sptr = std::make_shared<SPMRegistration<float> >();
        algo = SPM;
#else
        throw std::runtime_error("sirf_registration: SIRF not built with spm\n");
#endif
    }
    else
        throw std::runtime_error("sirf_registration: unknown algorithm - " + algorithm + ".\n");

    std::cout << "\nUsing " << algorithm << " registration algorithm...\n";
}

/// print usage
void print_usage()
{
    std::cout << "\n*** sirf_registration usage ***\n";

    // Required flags
    std::cout << "\n  Required flags:\n";
    std::cout << "    --algo <algo>:\t\tregistration algorithm (aladin/f3d/spm)\n";
    std::cout << "    --ref <fname> <eng>:\treference image (eng: Reg|STIR|Gadgetron)\n";
    std::cout << "    --flo <fname> <eng>:\tfloating image (eng: Reg|STIR|Gadgetron). (Can be used mulitple times for spm.)\n";

    // Optional flags
    std::cout << "\n  Optional flags:\n";
    std::cout << "    --warped_prefix <fname>:\twarped image filename\n";
    std::cout << "    --disp_fwd_prefix <fname>:\tforward displacement field image\n";
    std::cout << "    --disp_inv_prefix <fname>:\tinverse displacement field image\n";
    std::cout << "    --def_fwd_prefix <fname>:\tforward deformation field image\n";
    std::cout << "    --def_inv_prefix <fname>:\tinverse deformation field image\n";

    // Optional rigid/affine flags
    std::cout << "\n  Optional flags for rigid/affine algorithms (aladin/spm):\n";
    std::cout << "    --TM_fwd_prefix <fname>:\tforward transformation matrix\n";
    std::cout << "    --TM_inv_prefix <fname>:\tinverse transformation matrix\n";

    // Optional NiftyReg flags
    std::cout << "\n  Optional flags for NiftyReg (aladin/f3d):\n";
    std::cout << "    --rmask <fname> <eng>:\tmask of reference image (eng: Reg|STIR|Gadgetron)\n";
    std::cout << "    --fmask <fname> <eng>:\tmask of floating image (eng: Reg|STIR|Gadgetron)\n";
    std::cout << "    --print:\t\t\tprint all possible wrapped parameters and exit.\n";
    std::cout << "    --par_file <fname>:\t\tset parameter file\n";
    std::cout << "    --par \"<string>\":\t\tset wrapped parameter. Some examples (and note quotation marks):\n";
    std::cout << "            --par \"SetPerformRigid 1\"\n";
    std::cout << "            --par \"SetPerformAffine 0\"\n";
    std::cout << "            --par \"SetInterpolationToCubic\"\n";
    std::cout << "            --par \"SetFloatingThresholdUp 1 2\"\n";

    // Optional SPM flags
    std::cout << "\n  Optional flags for spm:\n";
    std::cout << "    --working_folder <fname>:\tfolder in which to save temporary files (default: cwd/spm_working_folder)\n";
    std::cout << "    --overwrite <bool>:\t\tshould I overwrite files if already present? (default: 1)\n";
    std::cout << "    --delete <bool>:\t\tshould I delete temporary files? (default: 1)\n";
}

/// main
int main(int argc, char* argv[])
{
    try {

        // Ignore executable name
        argc--; argv++;

        if (argc < 2) {
            print_usage();
            exit(EXIT_FAILURE);
        }

        // Check for help
        for (unsigned i=0; i<unsigned(argc); ++i) {
            if (strcmp(argv[i], "-h") == 0 || strcmp(argv[i], "--help") == 0) {
                print_usage();
                exit(EXIT_SUCCESS);
            }
        }

        // Parse
        std::string algo_str, ref_str, ref_eng_str,
                warped_str, disp_fwd_str, def_fwd_str, disp_inv_str, def_inv_str,
                TM_fwd_str, TM_inv_str,
                rmask_str, fmask_str, rmask_eng_str, fmask_eng_str,par_file_str,
                working_folder_str;
        std::vector<std::string> pars;
        std::vector<std::pair<std::string,std::string> > flo_strs;
        int overwrite(-1), delete_temp_file(-1);
        bool print(false);

        while (argc>0) {
            // Algo
            if (strcmp(argv[0],"--algo")==0) {
                if (argc<2) throw std::runtime_error("--algo requires an argument");
                algo_str = argv[1];
                argc-=2; argv+=2;
            }
            // ref
            else if (strcmp(argv[0],"--ref")==0) {
                if (argc<3) throw std::runtime_error("--ref requires two arguments");
                ref_str = argv[1];
                ref_eng_str = argv[2];
                argc-=3; argv+=3;
            }
            // flo
            else if (strcmp(argv[0],"--flo")==0) {
                if (argc<3) throw std::runtime_error("--flo requires two arguments");
                flo_strs.push_back(std::make_pair(argv[1],argv[2]));
                argc-=3; argv+=3;
            }
            // warped
            else if (strcmp(argv[0],"--warped_prefix")==0) {
                if (argc<2) throw std::runtime_error("--warped requires an argument");
                warped_str = argv[1];
                argc-=2; argv+=2;
            }
            // disp_fwd
            else if (strcmp(argv[0],"--disp_fwd_prefix")==0) {
                if (argc<2) throw std::runtime_error("--disp_fwd requires an argument");
                disp_fwd_str = argv[1];
                argc-=2; argv+=2;
            }
            // disp_inv
            else if (strcmp(argv[0],"--disp_inv_prefix")==0) {
                if (argc<2) throw std::runtime_error("--disp_inv requires an argument");
                disp_inv_str = argv[1];
                argc-=2; argv+=2;
            }
            // def_fwd
            else if (strcmp(argv[0],"--def_fwd_prefix")==0) {
                if (argc<2) throw std::runtime_error("--def_fwd requires an argument");
                def_fwd_str = argv[1];
                argc-=2; argv+=2;
            }
            // def_inv
            else if (strcmp(argv[0],"--def_inv_prefix")==0) {
                if (argc<2) throw std::runtime_error("--def_inv requires an argument");
                def_inv_str = argv[1];
                argc-=2; argv+=2;
            }
            // TM_fwd
            else if (strcmp(argv[0],"--TM_fwd_prefix")==0) {
                if (argc<2) throw std::runtime_error("--TM_fwd requires an argument");
                TM_fwd_str = argv[1];
                argc-=2; argv+=2;
            }
            // TM_inv
            else if (strcmp(argv[0],"--TM_inv_prefix")==0) {
                if (argc<2) throw std::runtime_error("--TM_inv requires an argument");
                TM_inv_str = argv[1];
                argc-=2; argv+=2;
            }
            // rmask
            else if (strcmp(argv[0],"--rmask")==0) {
                if (argc<3) throw std::runtime_error("--rmask requires two arguments");
                rmask_str = argv[1];
                rmask_eng_str = argv[2];
                argc-=3; argv+=3;
            }
            // fmask
            else if (strcmp(argv[0],"--fmask")==0) {
                if (argc<3) throw std::runtime_error("--fmask requires two arguments");
                fmask_str = argv[1];
                fmask_eng_str = argv[2];
                argc-=3; argv+=3;
            }
            // print
            else if (strcmp(argv[0],"--print")==0) {
                print = true;
                argc--; argv++;
            }
            // par_file
            else if (strcmp(argv[0],"--par_file")==0) {
                if (argc<2) throw std::runtime_error("--par_file requires an argument");
                par_file_str = argv[1];
                argc-=2; argv+=2;
            }
            // par
            else if (strcmp(argv[0],"--par")==0) {
                if (argc<2) throw std::runtime_error("--par_file requires an argument");
                pars.push_back(argv[1]);
                argc-=2; argv+=2;
            }
            // working_folder
            else if (strcmp(argv[0],"--working_folder")==0) {
                if (argc<2) throw std::runtime_error("--working_folder requires an argument");
                working_folder_str = argv[1];
                argc-=2; argv+=2;
            }
            // overwrite
            else if (strcmp(argv[0],"--overwrite")==0) {
                if (argc<2) throw std::runtime_error("--overwrite requires an argument");
                overwrite = std::stoi(argv[1]);
                if (!(overwrite==0 || overwrite==1))
                    throw std::runtime_error("--overwrite should be 0 or 1");
                argc-=2; argv+=2;
            }
            // delete
            else if (strcmp(argv[0],"--delete")==0) {
                if (argc<2) throw std::runtime_error("--delete requires an argument");
                delete_temp_file = std::stoi(argv[1]);
                if (!(delete_temp_file==0 || delete_temp_file==1))
                    throw std::runtime_error("--delete_temp_file should be 0 or 1");
                argc-=2; argv+=2;
            }
            else
                throw std::runtime_error("Unkown option " + std::string(argv[0]));
        }

        // algo
        if (algo_str.empty()) throw std::runtime_error("--algo not set");
        std::shared_ptr<Registration<float> > reg;
        Algorithm algo;
        algo_as_sptr(reg,algo,algo_str);

        // Ref
        if (ref_str.empty()) throw std::runtime_error("--ref not set");
        reg->set_reference_image(image_as_sptr(ref_str,ref_eng_str));
        // Flo
        if (flo_strs.empty()) throw std::runtime_error("--flo not set");
        for (unsigned i=0; i<flo_strs.size(); ++i)
            reg->add_floating_image(image_as_sptr(flo_strs.at(i).first,flo_strs.at(i).second));

        // rmask
        if (!rmask_str.empty()) {
            if (algo == SPM) throw std::runtime_error("--rmask not available for spm");
            std::dynamic_pointer_cast<NiftyRegistration<float> >(reg)->set_reference_mask(image_as_sptr(rmask_str,rmask_eng_str));
        }
        // fmask
        if (!fmask_str.empty()) {
            if (algo == SPM) throw std::runtime_error("--fmask not available for spm");
            std::dynamic_pointer_cast<NiftyRegistration<float> >(reg)->set_reference_mask(image_as_sptr(fmask_str,fmask_eng_str));
        }
        // print
        if (print) {
            if (algo == Aladin)
                std::dynamic_pointer_cast<NiftyAladinSym<float> >(reg)->print_all_wrapped_methods();
            else if (algo == F3d)
                std::dynamic_pointer_cast<NiftyF3dSym<float> >(reg)->print_all_wrapped_methods();
            else throw std::runtime_error("--print not available for spm");
            exit(EXIT_SUCCESS);
        }
        // par_file
        if (!par_file_str.empty()) {
            if (algo == SPM) throw std::runtime_error("--par_file not available for spm");
            std::dynamic_pointer_cast<NiftyRegistration<float> >(reg)->set_parameter_file(par_file_str);
        }
        // pars
        if (pars.size()>0) {
            if (algo == SPM) throw std::runtime_error("--par not available for spm");

            for (unsigned i=0; i<pars.size(); ++i) {
                std::istringstream iss(pars[i]);
                std::vector<std::string> tokens;
                std::copy(std::istream_iterator<std::string>(iss),
                     std::istream_iterator<std::string>(),
                     std::back_inserter(tokens));
                if (tokens.size() == 1)
                    std::dynamic_pointer_cast<NiftyRegistration<float> >(reg)->set_parameter(tokens[0]);
                else if (tokens.size() == 2)
                    std::dynamic_pointer_cast<NiftyRegistration<float> >(reg)->set_parameter(tokens[0],tokens[1]);
                else if (tokens.size() == 3)
                    std::dynamic_pointer_cast<NiftyRegistration<float> >(reg)->set_parameter(tokens[0],tokens[1],tokens[2]);
                else
                    throw std::runtime_error("Max number of NiftyReg args is 2.");
            }
        }
        // working folder
        if (algo == SPM && working_folder_str.empty())
            working_folder_str = boost::filesystem::current_path().append("spm_working_folder").string();
        if (!working_folder_str.empty()) {
            if (algo != SPM) throw std::runtime_error("--working_folder only available for spm");
#ifdef SIRF_SPM
            std::dynamic_pointer_cast<SPMRegistration<float> >(reg)->set_working_folder(working_folder_str);
#endif
        }
        // overwrite
        if (algo == SPM && overwrite == -1)
            overwrite = 1;
        if (overwrite!=-1) {
            if (algo != SPM) throw std::runtime_error("--overwrite only available for spm");
#ifdef SIRF_SPM
            std::dynamic_pointer_cast<SPMRegistration<float> >(reg)->set_working_folder_file_overwrite(bool(overwrite));
#endif
        }
        // delete temp files
        if (algo == SPM && delete_temp_file == -1)
            delete_temp_file = 1;
        if (delete_temp_file!=-1) {
            if (algo != SPM) throw std::runtime_error("--delete only available for spm");
#ifdef SIRF_SPM
            std::dynamic_pointer_cast<SPMRegistration<float> >(reg)->set_delete_temp_files(bool(delete_temp_file));
#endif
        }


        // ------------------------------------------------ //
        // Process
        // ------------------------------------------------ //

        reg->process();


        // ------------------------------------------------ //
        // Save results
        // ------------------------------------------------ //

        for (unsigned i=0; i<flo_strs.size(); ++i) {
            // warped
            if (!warped_str.empty())
                reg->get_output_sptr(i)->write(warped_str + std::to_string(i));

            // disp fwd
            if (!disp_fwd_str.empty())
                reg->get_displacement_field_forward_sptr(i)->write(disp_fwd_str + std::to_string(i));
            // disp inv
            if (!disp_inv_str.empty())
                reg->get_displacement_field_forward_sptr(i)->write(disp_inv_str + std::to_string(i));

            // def fwd
            if (!def_fwd_str.empty())
                reg->get_deformation_field_forward_sptr(i)->write(def_fwd_str + std::to_string(i));
            // def inv
            if (!def_inv_str.empty())
                reg->get_deformation_field_forward_sptr(i)->write(def_inv_str + std::to_string(i));

            // TM_fwd
            if (!TM_fwd_str.empty()) {
                if (algo == Aladin)
                    std::dynamic_pointer_cast<NiftyAladinSym<float> >(reg)->get_transformation_matrix_forward_sptr()->write(TM_fwd_str + std::to_string(i));
#ifdef SIRF_SPM
                else if (algo == SPM)
                    std::dynamic_pointer_cast<SPMRegistration<float> >(reg)->get_transformation_matrix_forward_sptr(i)->write(TM_fwd_str + std::to_string(i));
#endif
                else throw std::runtime_error("--TM_fwd only available for rigid/affine");
            }
            // TM_inv
            if (!TM_inv_str.empty()) {
                if (algo == Aladin)
                    std::dynamic_pointer_cast<NiftyAladinSym<float> >(reg)->get_transformation_matrix_inverse_sptr()->write(TM_inv_str + std::to_string(i));
#ifdef SIRF_SPM
                else if (algo == SPM)
                    std::dynamic_pointer_cast<SPMRegistration<float> >(reg)->get_transformation_matrix_inverse_sptr(i)->write(TM_inv_str + std::to_string(i));
#endif
                else throw std::runtime_error("--TM_inv only available for rigid/affine");
            }
        }
    }

    // If there was an error
    catch(const std::exception &error) {
        std::cerr << "\nError encountered:\n\t" << error.what() << "\n\n";
        return EXIT_FAILURE;
    }

    return EXIT_SUCCESS;
}
