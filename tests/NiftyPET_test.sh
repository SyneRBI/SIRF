#!/bin/bash

function print_usage() {
	echo "Usage: $0 [-h|--help] RAW_DATA_FOLDER TSTART TSTOP"
}

# Print full help
if [[ $1 == "-h" || $1 == "--help" ]]; then
	print_usage
	echo -e "\nThis bash script compares the STIR-wrapped NiftyPET functionality to NiftyPET and, where possible, STIR.\n"
	echo -e "### Test content\n"
	echo -e "From raw dicom data (RAW_DATA_FOLDER) and a given time window (TSTART and TSTOP), prompt sinograms are extracted with NiftyPET, STIR-NiftyPET and STIR. Randoms are estimated, and for NiftyPET and STIR-NiftyPET, a norm sinogram is extracted.\n"
	echo -e "Forward and back projections are performed with NiftyPET and STIR-NiftyPET, first with no corrections, then with norm, randoms and attenuation taken into account.\n"
	echo -e "Lastly, a single MLEM iteration is performed with NiftyPET and STIR-NiftyPET.\n"
	echo -e "At each step, the results are compared to ensure consistency.\n"
	echo -e "At the time of writing, NiftyPET requires python2, so two different python executables should be set as environmental variables - NIFTYPET_PYTHON_EXECUTABLE and SIRF_PYTHON_EXECUTABLE.\n"
	echo -e "### Data\n"
	echo -e "An example of test data can be found here - [https://doi.org/10.5281/zenodo.1472951](https://doi.org/10.5281/zenodo.1472951).\n"
	echo -e "Since the data needs to be read by NiftyPET, it should be in the raw form of .dcm/.bf, and not .l.hdr/.l.\n"
	exit 0
fi

# Check input arguments
if [[ "$#" -ne 3 ]]; then
	print_usage
	exit 1
fi
NTYP_RAW_DATA_FOLDER=$1
SINO_TSTART=$2
SINO_TSTOP=$3

# Check for python versions
if [[ $NIFTYPET_PYTHON_EXECUTABLE == "" ]]; then
	echo "Please set environmental variable: NIFTYPET_PYTHON_EXECUTABLE"
	exit 1
elif [[ $SIRF_PYTHON_EXECUTABLE == "" ]]; then
	echo "Please set environmental variable: SIRF_PYTHON_EXECUTABLE"
	exit 1
fi
NP_PY_EXE=$NIFTYPET_PYTHON_EXECUTABLE
SIRF_PY_EXE=$SIRF_PYTHON_EXECUTABLE

fail_count=0

function check_result() {
	res=$2
	if [ $res -eq 99 ]; then 
		fail_count=$((fail_count+1)); echo "Test failed. Fail count=$fail_count"
	elif [ $res -ne 0 ]; then 
		echo "Exiting at line $1"; exit $res
	fi
}

########################################################################################
#
#    COMPARE SINOGAMS
#
########################################################################################

compare_sinos_stir_space() {
	echo -e "\nComparing sinograms in STIR space..."
	sino_stir=$1
	sino_np=$2
	sino_np_2_stir=STIR_${sino_np%.*}.hs
	conv_niftypet_stir $sino_np_2_stir $sino_np sinogram toSTIR
	res=$?; if [ $res -ne 0 ]; then echo "Exiting at line $LINENO"; exit $res; fi
	$SIRF_PY_EXE -c \
"import sys
import numpy as np
f_sino1=sys.argv[1]
f_sino2=sys.argv[2]
from sirf.STIR import *
sino1 = AcquisitionData(f_sino1)
sino2 = AcquisitionData(f_sino2)
diff = sino1 - sino2
print('max sino STIR = ' + str(np.abs(sino1.as_array()).max()))
print('max sino NP->STIR = ' + str(np.abs(sino2.as_array()).max()))
print('max sino difference = ' + str(np.abs(diff.as_array()).max()))
if np.abs(diff.as_array()).max() > 1:
	sys.exit(99)" \
		$sino_stir $sino_np_2_stir
	check_result $LINENO $?
}

compare_sinos_np_space() {
	echo -e "\nComparing sinograms in NiftyPET space..."
	sino_np=$1
	sino_stir=$2
	sino_stir_2_np=NP_${sino_stir%.*}.dat
	conv_niftypet_stir $sino_stir_2_np $sino_stir sinogram toNP
	res=$?; if [ $res -ne 0 ]; then echo "Exiting at line $LINENO"; exit $res; fi
	$SIRF_PY_EXE -c \
"import sys
import numpy as np
f_sino1=sys.argv[1]
f_sino2=sys.argv[2]
sino1 = np.fromfile(f_sino1, dtype='float32')
sino2 = np.fromfile(f_sino2, dtype='float32')
diff = sino1 - sino2
norm = np.linalg.norm(diff.flatten()) / diff.size
print('max sino STIR->NP = ' + str(np.abs(sino1).max()))
print('max sino NP = ' + str(np.abs(sino2).max()))
print('max sino difference = ' + str(np.abs(diff).max()))
print('difference norm = ' + str(norm))
if norm > 1:
	sys.exit(99)" \
		$sino_stir_2_np $sino_np
	check_result $LINENO $?
}

########################################################################################
#
#    COMPARE IMAGES
#
########################################################################################

compare_ims_stir_space() {
	echo -e "\nComparing images in STIR space..."
	im_stir=$1
	im_np=$2
	im_np_2_stir_base=STIR_${im_np%.*}
	im_np_2_stir=STIR_${im_np%.*}.nii
	conv_niftypet_stir $im_np_2_stir_base $im_np image toSTIR && \
		sirf_convert_image_type $im_np_2_stir nii ${im_np_2_stir_base}.hv STIR && \
		rm ${im_np_2_stir_base}.*v
	res=$?; if [ $res -ne 0 ]; then echo "Exiting at line $LINENO"; exit $res; fi
	$SIRF_PY_EXE -c \
"import sys
import numpy as np
f_im1=sys.argv[1]
f_im2=sys.argv[2]
from sirf.STIR import *
im1 = ImageData(f_im1)
im2 = ImageData(f_im2)
diff = im1 - im2
norm = diff.norm() / diff.as_array().size
print('max im STIR = ' + str(np.abs(im1.as_array()).max()))
print('max im NP->STIR = ' + str(np.abs(im2.as_array()).max()))
print('max im difference = ' + str(np.abs(diff.as_array()).max()))
print('difference norm = ' + str(norm))
if norm > 1:
	sys.exit(99)" \
		$im_stir $im_np_2_stir
	check_result $LINENO $?
}

compare_ims_np_space() {
	echo -e "\nComparing images in NiftyPET space..."
	im_np=$1
	im_stir=$2
	im_stir_2_np=NP_${im_stir%.*}.dat
	conv_niftypet_stir $im_stir_2_np $im_stir image toNP
	res=$?; if [ $res -ne 0 ]; then echo "Exiting at line $LINENO"; exit $res; fi
	$SIRF_PY_EXE -c \
"import sys
import numpy as np
f_im1=sys.argv[1]
f_im2=sys.argv[2]
im1 = np.fromfile(f_im1, dtype='float32')
im2 = np.fromfile(f_im2, dtype='float32')
diff = im1 - im2
norm = np.linalg.norm(diff.flatten()) / diff.size
print('max im STIR->NP = ' + str(np.abs(im1).max()))
print('max im NP = ' + str(np.abs(im1).max()))
print('max im difference = ' + str(np.abs(diff).max()))
print('difference norm = ' + str(norm))
if norm > 1:
	sys.exit(99)" \
		$im_stir_2_np $im_np
	check_result $LINENO $?
}



########################################################################################
#
#    PROJECTION WITH NP
#
########################################################################################

project_np() {
	$NP_PY_EXE -c \
"import sys
from niftypet import nipet
import numpy as np
mMRpars = nipet.get_mmrparams()
f_out=sys.argv[1]
f_in=sys.argv[2]
type=sys.argv[3]
f_norm = sys.argv[4]
f_rands = sys.argv[5]
f_attn = sys.argv[6]

def read_sino(f_name):
	if f_name=='': return None
	sino = np.fromfile(f_name, dtype='float32')
	return np.reshape(sino, (837, 252, 344))

additive = read_sino(f_rands)
norm = read_sino(f_norm)
attn = read_sino(f_attn)
multiplicative = None
if norm is not None and attn is not None:
	multiplicative = norm * attn
elif norm is not None:
	multiplicative = norm
elif attn is not None:
	multiplicative = attn

def fwd(im):
	out = nipet.frwd_prj(im, mMRpars)
	if multiplicative is not None: out *= multiplicative
	if additive is not None: out += additive
	return out
def bck(sino):
	sino_to_proj = np.copy(sino)
	if multiplicative is not None: sino_to_proj *= multiplicative
	out = nipet.back_prj(sino_to_proj, mMRpars)
	return out
def non_zero(array):
	array[array==0] = array.max()*1e-8
	return array

if type=='fwd':
	im = np.fromfile(f_in, dtype='float32')
	im = np.reshape(im, (127, 320, 320))
	im = np.transpose(im, (1, 2, 0))
	out = fwd(im)
else:
	sino = read_sino(f_in)
	if type=='bck':
		out = bck(sino)
	else:
		im = np.ones((127,344,344), dtype=np.float32)
		bck_sino_ones = non_zero(bck(np.ones_like(sino)))
		if multiplicative is not None: sino *= multiplicative
		out = (im/bck_sino_ones) * bck(sino/non_zero(fwd(im)))
	out = nipet.img.mmrimg.convert2dev(out, mMRpars['Cnt'])
	out = np.transpose(out, (2, 0, 1))

out.astype('float32').tofile(f_out)" \
		"$@"
	res=$?; if [ $res -ne 0 ]; then echo "Exiting at line $LINENO"; exit $res; fi
}

########################################################################################
#
#    PROJECTION WITH STIR
#
########################################################################################

function project_stir() {
	echo -e "\nMLEM recon with STIR's NiftyPET wrapper"
	$SIRF_PY_EXE -c \
"import sys
import os
from sirf.STIR import *
import sirf.Reg
set_verbosity(0)
f_out = sys.argv[1]
f_in = sys.argv[2]
type = sys.argv[3]
f_norm = sys.argv[4]
f_rands = sys.argv[5]
f_attn = sys.argv[6]

acq_model = AcquisitionModelUsingNiftyPET()

asm_norm = None
asm_attn = None
if f_norm != '':
	asm_norm = AcquisitionSensitivityModel(AcquisitionData(f_norm))
if f_attn != '':
	asm_attn = AcquisitionSensitivityModel(AcquisitionData(f_attn))
asm = asm_norm or asm_attn
if asm_norm and asm_attn:
	asm = AcquisitionSensitivityModel(asm_norm, asm_attn)
if asm is not None:
	acq_model.set_acquisition_sensitivity(asm)
if f_rands != '':
	acq_model.set_background_term(AcquisitionData(f_rands))

if type=='fwd':
	im = ImageData(f_in)
	data_path = examples_data_path('PET')
	raw_data_file = existing_filepath(data_path, 'mMR/mMR_template_span11.hs')
	template_acq_data = AcquisitionData(raw_data_file)
	acq_model.set_up(template_acq_data, im)
	out = acq_model.forward(im)
else:
	sino = AcquisitionData(f_in)
	im = ImageData()
	im.initialise(dim=(127,320,320), vsize=(2.03125,2.08626,2.08626))
	im.fill(1.0)
	acq_model.set_up(sino, im)
	if type=='bck':
		out = acq_model.backward(sino)
	else:
		obj_fun = make_Poisson_loglikelihood(sino)
		obj_fun.set_acquisition_model(acq_model)
		recon = OSMAPOSLReconstructor()
		recon.set_objective_function(obj_fun)
		recon.set_num_subsets(1)
		recon.set_num_subiterations(1)
		recon.set_input(sino)
		recon.set_up(im)
		recon.set_current_estimate(im)
		recon.process()
		out = recon.get_output()
	out = sirf.Reg.ImageData(out)
out.write(f_out)" \
		"$@"
	res=$?; if [ $res -ne 0 ]; then echo "Exiting at line $LINENO"; exit $res; fi
}



echo "########################################################################################"
echo "#                                                                                      #"
echo "#    GETTING FILENAMES...                                                              #"
echo "#                                                                                      #"
echo "########################################################################################"

$NP_PY_EXE -c \
"import sys
from niftypet import nipet
mMRpars = nipet.get_mmrparams()
datain = nipet.classify_input(sys.argv[1], mMRpars)
if not all(k in datain for k in ('lm_dcm','lm_bf','nrm_dcm')):
	raise AssertionError('Missing some input data. Example data: https://doi.org/10.5281/zenodo.1472951.')
f = open('LM_DCM.txt','w+')
f.write('LM_DCM=' + datain['lm_dcm'] + '\n')
f.write('LM_BF=' + datain['lm_bf'] + '\n')
f.write('NORM_DCM=' + datain['nrm_dcm'] + '\n')
f.write('NORM_BF=' + datain['nrm_bf'] + '\n')
f.close()" \
	$NTYP_RAW_DATA_FOLDER
res=$?; if [ $res -ne 0 ]; then echo "Exiting at line $LINENO"; exit $res; fi
source LM_DCM.txt
rm LM_DCM.txt
if [ ! -f STIR_lm.l.hdr ]; then
	echo "Converting DICOM listmode to interfile..."
	nm_extract -i $LM_DCM -o . -p STIR_lm && \
		sed -i "s#STIR_lm.l#${LM_BF}#g" "STIR_lm.l.hdr" && \
		rm STIR_lm.l
	res=$?; if [ $res -ne 0 ]; then echo "Exiting at line $LINENO"; exit $res; fi
fi
if [ ! -f STIR_norm.n.hdr ]; then
	echo "Converting DICOM listmode to interfile..."
	nm_extract -i $NORM_DCM -o . -p STIR_norm
	res=$?; if [ $res -ne 0 ]; then echo "Exiting at line $LINENO"; exit $res; fi
fi




echo "########################################################################################"
echo "#                                                                                      #"
echo "#    LISTMODE EXTRACTION                                                               #"
echo "#                                                                                      #"
echo "########################################################################################"

# Extract NP sinogram
if [ ! -f NP_sino.dat ]; then
	echo "Extracting sinogram with NiftyPET..."
	$NP_PY_EXE -c \
"import sys
import shutil
folderin=sys.argv[1]
tstart=int(sys.argv[2])
tstop=int(sys.argv[3])
from niftypet import nipet
mMRpars = nipet.get_mmrparams()
datain = nipet.classify_input(folderin, mMRpars)
hst = nipet.mmrhist(datain, mMRpars, t0=tstart, t1=tstop)
hst['psino'].astype('float32').tofile('NP_sino.dat')
rands = nipet.randoms(hst, mMRpars)[0]
rands.astype('float32').tofile('NP_rands.dat')
norm = nipet.mmrnorm.get_norm_sino(datain, mMRpars, hst)
norm.astype('float32').tofile('NP_norm.dat')
mu_o = nipet.obj_mumap(datain, mMRpars, outpath='.', store=False)['im']
attn = nipet.frwd_prj(mu_o, mMRpars, attenuation=True)
attn.astype('float32').tofile('NP_attn.dat')
shutil.rmtree('mumap-obj')
" \
	$NTYP_RAW_DATA_FOLDER $SINO_TSTART $SINO_TSTOP
	res=$?; if [ $res -ne 0 ]; then echo "Exiting at line $LINENO"; exit $res; fi
	conv_niftypet_stir NP_attn_as_STIR NP_attn.dat sinogram toSTIR 
	res=$?; if [ $res -ne 0 ]; then echo "Exiting at line $LINENO"; exit $res; fi
fi

# Extract STIR sinogram and randoms
if [ ! -f STIR_sino.hs ]; then
	echo "Extracting sinogram with STIR..."
	$SIRF_PY_EXE -c \
"import sys
list_file=sys.argv[1]
tstart=int(sys.argv[2])
tstop=int(sys.argv[3])
from sirf.STIR import *
lm2sino = ListmodeToSinograms()
lm2sino.set_input(list_file)
lm2sino.set_output_prefix('STIR_sino')
tmpl_file = examples_data_path('PET') + '/mMR/mMR_template_span11.hs'
lm2sino.set_template(tmpl_file)
lm2sino.set_time_interval(tstart, tstop)
lm2sino.set_up()
lm2sino.process()
sino = lm2sino.get_output()
sino.write('STIR_sino')
randoms = lm2sino.estimate_randoms()
randoms.write('STIR_rands')" \
	STIR_lm.l.hdr $SINO_TSTART $SINO_TSTOP
	res=$?; if [ $res -ne 0 ]; then echo "Exiting at line $LINENO"; exit $res; fi
	rm *f1g1d0b0* STIR_lm.l*
fi

# Extract with STIR's NiftyPET wrapper
if [ ! -f STIR_sino2.hs ]; then
	echo "Extracting sinogram with STIR's NiftyPET wrapper..."
	lm_to_projdata_niftypet $LM_BF $SINO_TSTART $SINO_TSTOP -N $NORM_BF -p STIR_sino2 -r STIR_rands2 -n STIR_norm2
	res=$?; if [ $res -ne 0 ]; then echo "Exiting at line $LINENO"; exit $res; fi
fi



echo "########################################################################################"
echo "#                                                                                      #"
echo "#    COMPARE LISTMODE EXTRACTION AND SINOGRAM CONVERSION                               #"
echo "#                                                                                      #"
echo "########################################################################################"

echo -e "\nComparing listmode extraction. STIR vs. NP in STIR space..."
compare_sinos_stir_space STIR_sino.hs NP_sino.dat
echo -e "\nComparing listmode extraction. NP vs. STIR in NP space..."
compare_sinos_np_space NP_sino.dat STIR_sino.hs
echo -e "\nComparing listmode extraction. NP vs. NP-STIR in NP space..."
compare_sinos_np_space NP_sino.dat STIR_sino2.hs




echo "########################################################################################"
echo "#                                                                                      #"
echo "#    RANDOMS COMPARISON                                                                #"
echo "#                                                                                      #"
echo "########################################################################################"

compare_sinos_np_space NP_rands.dat STIR_rands.hs
compare_sinos_np_space NP_rands.dat STIR_rands2.hs




echo "########################################################################################"
echo "#                                                                                      #"
echo "#    NORM COMPARISON                                                                   #"
echo "#                                                                                      #"
echo "########################################################################################"

compare_sinos_np_space NP_norm.dat STIR_norm2.hs




echo "########################################################################################"
echo "#                                                                                      #"
echo "#    COMPARE BACK PROJECTION (NO CORRECTIONS)                                          #"
echo "#                                                                                      #"
echo "########################################################################################"

conv_niftypet_stir NP_sino_as_STIR NP_sino.dat sinogram toSTIR

project_np NP_im.dat NP_sino.dat bck "" "" ""
project_stir STIR_im NP_sino_as_STIR.hs bck "" "" ""
compare_ims_stir_space STIR_im.nii NP_im.dat
compare_ims_np_space NP_im.dat STIR_im.nii




echo "########################################################################################"
echo "#                                                                                      #"
echo "#    COMPARE BACK PROJECTION WITH NORM, RANDOMS AND ATTENUATION                        #"
echo "#                                                                                      #"
echo "########################################################################################"

project_np NP_im_norm_rands_attn.dat NP_sino.dat bck NP_norm.dat NP_rands.dat NP_attn.dat
project_stir STIR_im_norm_rands_attn NP_sino_as_STIR.hs bck STIR_norm2.hs STIR_rands.hs NP_attn_as_STIR.hs
compare_ims_np_space NP_im_norm_rands_attn.dat STIR_im_norm_rands_attn.nii





echo "########################################################################################"
echo "#                                                                                      #"
echo "#    COMPARE FORWARD PROJECTION (NO CORRECTIONS)                                       #"
echo "#                                                                                      #"
echo "########################################################################################"

project_np NP_fwd.dat NP_im.dat fwd "" "" ""
project_stir STIR_fwd.hs STIR_NP_im.nii fwd "" "" ""
compare_sinos_np_space NP_fwd.dat STIR_fwd.hs




echo "########################################################################################"
echo "#                                                                                      #"
echo "#    COMPARE FORWARD PROJECTION WITH NORM, RANDOMS AND ATTENUATION                     #"
echo "#                                                                                      #"
echo "########################################################################################"

project_np NP_fwd_norm_rands_attn.dat NP_im.dat fwd NP_norm.dat NP_rands.dat NP_attn.dat
project_stir STIR_fwd_norm_rands_attn.hs STIR_NP_im.nii fwd STIR_norm2.hs STIR_rands.hs NP_attn_as_STIR.hs
compare_sinos_np_space NP_fwd_norm_rands_attn.dat STIR_fwd_norm_rands_attn.hs




echo "########################################################################################"
echo "#                                                                                      #"
echo "#    COMPARE OSEM ITERATION (NO CORRECTIONS)                                           #"
echo "#                                                                                      #"
echo "########################################################################################"

project_np NP_mlem.dat NP_sino.dat mlem "" "" ""
project_stir STIR_mlem NP_sino_as_STIR.hs mlem "" "" ""
compare_ims_np_space NP_mlem.dat STIR_mlem.nii



echo "########################################################################################"
echo "#                                                                                      #"
echo "#    COMPARE OSEM ITERATION WITH NORM, RANDOMS AND ATTENUATION                         #"
echo "#                                                                                      #"
echo "########################################################################################"

project_np NP_mlem_norm_rands_attn.dat NP_sino.dat mlem NP_norm.dat NP_rands.dat NP_attn.dat 1
project_stir STIR_mlem_norm_rands_attn NP_sino_as_STIR.hs mlem STIR_norm2.hs STIR_rands.hs NP_attn_as_STIR.hs 1
compare_ims_np_space NP_mlem_norm_rands_attn.dat STIR_mlem_norm_rands_attn.nii




echo "########################################################################################"
echo "#                                                                                      #"
echo "#    RESULTS                                                                           #"
echo "#                                                                                      #"
echo "########################################################################################"

rm tmp_* eff_out.txt geom_out.txt inter_out.txt 2> /dev/null

if [ $fail_count -eq 0 ]; then
	echo -e "\nFinished with no errors!"
	exit 0
else
	echo -e "\nFinished with $fail_count errors..."
	exit 1
fi