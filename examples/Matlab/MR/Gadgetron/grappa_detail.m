function grappa_detail
% GRAPPA_DETAIL Annotated demo for reconstruction of GRAPPA acquired data.
% 
% Demonstrates use of the EPSRC-funded CCP-PETMR code (SIRF). 
% See function grappa_basic.m for a simpler example.
%
% Pre-requisites:
% 1) This MATLAB code needs to be able to access a listening gadgetron.
%    On the Virtual Machine, gadgetron is installed and the user just needs
%    to type 'gadgetron' in a terminal window.
%    On standalone systems, the user will need to have installed ISMRMRD
%    and gadgetron code.
%
% 2) An input data file from a GRAPPA MRI acquisition in the ISMRMRD format.
%    Example GRAPPA datasets:
%    a) 'meas_MID00108_FID57249_test_2D_2x.dat' is 
%       available from https://www.ccppetmr.ac.uk/downloads
%       This is in the manufacturer's raw data format and needs to be
%       converted to ISMRMRD format using 'siemens_to_ismrmrd'.
%       This executable is installed on the Virtual Machine.
%
%    b) An ISMRMRD h5 file can be simulated using the MATLAB function
%    gen_us_data.m . This will reconstruct faster than the real data in a).
%
% Usage:
%  grappa_detail
%
%
% See also GRAPPA_BASIC GEN_US_DATA

% CCP PETMR Synergistic Image Reconstruction Framework (SIRF).
% Copyright 2015 - 2017 Rutherford Appleton Laboratory STFC.
% Copyright 2015 - 2017 University College London.
% 
% This is software developed for the Collaborative Computational
% Project in Positron Emission Tomography and Magnetic Resonance imaging
% (http://www.ccppetmr.ac.uk/).
% 
% Licensed under the Apache License, Version 2.0 (the "License");
% you may not use this file except in compliance with the License.
% You may obtain a copy of the License at
% http://www.apache.org/licenses/LICENSE-2.0
% Unless required by applicable law or agreed to in writing, software
% distributed under the License is distributed on an "AS IS" BASIS,
% WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
% See the License for the specific language governing permissions and
% limitations under the License.

% default engine to be used if none given
if nargin < 1
    engine = [];
end
import_str = set_up_MR(engine);
eval(import_str)

% Get the filename of the input ISMRMRD h5 file
[fn,pn] = uigetfile('*.h5','Select ISMRMRD H5 file', mr_data_path) ;
filein = fullfile(pn,fn) ;

% Initially we create a container that points to the h5 file. Data is
% not read from file until the gadgetron is called using
% the 'process' method.

% Create an Acquisition Container. Here because of the previous
% 'import mGadgetron.*', this will be of type mGadgetron.AcquisitionData
acq_data = AcquisitionData(filein);

% Pre-process this input data using three preparation gadgets
% from gadgetron.
% List gadgets to use (not all may be required for this test data).
prep_gadgets = [{'NoiseAdjustGadget'}, ...
    {'AsymmetricEchoAdjustROGadget'} ...
    {'RemoveROOversamplingGadget'} ];

% Call gadgetron by using the 'process' method. This runs the gadgets
% specified in prep_gadgets, returning an instance
% of an mGadgetron.AcquisitionsContainer
preprocessed_data = acq_data.process(prep_gadgets);

% Extract sorted k-space, permute dimensions and display
preprocessed_array = preprocessed_data.as_array ;
[ns,nc,nro] = preprocessed_data.dimensions ; % [nx ncoil ny]
if exist('imshow','file') && exist('imadjust','file') && exist('mat2gray','file')
    preprocessed_array = permute(preprocessed_array,[3 1 2]) ; %  [ny nx ncoil]
    preprocessed_coil_as_col = reshape(preprocessed_array,[nro nc*ns]) ; 
    figure('Name','Input k-space')
    fprintf('Displaying k-space\n')
    imshow(imadjust(mat2gray(abs(preprocessed_coil_as_col)),[0 0.7],[],0.2))
else
    preprocessed_array = permute(preprocessed_array,[1 3 2]) ; %  [nx ny ncoil]
    title = 'Acquisition data (magnitude)';
    mUtilities.show_3D_array...
        (abs(preprocessed_array).^0.2, title, 'samples', 'readouts', 'coil');
    mUtilities.set_window(0.1, 0.1, 0.8, 0.8)
end

% Perform reconstruction of the preprocessed data.

% 1) Create a recon object for the desired reconstruction.

% In this demo, the recon object is created using the class
% Reconstructor(). A simpler class is available in the SIRF code
% for a GRAPPA reconstruction:
%   recon = CartesianGRAPPAReconstructor()
%
%    To find what this does behind the scenes:
%     type edit mGadgetron.CartesianGRAPPAReconstructor
%     and note the name assigned in the self function, here
%       'SimpleGRAPPAReconstructionProcessor'.
%     Then find the gadget chain defined by the class with the same
%     name in the file xGadgetron/cGadgetron/chain_lib.h
%

recon_gadgets =  [...
    {'AcquisitionAccumulateTriggerGadget'}, ...
    {'BucketToBufferGadget'}, ...
    {'GenericReconCartesianReferencePrepGadget'}, ...
    {'GRAPPA:GenericReconCartesianGrappaGadget'}, ...
    {'GenericReconFieldOfViewAdjustmentGadget'}, ...
    {'GenericReconImageArrayScalingGadget'}, ...
    {'ImageArraySplitGadget'} ...
    ];

recon = Reconstructor(recon_gadgets) ;


% 2) The GRAPPA gadget can compute G-factors in addition to
% reconstructed images. We can set a gadget property as below if the gadget
% has been identified with a label. In the above list of recon_gadgets,
% the 4th is labelled 'GRAPPA' and we can use this label as below:
recon.set_gadget_property('GRAPPA', 'send_out_gfactor', true)

% If the chain had been set using
% recon = CartesianGRAPPAReconstructor(), an alternative method
% would be available:
%  recon.compute_gfactors(true);


% 3) set the reconstruction input to be the data we just preprocessed.
recon.set_input(preprocessed_data);

% 4) Run the reconstruction using 'process' to call gadgetron.
fprintf('---\n reconstructing...\n');
recon.process();

% Output

% Reconstructed data sits in memory. We need to first get containers
% for both the reconstructed images and g-factors, before extracting the
% data as MATLAB arrays. Containers in effect point to the data.

% Get images and gfactors as containers with type mGadgetron.ImagesContainer
% (Note this syntax may change in the future with the addition of a
%  method '.get_gfactor'.)
image_data = recon.get_output('image');
gfact_data = recon.get_output('gfactor');

% Return as MATLAB matrices the data pointed to by the containers.
% Note the image data is complex.
image_array = image_data.as_array();
gfact_array = gfact_data.as_array();

% Display figures of modulus image, phase image and gfactor map
% Requires ImageProcessing Toolbox
if exist('mat2gray','file') && exist('montage','file')
    figure('Name',['recon modulus images of file: ',fn])
    idisp = mat2gray(abs(image_array));
    montage(reshape(idisp,[size(idisp,1) size(idisp,2) 1 size(idisp,3)])) ;
    
    figure('Name',['recon phase images of file: ',fn])
    idisp = mat2gray(angle(image_array),[-pi pi]);
    montage(reshape(idisp,[size(idisp,1) size(idisp,2) 1 size(idisp,3)])) ;
    
    figure('Name','Gfactor')
    idisp = mat2gray(abs(gfact_array));
    montage(reshape(idisp,[size(idisp,1) size(idisp,2) 1 size(idisp,3)])) ;
else
    image_array = abs(image_array);
    gfact_array = abs(gfact_array);
    title = 'Reconstructed image data (magnitude)';
    mUtilities.show_3D_array(image_array, title, 'samples', 'readouts', 'slice');
    title = 'G-factor data (magnitude)';
    mUtilities.show_3D_array(gfact_array, title, 'samples', 'readouts', 'slice');
end




