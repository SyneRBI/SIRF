function fully_sampled_recon_single_chain(engine)
% FULLY_SAMPLED_RECON_SINGLE_CHAIN Complete recon in one Gadgetron chain
%
% See FULLY_SAMPLED_RECON for example data files.

% SyneRBI Synergistic Image Reconstruction Framework (SIRF).
% Copyright 2015 - 2019 Rutherford Appleton Laboratory STFC.
% Copyright 2015 - 2019 University College London.
% 
% This is software developed for the Collaborative Computational
% Project in Synergistic Reconstruction for Biomedical Imaging (formerly CCP PETMR)
% (http://www.ccpsynerbi.ac.uk/).
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
% import_str = set_up_MR(engine);
% eval(import_str)
MR = set_up_MR(engine);
mr_data_path = sirf.Utilities.examples_data_path('MR');

% In this demo, rather than using a predefined image reconstruction 
% object, here an image reconstruction object is created by concatinating 
% multiple gadgets - for more information on Gadgetron and its gadgets,
% please see: 
%  https://github.com/gadgetron/.

% Parameters for individual gadgets can be defined either during the 
% creation of the reconstruction object:
% e.g. AcquisitionAccumulateTriggerGadget(trigger_dimension=repetition)
% or by giving a gadget a label (cf. label ex: for the last gadget)
% and using set_gadget_property(label, propery, value).
% The gadgets will be concatenated and will be executed as soon as 
% process() is called.
gadgets = [...
    {'RemoveROOversamplingGadget'}, ...
    {'AcquisitionAccumulateTriggerGadget'}, ...
    {'BucketToBufferGadget'}, ...
    {'SimpleReconGadget'}, ...
    {'ImageArraySplitGadget'}, ...
    {'ex:ExtractGadget'} ...
    ];

% create reconstructor
recon = MR.Reconstructor(gadgets);

% change a property of the gadget labelled 'ex'
% ExtractGadget defines which type of image should be returned:
% none      0
% magnitude 1
% real      2
% imag      4
% phase     8
% max       16  
% in this example '5' returns both magnitude and imag 
%recon.set_gadget_property('ex', 'extract_mask', 5);
recon.set_gadget_property('ex', 'extract_magnitude', true);
recon.set_gadget_property('ex', 'extract_imag', true);

% define raw data source
[filename, pathname] = uigetfile('*.h5', 'Select raw data file', mr_data_path);
acq_data = MR.AcquisitionData(fullfile(pathname, filename));
recon.set_input(acq_data)

% optionally set Gadgetron server host and port
recon.set_host('localhost')
% On VM you can try a port other than the default 9002, e.g. 9003, by taking
% the following steps:
% 1) in ~/devel/install/share/gadgetron/config/gadgetron.xml replace
%    <port>9002</port> with <port>9003</port>
% 2) go to Settings->Network->Advanced->Port Forwarding and add new rule
%    (click on green + in the upper right corner) with Host and Guest ports
%    set to 9003
% 3) uncomment the next line
%recon.set_port('9003')
% Note: each gadget chain can run on a different VM - to try, start two VMs
% and do the above steps 1 and 2 on one of them, then add recon.set_port('9003')
% before recon.process in grappa_detail.py (where preprocessing runs on
% port 9002).

% perform reconstruction
recon.process()

% get reconstructed image object
image_data = recon.get_output();

% show reconstructed image data
if exist('montage','file') && exist('mat2gray','file')
    figure('Name','image data')
    idisp = mat2gray(abs(image_data.as_array()));
    montage(reshape(idisp,[size(idisp,1) size(idisp,2) 1 size(idisp,3)])) ;
else
    image_array = image_data.as_array();
    title = 'Reconstructed image data (magnitude)';
    sirf.Utilities.show_3D_array(abs(image_array(:,:,1:2:end)), title, ...
        'samples', 'readouts', 'slice');
    title = 'Reconstructed image data (imaginary part)';
    sirf.Utilities.show_3D_array(imag(image_array(:,:,2:2:end)), title, ...
        'samples', 'readouts', 'slice');
end
