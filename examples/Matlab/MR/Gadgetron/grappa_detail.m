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

% load mutilities and mgadgetron libraries
ccp_libload

% import mGadgetron MATLAB classes so that they can be called in this
% function without using the prefix 'mGadgetron.'
import mGadgetron.*


% Get the filename of the input ISMRMRD h5 file
disp('Select ISMRMRD H5 file')
[fn,pn] = uigetfile('*.h5','Select ISMRMRD H5 file', mr_data_path) ;
filein = fullfile(pn,fn) ;

% Initially we create a container that points to the h5 file. Data is
% not read from file until the gadgetron is called using
% the 'process' method.

% Create an Acquisition Container. Here because of the previous
% 'import mGadgetron.*', this will be of type mGadgetron.AcquisitionData
input_Cont = AcquisitionData(filein);

% Pre-process this input data using three preparation gadgets
% from gadgetron.
% List gadgets to use (not all may be required for this test data).
prep_gadgets = [{'NoiseAdjustGadget'}, ...
    {'AsymmetricEchoAdjustROGadget'} ...
    {'RemoveROOversamplingGadget'} ];

% Call gadgetron by using the 'process' method. This runs the gadgets
% specified in prep_gadgets, returning an instance
% of an mGadgetron.AcquisitionsContainer
preprocessed_AcCont = input_Cont.process(prep_gadgets);

% Extract sorted k-space, permute dimensions and display
ksp = preprocessed_AcCont.as_array ;
[ns,nc,nro] = preprocessed_AcCont.dimensions ; % [nx ncoil ny]
ksp = permute(ksp,[3 1 2]) ; %  [ny nx ncoil]
ksp_coil_as_col = reshape(ksp,[nro nc*ns]) ; 
if exist('imshow','file') && exist('imadjust','file') && exist('mat2gray','file')
    figure('Name','Input k-space')
    disp(['Displaying k-space'])
    imshow(imadjust(mat2gray(abs(ksp_coil_as_col)),[0 0.7],[],0.2))
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
recon.set_input(preprocessed_AcCont);

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
images_imcont = recon.get_output('image');
gfacts_imcont = recon.get_output('gfactor');

% Return as MATLAB matrices the data pointed to by the containers.
% Note the image data is complex.
idata = images_imcont.as_array();
gdata = gfacts_imcont.as_array();

% Display figures of modulus image, phase image and gfactor map
% Requires ImageProcessing Toolbox
if exist('mat2gray','file') && exist('montage','file')
    figure('Name',['recon modulus images of file: ',fn])
    idisp = mat2gray(abs(idata));
    montage(reshape(idisp,[size(idisp,1) size(idisp,2) 1 size(idisp,3)])) ;
    
    figure('Name',['recon phase images of file: ',fn])
    idisp = mat2gray(angle(idata),[-pi pi]);
    montage(reshape(idisp,[size(idisp,1) size(idisp,2) 1 size(idisp,3)])) ;
    
    figure('Name','Gfactor')
    idisp = mat2gray(abs(gdata));
    montage(reshape(idisp,[size(idisp,1) size(idisp,2) 1 size(idisp,3)])) ;
else
    idata = abs(idata);
    gdata = abs(gdata);
    idata = idata/max(max(max(idata)));
    gdata = gdata/max(max(max(gdata)));
    n = images_imcont.number();
    fprintf('Enter slice number to view its data\n')
    fprintf('(a value outside the range [1 : %d] will stop this loop).\n', n)
    while (true)
        i = input('slice: ');
        if i < 1 || i > n
            break
        end
        figure(i)
        imshow(idata(:,:,i));
        title(['image ' num2str(i)])
        figure(i + n)
        imshow(gdata(:,:,i));
        title(['G-factor ' num2str(i)])
    end
end




