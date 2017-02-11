function fully_sampled_recon
% FULLY_SAMPLED_RECON  Recon of fully sampled data
% See FULLY_SAMP_BASIC for information about example data files.
%
% See also FULLY_SAMP_BASIC FULLY_SAMPLED_RECON_SINGLE_CHAIN
% FULLY_SAMPLED_RECON_THREE_CHAINS


% load mutilities and mgadgetron libraries
ccp_libload

% import mGadgetron MATLAB classes so that they can be called in this
% function without using the prefix 'mGadgetron.'
import mGadgetron.*


% MR raw data formats from different vendors can be transformed to 
% HDF file format using siemens_to_ismrmrd, philips_to_ismrmrd or
% bruker_to_ismrmrd on https://github.com/ismrmrd/.
% acquisitions will be read from this HDF file
[filename, pathname] = uigetfile('*.h5', 'Select raw data file');
input_data = AcquisitionData(fullfile(pathname, filename));

% pre-process acquisition data
% Prior to image reconstruction several pre-processing steps such as 
% asymetric echo compensation, noise decorrelation for multi-coil data or 
% removal of oversampling along frequency encoding (i.e. readout or kx)
% direction. So far only the removal of readout oversampling and noise and
% asymmetric echo adjusting is implemented
fprintf('processing acquisitions...\n')
processed_data = preprocess_acquisitions(input_data);

% perform reconstruction
% Create a reconstruction object (in this case simple 2D Cartesian FFT) and
% provide pre-processed k-space data as input
recon = SimpleReconstruction();
recon.set_input(processed_data)

% perform reconstruction
fprintf('reconstructing...\n')
recon.process()

% retrieve reconstruction as object
images = recon.get_output();

% use as_array method and plot obtained images

idisp = mat2gray(abs(images.as_array()));
montage(reshape(idisp,[size(idisp,1) size(idisp,2) 1 size(idisp,3)])) ;



