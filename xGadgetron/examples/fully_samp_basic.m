function fully_samp_basic
% FULLY_SAMP_BASIC Demo for reconstruction of fully sampled MR data.
% 
% Demonstrates use of the EPSRC-funded CCP-PETMR code (SIRF). 
%
% Pre-requisites:
% 1) This MATLAB code needs to be able to access a listening gadgetron.
%    On the Virtual Machine, gadgetron is installed and the user just needs
%    to type 'gadgetron' in a terminal window.
%    On standalone systems, the user will need to have installed ISMRMRD
%    and gadgetron code.
%
% 2) An input data file in the ISMRMRD format.
%    Example datasets:
%    a) 'meas_MID00103_FID57244_test.dat' is 
%       available from https://www.ccppetmr.ac.uk/downloads
%       This is in the manufacturer's raw data format and needs to be
%       converted to ISMRMRD format using 'siemens_to_ismrmrd'.
%       This executable is installed on the Virtual Machine.
%       This is a large dataset from a fully sampled, 3D acquisition.
% 
%    b) 'meas_MID00107_FID57248_test_2D.dat' similar to above but is
%    multi-slce 2D (currently not avaialable for download). Slices not
%    acquired in sequential order.
%
%    c) An ISMRMRD h5 file can be simulated using test_create_dataset.m
%    available in the ISMRMRD/examples/matlab folder at:
%    https://github.com/ismrmrd/ismrmrd
%    This gives 5 repetitions of a single slice.
%
% Usage:
%  fully_samp_basic
%
%
% Adapted by David Atkinson (D.Atkinson@ucl.ac.uk) from original 
% code by Evgueni Ovtchinnikov
%
% See also GRAPPA_BASIC


% load the SIRF mutilities library
if ~libisloaded('mutilities')
    fprintf('loading mutilities library...\n')
    [notfound, warnings] = loadlibrary('mutilities');
end

% This demo uses gadgetron as the MR reconstruction engine. 
% Load the SIRF library and import the SIRF class.
if ~libisloaded('mgadgetron')
    fprintf('loading mgadgetron library...\n')
    [notfound, warnings] = loadlibrary('mgadgetron');
end

import mGadgetron.*

% Get the filename of the input ISMRMRD h5 file
disp('Select ISMRMRD H5 file')
[fn,pn] = uigetfile('*','Select ISMRMRD H5 file') ;
filein = fullfile(pn,fn) ;

% Load this ISMRMRD h5 file
input_data = AcquisitionData(filein);

% Pre-process to remove oversampling
preprocessed = MR_remove_x_oversampling(input_data);


% Perform reconstruction of the preprocessed data.
% 1. set the reconstruction to be for Cartesian fully sampled data.
%    SimpleReconstruction() sets up a default gadget train.
recon = SimpleReconstruction();

% 2. set the reconstruction input to be the data we just preprocessed.
recon.set_input(preprocessed);

% 3. run (i.e. 'process') the reconstruction.
fprintf('---\n reconstructing...\n');
recon.process();

% Extract an image object from the reconstruction and convert this
% to a MATLAB array.
image_obj = recon.get_output();
idata = image_obj.as_array();  

% For these demos, individual coil images are combined in the gadgetron 
% by square root sum of squres and so the next line is not necessary:
idata = abs(idata) ;

idata = idata./max(idata(:)) ;
idata = reshape(idata,[size(idata,1) size(idata,2) 1 size(idata,3)]) ;
figure('Name',['recon of file: ',fn])
montage(idata)


    



