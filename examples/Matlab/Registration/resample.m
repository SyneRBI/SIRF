function resample(varargin)
% Resampling of SIRF images. 
%   --eng_ref <eng>              engine for reference image [default: Reg]
%   --eng_flo <eng>              engine for floating image [default: Reg]
%   --ref <file>                 reference image (default: test.nii.gz)
%   --flo <file>                 floating image (default: test2.nii.gz)
%   --algo <algo>                resampling algorithm [default: NiftyResample]
%   --output <file>              output image filename [default: output]
%   --intrp <intrp>              interpolation order, defaults to cubic [default: 3]
%   --trans_filenames ...        transformation filenames, (with quotations): "filename1,filename2,filename3"
%   --trans_types ...            transformation types, e.g. (with quotations): "AffineTransformation,NiftiImageData3DDeformation,NiftiImageData3DDisplacement"

% CCP PETMR Synergistic Image Reconstruction Framework (SIRF).
% Copyright 2018 - 2019 University College London.
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

trans_filenames={};
trans_types={};

%% Parse the input
i=1;
while i <= length(varargin)
    if     strcmp(varargin{i},'eng_ref')
        eng_ref = get_arg(varargin,i,1);
        i=i+2;
    elseif strcmp(varargin{i},'eng_flo')
        eng_flo = get_arg(varargin,i,1);
        i=i+2;
    elseif strcmp(varargin{i},'ref')
        ref_file = get_arg(varargin,i,1);
        i=i+2;
    elseif strcmp(varargin{i},'flo')
        flo_file = get_arg(varargin,i,1);
        i=i+2;
   elseif strcmp(varargin{i},'algo')
        algo = get_arg(varargin,i,1);
        i=i+2;
   elseif strcmp(varargin{i},'output')
        output = get_arg(varargin,i,1);
        i=i+2;
    elseif strcmp(varargin{i},'intrp')
        intrp = str2num(get_arg(varargin,i,1));
        i=i+2;
    elseif strcmp(varargin{i},'trans')
        trans_filenames = [trans_filenames; get_arg(varargin,i,1)];
        trans_types     = [trans_types;     get_arg(varargin,i,2)];
        i=i+3;
    else
        error(['Unknown argument: ' varargin{i} '. Use help(function) for help.']);  
    end
end

%% Default values
% If using default ref or flo images, need SIRF data
if ~exist('ref_file','var') || ~exist('flo_file','var')
    SIRF_PATH     = getenv('SIRF_PATH');
    examples_path = fullfile(SIRF_PATH, '/data/examples/Registration');
end

if ~exist('eng_ref','var')  eng_ref  = 'Reg'; end
if ~exist('eng_flo','var')  eng_flo  = 'Reg'; end
if ~exist('ref_file','var') ref_file = fullfile(examples_path, 'test.nii.gz');  end
if ~exist('flo_file','var') flo_file = fullfile(examples_path, 'test2.nii.gz'); end
if ~exist('algo','var')     algo     = 'NiftyResample'; end
if ~exist('output','var')   output   = 'output'; end
if ~exist('intrp','var')    intrp    = 3; end

%% Resampling
% Open reference and floating images
disp(['Engine for reference image: ', eng_ref])
disp(['Engine for floating image: ', eng_flo])
disp(['Reference image: ', ref_file])
disp(['Floating image: ', flo_file])

% Dynamically set up the engines required
set_up_Reg();
eng_ref = set_up_engine(eng_ref);
eng_flo = set_up_engine(eng_flo);

ref = eng_ref.ImageData(ref_file);
flo = eng_flo.ImageData(flo_file);

% Dynamically create resample algorithm
res = eval(['mReg.' algo]);
res.set_reference_image(ref)
res.set_floating_image(flo)
res.set_interpolation_type(intrp)

% create and add each transformation
for i=1:size(trans_filenames)
  disp(['Transformation ' i ' filename: ' trans_filenames(i)])
  disp(['Transformation ' i ' type: ' trans_types(i)])
  trans = eval(['mReg.' trans_types(i) '(' trans_filenames(i) ');']);
  res.add_transformation(trans);
end
 
% Resample
res.process()
 
% Output
res.get_output().write(output)

end

function arg = get_arg(list,index,increment)
    if (index+increment > length(list))
        error(['Not sufficient arguments following: ' list{index}])
    end
    arg = list{index+increment};
end