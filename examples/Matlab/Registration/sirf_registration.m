function registration(varargin)
% Registration of SIRF images. 
%   --eng_ref <eng>              engine for reference image [default: Reg]
%   --eng_flo <eng>              engine for floating image [default: Reg]
%   --ref <file>                 reference image (default: test.nii.gz)
%   --flo <file>                 floating image (default: test2.nii.gz)
%   --par <file>                 parameter file (default: niftyreg_aladin.par)
%   --algo <algo>                registration algorithm [default: NiftyAladinSym]
%   --rmask                      mask of reference image
%   --fmask                      mask of floating image
%   --warped <file>              warped image filename [default: warped]
%   --TM_forward                 forward transformation matrix (if rigid/affine)
%   --TM_inverse                 inverse transformation matrix (if rigid/affine)
%   --disp_fwd_4D                4D forward displacement field image
%   --def_fwd_4D                 4D forward deformation field image
%   --disp_inv_4D                4D inverse displacement field image
%   --def_inv_4D                 4D inverse deformation field image

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
   elseif strcmp(varargin{i},'warped')
        warped = get_arg(varargin,i,1);
        i=i+2;
    elseif strcmp(varargin{i},'par')
        par_file = get_arg(varargin,i,1);
        i=i+2;
        
    elseif strcmp(varargin{i},'rmask')
        rmask = get_arg(varargin,i,1);
        i=i+2;        
    elseif strcmp(varargin{i},'fmask')
        fmask = get_arg(varargin,i,1);
        i=i+2;        
        
    elseif strcmp(varargin{i},'TM_forward')
        TM_forward = get_arg(varargin,i,1);
        i=i+2;        
    elseif strcmp(varargin{i},'TM_inverse')
        TM_inverse = get_arg(varargin,i,1);
        i=i+2;        
    elseif strcmp(varargin{i},'disp_fwd_4D')
        disp_fwd_4D = get_arg(varargin,i,1);
        i=i+2;        
    elseif strcmp(varargin{i},'disp_inv_4D')
        disp_inv_4D = get_arg(varargin,i,1);
        i=i+2;        
    elseif strcmp(varargin{i},'def_fwd_4D')
        def_fwd_4D = get_arg(varargin,i,1);
        i=i+2;        
    elseif strcmp(varargin{i},'def_inv_4D')
        def_inv_4D = get_arg(varargin,i,1);
        i=i+2;
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
if ~exist('algo','var')     algo     = 'NiftyAladinSym'; end
if ~exist('warped','var')   warped   = 'warped'; end
if ~exist('par_file','var') par_file = fullfile(examples_path, 'paramFiles/niftyreg_aladin.par'); end

%% Registration
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
reg = eval(['mReg.' algo]);
reg.set_reference_image(ref)
reg.set_floating_image(flo)
if exist('par_file','var')
    reg.set_parameter_file(par_file);
end
if exist('rmask_file','var')
    rmask = eng_ref.ImageData(rmask_file);
    reg.set_reference_mask(rmask)
end
if exist('fmask_file','var')
    fmask = eng_ref.ImageData(fmask_file);
    reg.set_floating_mask(fmask)
end

reg.process();

%% Output 
reg.get_output().write(warped);

% TMs
if exist('TM_forward','var')
    reg.get_transformation_matrix_forward().write(TM_forward)
end
if exist('TM_inverse','var')
    reg.get_transformation_matrix_forward().write(TM_inverse)
end
 
% Disp fields
if exist('disp_fwd_4D','var')
    reg.get_displacement_field_forward().write(disp_fwd_4D)
end
if exist('disp_inv_4D','var')
    reg.get_displacement_field_inverse().write(disp_inv_4D)
end

% Def fields
if exist('def_fwd_4D','var')
    reg.get_deformation_field_forward().write(def_fwd_4D)
end
if exist('disp_inv_4D','var')
    reg.get_deformation_field_inverse().write(def_inv_4D)
end

end

function arg = get_arg(list,index,increment)
    if (index+increment > length(list))
        error(['Not sufficient arguments following: ' list{index}])
    end
    arg = list{index+increment};
end