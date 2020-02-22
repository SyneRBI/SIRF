function sirf_registration(varargin)
%SIRF_REGISTRATION Registration of SIRF images.
%
%   Use name-pair values to set required and optional arguments
%
%   Required arguments
%   -------------------------------------------------------------
%
%   --ref <file> <eng>           reference image (eng: Reg|STIR|Gadgetron)
%       default: --ref test.nii.gz Reg
%   --flo <file> <eng>           floating image (eng: Reg|STIR|Gadgetron)
%       default: --flo test2.nii.gz Reg
%   --algo <algo>                registration algorithm (aladin/f3d/spm)
%       default: --algo aladin
%
%   Optional arguments
%   -------------------------------------------------------------
%
%   --warped_prefix <fname>      warped image filename
%   --disp_fwd_prefix <fname>    forward displacement field image
%   --disp_inv_prefix <fname>    inverse displacement field image
%   --def_fwd_prefix <fname>     forward deformation field image
%   --def_inv_prefix <fname>     inverse deformation field image
%
%   Optional arguments for rigid/affine algorithms (aladin/spm)
%   -------------------------------------------------------------
%
%   --TM_fwd_prefix <fname>      forward transformation matrix
%   --TM_inv_prefix <fname>      inverse transformation matrix
%
%   Optional flags for NiftyReg (aladin/f3d)
%   -------------------------------------------------------------
%   --rmask <fname> <eng>        mask of reference image (eng: Reg|STIR|Gadgetron)
%   --fmask <fname> <eng>        mask of floating image (eng: Reg|STIR|Gadgetron)
%   --print                      print all possible wrapped parameters and exit.
%   --par_file <fname>           set parameter file
%   --par <string>               set wrapped parameter. Some examples (and note quotation marks):
%           '--par', 'SetPerformRigid 1'
%           '--par', 'SetPerformAffine 0'
%           '--par', 'SetInterpolationToCubic'
%           '--par', 'SetFloatingThresholdUp 1 2'
%
%   Optional flags for spm
%   -------------------------------------------------------------

%   --working_folder <fname>     folder in which to save temporary files (default: cwd/spm_working_folder)
%   --overwrite <bool>           should I overwrite files if already present? (default: true)
%   --delete <bool>              should I delete temporary files? (default: true)
%
%   If any are missing, these defaults are used
%   ------------------------------------------------------------- 
%   --ref test.nii.gz Reg
%   --flo test2.nii.gz Reg
%   --algo aladin

% CCP PETMR Synergistic Image Reconstruction Framework (SIRF).
% Copyright 2018 - 2020 University College London.
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
flo_files={};
flo_engs={};
pars={};
while i <= length(varargin)
    % Required
    if strcmp(varargin{i},'--ref')
        ref_file = get_arg(varargin,i,1);
        eng_ref = get_arg(varargin,i,2);
        i=i+3;
    elseif strcmp(varargin{i},'--flo')
        flo_files{end+1} = get_arg(varargin,i,1);
        flo_engs{end+1} = get_arg(varargin,i,2);
        i=i+3;
    elseif strcmp(varargin{i},'--algo')
        algo = get_arg(varargin,i,1);
        i=i+2;
        
    % Optional
    elseif strcmp(varargin{i},'--warped_prefix')
        warped = get_arg(varargin,i,1);
        i=i+2;
    elseif strcmp(varargin{i},'--disp_fwd_prefix')
        disp_fwd = get_arg(varargin,i,1);
        i=i+2;        
    elseif strcmp(varargin{i},'--disp_inv_prefix')
        disp_inv = get_arg(varargin,i,1);
        i=i+2;        
    elseif strcmp(varargin{i},'--def_fwd_prefix')
        def_fwd = get_arg(varargin,i,1);
        i=i+2;        
    elseif strcmp(varargin{i},'--def_inv_prefix')
        def_inv = get_arg(varargin,i,1);
        i=i+2;
        
    % Optional NiftyReg
    elseif strcmp(varargin{i},'--rmask')
        rmask = get_arg(varargin,i,1);
        rmask_eng = get_arg(varargin,i,2);
        i=i+3;
    elseif strcmp(varargin{i},'--fmask')
        fmask = get_arg(varargin,i,1);
        fmask_eng = get_arg(varargin,i,2);
        i=i+3;
    elseif strcmp(varargin{i},'--print')
        print = true;
        i=i+1;
    elseif strcmp(varargin{i},'--par_file')
        par_file = get_arg(varargin,i,1);
        i=i+2;
    elseif strcmp(varargin{i},'--par')
        pars{end+1} = get_arg(varargin,i,1);
        i=i+2;
        
    % Optional rigid/affine
    elseif strcmp(varargin{i},'--TM_fwd_prefix')
        TM_fwd = get_arg(varargin,i,1);
        i=i+2;        
    elseif strcmp(varargin{i},'--TM_inv_prefix')
        TM_inv = get_arg(varargin,i,1);
        i=i+2; 
        
    % Optional spm
    elseif strcmp(varargin{i},'--working_folder')
        working_folder = get_arg(varargin,i,1);
        i=i+2;
    elseif strcmp(varargin{i},'--overwrite')
        overwrite = get_arg(varargin,i,1);
        i=i+2;
    elseif strcmp(varargin{i},'--delete')
        delete_temp_files = get_arg(varargin,i,1);
        i=i+2;

    else
        error(['Unknown argument: ' varargin{i} '. Use help(function) for help.']);  
    end    
end

%% Get defaults
if ~exist('algo','var'); algo='aladin'; end
if ~exist('ref_file','var') || ~exist('flo_file','var')
    SIRF_PATH     = getenv('SIRF_PATH');
    examples_path = fullfile(SIRF_PATH,'/data/examples/Registration');
end
if ~exist('ref_file','var')
    ref_file=fullfile(examples_path,'test.nii.gz');
    eng_ref = 'Reg';
end
if isempty(flo_files)
    flo_files{1} = fullfile(examples_path,'test2.nii.gz');
    flo_engs{1}  = 'Reg';
end

%% Get algorithm
set_up_Reg();
reg = get_algorithm(algo);

%% Print wrapped methods
if exist('print','var')
    assert(strcmp(algo,'spm')~=1, '--print not available for spm')
    reg.print_all_wrapped_methods();
    return;
end

%% Input images
reg.set_reference_image(get_im(ref_file,eng_ref));
for i=1:numel(flo_files)
    reg.add_floating_image(get_im(flo_files{i},flo_engs{i}));
end

%% Set parameters
% rmask
if exist('rmask','var')
    assert(strcmp(algo,'spm')~=1, '--rmask not available for spm')
    reg.set_reference_mask(get_im(rmask,rmask_eng));
end
% fmask
if exist('fmask','var')
    assert(strcmp(algo,'spm')~=1, '--fmask not available for spm')
    reg.set_floating_mask(get_im(fmask,fmask_eng));
end
% par_file
if exist('par_file','var')
    assert(strcmp(algo,'spm')~=1, '--par_file not available for spm')
    reg.set_parameter_file(par_file);
end
% pars
if (strcmp(algo,'spm')); assert(isempty(pars), '--pars not available for spm'); end
for i=1:numel(pars)
    disp(pars)
    pars_split = split(pars);
    if numel(pars_split)==1
        reg.set_parameter(pars_split{1});
    elseif numel(pars_split)==2
        reg.set_parameter(pars_split{1},pars_split{2});
    elseif numel(pars_split)==3
        reg.set_parameter(pars_split{1},pars_split{2},pars_split{3});
    else
        error('Max number of NiftyReg args is 2.')
    end
end

% working_folder
if strcmp(algo,'spm') && ~exist('working_folder','var')
    working_folder = [pwd filesep 'spm_working_folder'];
end
if exist('working_folder','var')
    assert(strcmp(algo,'spm'), '--working_folder only available for spm')
    reg.set_working_folder(working_folder);
end
% overwrite
if strcmp(algo,'spm') && ~exist('overwrite','var')
    overwrite = true;
end
if exist('overwrite','var')
    assert(strcmp(algo,'spm'), '--overwrite only available for spm')
    reg.set_working_folder_file_overwrite(overwrite);
end
% delete_temp_files
if strcmp(algo,'spm') && ~exist('delete_temp_files','var')
    delete_temp_files = true;
end
if exist('delete_temp_files','var')
    assert(strcmp(algo,'spm'), '--delete only available for spm')
    reg.set_delete_temp_files(delete_temp_files);
end

%% Register
reg.process();

%% Save results
for i = 1 : numel(flo_files)
    
    if exist('warped','var')
        reg.get_output(i).write([warped num2str(i)]);
    end
    if exist('disp_fwd','var')
        reg.get_displacement_field_forward(i).write([disp_fwd num2str(i)]);
    end
    if exist('disp_inv','var')
        reg.get_displacement_field_inverse(i).write([disp_inv num2str(i)]);
    end
    if exist('def_fwd','var')
        reg.get_deformation_field_forward(i).write([def_fwd num2str(i)]);
    end
    if exist('def_inv','var')
        reg.get_deformation_field_inverse(i).write([def_inv num2str(i)]);
    end
    if exist('TM_fwd','var')
        assert(strcmp(algo,'f3d'), '--TM_fwd only available for rigid/affine algorithms')
        reg.get_transformation_matrix_forward(i).write([TM_fwd num2str(i)]);
    end
    if exist('TM_inv','var')
        assert(strcmp(algo,'f3d'), '--TM_inv only available for rigid/affine algorithms')
        reg.get_transformation_matrix_inverse(i).write([TM_inv num2str(i)]);
    end
end
end

function im = get_im(filename, eng)
    if strcmp(eng, 'Reg')
        im = sirf.Reg.ImageData(filename);
    elseif strcmp(eng, 'STIR')
        im = sirf.STIR.ImageData(filename);
    elseif strcmp(eng, 'Gadgetron')
        im = sirf.Gadgetron.ImageData(filename);
    else
        error(['Unknown engine: ' engine]);
    end
end

function reg = get_algorithm(algo_str)
    if strcmp(algo_str,'aladin')
        reg = sirf.Reg.NiftyAladinSym();
    elseif strcmp(algo_str,'f3d')
        reg = sirf.Reg.NiftyF3DSym();
    elseif strcmp(algo_str,'spm')
        reg = sirf.Reg.SPMRegistration();
    else
        error(['Unknown algorithm: ' algo_str]);
    end
end

function arg = get_arg(list,index,increment)
    if (index+increment > length(list))
        error(['Not sufficient arguments following: ' list{index}])
    end
    arg = list{index+increment};
end
