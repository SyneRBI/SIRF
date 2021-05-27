% SyneRBI Synergistic Image Reconstruction Framework (SIRF).
% Copyright 2018 - 2020 University College London
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

set_up_Reg();

% Paths
SIRF_PATH     = getenv('SIRF_PATH');
examples_path = fullfile(SIRF_PATH, '/data/examples/Registration');
output_prefix = fullfile(pwd, '/results/');

% Input filenames
g.ref_aladin_filename                        = fullfile(examples_path, '/test.nii.gz');
g.flo_aladin_filename                        = fullfile(examples_path, '/test2.nii.gz');
g.ref_f3d_filename                           = fullfile(examples_path, '/mouseFixed.nii.gz');
g.flo_f3d_filename                           = fullfile(examples_path, '/mouseMoving.nii.gz');
g.parameter_file_aladin                      = fullfile(examples_path, '/paramFiles/niftyreg_aladin.par');
g.parameter_file_f3d                         = fullfile(examples_path, '/paramFiles/niftyreg_f3d.par');

% Output filenames
g.save_nifti_image                           = fullfile(output_prefix, 'matlab_save_NiftiImageData.nii');
g.save_nifti_image_3d                        = fullfile(output_prefix, 'matlab_save_NiftiImageData3D.nii');
g.save_nifti_image_3d_tensor_not_split       = fullfile(output_prefix, 'matlab_save_NiftiImageData3DTensor_not_split.nii');
g.save_nifti_image_3d_tensor_split           = fullfile(output_prefix, 'matlab_save_NiftiImageData3DTensor_split_%s.nii');
g.save_nifti_image_3d_deformation_not_split  = fullfile(output_prefix, 'matlab_save_NiftiImageData3DDeformation_not_split.nii');
g.save_nifti_image_3d_deformation_split      = fullfile(output_prefix, 'matlab_save_NiftiImageData3DDeformation_split_%s.nii');
g.save_nifti_image_3d_displacement_not_split = fullfile(output_prefix, 'matlab_save_NiftiImageData3DDisplacement_not_split.nii');
g.save_nifti_image_3d_displacement_split     = fullfile(output_prefix, 'matlab_save_NiftiImageData3DDisplacement_split_%s.nii');
g.save_nifti_image_upsample                  = fullfile(output_prefix, 'matlab_save_NiftiImageData_upsample.nii');
g.save_nifti_image_downsample                = fullfile(output_prefix, 'matlab_save_NiftiImageData_downsample.nii');
g.save_nifti_image_up_downsample             = fullfile(output_prefix, 'matlab_save_NiftiImageData_upsample_downsample.nii');
g.aladin_warped                              = fullfile(output_prefix, 'matlab_aladin_warped.nii');
g.f3d_warped                                 = fullfile(output_prefix, 'matlab_f3d_warped.nii');
g.TM_forward		                     = fullfile(output_prefix, 'matlab_TM_forward.txt');
g.TM_inverse		                     = fullfile(output_prefix, 'matlab_TM_inverse.txt');
g.aladin_def_forward                         = fullfile(output_prefix, 'matlab_aladin_def_forward.nii');
g.aladin_def_inverse_xyz                     = fullfile(output_prefix, 'matlab_aladin_def_inverse_%s.nii');
g.aladin_def_inverse                         = fullfile(output_prefix, 'matlab_aladin_def_inverse.nii');
g.aladin_def_fwd_inv                         = fullfile(output_prefix, 'matlab_aladin_def_fwd_then_inv.nii');
g.aladin_disp_forward                        = fullfile(output_prefix, 'matlab_aladin_disp_forward.nii');
g.aladin_disp_inverse                        = fullfile(output_prefix, 'matlab_aladin_disp_inverse_%s.nii');
g.f3d_def_forward                            = fullfile(output_prefix, 'matlab_f3d_def_forward.nii');
g.f3d_def_inverse                            = fullfile(output_prefix, 'matlab_f3d_def_inverse_%s.nii');
g.f3d_disp_forward                           = fullfile(output_prefix, 'matlab_f3d_disp_forward.nii');
g.f3d_disp_inverse                           = fullfile(output_prefix, 'matlab_f3d_disp_inverse_%s.nii');

g.rigid_resample                             = fullfile(output_prefix, 'matlab_rigid_resample.nii');
g.nonrigid_resample_disp                     = fullfile(output_prefix, 'matlab_nonrigid_resample_disp.nii');
g.nonrigid_resample_def                      = fullfile(output_prefix, 'matlab_nonrigid_resample_def.nii');
g.niftymomo_resample_adj                     = fullfile(output_prefix, 'matlab_niftymomo_resample_adj.nii');
g.output_weighted_mean                       = fullfile(output_prefix, 'matlab_weighted_mean.nii');
g.output_weighted_mean_def                   = fullfile(output_prefix, 'matlab_weighted_mean_def.nii');
g.output_float                               = fullfile(output_prefix, 'matlab_reg_aladin_float.nii');

g.ref_aladin                                 = sirf.Reg.NiftiImageData3D( g.ref_aladin_filename );
g.flo_aladin                                 = sirf.Reg.NiftiImageData3D( g.flo_aladin_filename );
g.ref_f3d                                    = sirf.Reg.NiftiImageData3D(   g.ref_f3d_filename  );
g.flo_f3d                                    = sirf.Reg.NiftiImageData3D(   g.flo_f3d_filename  );

% Check if we have niftiread
toolkits=ver;
have_niftiread = any(strcmp({toolkits.Name},'Image Processing Toolbox'));

% You can change these when debugging
try_niftiimage = true;
try_niftiimage3d = true;
try_niftiimage3dtensor = true;
try_niftiimage3ddisplacement = true;
try_niftiimage3ddeformation = true;
try_niftyaladin = true;
try_niftyf3d = true;
try_transformations = true;
try_resample = true;
try_niftymomo = true;
try_weighted_mean = true;
try_affinetransformation = true;
try_quaternion = true;

if try_niftiimage
	disp('% ----------------------------------------------------------------------- %')
	disp('%                  Starting NiftiImageData test...')
	disp('%------------------------------------------------------------------------ %')

    % default constructor
    a = sirf.Reg.NiftiImageData();

    % Read from file
    b = sirf.Reg.NiftiImageData(g.ref_aladin_filename);

    % Save to file
    b.write(g.save_nifti_image);

    % Fill
    b.fill(100);

    % Get max
    assert(b.get_max() == 100, 'NiftiImageData fill()/get_max() failed.');

    % Get min
    assert(b.get_min() == 100, 'NiftiImageData fill()/get_min() failed.');

    % Deep copy
    d = b.deep_copy();
    assert(d.handle_ ~= b.handle_, 'NiftiImageData deep_copy failed.');
    assert(d == b, 'NiftiImageData deep_copy failed.');

    % Addition
    e = d + d;
    assert(abs(e.get_max() - 2 * d.get_max()) < 0.0001, 'NiftiImageData __add__/get_max() failed.')

    % Subtraction
    e = d - d;
    assert(abs(e.get_max()) < 0.0001, 'NiftiImageData __sub__ failed.')

    % Sum
    assert(abs(e.get_sum()) < 0.0001, 'NiftiImageData get_sum() failed.')

    % Add num to image
    q = e + 1;
    assert(q.get_max() == e.get_max() + 1, 'NiftiImageData __add__ val failed.');

    % Subtract num from image
    r = e - 1;
    assert(r.get_max() == e.get_max() - 1, 'NiftiImageData __sub__ val failed.');

    % Multiply image by num
    s = e * 10;
    assert(s.get_max() == e.get_max() * 10, 'NiftiImageData __mul__ val failed.');

    % Dimensions
    f = e.get_dimensions();
    assert(all(f == [3, 64, 64, 64, 1, 1, 1, 1]), 'NiftiImageData get_dimensions() failed.')

    % Get as array
    arr = d.as_array();
    assert(max(arr(:)) == 100, 'NiftiImageData as_array().max() failed.')
    assert(ndims(arr) == 3, 'NiftiImageData as_array() ndims failed.')
    assert(all(size(arr) == [64, 64, 64]), 'NiftiImageData as_array().shape failed.')

    % Test saving to datatype
    g.ref_aladin.write(g.output_float,16); % save to float
    ref_aladin_float = sirf.Reg.NiftiImageData3D(g.output_float);
    arr1 = g.ref_aladin.as_array();
    arr2 = ref_aladin_float.as_array();
    assert(all(arr1(:)==arr2(:)), 'NiftiImageData::write()/change_datatype() failed.');

    % Test print methods
    q.print_header();
    sirf.Reg.NiftiImageData.print_headers([q s]);

    % Crop image
    min_ = [];
    max_ = [];
    for i=1:7
        min_(i) = 0;
        max_(i) = f(i+1) - 1;
    end
    max_(3) = 62;
    e = e;
    s.crop(min_,max_);
    assert(all(size(s.as_array()) == [64, 64, 63]), 'NiftiImageData crop() failed.')

    % Get voxel sizes
    s = b.get_voxel_sizes();
    assert(all(s == [0, 4.0625, 4.0625, 4.0625, 0, 0, 0, 0]), 'NiftiImageData get_voxel_sizes() failed.')

    % Check upsampling/downsampling
    u = sirf.Reg.NiftiImageData(g.ref_aladin_filename);
    original_spacing    = u.get_voxel_sizes();
    original_spacing    = original_spacing(2:4);
    upsampled_spacing   = [original_spacing(1)/2, original_spacing(2)/4, original_spacing(3)];
    downsampled_spacing = [original_spacing(1)*2, original_spacing(2)*4, original_spacing(3)];
    % Downsample
    v = u.deep_copy();
    v.set_voxel_spacing(downsampled_spacing,3);
    v.write(g.save_nifti_image_downsample);
    % Upsample then downsample, check nothing has changed
    w = u.deep_copy();
    w.set_voxel_spacing(upsampled_spacing,0);
    w.write(g.save_nifti_image_upsample);
    x = w.deep_copy();
    x.set_voxel_spacing(original_spacing,0);
    x.write(g.save_nifti_image_up_downsample);
    sirf.Reg.NiftiImageData.print_headers([u v w x]);
    assert(x == u, 'NiftiImageData::upsample()/downsample() failed.')

    % Check get_contains_nans
    x_arr = x.as_array();
    x_arr(:)=0;
    x.fill(x_arr);
    assert(~x.get_contains_nans(),'NiftiImageData::get_contains_nans() 1 failed.')
    x_arr(1) = nan;
    x.fill(x_arr);
    assert(x.get_contains_nans(),'NiftiImageData::get_contains_nans() 2 failed.')

    if have_niftiread
        arr1 = sirf.Reg.NiftiImageData(g.ref_aladin_filename).as_array();
        arr2 = niftiread(g.ref_aladin_filename);
        assert(all(all(all(arr1 == arr2))), 'NiftiImageData as_array() failed.')
    end

    % Test geom info
    im = sirf.Reg.NiftiImageData(g.ref_aladin_filename);
    geom_info = im.get_geometrical_info();
    geom_info.print_info();
    % Get voxel sizes
    assert(all(geom_info.get_size() == [64, 64, 64]), 'SIRF get_geometrical_info().get_size() failed.');
    assert(all(geom_info.get_spacing() == [4.0625, 4.0625, 4.0625]), 'SIRF get_geometrical_info().get_spacing() failed.');

    im.standardise();
    assert(abs(im.get_standard_deviation() - 1) < 0.01, 'NiftiImageData standardise() or get_standard_deviation() failed.')
    assert(abs(im.get_variance() - 1) < 0.01, 'NiftiImageData standardise() or get_variance() failed.')
    assert(abs(im.get_mean()) < 0.0001, 'NiftiImageData standardise() or get_mean() failed.')

    % Check normalise 
    im.normalise_zero_and_one();
    assert(abs(im.get_min()) < 0.0001 && abs(im.get_max() - 1) < 0.0001, 'NiftiImageData normalise_between_zero_and_one() failed.')

    % Test inner product
    in1 = x.deep_copy();
    in2 = x.deep_copy();
    in1_arr = in1.as_array();
    in2_arr = in2.as_array();
    dims = in1.get_dimensions();
    for idx_x = 1 : dims(2)
        for idx_y = 1 : dims(3)
            for idx_z = 1 : dims(4)
                in1_arr(idx_x, idx_y, idx_z) = single(i-1);
                in2_arr(idx_x, idx_y, idx_z) = single(3*(i-1)-1);
            end
        end
    end
    inner_product = sum(in1_arr(:) .* in2_arr(:));
    in1.fill(in1_arr);
    in2.fill(in2_arr);
    disp(inner_product)
    disp(in1.get_inner_product(in2))
    assert(abs(inner_product - in1.get_inner_product(in2)) < 1e-4, 'NiftiImageData::get_inner_product() failed.');

    % Pad then crop, should be the same
    aa = g.ref_aladin;
    cc = aa.clone();
    original_dims = aa.get_dimensions();

    pad_in_min_dir = [1, 2, 3, 0, 0, 0, 0];
    pad_in_max_dir = [4, 5, 6, 0, 0, 0, 0];
    cc.pad(pad_in_min_dir, pad_in_max_dir, 100.);

    padded_dims = cc.get_dimensions();
    for i = 1:7
        assert(padded_dims(i+1) == original_dims(i+1) + pad_in_min_dir(i) + pad_in_max_dir(i), ...
            'NiftiImageData::pad failed')
    end

    % Crop back to beginning
    cropped_min_dir = pad_in_min_dir;
    for i = 1:7
        cropped_max_dir(i) = original_dims(i+1) + cropped_min_dir(i) - 1;
    end
    cc.crop(cropped_min_dir, cropped_max_dir);
    assert(aa == cc, 'NiftiImageData::pad/crop failed')

    disp('% ----------------------------------------------------------------------- %')
    disp('%                  Finished NiftiImageData test.')
    disp('%------------------------------------------------------------------------ %')
end

if try_niftiimage3d
    disp('% ----------------------------------------------------------------------- %')
    disp('%                  Starting NiftiImageData3D test...')
    disp('%------------------------------------------------------------------------ %')

    % default constructor
    a = sirf.Reg.NiftiImageData3D();

    % Read from file
    b = sirf.Reg.NiftiImageData3D(g.ref_aladin_filename);

    % Save to file
    b.write(g.save_nifti_image_3d);

    % Fill
    b.fill(100);

    % Get max
    assert(b.get_max() == 100, 'NiftiImageData3D fill()/get_max() failed.');

    % Get min
    assert(b.get_min() == 100, 'NiftiImageData3D fill()/get_min() failed.');

    % Deep copy
    d = b.deep_copy();
    assert(d.handle_ ~= b.handle_, 'NiftiImageData3D deep_copy failed.');
    assert(d == b, 'NiftiImageData3D deep_copy failed.');

    % Addition
    e = d + d;
    assert(abs(e.get_max() - 2 * d.get_max()) < 0.0001, 'NiftiImageData3D __add__/get_max() failed.')

    % Subtraction
    e = d - d;
    assert(abs(e.get_max()) < 0.0001, 'NiftiImageData3D __sub__ failed.')

    % Sum
    assert(abs(e.get_sum()) < 0.0001, 'NiftiImageData3D get_sum() failed.')

    % Dimensions
    f = e.get_dimensions();
    assert(all(f == [3, 64, 64, 64, 1, 1, 1, 1]), 'NiftiImageData3D get_dimensions() failed.')

    % Get as array
    arr = d.as_array();
    assert(max(arr(:)) == 100, 'NiftiImageData3D as_array().max() failed.')
    assert(ndims(arr) == 3, 'NiftiImageData3D as_array() ndims failed.')
    assert(all(size(arr) == [64, 64, 64]), 'NiftiImageData3D as_array().shape failed.')

    % try linear algebra
    h = d/10000;
    assert(abs(h.get_max()-d.get_max()/10000) < 1e-4,'NiftiImageData3D linear algebra failed.')

    % Check as_array and fill is symmetric
    ref_aladin_arr = g.ref_aladin.as_array();
    ref_aladin2 = g.ref_aladin.clone();
    ref_aladin2.fill(ref_aladin_arr);
    assert(ref_aladin2 == g.ref_aladin, 'NiftiImageData3D::as_array()/fill() failed.')

    disp('% ----------------------------------------------------------------------- %')
    disp('%                  Finished NiftiImageData3D test.')
    disp('%------------------------------------------------------------------------ %')
end

if try_niftiimage3dtensor
    disp('% ----------------------------------------------------------------------- %')
    disp('%                  Starting NiftiImageData3DTensor test...')
    disp('%------------------------------------------------------------------------ %')

    % Create NiftiImageData3DTensor from NiftiImageData3D
    b = sirf.Reg.NiftiImageData3DTensor();
    b.create_from_3D_image(g.ref_aladin);

    % Save to file
    b.write(g.save_nifti_image_3d_tensor_not_split);
    b.write_split_xyz_components(g.save_nifti_image_3d_tensor_split);

    % Constructor from file
    c = sirf.Reg.NiftiImageData3DTensor(g.save_nifti_image_3d_tensor_not_split);

    % Fill
    c.fill(100)

    % Get max
    assert(c.get_max() == 100, 'NiftiImageData3DTensor fill()/get_max() failed.');

    % Get min
    assert(c.get_min() == 100, 'NiftiImageData3DTensor fill()/get_min() failed.');

    % Deep copy
    d = c.deep_copy();
    assert(d.handle_ ~= c.handle_, 'NiftiImageData3DTensor deep_copy failed (they have the same handle).');
    assert(d == c, 'NiftiImageData3DTensor deep_copy failed (values do not match).');

    % Addition
    e = d + d;
    assert(abs(e.get_max() - 2 * d.get_max()) < 0.0001, 'NiftiImageData3DTensor __add__/get_max() failed.')

    % Subtraction
    e = d - d;
    assert(abs(e.get_max()) < 0.0001, 'NiftiImageData3DTensor __sub__ failed.')

    % Sum
    assert(abs(e.get_sum()) < 0.0001, 'NiftiImageData3DTensor get_sum() failed.')

    % Dimensions
    f = e.get_dimensions();
    assert(all(f == [5, 64, 64, 64, 1, 3, 1, 1]), 'NiftiImageData3DTensor get_dimensions() failed.')

    % Get as array
    arr = d.as_array();
    assert(max(arr(:)) == 100, 'NiftiImageData3DTensor as_array().max() failed.')
    assert(ndims(arr) == 5, 'NiftiImageData3DTensor as_array() ndims failed.')
    assert(all(size(arr) == [64, 64, 64, 1, 3]), 'NiftiImageData3DTensor as_array().shape failed.')

    % Constructor from single components
    im1 = g.ref_aladin.deep_copy();
    im2 = g.ref_aladin.deep_copy();
    im3 = g.ref_aladin.deep_copy();
    im1.fill(30);
    im2.fill(20);
    im3.fill(-10);
    h = sirf.Reg.NiftiImageData3DTensor(im1, im2, im3);

    % Test flip components
    h.flip_component(0);
    assert(h.get_max() ==  20, 'NiftiImageData3DTensor flip_component() failed.');
    assert(h.get_min() == -30, 'NiftiImageData3DTensor flip_component() failed.');


    disp('% ----------------------------------------------------------------------- %')
    disp('%                  Finished NiftiImageData3DTensor test.')
    disp('%------------------------------------------------------------------------ %')
end

if try_niftiimage3ddisplacement
    disp('% ----------------------------------------------------------------------- %')
    disp('%                  Starting NiftiImageData3DDisplacement test...')
    disp('%------------------------------------------------------------------------ %')

    % Create NiftiImageData3DDisplacement from NiftiImageData3D
    b = sirf.Reg.NiftiImageData3DDisplacement();
    b.create_from_3D_image(g.ref_aladin);

    % Save to file
    b.write(g.save_nifti_image_3d_displacement_not_split);
    b.write_split_xyz_components(g.save_nifti_image_3d_displacement_split);

    % Constructor from file
    c = sirf.Reg.NiftiImageData3DDisplacement(g.save_nifti_image_3d_displacement_not_split);

    % Constructor from 3x3D
    h = sirf.Reg.NiftiImageData3DDisplacement(g.ref_aladin, g.ref_aladin, g.ref_aladin);

    % Fill
    c.fill(100)

    % Get max
    assert(c.get_max() == 100, 'NiftiImageData3DDisplacement fill()/get_max() failed.');

    % Get min
    assert(c.get_min() == 100, 'NiftiImageData3DDisplacement fill()/get_min() failed.');

    % Deep copy
    d = c.deep_copy();
    assert(d.handle_ ~= c.handle_, 'NiftiImageData3DDisplacement deep_copy failed (they have the same handle).');
    assert(d == c, 'NiftiImageData3DDisplacement deep_copy failed (values do not match).');

    % Addition
    e = d + d;
    assert(abs(e.get_max() - 2 * d.get_max()) < 0.0001, 'NiftiImageData3DDisplacement __add__/get_max() failed.')

    % Subtraction
    e = d - d;
    assert(abs(e.get_max()) < 0.0001, 'NiftiImageData3DDisplacement __sub__ failed.')

    % Sum
    assert(abs(e.get_sum()) < 0.0001, 'NiftiImageData3DDisplacement get_sum() failed.')

    % Dimensions
    f = e.get_dimensions();
    assert(all(f == [5, 64, 64, 64, 1, 3, 1, 1]), 'NiftiImageData3DDisplacement get_dimensions() failed.')

    % Get as array
    arr = d.as_array();
    assert(max(arr(:)) == 100, 'NiftiImageData3DDisplacement as_array().max() failed.')
    assert(ndims(arr) == 5, 'NiftiImageData3DDisplacement as_array() ndims failed.')
    assert(all(size(arr) == [64, 64, 64, 1, 3]), 'NiftiImageData3DDisplacement as_array().shape failed.')

    % Check upsampling/downsampling
    u = sirf.Reg.NiftiImageData3DDisplacement(g.save_nifti_image_3d_displacement_not_split);
    original_spacing    = u.get_voxel_sizes();
    original_spacing    = original_spacing(2:4);
    upsampled_spacing   = [original_spacing(1)/2, original_spacing(2)/4, original_spacing(3)];
    downsampled_spacing = [original_spacing(1)*2, original_spacing(2)*4, original_spacing(3)];
    % Downsample
    v = u.deep_copy();
    v.set_voxel_spacing(downsampled_spacing,3);
    v.write(g.save_nifti_image_downsample);
    % Upsample then downsample, check nothing has changed
    w = u.deep_copy();
    w.set_voxel_spacing(upsampled_spacing,0);
    w.write(g.save_nifti_image_upsample);
    x = w.deep_copy();
    x.set_voxel_spacing(original_spacing,0);
    x.write(g.save_nifti_image_up_downsample);
    sirf.Reg.NiftiImageData.print_headers([u v w x]);
    assert(x == u, 'NiftiImageData3DDisplacement::upsample()/downsample() failed.')


    disp('% ----------------------------------------------------------------------- %')
    disp('%                  Finished NiftiImageData3DDisplacement test.')
    disp('%------------------------------------------------------------------------ %')
end

if try_niftiimage3ddeformation
    disp('% ----------------------------------------------------------------------- %')
    disp('%                  Starting NiftiImageData3DDeformation test...')
    disp('%------------------------------------------------------------------------ %')

    % Create NiftiImageData3DDeformation from NiftiImageData3D
    b = sirf.Reg.NiftiImageData3DDeformation();
    b.create_from_3D_image(g.ref_aladin);

    % Save to file
    b.write(g.save_nifti_image_3d_deformation_not_split);
    b.write_split_xyz_components(g.save_nifti_image_3d_deformation_split);

    % Constructor from file
    c = sirf.Reg.NiftiImageData3DDeformation(g.save_nifti_image_3d_deformation_not_split);

    % Constructor from 3x3D
    h = sirf.Reg.NiftiImageData3DDeformation(g.ref_aladin, g.ref_aladin, g.ref_aladin);

    % Fill
    c.fill(100)

    % Get max
    assert(c.get_max() == 100, 'NiftiImageData3DDeformation fill()/get_max() failed.');

    % Get min
    assert(c.get_min() == 100, 'NiftiImageData3DDeformation fill()/get_min() failed.');

    % Deep copy
    d = c.deep_copy();
    assert(d.handle_ ~= c.handle_, 'NiftiImageData3DDeformation deep_copy failed (they have the same handle).');
    assert(d == c, 'NiftiImageData3DDeformation deep_copy failed (values do not match).');

    % Addition
    e = d + d;
    assert(abs(e.get_max() - 2 * d.get_max()) < 0.0001, 'NiftiImageData3DDeformation __add__/get_max() failed.')

    % Subtraction
    e = d - d;
    assert(abs(e.get_max()) < 0.0001, 'NiftiImageData3DDeformation __sub__ failed.')

    % Sum
    assert(abs(e.get_sum()) < 0.0001, 'NiftiImageData3DDeformation get_sum() failed.')

    % Dimensions
    f = e.get_dimensions();
    assert(all(f == [5, 64, 64, 64, 1, 3, 1, 1]), 'NiftiImageData3DDeformation get_dimensions() failed.')

    % Get as array
    arr = d.as_array();
    assert(max(arr(:)) == 100, 'NiftiImageData3DDeformation as_array().max() failed.')
    assert(ndims(arr) == 5, 'NiftiImageData3DDeformation as_array() ndims failed.')
    assert(all(size(arr) == [64, 64, 64, 1, 3]), 'NiftiImageData3DDeformation as_array().shape failed.')

    disp('% ----------------------------------------------------------------------- %')
    disp('%                  Finished NiftiImageData3DDeformation test.')
    disp('%------------------------------------------------------------------------ %')
end

if try_niftyaladin
	disp('% ----------------------------------------------------------------------- %')
	disp('%                  Starting Nifty aladin test...')
	disp('%------------------------------------------------------------------------ %')

    % Print all wrapped methods.
    sirf.Reg.NiftyAladinSym.print_all_wrapped_methods();

    % First set up some masks
    ref_mask = g.ref_aladin.deep_copy();
    flo_mask = g.flo_aladin.deep_copy();
    ref_mask.fill(1);
    flo_mask.fill(1);

	% default constructor
    na = sirf.Reg.NiftyAladinSym();
    na.set_reference_image(g.ref_aladin);
    na.set_floating_image(g.flo_aladin);
    na.set_parameter_file(g.parameter_file_aladin);
    na.set_parameter('SetInterpolationToCubic');
    na.set_parameter('SetLevelsToPerform', '1');
    na.set_parameter('SetMaxIterations', '5');
    na.set_parameter('SetPerformRigid', '1');
    na.set_parameter('SetPerformAffine', '0');
    na.set_reference_mask(ref_mask);
    na.set_floating_mask(flo_mask);
    na.process();

    % Get outputs
    warped = na.get_output().deep_copy();
    def_forward = na.get_deformation_field_forward().deep_copy();
    def_inverse = na.get_deformation_field_inverse().deep_copy();
    disp_forward = na.get_displacement_field_forward().deep_copy();
    disp_inverse = na.get_displacement_field_inverse().deep_copy();
    TM_forward_ = na.get_transformation_matrix_forward().deep_copy();
    TM_inverse_ = na.get_transformation_matrix_inverse().deep_copy();

    % Test via filenames
    na.set_reference_image_filename(g.ref_aladin_filename);
    na.set_floating_image_filename(g.flo_aladin_filename);
    na.process();

    assert(warped == na.get_output() && ...
        def_forward == na.get_deformation_field_forward() && ...
        def_inverse == na.get_deformation_field_inverse() && ...
        disp_forward == na.get_displacement_field_forward() && ...
        disp_inverse == na.get_displacement_field_inverse() && ...
        TM_forward_ == na.get_transformation_matrix_forward() && ...
        TM_inverse_ == na.get_transformation_matrix_inverse(),...
        'Registration via filenames failed')

    warped.write(g.aladin_warped);
    na.get_transformation_matrix_forward().write(g.TM_forward);
    na.get_transformation_matrix_inverse().write(g.TM_inverse);
    def_forward.write(g.aladin_def_forward);
    def_inverse.write_split_xyz_components(g.aladin_def_inverse_xyz);
    def_inverse.write(g.aladin_def_inverse);
    disp_forward.write(g.aladin_disp_forward);
    disp_inverse.write_split_xyz_components(g.aladin_disp_inverse);

    % forward TM
    forward_tm = na.get_transformation_matrix_forward().as_array()

    % Inverse TM
    inverse_tm = na.get_transformation_matrix_inverse().as_array()

    % Test converting disp to def
    a = sirf.Reg.NiftiImageData3DDeformation(disp_forward);
    assert(a == def_forward, 'NiftiImageData3DDeformation::create_from_disp() failed.');

    % Test converting def to disp
    b = sirf.Reg.NiftiImageData3DDisplacement(def_forward);
    assert(b == disp_forward, 'NiftiImageData3DDisplacement::create_from_def() failed.');

    % Check NiftiImageData3DDeformation::get_inverse()
    def_fwd_then_inv = def_forward.get_inverse(g.flo_aladin);
    def_fwd_then_inv.write(g.aladin_def_fwd_inv);
    sirf.Reg.NiftiImageData.print_headers([g.ref_aladin, g.flo_aladin, def_inverse, def_fwd_then_inv]);

    % Reference forward with def_inv
    resample = sirf.Reg.NiftyResampler();
    resample.set_reference_image(g.flo_aladin);
    resample.set_floating_image(g.ref_aladin);
    resample.set_padding_value(0.);
    resample.set_interpolation_type_to_linear();
    resample.add_transformation(def_inverse);
    out1 = resample.forward(g.ref_aladin);

    % Reference forward with def_fwd_then_inv
    resample.clear_transformations();
    resample.add_transformation(def_fwd_then_inv);
    out2 = resample.forward(g.ref_aladin);

    sirf.Reg.NiftiImageData.print_headers([out1, out2]);
    assert(out1 == out2, 'NiftiImageData3DDeformation::get_inverse() failed.')

	disp('% ----------------------------------------------------------------------- %')
	disp('%                  Finished Nifty aladin test.')
	disp('%------------------------------------------------------------------------ %')
end

if try_niftyf3d
	disp('% ----------------------------------------------------------------------- %')
	disp('%                  Starting Nifty f3d test...')
	disp('%------------------------------------------------------------------------ %')

    % Crop input to increase speed
    dim = g.ref_f3d.get_dimensions();
    mid = dim(2:4)/2;
    min_idx = [ mid(1)-5,mid(2)-5,mid(3)-5,0,0,0,0 ];
    max_idx = [ mid(1)+5,mid(2)+4,mid(3)+3,0,0,0,0 ];
    ref_f3d_crop = g.ref_f3d.clone();
    ref_f3d_crop.crop(min_idx, max_idx);
    flo_f3d_crop = g.flo_f3d.clone();
    flo_f3d_crop.crop(min_idx, max_idx);

	% Print all wrapped methods.
	sirf.Reg.NiftyF3dSym.print_all_wrapped_methods();

    % Get initial transformation
    tm_init = sirf.Reg.AffineTransformation(g.TM_forward);

	% default constructor
    nf = sirf.Reg.NiftyF3dSym();
    nf.set_reference_image(ref_f3d_crop);
    nf.set_floating_image(flo_f3d_crop);
    nf.set_parameter_file(g.parameter_file_f3d);
    % nf.set_reference_time_point(1);
    % nf.set_floating_time_point(1);
    nf.set_initial_affine_transformation(tm_init);
    nf.process();

    % Get outputs
    warped = nf.get_output();
    def_forward = nf.get_deformation_field_forward();
    def_inverse = nf.get_deformation_field_inverse();
    disp_forward = nf.get_displacement_field_forward();
    disp_inverse = nf.get_displacement_field_inverse();

    warped.write(g.f3d_warped);
    def_forward.write(g.f3d_def_forward);
    def_inverse.write_split_xyz_components(g.f3d_def_inverse);
    disp_forward.write(g.f3d_disp_forward);
    disp_inverse.write_split_xyz_components(g.f3d_disp_inverse);

    % Compare between sirf.Reg.NiftiImageData3DDefofmation::as_array() and niftiread
    if have_niftiread
        deff_arr = def_forward.as_array();
        deff_matlab_arr = niftiread(g.f3d_def_forward);
        assert(all(deff_arr(:) == deff_matlab_arr(:)), 'NiftiImageData3DDeformation as_array() failed.')
    end

    % Check as_array and fill for deformation fields
    deff2 = def_forward.clone();
    deff2.fill(deff_arr);
    assert(def_forward == deff2, 'NiftiImageData3DDeformation::as_array()/fill() failed.')

	disp('% ----------------------------------------------------------------------- %')
	disp('%                  Finished Nifty f3d test.')
	disp('%------------------------------------------------------------------------ %')
end

if try_transformations
	disp('% ----------------------------------------------------------------------- %')
	disp('%                  Starting Transformation test...')
	disp('%------------------------------------------------------------------------ %')


    % Get transformations
    a3 = na.get_transformation_matrix_forward();
    b3 = na.get_displacement_field_forward();
    c3 = na.get_deformation_field_forward();

    % Get as deformations
    a_def = a3.get_as_deformation_field(g.ref_aladin);
    b_def = b3.get_as_deformation_field(g.ref_aladin);
    c_def = c3.get_as_deformation_field(g.ref_aladin);
    assert(a_def == na.get_deformation_field_forward(), 'TransformationAffine get_as_deformation_field() failed.')
    assert(b_def == na.get_deformation_field_forward(), 'TransformationDisplacement get_as_deformation_field() failed.')
    assert(c_def == na.get_deformation_field_forward(), 'TransformationDeformation get_as_deformation_field() failed.')

    % Compose into single deformation. Use two identity matrices and the disp field. Get as def and should be the same.
    tm_iden = sirf.Reg.AffineTransformation.get_identity();
    trans = [tm_iden, tm_iden, c3];
    composed = sirf.Reg.NiftiImageData3DDeformation.compose_single_deformation(trans, g.ref_aladin);
    assert(composed == na.get_deformation_field_forward(), 'compose_single_deformation failed.')

    % Test get_inverse
    tm_inv = tm_iden.get_inverse();


	disp('% ----------------------------------------------------------------------- %')
	disp('%                  Finished Transformation test.')
	disp('%------------------------------------------------------------------------ %')
end

if try_resample
    disp('% ----------------------------------------------------------------------- %')
    disp('%                  Starting Nifty resample test...')
    disp('%------------------------------------------------------------------------ %')

	tm_iden = sirf.Reg.AffineTransformation.get_identity();
    tm      = na.get_transformation_matrix_forward();
    displ   = na.get_displacement_field_forward();
    deff    = na.get_deformation_field_forward();
    padding_value = -20;

    disp('Testing rigid resample...')
    nr1 = sirf.Reg.NiftyResampler();
    nr1.set_reference_image(g.ref_aladin);
    nr1.set_floating_image(g.flo_aladin);
    nr1.set_interpolation_type_to_cubic_spline();  % try different interpolations
    nr1.set_interpolation_type(3);  % try different interpolations (cubic)
    nr1.add_transformation(tm_iden);
    nr1.clear_transformations();
    nr1.add_transformation(tm_iden);
    nr1.add_transformation(tm);
    nr1.process();
    nr1.get_output().write(g.rigid_resample);

    disp('Testing non-rigid displacement...')
    nr2 = sirf.Reg.NiftyResampler();
    nr2.set_reference_image(g.ref_aladin);
    nr2.set_floating_image(g.flo_aladin);
    nr2.set_interpolation_type_to_sinc();  % try different interpolations
    nr2.set_interpolation_type_to_nearest_neighbour();  % try different interpolations
    nr2.add_transformation(displ);
    nr2.set_padding_value(padding_value);
    nr2.process();
    nr2.get_output().write(g.nonrigid_resample_disp);

    assert(nr2.get_output().get_min() == padding_value, 'NiftyResampler:set_padding_value failed.')

    disp('Testing non-rigid deformation...')
    nr3 = sirf.Reg.NiftyResampler();
    nr3.set_reference_image(g.ref_aladin)
    nr3.set_floating_image(g.flo_aladin)
    nr3.set_interpolation_type_to_linear()  % try different interpolations
    nr3.add_transformation(deff);
    nr3.set_interpolation_type_to_linear()
    nr3.process()
    nr3.get_output().write(g.nonrigid_resample_def)

    % Check that the following give the same result
    %       out = resample.forward(in)
    %       resample.forward(out, in)
    out1 = nr3.forward(g.flo_aladin);
    out2 = g.ref_aladin.deep_copy();
    nr3.forward(out2, g.flo_aladin);
    assert(out1 == out2, 'out = NiftyResampler::forward(in) and NiftyResampler::forward(out, in) do not give same result.')

    % TODO this doesn't work. For some reason (even with NiftyReg directly), resampling with the TM from the registration
    % doesn't give the same result as the output from the registration itself (even with same interpolations). Even though 
    % ref and flo images are positive, the output of the registration can be negative. This implies that linear interpolation 
    % is not used to generate final image. You would hope it's used throughout the registration process, otherwise why is it there?
    % assert(na.get_output() == nr1.get_output(), 'Rigid resampled output should match registration (aladin) output.')

    disp('% ----------------------------------------------------------------------- %')
    disp('%                  Finished Nifty resample test.')
    disp('%------------------------------------------------------------------------ %')
end

if try_niftymomo
    disp('% ----------------------------------------------------------------------- %')
    disp('%                  Starting NiftyMomMo test...')
    disp('%------------------------------------------------------------------------ %')

    % The forward and the adjoint should meet the following criterion:
    % | < x, Ty > - < y, Tsx > | / 0.5 * (| < x, Ty > | + | < y, Tsx > |) < epsilon
    % for all images x and y, where T is the transform and Ts is the adjoint.

    x = g.ref_aladin;
    T = na.get_transformation_matrix_forward();
    y = g.flo_aladin;

    % Add in a magnification to make things interesting
    t = T.as_array();
    t(1,1) = 1.5;
    T = sirf.Reg.AffineTransformation(t);

    % make it slightly unsquare to spice things up
    min_idx = [1,2,3];
    y_dims = y.get_dimensions();
    max_idx = [y_dims(2) - 3, y_dims(3) - 1, y_dims(4) - 5];
    y.crop(min_idx, max_idx);

    disp('Testing adjoint resample..')
    nr = sirf.Reg.NiftyResampler();
    nr.set_reference_image(x);
    nr.set_floating_image(y);
    nr.set_interpolation_type_to_linear();
    nr.add_transformation(T);

    % Do the forward
    Ty = nr.forward(y);

    % Do the adjoint
    Tsx = nr.adjoint(x);

    % Check the adjoint is truly the adjoint with: |<x, Ty> - <y, Tsx>| / 0.5*(|<x, Ty>|+|<y, Tsx>|) < epsilon
    inner_x_Ty = x.get_inner_product(Ty);
    inner_y_Tsx = y.get_inner_product(Tsx);
    adjoint_test = abs(inner_x_Ty - inner_y_Tsx) / (0.5 * (abs(inner_x_Ty) + abs(inner_y_Tsx)));
    disp(['<x, Ty>  = ' num2str(inner_x_Ty)])
    disp(['<y, Tsx> = ' num2str(inner_y_Tsx)])
    disp(['|<x, Ty> - <y, Tsx>| / 0.5*(|<x, Ty>|+|<y, Tsx>|) = ' num2str(adjoint_test)])
    assert(adjoint_test < 1e-4, 'NiftyResampler::adjoint() failed')

    % Check that the following give the same result
    %       out = resample.adjoint(in)
    %       resample.adjoint(out, in)
    out1 = nr.adjoint(x);
    out2 = y.deep_copy();
    nr.backward(out2, x);
    assert(out1 == out2, 'out = NiftyResampler::adjoint(in) and NiftyResampler::adjoint(out, in) do not give same result.')

    disp('% ----------------------------------------------------------------------- %')
    disp('%                  Finished NiftyMomMo test.')
    disp('%------------------------------------------------------------------------ %')
end

if try_weighted_mean
    disp('% ----------------------------------------------------------------------- %')
    disp('%                  Starting weighted mean test...')
    disp('%------------------------------------------------------------------------ %')

		% Do 3D
		wm1 = sirf.Reg.ImageWeightedMean();
        % Change to float to avoid rounding errors
        im1 = g.ref_aladin.deep_copy();
        im2 = g.ref_aladin.deep_copy();
        im3 = g.ref_aladin.deep_copy();
        im4 = g.ref_aladin.deep_copy();
        im1.fill(1);
		im2.fill(4);
		im3.fill(7);
		im4.fill(6);
		wm1.add_image(im1, 2);
		wm1.add_image(im2, 4);
		wm1.add_image(im3, 3);
		wm1.add_image(im4, 1);
                wm1.process();
		wm1.get_output().write(g.output_weighted_mean);
		% Answer should be 4.5, so compare it to that!
        res = g.ref_aladin.deep_copy();
		res.fill(4.5);
		assert(wm1.get_output() == res, '3D weighted mean test failed.')

		% Do 4D
		wm2 = sirf.Reg.ImageWeightedMean();
		im1 = na.get_deformation_field_forward().deep_copy();
		im2 = na.get_deformation_field_forward().deep_copy();
		im3 = na.get_deformation_field_forward().deep_copy();
		im4 = na.get_deformation_field_forward().deep_copy();
		im1.fill(1);
		im2.fill(4);
		im3.fill(7);
		im4.fill(6);
		wm2.add_image(im1, 2);
		wm2.add_image(im2, 4);
		wm2.add_image(im3, 3);
		wm2.add_image(im4, 1);
                wm2.process();
		wm2.get_output().write(g.output_weighted_mean_def);
		% Answer should be 4.5, so compare it to that!
		res = na.get_deformation_field_forward().deep_copy();
		res.fill(4.5);
		assert(wm2.get_output() == res, '4D weighted mean test failed.')


    disp('% ----------------------------------------------------------------------- %')
    disp('%                  Finished weighted mean test.')
    disp('%------------------------------------------------------------------------ %')
end

if try_affinetransformation
    disp('% ----------------------------------------------------------------------- %')
    disp('%                  Starting AffineTransformation test...')
    disp('%------------------------------------------------------------------------ %')

    % Construct from file
    a = sirf.Reg.AffineTransformation(g.TM_forward);

    % Multiply forward and inverse, should equal identity
    b = na.get_transformation_matrix_forward();
    c = na.get_transformation_matrix_inverse();
    d = b * c;
    e = sirf.Reg.AffineTransformation.get_identity();
    assert(d == e, 'AffineTransformation::mult/comparison failed.');

    assert(e.get_determinant() - 1. < 1.e-7, 'AffineTransformation::get_determinant failed.');

    % Test get_Euler_angles
    array(4,4) =  0;
    array(1,3) =  1;
    array(2,2) = -1;
    array(3,1) = -1;
    array(4,4) =  1;
    test_Eul = sirf.Reg.AffineTransformation(array);
    % Example given by rotm2eul for MATLAB is [0 0 1; 0 -1 0; -1 0 0] -> XYZ = [-3.1416 1.5708 0]
    Eul = test_Eul.get_Euler_angles();
    Eul_expected = [-3.1416, 1.5708, 0];
    assert(all(abs(Eul-Eul_expected) < 1e-4), 'AffineTransformation get_Euler_angles() failed.')

    % Check as_array
    f = b.as_array();
    gg = sirf.Reg.AffineTransformation(f);
    h = gg.as_array();
    assert(all(all(abs(f-h) < 1e-4)), 'AffineTransformation as_array() failed.')

    % Average!
    trans = [0., 0., 0.];
    quat_1_array = [0.92707,  0.02149,   0.19191,  0.32132];
    quat_2_array = [0.90361,  0.0025836, 0.097279, 0.41716];
    quat_3_array = [0.75868, -0.21289,   0.53263,  0.30884];
    quat_1 = sirf.Reg.Quaternion(quat_1_array);
    quat_2 = sirf.Reg.Quaternion(quat_2_array);
    quat_3 = sirf.Reg.Quaternion(quat_3_array);
    tm_1   = sirf.Reg.AffineTransformation(trans,quat_1);
    tm_2   = sirf.Reg.AffineTransformation(trans,quat_2);
    tm_3   = sirf.Reg.AffineTransformation(trans,quat_3);
    average = sirf.Reg.AffineTransformation.get_average([tm_1, tm_2, tm_3]);
    exptd_avg_array = [ 0.5836, -0.6736, 0.4535, 0;,...
                        0.6007,  0.7339, 0.3171, 0;,...
                       -0.5464,  0.0874, 0.8329, 0;,...
                        0,       0,      0,      1];
    exptd_average = sirf.Reg.AffineTransformation(exptd_avg_array);
    average_array = average.as_array();
    assert(all(all(abs(exptd_avg_array-average_array) < 1e-4)), 'AffineTransformation average failed.')
    disp(average_array)


    disp('% ----------------------------------------------------------------------- %')
    disp('%                  Finished AffineTransformation test.')
    disp('%------------------------------------------------------------------------ %')
end

if try_quaternion
    disp('% ----------------------------------------------------------------------- %')
    disp('%                  Starting Quaternion test...')
    disp('%------------------------------------------------------------------------ %')

    % Construct TM
    array(4,4) =  0;
    array(1,3) =  1;
    array(2,2) =  1;
    array(3,1) = -1;
    array(4,4) =  1;
    rotm = sirf.Reg.AffineTransformation(array);

    % Convert to quaternion
    quat = sirf.Reg.Quaternion(rotm);
    a = quat.as_array();

    % Construct from numpy array
    expt_array = [0.707107, 0., 0.707107, 0.];
    expt = sirf.Reg.Quaternion(expt_array);

    % Compare to expected values
    quat_array = quat.as_array();
    assert(all(abs(quat_array-expt_array)) < 1e-4, 'Quaternion from TM failed.')
    
    % Convert back to TM
    trans_array = [0., 0., 0.];
    affine = sirf.Reg.AffineTransformation(trans_array,quat);
    assert(affine == rotm, 'TM to quaternion failed.');

    % Convert TM to quaternion
    quat2 = affine.get_quaternion();
    quat2_array = quat2.as_array();
    assert(all(abs(quat_array-quat2_array)) < 1e-4, 'AffineTransformation:get_quaternion() failed.')

    % Average!
    quat_1_array = [0.92707,  0.02149,   0.19191,  0.32132];
    quat_2_array = [0.90361,  0.0025836, 0.097279, 0.41716];
    quat_3_array = [0.75868, -0.21289,   0.53263,  0.30884];
    quat_1 = sirf.Reg.Quaternion(quat_1_array);
    quat_2 = sirf.Reg.Quaternion(quat_2_array);
    quat_3 = sirf.Reg.Quaternion(quat_3_array);
    exptd_avg_array = [0.88748, -0.0647152, 0.281671, 0.35896];
    exptd_average = sirf.Reg.Quaternion(exptd_avg_array);
    average = sirf.Reg.Quaternion.get_average([quat_1, quat_2, quat_3]);
    average_array = average.as_array();
    assert(all(abs(exptd_avg_array-average_array) < 1e-4), 'Quaternion average failed.')
    disp(average_array)


    disp('% ----------------------------------------------------------------------- %')
    disp('%                  Finished Quaternion test.')
    disp('%------------------------------------------------------------------------ %')
end