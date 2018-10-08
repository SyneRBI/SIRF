[~] = set_up_Reg([]);
[~] = set_up_PET([]);

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
g.save_nifti_image                           = fullfile(output_prefix, 'matlab_save_NiftiImage.nii');
g.save_nifti_image_3d                        = fullfile(output_prefix, 'matlab_save_NiftiImage3D.nii');
g.save_nifti_image_3d_tensor_not_split       = fullfile(output_prefix, 'matlab_save_NiftiImage3DTensor_not_split.nii');
g.save_nifti_image_3d_tensor_split           = fullfile(output_prefix, 'matlab_save_NiftiImage3DTensor_split_%s.nii');
g.save_nifti_image_3d_deformation_not_split  = fullfile(output_prefix, 'matlab_save_NiftiImage3DDeformation_not_split.nii');
g.save_nifti_image_3d_deformation_split      = fullfile(output_prefix, 'matlab_save_NiftiImage3DDeformation_split_%s.nii');
g.save_nifti_image_3d_displacement_not_split = fullfile(output_prefix, 'matlab_save_NiftiImage3DDisplacement_not_split.nii');
g.save_nifti_image_3d_displacement_split     = fullfile(output_prefix, 'matlab_save_NiftiImage3DDisplacement_split_%s.nii');
g.aladin_warped                              = fullfile(output_prefix, 'matlab_aladin_warped.nii');
g.f3d_warped                                 = fullfile(output_prefix, 'matlab_f3d_warped.nii');
g.TM_fwrd				                     = fullfile(output_prefix, 'matlab_TM_fwrd.txt');
g.TM_back				                     = fullfile(output_prefix, 'matlab_TM_back.txt');
g.aladin_def_fwrd                            = fullfile(output_prefix, 'matlab_aladin_def_fwrd.nii');
g.aladin_def_back                            = fullfile(output_prefix, 'matlab_aladin_def_back_%s.nii');
g.aladin_disp_fwrd                           = fullfile(output_prefix, 'matlab_aladin_disp_fwrd.nii');
g.aladin_disp_back                           = fullfile(output_prefix, 'matlab_aladin_disp_back_%s.nii');
g.f3d_def_fwrd                               = fullfile(output_prefix, 'matlab_f3d_disp_fwrd.nii');
g.f3d_def_back                               = fullfile(output_prefix, 'matlab_f3d_disp_back_%s.nii');
g.f3d_disp_fwrd                              = fullfile(output_prefix, 'matlab_f3d_disp_fwrd.nii');
g.f3d_disp_back                              = fullfile(output_prefix, 'matlab_f3d_disp_back_%s.nii');

g.rigid_resample                             = fullfile(output_prefix, 'matlab_rigid_resample.nii');
g.nonrigid_resample_disp                     = fullfile(output_prefix, 'matlab_nonrigid_resample_disp.nii');
g.nonrigid_resample_def                      = fullfile(output_prefix, 'matlab_nonrigid_resample_def.nii');
g.output_weighted_mean                       = fullfile(output_prefix, 'matlab_weighted_mean.nii');
g.output_weighted_mean_def                   = fullfile(output_prefix, 'matlab_weighted_mean_def.nii');
g.output_float                               = fullfile(output_prefix, 'matlab_reg_aladin_float.nii');

g.ref_aladin                                 = mSIRFReg.NiftiImage3D( g.ref_aladin_filename );
g.flo_aladin                                 = mSIRFReg.NiftiImage3D( g.flo_aladin_filename );
g.ref_f3d                                    = mSIRFReg.NiftiImage3D(   g.ref_f3d_filename  );
g.flo_f3d                                    = mSIRFReg.NiftiImage3D(   g.flo_f3d_filename  );

try_niftiimage(g);
try_niftiimage3d(g);
try_niftiimage3dtensor(g);
try_niftiimage3ddisplacement(g);
try_niftiimage3ddeformation(g);
na = try_niftyaladin(g);
try_niftyf3d(g);
try_transformations(g,na);
try_resample(g,na);
try_weighted_mean(g,na);
try_stir_to_sirfreg(g);

function try_niftiimage(g)
	disp('% ----------------------------------------------------------------------- %')
	disp('%                  Starting NiftiImage test...')
	disp('%------------------------------------------------------------------------ %')

    % default constructor
    a = mSIRFReg.NiftiImage();

    % Read from file
    b = mSIRFReg.NiftiImage(g.ref_aladin_filename);

    % Save to file
    b.save_to_file(g.save_nifti_image);

    % Fill
    b.fill(100);

    % Get max
    assert(b.get_max() == 100, 'NiftiImage fill()/get_max() failed.');

    % Get min
    assert(b.get_min() == 100, 'NiftiImage fill()/get_min() failed.');

    % Deep copy
    d = b.deep_copy();
    assert(d.handle_ ~= b.handle_, 'NiftiImage deep_copy failed.');
    assert(d == b, 'NiftiImage deep_copy failed.');

    % Addition
    e = d + d;
    assert(abs(e.get_max() - 2 * d.get_max()) < 0.0001, 'NiftiImage __add__/get_max() failed.')

    % Subtraction
    e = d - d;
    assert(abs(e.get_max()) < 0.0001, 'NiftiImage __sub__ failed.')

    % Sum
    assert(abs(e.get_sum()) < 0.0001, 'NiftiImage get_sum() failed.')

    % Add num to image
    q = e + 1;
    assert(q.get_max() == e.get_max() + 1, 'NiftiImage __add__ val failed.');

    % Subtract num from image
    r = e - 1;
    assert(r.get_max() == e.get_max() - 1, 'NiftiImage __sub__ val failed.');

    % Multiply image by num
    s = e * 10;
    assert(s.get_max() == e.get_max() * 10, 'NiftiImage __mul__ val failed.');

    % Dimensions
    f = e.get_dimensions();
    assert(all(f == [3, 64, 64, 64, 1, 1, 1, 1]), 'NiftiImage get_dimensions() failed.')

    % Get as array
    arr = d.as_array();
    assert(max(arr(:)) == 100, 'NiftiImage as_array().max() failed.')
    assert(ndims(arr) == 3, 'NiftiImage as_array() ndims failed.')
    assert(all(size(arr) == [64, 64, 64]), 'NiftiImage as_array().shape failed.')

    % Test saving to datatype
    g.ref_aladin.save_to_file(g.output_float,16); % save to float
    ref_aladin_float = mSIRFReg.NiftiImage3D(g.output_float);
    arr1 = g.ref_aladin.as_array();
    arr2 = ref_aladin_float.as_array();
    assert(all(arr1(:)==arr2(:)), "SIRFRegMisc::save_to_file()/change_datatype() failed.");

    % Test print methods
    q.print_header();
    mSIRFReg.NiftiImage.print_headers([q s]);

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
    size(s.as_array())
    assert(all(size(s.as_array()) == [64, 64, 63]), 'NiftiImage crop() failed.')


    disp('% ----------------------------------------------------------------------- %')
    disp('%                  Finished NiftiImage test.')
    disp('%------------------------------------------------------------------------ %')
end

function try_niftiimage3d(g)
    disp('% ----------------------------------------------------------------------- %')
    disp('%                  Starting NiftiImage3D test...')
    disp('%------------------------------------------------------------------------ %')

    % default constructor
    a = mSIRFReg.NiftiImage3D();

    % Read from file
    b = mSIRFReg.NiftiImage3D(g.ref_aladin_filename);

    % Save to file
    b.save_to_file(g.save_nifti_image_3d);

    % Fill
    b.fill(100);

    % Get max
    assert(b.get_max() == 100, 'NiftiImage3D fill()/get_max() failed.');

    % Get min
    assert(b.get_min() == 100, 'NiftiImage3D fill()/get_min() failed.');

    % Deep copy
    d = b.deep_copy();
    assert(d.handle_ ~= b.handle_, 'NiftiImage3D deep_copy failed.');
    assert(d == b, 'NiftiImage3D deep_copy failed.');

    % Addition
    e = d + d;
    assert(abs(e.get_max() - 2 * d.get_max()) < 0.0001, 'NiftiImage3D __add__/get_max() failed.')

    % Subtraction
    e = d - d;
    assert(abs(e.get_max()) < 0.0001, 'NiftiImage3D __sub__ failed.')

    % Sum
    assert(abs(e.get_sum()) < 0.0001, 'NiftiImage3D get_sum() failed.')

    % Dimensions
    f = e.get_dimensions();
    assert(all(f == [3, 64, 64, 64, 1, 1, 1, 1]), 'NiftiImage3D get_dimensions() failed.')

    % Get as array
    arr = d.as_array();
    assert(max(arr(:)) == 100, 'NiftiImage3D as_array().max() failed.')
    assert(ndims(arr) == 3, 'NiftiImage3D as_array() ndims failed.')
    assert(all(size(arr) == [64, 64, 64]), 'NiftiImage3D as_array().shape failed.')

    disp('% ----------------------------------------------------------------------- %')
    disp('%                  Finished NiftiImage3D test.')
    disp('%------------------------------------------------------------------------ %')
end

function try_niftiimage3dtensor(g)
    disp('% ----------------------------------------------------------------------- %')
    disp('%                  Starting NiftiImage3DTensor test...')
    disp('%------------------------------------------------------------------------ %')

    % Create NiftiImage3DTensor from NiftiImage3D
    b = mSIRFReg.NiftiImage3DTensor();
    b.create_from_3D_image(g.ref_aladin);

    % Save to file
    b.save_to_file(g.save_nifti_image_3d_tensor_not_split);
    b.save_to_file_split_xyz_components(g.save_nifti_image_3d_tensor_split);

    % Constructor from file
    c = mSIRFReg.NiftiImage3DTensor(g.save_nifti_image_3d_tensor_not_split);

    % Fill
    c.fill(100)

    % Get max
    assert(c.get_max() == 100, 'NiftiImage3DTensor fill()/get_max() failed.');

    % Get min
    assert(c.get_min() == 100, 'NiftiImage3DTensor fill()/get_min() failed.');

    % Deep copy
    d = c.deep_copy();
    assert(d.handle_ ~= c.handle_, 'NiftiImage3DTensor deep_copy failed (they have the same handle).');
    assert(d == c, 'NiftiImage3DTensor deep_copy failed (values do not match).');

    % Addition
    e = d + d;
    assert(abs(e.get_max() - 2 * d.get_max()) < 0.0001, 'NiftiImage3DTensor __add__/get_max() failed.')

    % Subtraction
    e = d - d;
    assert(abs(e.get_max()) < 0.0001, 'NiftiImage3DTensor __sub__ failed.')

    % Sum
    assert(abs(e.get_sum()) < 0.0001, 'NiftiImage3DTensor get_sum() failed.')

    % Dimensions
    f = e.get_dimensions();
    assert(all(f == [5, 64, 64, 64, 1, 3, 1, 1]), 'NiftiImage3DTensor get_dimensions() failed.')

    % Get as array
    arr = d.as_array();
    assert(max(arr(:)) == 100, 'NiftiImage3DTensor as_array().max() failed.')
    assert(ndims(arr) == 5, 'NiftiImage3DTensor as_array() ndims failed.')
    assert(all(size(arr) == [64, 64, 64, 1, 3]), 'NiftiImage3DTensor as_array().shape failed.')

    % Constructor from single components
    im1 = g.ref_aladin.deep_copy();
    im2 = g.ref_aladin.deep_copy();
    im3 = g.ref_aladin.deep_copy();
    im1.fill(30);
    im2.fill(20);
    im3.fill(-10);
    h = mSIRFReg.NiftiImage3DTensor(im1, im2, im3);

    % Test flip components
    h.flip_component(0);
    assert(h.get_max() ==  20, "NiftiImage3DTensor flip_component() failed.");
    assert(h.get_min() == -30, "NiftiImage3DTensor flip_component() failed.");


    disp('% ----------------------------------------------------------------------- %')
    disp('%                  Finished NiftiImage3DTensor test.')
    disp('%------------------------------------------------------------------------ %')
end

function try_niftiimage3ddisplacement(g)
    disp('% ----------------------------------------------------------------------- %')
    disp('%                  Starting NiftiImage3DDisplacement test...')
    disp('%------------------------------------------------------------------------ %')

    % Create NiftiImage3DDisplacement from NiftiImage3D
    b = mSIRFReg.NiftiImage3DDisplacement();
    b.create_from_3D_image(g.ref_aladin);

    % Save to file
    b.save_to_file(g.save_nifti_image_3d_displacement_not_split);
    b.save_to_file_split_xyz_components(g.save_nifti_image_3d_displacement_split);

    % Constructor from file
    c = mSIRFReg.NiftiImage3DDisplacement(g.save_nifti_image_3d_displacement_not_split);

    % Constructor from 3x3D
    h = mSIRFReg.NiftiImage3DDisplacement(g.ref_aladin, g.ref_aladin, g.ref_aladin);

    % Fill
    c.fill(100)

    % Get max
    assert(c.get_max() == 100, 'NiftiImage3DDisplacement fill()/get_max() failed.');

    % Get min
    assert(c.get_min() == 100, 'NiftiImage3DDisplacement fill()/get_min() failed.');

    % Deep copy
    d = c.deep_copy();
    assert(d.handle_ ~= c.handle_, 'NiftiImage3DDisplacement deep_copy failed (they have the same handle).');
    assert(d == c, 'NiftiImage3DDisplacement deep_copy failed (values do not match).');

    % Addition
    e = d + d;
    assert(abs(e.get_max() - 2 * d.get_max()) < 0.0001, 'NiftiImage3DDisplacement __add__/get_max() failed.')

    % Subtraction
    e = d - d;
    assert(abs(e.get_max()) < 0.0001, 'NiftiImage3DDisplacement __sub__ failed.')

    % Sum
    assert(abs(e.get_sum()) < 0.0001, 'NiftiImage3DDisplacement get_sum() failed.')

    % Dimensions
    f = e.get_dimensions();
    assert(all(f == [5, 64, 64, 64, 1, 3, 1, 1]), 'NiftiImage3DDisplacement get_dimensions() failed.')

    % Get as array
    arr = d.as_array();
    assert(max(arr(:)) == 100, 'NiftiImage3DDisplacement as_array().max() failed.')
    assert(ndims(arr) == 5, 'NiftiImage3DDisplacement as_array() ndims failed.')
    assert(all(size(arr) == [64, 64, 64, 1, 3]), 'NiftiImage3DDisplacement as_array().shape failed.')

    disp('% ----------------------------------------------------------------------- %')
    disp('%                  Finished NiftiImage3DDisplacement test.')
    disp('%------------------------------------------------------------------------ %')
end

function try_niftiimage3ddeformation(g)
    disp('% ----------------------------------------------------------------------- %')
    disp('%                  Starting NiftiImage3DDeformation test...')
    disp('%------------------------------------------------------------------------ %')

    % Create NiftiImage3DDeformation from NiftiImage3D
    b = mSIRFReg.NiftiImage3DDeformation();
    b.create_from_3D_image(g.ref_aladin);

    % Save to file
    b.save_to_file(g.save_nifti_image_3d_deformation_not_split);
    b.save_to_file_split_xyz_components(g.save_nifti_image_3d_deformation_split);

    % Constructor from file
    c = mSIRFReg.NiftiImage3DDeformation(g.save_nifti_image_3d_deformation_not_split);

    % Constructor from 3x3D
    h = mSIRFReg.NiftiImage3DDeformation(g.ref_aladin, g.ref_aladin, g.ref_aladin);

    % Fill
    c.fill(100)

    % Get max
    assert(c.get_max() == 100, 'NiftiImage3DDeformation fill()/get_max() failed.');

    % Get min
    assert(c.get_min() == 100, 'NiftiImage3DDeformation fill()/get_min() failed.');

    % Deep copy
    d = c.deep_copy();
    assert(d.handle_ ~= c.handle_, 'NiftiImage3DDeformation deep_copy failed (they have the same handle).');
    assert(d == c, 'NiftiImage3DDeformation deep_copy failed (values do not match).');

    % Addition
    e = d + d;
    assert(abs(e.get_max() - 2 * d.get_max()) < 0.0001, 'NiftiImage3DDeformation __add__/get_max() failed.')

    % Subtraction
    e = d - d;
    assert(abs(e.get_max()) < 0.0001, 'NiftiImage3DDeformation __sub__ failed.')

    % Sum
    assert(abs(e.get_sum()) < 0.0001, 'NiftiImage3DDeformation get_sum() failed.')

    % Dimensions
    f = e.get_dimensions();
    assert(all(f == [5, 64, 64, 64, 1, 3, 1, 1]), 'NiftiImage3DDeformation get_dimensions() failed.')

    % Get as array
    arr = d.as_array();
    assert(max(arr(:)) == 100, 'NiftiImage3DDeformation as_array().max() failed.')
    assert(ndims(arr) == 5, 'NiftiImage3DDeformation as_array() ndims failed.')
    assert(all(size(arr) == [64, 64, 64, 1, 3]), 'NiftiImage3DDeformation as_array().shape failed.')

    disp('% ----------------------------------------------------------------------- %')
    disp('%                  Finished NiftiImage3DDeformation test.')
    disp('%------------------------------------------------------------------------ %')
end

function na =try_niftyaladin(g)
	disp('% ----------------------------------------------------------------------- %')
	disp('%                  Starting Nifty aladin test...')
	disp('%------------------------------------------------------------------------ %')

	% default constructor
    na = mSIRFReg.NiftyAladinSym();
    na.set_reference_image(g.ref_aladin);
    na.set_floating_image(g.flo_aladin);
    na.set_parameter_file(g.parameter_file_aladin);
    na.set_parameter('SetInterpolationToCubic');
    na.set_parameter('SetLevelsToPerform', '1');
    na.set_parameter('SetMaxIterations', '5');
    na.update();

    % Get outputs
    warped = na.get_output();
    def_fwrd = na.get_deformation_field_fwrd();
    def_back = na.get_deformation_field_back();
    disp_fwrd = na.get_displacement_field_fwrd();
    disp_back = na.get_displacement_field_back();

    warped.save_to_file(g.aladin_warped);
    na.get_transformation_matrix_fwrd().save_to_file(g.TM_fwrd);
    na.get_transformation_matrix_back().save_to_file(g.TM_back);
    def_fwrd.save_to_file(g.aladin_def_fwrd);
    def_back.save_to_file_split_xyz_components(g.aladin_def_back);
    disp_fwrd.save_to_file(g.aladin_disp_fwrd);
    disp_back.save_to_file_split_xyz_components(g.aladin_disp_back);

    % Fwrd TM
    fwrd_tm = na.get_transformation_matrix_fwrd().as_array()

    % Back TM
    back_tm = na.get_transformation_matrix_back().as_array()

    % Test converting disp to def
    a = mSIRFReg.NiftiImage3DDeformation();
    a.create_from_disp(disp_fwrd);
    assert(a == def_fwrd, "NiftiImage3DDeformation::create_from_disp() failed.");

    % Test converting def to disp
    b = mSIRFReg.NiftiImage3DDisplacement();
    b.create_from_def(def_fwrd);
    assert(b == disp_fwrd, "NiftiImage3DDisplacement::create_from_def() failed.");

	disp('% ----------------------------------------------------------------------- %')
	disp('%                  Finished Nifty aladin test.')
	disp('%------------------------------------------------------------------------ %')
end

function try_niftyf3d(g)
	disp('% ----------------------------------------------------------------------- %')
	disp('%                  Starting Nifty f3d test...')
	disp('%------------------------------------------------------------------------ %')

    % Get initial transformation
    tm_init = mSIRFReg.Mat44(g.TM_fwrd);

	% default constructor
    nf = mSIRFReg.NiftyF3dSym();
    nf.set_reference_image(g.ref_f3d);
    nf.set_floating_image(g.flo_f3d);
    nf.set_parameter_file(g.parameter_file_f3d);
    nf.set_reference_time_point(1);
    nf.set_floating_time_point(1);
    nf.set_initial_affine_transformation(tm_init);
    nf.update();

    % Get outputs
    warped = nf.get_output();
    def_fwrd = nf.get_deformation_field_fwrd();
    def_back = nf.get_deformation_field_back();
    disp_fwrd = nf.get_displacement_field_fwrd();
    disp_back = nf.get_displacement_field_back();

    warped.save_to_file(g.f3d_warped);
    def_fwrd.save_to_file(g.f3d_def_fwrd);
    def_back.save_to_file_split_xyz_components(g.f3d_def_back);
    disp_fwrd.save_to_file(g.f3d_disp_fwrd);
    disp_back.save_to_file_split_xyz_components(g.f3d_disp_back);

	disp('% ----------------------------------------------------------------------- %')
	disp('%                  Finished Nifty f3d test.')
	disp('%------------------------------------------------------------------------ %')
end

function try_transformations(g,na)
	disp('% ----------------------------------------------------------------------- %')
	disp('%                  Starting Transformation test...')
	disp('%------------------------------------------------------------------------ %')


    % Get transformations
    a3 = na.get_transformation_matrix_fwrd();
    b3 = na.get_displacement_field_fwrd();
    c3 = na.get_deformation_field_fwrd();

    % Get as deformations
    a_def = a3.get_as_deformation_field(g.ref_aladin);
    b_def = b3.get_as_deformation_field(g.ref_aladin);
    c_def = c3.get_as_deformation_field(g.ref_aladin);
    assert(a_def == na.get_deformation_field_fwrd(), 'SIRFRegTransformationAffine get_as_deformation_field() failed.')
    assert(b_def == na.get_deformation_field_fwrd(), 'SIRFRegTransformationDisplacement get_as_deformation_field() failed.')
    assert(c_def == na.get_deformation_field_fwrd(), 'SIRFRegTransformationDeformation get_as_deformation_field() failed.')

    % Compose into single deformation. Use two identity matrices and the disp field. Get as def and should be the same.
    tm_iden = mSIRFReg.Mat44.get_identity();
    trans = [tm_iden, tm_iden, c3];
    composed = mSIRFReg.NiftiImage3DDeformation.compose_single_deformation(trans, g.ref_aladin);
    assert(composed == na.get_deformation_field_fwrd(), 'compose_single_deformation failed.')


	disp('% ----------------------------------------------------------------------- %')
	disp('%                  Finished Transformation test.')
	disp('%------------------------------------------------------------------------ %')
end

function try_resample(g,na)
    disp('% ----------------------------------------------------------------------- %')
    disp('%                  Starting Nifty resample test...')
    disp('%------------------------------------------------------------------------ %')

	tm_iden = mSIRFReg.Mat44.get_identity();
    tm      = na.get_transformation_matrix_fwrd();
    displ   = na.get_displacement_field_fwrd();
    deff    = na.get_deformation_field_fwrd();

    disp('Testing rigid resample...')
    nr1 = mSIRFReg.NiftyResample();
    nr1.set_reference_image(g.ref_aladin);
    nr1.set_floating_image(g.flo_aladin);
    nr1.set_interpolation_type_to_cubic_spline();  % try different interpolations
    nr1.set_interpolation_type(3);  % try different interpolations (cubic)
    nr1.add_transformation_affine(tm_iden);
		nr1.add_transformation_affine(tm);
    nr1.update();
    nr1.get_output().save_to_file(g.rigid_resample);

    disp('Testing non-rigid displacement...')
    nr2 = mSIRFReg.NiftyResample();
    nr2.set_reference_image(g.ref_aladin);
    nr2.set_floating_image(g.flo_aladin);
    nr2.set_interpolation_type_to_sinc();  % try different interpolations
    nr2.set_interpolation_type_to_linear();  % try different interpolations
    nr2.add_transformation_disp(displ);
    nr2.update();
    nr2.get_output().save_to_file(g.nonrigid_resample_disp);

    disp('Testing non-rigid deformation...')
    nr3 = mSIRFReg.NiftyResample();
    nr3.set_reference_image(g.ref_aladin)
    nr3.set_floating_image(g.flo_aladin)
    nr3.set_interpolation_type_to_nearest_neighbour()  % try different interpolations
    nr3.add_transformation_def(deff);
    nr3.set_interpolation_type_to_linear()
    nr3.update()
    nr3.get_output().save_to_file(g.nonrigid_resample_def)

    assert(na.get_output() == nr1.get_output(), 'Rigid resampled output should match registration (aladin) output.')

    disp('% ----------------------------------------------------------------------- %')
    disp('%                  Finished Nifty resample test.')
    disp('%------------------------------------------------------------------------ %')
end

function try_weighted_mean(g,na)
    disp('% ----------------------------------------------------------------------- %')
    disp('%                  Starting weighted mean test...')
    disp('%------------------------------------------------------------------------ %')

		% Do 3D
		wm1 = mSIRFReg.ImageWeightedMean();
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
		wm1.update();
		wm1.get_output().save_to_file(g.output_weighted_mean);
		% Answer should be 4.5, so compare it to that!
        res = g.ref_aladin.deep_copy();
		res.fill(4.5);
		assert(wm1.get_output() == res, '3D weighted mean test failed.')

		% Do 4D
		wm2 = mSIRFReg.ImageWeightedMean();
		im1 = na.get_deformation_field_fwrd().deep_copy();
		im2 = na.get_deformation_field_fwrd().deep_copy();
		im3 = na.get_deformation_field_fwrd().deep_copy();
		im4 = na.get_deformation_field_fwrd().deep_copy();
		im1.fill(1);
		im2.fill(4);
		im3.fill(7);
		im4.fill(6);
		wm2.add_image(im1, 2);
		wm2.add_image(im2, 4);
		wm2.add_image(im3, 3);
		wm2.add_image(im4, 1);
		wm2.update();
		wm2.get_output().save_to_file(g.output_weighted_mean_def);
		% Answer should be 4.5, so compare it to that!
		res = na.get_deformation_field_fwrd().deep_copy();
		res.fill(4.5);
		assert(wm2.get_output() == res, '4D weighted mean test failed.')


    disp('% ----------------------------------------------------------------------- %')
    disp('%                  Finished weighted mean test.')
    disp('%------------------------------------------------------------------------ %')
end

function try_stir_to_sirfreg(g)
    disp('% ----------------------------------------------------------------------- %')
    disp('%                  Starting STIR to SIRFReg test...')
    disp('%------------------------------------------------------------------------ %')

		% Open stir image
                pet_image_data = mSTIR.ImageData(g.ref_aladin_filename);
		image_data_from_stir = mSIRFReg.NiftiImage3D(pet_image_data);

		% Now fill the stir and sirfreg images with 1 and 100, respectively
		pet_image_data.fill(1.);
		image_data_from_stir.fill(100);
		arr_pet = pet_image_data.as_array();
		assert(max(arr_pet(:)) ~= image_data_from_stir.get_max(), 'Maxes of STIR and SIRFReg images should be different.');

		% Fill the stir image with the sirfreg
		image_data_from_stir.copy_data_to(pet_image_data);
		arr_pet = pet_image_data.as_array();
		assert(max(arr_pet(:)) == image_data_from_stir.get_max(), 'Maxes of STIR and SIRFReg images should match.');


    disp('% ----------------------------------------------------------------------- %')
    disp('%                  Finished STIR to SIRFReg test.')
    disp('%------------------------------------------------------------------------ %')
end

function try_sirfregmat44(g,na)
    disp('% ----------------------------------------------------------------------- %')
    disp('%                  Starting SIRFRegMat44 test...')
    disp('%------------------------------------------------------------------------ %')

    % Construct from file
    a = mSIRFReg.Mat44(TM_fwrd);

    % Multiply fwrd and inverse, should equal identity
    b = na.get_transformation_matrix_fwrd();
    c = na.get_transformation_matrix_back();
    d = b * c;
    e = mSIRFReg.Mat44.get_identity();
    assert(d == e, 'SIRFRegMat44::mult/comparison failed.');

    d.fill(3);
    f = d.as_array();
    assert(np.all(f == 3), 'SIRFRegMat44::fill/operator[] failed.');

    assert(d.get_determinant() < 1.e-7, 'SIRFRegMat44::get_determinant failed.');
    assert(e.get_determinant() - 1. < 1.e-7, 'SIRFRegMat44::get_determinant failed.');

    disp('% ----------------------------------------------------------------------- %')
    disp('%                  Finished SIRFRegMat44 test.')
    disp('%------------------------------------------------------------------------ %')
end
