[~] = set_up_Reg([]);
[~] = set_up_PET([]);

% Paths
SIRF_PATH     = getenv('SIRF_PATH');
examples_path = [SIRF_PATH  '/data/examples/Registration'];
output_path   = [pwd  '/results/matlab_'];

% Input filenames
g.ref_aladin_filename      = [examples_path  '/test.nii.gz'];
g.flo_aladin_filename      = [examples_path  '/test2.nii.gz'];
g.ref_f3d_filename         = [examples_path  '/mouseFixed.nii.gz'];
g.flo_f3d_filename         = [examples_path  '/mouseMoving.nii.gz'];
g.parameter_file_aladin    = [examples_path  '/paramFiles/aladin.par'];
g.parameter_file_f3d       = [examples_path  '/paramFiles/f3d.par'];
g.matrix                   = [examples_path  '/transformation_matrix.txt'];
g.stir_nifti               = [examples_path  '/nifti_created_by_stir.nii'];

% Output filenames
g.img_data                 = [output_path    'save_SIRFImageData'];
g.img_data_def_not_split   = [output_path    'save_SIRFImageDataDeformation_not_split'];
g.img_data_def_split       = [output_path    'save_SIRFImageDataDeformation_split'];
g.aladin_warped            = [output_path    'aladin_warped'];
g.f3d_warped               = [output_path    'f3d_warped'];
g.TM_fwrd				   = [output_path    'TM_fwrd.txt'];
g.TM_back				   = [output_path    'TM_back.txt'];
g.aladin_def_fwrd          = [output_path    'aladin_def_fwrd'];
g.aladin_def_back          = [output_path    'aladin_def_back'];
g.aladin_disp_fwrd         = [output_path    'aladin_disp_fwrd'];
g.aladin_disp_back         = [output_path    'aladin_disp_back'];
g.f3d_def_fwrd             = [output_path    'f3d_disp_fwrd'];
g.f3d_def_back             = [output_path    'f3d_disp_back'];
g.f3d_disp_fwrd            = [output_path    'f3d_disp_fwrd'];
g.f3d_disp_back            = [output_path    'f3d_disp_back'];

g.rigid_resample           = [output_path    'rigid_resample'];
g.nonrigid_resample_disp   = [output_path    'nonrigid_resample_disp'];
g.nonrigid_resample_def    = [output_path    'nonrigid_resample_def'];
g.output_weighted_mean     = [output_path    'weighted_mean'];
g.output_weighted_mean_def = [output_path    'weighted_mean_def'];

g.output_stir_nifti        = [output_path    'stir_nifti.nii'];

g.ref_aladin = mSIRFReg.ImageData( g.ref_aladin_filename );
g.flo_aladin = mSIRFReg.ImageData( g.flo_aladin_filename );
g.ref_f3d    = mSIRFReg.ImageData(   g.ref_f3d_filename  );
g.flo_f3d    = mSIRFReg.ImageData(   g.flo_f3d_filename  );
g.nifti      = mSIRFReg.ImageData(        g.stir_nifti   );

g.required_percentage_accuracy = single(1);

try_misc_functions(g);
try_sirfimagedata(g);
try_sirfimagedatadeformation(g);
na = try_niftyaladin(g);
try_niftyf3d(g);
try_transformations(g,na);
try_resample(g,na);
try_weighted_mean(g,na);
try_stir_to_sirfreg(g);

function try_misc_functions(g)
	disp('% ----------------------------------------------------------------------- %')
	disp('%                  Starting misc functions test...')
	disp('%------------------------------------------------------------------------ %')

    % do nifti images match?
    assert(mSIRFReg.Misc.do_nifti_images_match(g.ref_aladin, g.ref_aladin, g.required_percentage_accuracy) == 1, "Images don\'t match, but they should.")
    assert(mSIRFReg.Misc.do_nifti_images_match(g.ref_aladin, g.flo_aladin, g.required_percentage_accuracy) == 0, "Images match, but they shouldn't.")

    % dump from filename
    mSIRFReg.Misc.dump_nifti_info(g.ref_aladin_filename)
    % dump from SIRFImageData
    mSIRFReg.Misc.dump_nifti_info(g.ref_aladin)
    % dump from multiple images
    mSIRFReg.Misc.dump_nifti_info([g.ref_aladin, g.flo_aladin, g.nifti])
    % dump from SIRFImageDataDeformation
    deform = mSIRFReg.ImageDataDeformation();
    deform.create_from_3D_image(g.ref_aladin);
    mSIRFReg.Misc.dump_nifti_info(deform)

    % identity matrix
    tm_iden = mSIRFReg.Misc.get_matrix();
    disp(tm_iden)

	disp('% ----------------------------------------------------------------------- %')
	disp('%                  Finished misc functions test.')
	disp('%------------------------------------------------------------------------ %')
end

function try_sirfimagedata(g)
	disp('% ----------------------------------------------------------------------- %')
	disp('%                  Starting SIRFImageData test...')
	disp('%------------------------------------------------------------------------ %')

    % default constructor
    a = mSIRFReg.ImageData();

    % Read from file
    b = mSIRFReg.ImageData(g.ref_aladin_filename);

    % Construct from stir mSTIR.ImageData
    stir = mSTIR.ImageData(g.stir_nifti);
    c = mSIRFReg.ImageData(stir);

    % Save to file
    c.save_to_file(g.img_data);

    % Fill
    c.fill(100);

    % Get max
    assert(c.get_max() == 100, 'SIRFImageData fill()/get_max() failed.');

    % Get min
    assert(c.get_min() == 100, 'SIRFImageData fill()/get_min() failed.');

    % Copy data to mSTIRImageData
    stir.fill(3.);
    c.copy_data_to(stir);
    arr = stir.as_array();
    assert(max(arr(:)) == 100, 'SIRFImageData copy_data_to stir ImageData failed.');

    % Deep copy
    d = c.deep_copy();
    assert(d.handle_ ~= c.handle_, 'SIRFImageData deep_copy failed.');
    assert(mSIRFReg.Misc.do_nifti_images_match(d, c, g.required_percentage_accuracy) == 1, 'SIRFImageData deep_copy failed.');

    % Addition
    e = d + d;
    assert(abs(e.get_max() - 2 * d.get_max()) < 0.0001, 'SIRFImageData __add__/get_max() failed.')

    % Subtraction
    e = d - d;
    assert(abs(e.get_max()) < 0.0001, 'SIRFImageData __sub__ failed.')

    % Sum
    assert(abs(e.get_sum()) < 0.0001, 'SIRFImageData get_sum() failed.')

    % Dimensions
    f = e.get_dimensions();
    assert(all(f == [3, 285, 285, 127, 1, 1, 1, 1]), 'SIRFImageData get_dimensions() failed.')

    % Get as array
    arr = d.as_array();
    assert(max(arr(:)) == 100, 'SIRFImageData as_array().max() failed.')
    assert(ndims(arr) == 3, 'SIRFImageData as_array() ndims failed.')
    assert(all(size(arr) == [285, 285, 127]), 'SIRFImageData as_array().shape failed.')

	disp('% ----------------------------------------------------------------------- %')
	disp('%                  Finished SIRFImageData test.')
	disp('%------------------------------------------------------------------------ %')
end

function try_sirfimagedatadeformation(g)
	disp('% ----------------------------------------------------------------------- %')
	disp('%                  Starting SIRFImageDataDeformation test...')
	disp('%------------------------------------------------------------------------ %')

    % Create SIRFImageDataDeformation from SIRFImageData
    b = mSIRFReg.ImageDataDeformation();
    b.create_from_3D_image(g.ref_aladin);

    % Save to file
    b.save_to_file(g.img_data_def_not_split);
    b.save_to_file_split_xyz_components(g.img_data_def_split);

    % Constructor from file
    c = mSIRFReg.ImageDataDeformation([g.img_data_def_not_split '.nii']);

    % Fill
    c.fill(100)

    % Get max
    assert(c.get_max() == 100, 'SIRFImageDataDeformation fill()/get_max() failed.');

    % Get min
    assert(c.get_min() == 100, 'SIRFImageDataDeformation fill()/get_min() failed.');

    % Deep copy
    d = c.deep_copy();
    assert(d.handle_ ~= c.handle_, 'SIRFImageDataDeformation deep_copy failed (they have the same handle).');
    assert(mSIRFReg.Misc.do_nifti_images_match(d, c, g.required_percentage_accuracy) == 1, 'SIRFImageDataDeformation deep_copy failed (values do not match).');

    % Addition
    e = d + d;
    assert(abs(e.get_max() - 2 * d.get_max()) < 0.0001, 'SIRFImageDataDeformation __add__/get_max() failed.')

    % Subtraction
    e = d - d;
    assert(abs(e.get_max()) < 0.0001, 'SIRFImageDataDeformation __sub__ failed.')

    % Sum
    assert(abs(e.get_sum()) < 0.0001, 'SIRFImageDataDeformation get_sum() failed.')

    % Dimensions
    f = e.get_dimensions();
    assert(all(f == [5, 64, 64, 64, 1, 3, 1, 1]), 'SIRFImageDataDeformation get_dimensions() failed.')

    % Get as array
    arr = d.as_array();
    assert(max(arr(:)) == 100, 'SIRFImageDataDeformation as_array().max() failed.')
    assert(ndims(arr) == 5, 'SIRFImageDataDeformation as_array() ndims failed.')
    assert(all(size(arr) == [64, 64, 64, 1, 3]), 'SIRFImageDataDeformation as_array().shape failed.')

	disp('% ----------------------------------------------------------------------- %')
	disp('%                  Finished SIRFImageDataDeformation test.')
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
    na.update();

    % Get outputs
    warped = na.get_output();
    def_fwrd = na.get_deformation_field_fwrd();
    def_back = na.get_deformation_field_back();
    disp_fwrd = na.get_displacement_field_fwrd();
    disp_back = na.get_displacement_field_back();

    warped.save_to_file(g.aladin_warped);
    na.save_transformation_matrix_fwrd(g.TM_fwrd);
    na.save_transformation_matrix_back(g.TM_back);
    def_fwrd.save_to_file(g.aladin_def_fwrd);
    def_back.save_to_file_split_xyz_components(g.aladin_def_back);
    disp_fwrd.save_to_file(g.aladin_disp_fwrd);
    disp_back.save_to_file_split_xyz_components(g.aladin_disp_back);

    % Fwrd TM
    fwrd_tm = na.get_transformation_matrix_fwrd()
    disp('Fwrd tm:')
    disp(fwrd_tm)

    % Back TM
    back_tm = na.get_transformation_matrix_back()
    disp('Back tm:')
    disp(back_tm)

	disp('% ----------------------------------------------------------------------- %')
	disp('%                  Finished Nifty aladin test.')
	disp('%------------------------------------------------------------------------ %')
end

function try_niftyf3d(g)
	disp('% ----------------------------------------------------------------------- %')
	disp('%                  Starting Nifty f3d test...')
	disp('%------------------------------------------------------------------------ %')

	% default constructor
    nf = mSIRFReg.NiftyF3dSym();
    nf.set_reference_image(g.ref_f3d);
    nf.set_floating_image(g.flo_f3d);
    nf.set_parameter_file(g.parameter_file_f3d);
    nf.set_reference_time_point(1);
    nf.set_floating_time_point(1);
    nf.set_initial_affine_transformation(g.TM_fwrd);
    nf.update();

    % Get outputs
    warped = nf.get_output();
    def_fwrd = nf.get_deformation_field_fwrd();
    def_back = nf.get_deformation_field_back();
    disp_fwrd = nf.get_displacement_field_fwrd();
    disp_back = nf.get_displacement_field_back();

    warped.save_to_file(g.f3d_warped);
    def_fwrd.save_to_file_split_xyz_components(g.f3d_def_fwrd);
    def_back.save_to_file(g.f3d_def_back);
    disp_fwrd.save_to_file_split_xyz_components(g.f3d_disp_fwrd);
    disp_back.save_to_file(g.f3d_disp_back);

	disp('% ----------------------------------------------------------------------- %')
	disp('%                  Finished Nifty f3d test.')
	disp('%------------------------------------------------------------------------ %')
end

function try_transformations(g,na)
	disp('% ----------------------------------------------------------------------- %')
	disp('%                  Starting Transformation test...')
	disp('%------------------------------------------------------------------------ %')


    % Affine
    disp('Testing affine...')
    a1 = mSIRFReg.TransformationAffine();
    a2 = mSIRFReg.TransformationAffine(g.TM_fwrd);
    a3 = mSIRFReg.TransformationAffine(na.get_transformation_matrix_fwrd());

    % Displacement
    disp('Testing displacement...')
    b1 = mSIRFReg.TransformationDisplacement();
    b2 = mSIRFReg.TransformationDisplacement([g.aladin_disp_fwrd '.nii']);
    b3 = mSIRFReg.TransformationDisplacement(na.get_displacement_field_fwrd());

    % Deformation
    disp('Testing deformation...')
    c1 = mSIRFReg.TransformationDeformation();
    c2 = mSIRFReg.TransformationDeformation([g.aladin_def_fwrd '.nii']);
    c3 = mSIRFReg.TransformationDeformation(na.get_deformation_field_fwrd());

    % Get as deformations
    a_def = a3.get_as_deformation_field(g.ref_aladin);
    b_def = b3.get_as_deformation_field(g.ref_aladin);
    c_def = c3.get_as_deformation_field(g.ref_aladin);
    assert(mSIRFReg.Misc.do_nifti_images_match(a_def, na.get_deformation_field_fwrd(), g.required_percentage_accuracy)==1, 'SIRFRegTransformationAffine get_as_deformation_field() failed.')
    assert(mSIRFReg.Misc.do_nifti_images_match(b_def, na.get_deformation_field_fwrd(), g.required_percentage_accuracy)==1, 'SIRFRegTransformationDisplacement get_as_deformation_field() failed.')
    assert(mSIRFReg.Misc.do_nifti_images_match(c_def, na.get_deformation_field_fwrd(), g.required_percentage_accuracy)==1, 'SIRFRegTransformationDeformation get_as_deformation_field() failed.')

    % Compose into single deformation. Use two identity matrices and the disp field. Get as def and should be the same.
    tm_iden = mSIRFReg.Misc.get_matrix();
    trans_aff_iden = mSIRFReg.TransformationAffine(tm_iden);
    trans = [trans_aff_iden, trans_aff_iden, c3];
    composed = mSIRFReg.Misc.compose_transformations_into_single_deformation(trans, g.ref_aladin);
    assert(mSIRFReg.Misc.do_nifti_images_match(composed.get_as_deformation_field(g.ref_aladin), na.get_deformation_field_fwrd(), g.required_percentage_accuracy) == 1, 'compose_transformations_into_single_deformation failed.')


	disp('% ----------------------------------------------------------------------- %')
	disp('%                  Finished Transformation test.')
	disp('%------------------------------------------------------------------------ %')
end

function try_resample(g,na)
    disp('% ----------------------------------------------------------------------- %')
    disp('%                  Starting Nifty resample test...')
    disp('%------------------------------------------------------------------------ %')

		tm_iden = mSIRFReg.Misc.get_matrix();
		tm_iden = mSIRFReg.TransformationAffine(tm_iden);
    tm   = mSIRFReg.TransformationAffine(na.get_transformation_matrix_fwrd());
    displ = mSIRFReg.TransformationDisplacement(na.get_displacement_field_fwrd());
    deff = mSIRFReg.TransformationDeformation(na.get_deformation_field_fwrd());

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
    nr3 = mSIRFReg.NiftyResample()
    nr3.set_reference_image(g.ref_aladin)
    nr3.set_floating_image(g.flo_aladin)
    nr3.set_interpolation_type_to_nearest_neighbour()  % try different interpolations
    nr3.add_transformation_def(deff);
    nr3.set_interpolation_type_to_linear()
    nr3.update()
    nr3.get_output().save_to_file(g.nonrigid_resample_def)

    assert(mSIRFReg.Misc.do_nifti_images_match(na.get_output(), nr1.get_output(), g.required_percentage_accuracy) == 1, 'Rigid resampled output should match registration (aladin) output.')

    disp('% ----------------------------------------------------------------------- %')
    disp('%                  Finished Nifty resample test.')
    disp('%------------------------------------------------------------------------ %')
end

function try_weighted_mean(g,na)
    disp('% ----------------------------------------------------------------------- %')
    disp('%                  Starting weighted mean test...')
    disp('%------------------------------------------------------------------------ %')

		% Do 3D
		wm1 = mSIRFReg.ImageWeightedMean3D();
		im1 = mSIRFReg.ImageData(g.stir_nifti);
		im2 = mSIRFReg.ImageData(g.stir_nifti);
		im3 = mSIRFReg.ImageData(g.stir_nifti);
		im4 = mSIRFReg.ImageData(g.stir_nifti);
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
		res = mSIRFReg.ImageData(g.stir_nifti);
		res.fill(4.5);
		assert(mSIRFReg.Misc.do_nifti_images_match(wm1.get_output(), res, g.required_percentage_accuracy) == 1, '3D weighted mean test failed.')

		% Do 4D
		wm2 = mSIRFReg.ImageWeightedMean4D()
		im1 = na.get_deformation_field_fwrd().deep_copy()
		im2 = na.get_deformation_field_fwrd().deep_copy()
		im3 = na.get_deformation_field_fwrd().deep_copy()
		im4 = na.get_deformation_field_fwrd().deep_copy()
		im1.fill(1)
		im2.fill(4)
		im3.fill(7)
		im4.fill(6)
		wm2.add_image(im1, 2)
		wm2.add_image(im2, 4)
		wm2.add_image(im3, 3)
		wm2.add_image(im4, 1)
		wm2.update()
		wm2.get_output().save_to_file(g.output_weighted_mean_def)
		% Answer should be 4.5, so compare it to that!
		res = na.get_deformation_field_fwrd().deep_copy()
		res.fill(4.5)
		assert(mSIRFReg.Misc.do_nifti_images_match(wm2.get_output(), res, g.required_percentage_accuracy) == 1, '4D weighted mean test failed.')


    disp('% ----------------------------------------------------------------------- %')
    disp('%                  Finished weighted mean test.')
    disp('%------------------------------------------------------------------------ %')
end

function try_stir_to_sirfreg(g)
    disp('% ----------------------------------------------------------------------- %')
    disp('%                  Starting STIR to SIRFReg test...')
    disp('%------------------------------------------------------------------------ %')

		% Open stir image
		pet_image_data = mSTIR.ImageData(g.stir_nifti);
		image_data_from_stir = mSIRFReg.ImageData(pet_image_data);

		% Compare to nifti IO (if they don't match, you'll see a message but don't throw an error for now)
		image_data_from_nifti = mSIRFReg.ImageData(g.stir_nifti);
		mSIRFReg.Misc.do_nifti_images_match(image_data_from_stir, image_data_from_nifti, g.required_percentage_accuracy);

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
