% equivalent of demo5 using python interpreter

PYPATH = py.sys.path;
PGPATH = [getenv('SRC_PATH') '/xGadgetron/pGadgetron'];
if count(PYPATH, PGPATH) == 0
    insert(PYPATH, int32(0), PGPATH);
end

try
    % acquisitions will be read from this HDF file
    file = input('raw data file: ', 's');
    input_data = py.pGadgetron.MR_Acquisitions(file);

    prep_gadgets = py.list({'NoiseAdjustGadget', 'AsymmetricEchoGadget', ...
         'RemoveROOversamplingGadget'});
    processed_data = input_data.process(prep_gadgets);

    fprintf('---\n processed acquisition data norm: %e\n', processed_data.norm());
    
    recon = py.pGadgetron.MR_BasicReconstruction();
    recon.set_input(processed_data);
    recon.process();
    complex_images = recon.get_output();

    fprintf('---\n reconstructed images norm: %e\n', complex_images.norm());

    csms = py.pGadgetron.MR_CoilSensitivityMaps();

    fprintf('---\n sorting acquisitions...');
    processed_data.sort();
    fprintf('ok\n');
    fprintf('---\n computing sensitivity maps...');
    csms.calculate(processed_data);
    fprintf('ok\n');

    am = py.pGadgetron.MR_AcquisitionModel(processed_data, complex_images);
    am.set_coil_sensitivity_maps(csms);
    fwd_acqs = am.forward(complex_images);
    fprintf('---\n reconstructed images forward projection norm %e\n', fwd_acqs.norm());

    acq_diff = fwd_acqs - processed_data;
    rr = acq_diff.norm()/fwd_acqs.norm();
    fprintf('---\n reconstruction residual norm (rel): %e\n', rr);

    bwd_images = am.backward(processed_data);
    img_diff = bwd_images - complex_images;
    fprintf(...
        '---\n difference between reconstructed and back-projected images: %e\n', ...
        (img_diff.norm()/complex_images.norm()));

    xFy = processed_data * fwd_acqs;
    Bxy = bwd_images * complex_images;
    fprintf('---\n (x, F y) = (%e, %e)\n', real(xFy), imag(xFy));
    fprintf('= (B x, y) = (%e, %e)\n', real(Bxy), imag(Bxy));

catch err
    % display error information
    fprintf('%s\n', err.message)
    fprintf('error id is %s\n', err.identifier)
end
