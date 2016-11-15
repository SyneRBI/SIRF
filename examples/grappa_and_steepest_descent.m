% GRAPPA reconstruction with the steepest descent step
% to illustrate the use of Acquisition Model projections.

if ~libisloaded('mutilities')
    fprintf('loading mutilities library...\n')
    [notfound, warnings] = loadlibrary('mutilities');
end
if ~libisloaded('mgadgetron')
    fprintf('loading mgadgetron library...\n')
    [notfound, warnings] = loadlibrary('mgadgetron');
end

try

    % define raw data source
    file = input('raw data file: ', 's');
    input_data = gadgetron.MR_Acquisitions(file);

    % pre-process acquisitions
    prep_gadgets = [{'NoiseAdjustGadget'} {'AsymmetricEchoAdjustROGadget'} ...
         {'RemoveROOversamplingGadget'}];
    preprocessed_data = input_data.process(prep_gadgets);
    pp_norm = preprocessed_data.norm();

    % perform reconstruction
    recon = gadgetron.MR_BasicGRAPPAReconstruction();
    recon.set_input(preprocessed_data);
    fprintf('---\n reconstructing...\n');
    recon.process();
    output = recon.get_output();
    complex_images = output.select(2);

    % compute coil sensitivity maps
    csms = gadgetron.MR_CoilSensitivityMaps();
    fprintf('---\n sorting acquisitions...\n')
    preprocessed_data.sort()
    fprintf('---\n calculating sensitivity maps...\n')
    csms.calculate(preprocessed_data)

    % create acquisition model based on the acquisition parameters
    % stored in preprocessed_data and image parameters stored in complex_images
    am = gadgetron.MR_AcquisitionModel(preprocessed_data, complex_images);
    am.set_coil_sensitivity_maps(csms)

    % use the acquisition model (forward projection) to simulate acquisitions
    fwd_data = am.forward(complex_images);
    fwd_norm = fwd_data.norm();
    % compute the difference between real and simulated acquisitions
    diff = fwd_data - preprocessed_data * (fwd_norm/pp_norm);
    rr = diff.norm()/fwd_norm;
    fprintf('---\n reconstruction residual norm (rel): %e\n', rr)

    % try to improve the reconstruction by the steepest descent step
    g = am.backward(diff);
    w = am.forward(g);
    alpha = (g*g)/(w*w);
    r_complex_imgs = complex_images - g*alpha;

    % get real-valued reconstructed and refined images
    images = complex_images.real();
    r_imgs = r_complex_imgs.real();

    % plot images
    n = images.number();
    fprintf('Enter slice number to view it\n')
    fprintf('(a value outside the range [1 : %d] will stop this loop)\n', n)
    while (true)
        i = input('slice: ');
        if i < 1 || i > n
            break
        end
        data = images.image_as_array(i);
        rdata = r_imgs.image_as_array(i);
        data = data/max(max(max(data)));
        rdata = rdata/max(max(max(rdata)));
        figure(i)
        imshow(data(:,:,1));
        title('GRAPPA image')
        figure(i + n)
        imshow(rdata(:,:,1));
        title('refined image')
    end
    
catch err
    % display error information
    fprintf('%s\n', err.message)
    fprintf('error id is %s\n', err.identifier)
end
