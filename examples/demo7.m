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
    prep_gadgets = [{'NoiseAdjustGadget'} {'AsymmetricEchoGadget'} ...
         {'RemoveROOversamplingGadget'}];
    acq_proc = gadgetron.AcquisitionsProcessor(prep_gadgets);
    fprintf('---\n pre-processing acquisitions...\n')
    preprocessed_data = acq_proc.process(input_data);

    % perform reconstruction
    recon = gadgetron.MR_BasicGRAPPAReconstruction();
    recon.set_input(preprocessed_data);
    fprintf('---\n reconstructing...\n');
    recon.process();
    % for undersampled acquisition data GRAPPA computes G-factors
    % in addition to reconstructed images
    output = recon.get_output();

    % get real-valued reconstructed images and G-factors
    images = gadgetron.MR_extract_real_images(output.select(2));
    gfacts = gadgetron.MR_extract_real_images(output.select(2,1));

    n = images.number();
    fprintf('Enter slice number to view its data\n')
    fprintf('(a value outside the range [1 : %d] will stop this loop).\n', n)
    while (true)
        i = input('slice: ');
        if i < 1 || i > n
            break
        end
        idata = images.image_as_array(i);
        gdata = gfacts.image_as_array(i);
        idata = idata/max(max(max(idata)));
        gdata = gdata/max(max(max(gdata)));
        figure(i)
        imshow(idata(:,:,1));
        title(['image ' num2str(i)])
        figure(i + n)
        imshow(gdata(:,:,1));
        title(['G-factor ' num2str(i)])
    end
    
catch err
    % display error information
    fprintf('%s\n', err.message)
    fprintf('error id is %s\n', err.identifier)
end
