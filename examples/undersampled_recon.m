% Upper-level demo, GRAPPA reconstruction of undersampled data.

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
    input_data = gadgetron.AcquisitionData(file);

    % pre-process acquisitions
    prep_gadgets = [{'NoiseAdjustGadget'} {'AsymmetricEchoAdjustROGadget'} ...
         {'RemoveROOversamplingGadget'}];
    preprocessed_data = input_data.process(prep_gadgets);

    % perform reconstruction
    recon = gadgetron.GenericCartesianGRAPPAReconstruction();
    % for undersampled acquisition data GRAPPA can compute G-factors
    % in addition to reconstructed images
    recon.compute_gfactors(true);
    recon.set_input(preprocessed_data);
    fprintf('---\n reconstructing...\n');
    recon.process();
    images = recon.get_output('image');
    gfacts = recon.get_output('gfactor');

    % get real-valued reconstructed images and G-factors
    images = images.real();
    gfacts = gfacts.real();

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
