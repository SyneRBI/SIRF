    data = image.density();
    data = data/max(max(max(data)));
    figure(1000000)
    imshow(data(:,:,1));

%opengl software

%     scale = 1.0/max(max(max(data)));
%     stir.show(data, scale, 10)
%     %drawnow

%         scale = 1.0/max(max(max(data)));
%         stir.show(data, scale, 10)

%    drawnow

    %image0.initialise(60, 60, 31, 4.44114, 4.44114, 3.375)    
    image0.fill(1.0)
    image1 = image0.clone();
    image = image0.get_empty_copy();

%    obj_fun.set_input_filename('Utahscat600k_ca_seg4.hs')

%    obj_fun.set_up()

%    recon.set_start_subset_num(0)

%    recon.set_start_subiteration_num(1)

    fprintf('subiteration range: %d to %d\n', start, stop)

    recon.set_subiteration_num(start)

    