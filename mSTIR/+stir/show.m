function show(imageData, scale, z)
    a = scale*imageData(:,:,z);
    imshow(a)
end
