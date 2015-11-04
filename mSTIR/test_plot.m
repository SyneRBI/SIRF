opengl software

load data
size(data, 1)
size(data, 2)
size(data, 3)
scale = 1.0/max(max(max(data)))
figure(1000000)
stir.show(data, scale, 1)
drawnow
