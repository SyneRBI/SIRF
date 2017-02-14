function K = add_noise(K, sd)
% ADD_NOISE Emulates add_noise function in ISMRMRD C++ test.
%
% K = add_noise(K, sd)
% K is complex k-space
% sd is Std Dev of noise with zero mean and Normal distribution in real and
% imaginary channels.
%
% David Atkinson  D.Atkinson@ucl.ac.uk
%


nreal = random('Normal',0, sd, size(K)) ;
nimag = random('Normal',0, sd, size(K)) ;

K = K + complex(nreal, nimag) ;
