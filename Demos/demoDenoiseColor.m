% Demo for denoising of a color image using the Mumford-Shah model

% prepare image
rng(123)
scale = 1;
f = imresize(double(imread('churchColor.jpg'))/255, scale);
sigma = 0.1;
noise = sigma * randn(size(f)) * max(f(:));
fNoisy = f + noise;
figure
imshow(fNoisy)


% compute MumfordShah solution
clear opts;
opts.verbose = 1;
opts.method = 'L2';
opts.isotropic = 2;
opts.groundTruth = f;
opts.muSeq = @(k) 2^k * 1e-6;

% model parameters
j = 0.1;
alpha = 10; % smoothing parameter
gamma = alpha*j^2; % jump parameter

proxHandle = makeProxL2w( fNoisy, ones(size(fNoisy)) );
u = mumfordShah2D(gamma, alpha, proxHandle, opts);

% show solution
subplot(1,2,1)
imshow(fNoisy)
title('Noisy image')
subplot(1,2,2)
imshow(u)
title('Mumford-Shah result')
