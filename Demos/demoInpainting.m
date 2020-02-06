% Demo for inpainting of a color image using the Mumford-Shah model

% load image
img = double(imread('fruitsColor.jpg'))/255;

% add Gaussian noise
sigma = 0.02;
rng(123)
imgNoisy = img + sigma * randn(size(img));

% set weights and destroy image
[m,n,l] = size(img);
missingFraction = 0.6;
weights = rand(m,n) > missingFraction;
imgNoisy(~cat(3, weights, weights, weights)) = 0; 

% MumfordShah restoration
% define restoration options
s = 0.07; % jump size
alpha = 0.8; % smoothing parameter
gamma = s^2 * alpha;
opts.method = 'L2';
opts.isotropic = 2; 
opts.muSeq = @(k) k^2.01 * 1e-6; % 
%opts.muSeq = @(k) 1.5^k * 1e-6; % 
opts.verbose = true;
opts.tol = 0.01;

prox = makeProxL2w(imgNoisy, weights);
u = mumfordShah2D(gamma, alpha, prox, opts);
imshow(u)

% Show result
subplot(2,2,1)
imshow(img)
title('Original')
subplot(2,2,2)
imshow(imgNoisy)
title('Noisy image and black pixels are missing')
subplot(2,2,3)
imshow(u)
title('Mumford Shah restoration')

