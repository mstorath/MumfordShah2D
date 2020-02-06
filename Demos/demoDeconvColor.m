% Demo for deconvolution of a color image using the Mumford-Shah model

% prepare convolved image
rng(123)
scale = 1;
img = imresize(double(imread('churchColor.jpg'))/255, scale);
[m,n, ~] = size(img);
stdDevKernel = scale*4;
K = fspecial('Gaussian', [m,n], stdDevKernel);
fftK = fft2(fftshift(K));
A = convop(cat(3,fftK,fftK, fftK));
fBlurry = A * img;
sigma = 0.025;
noise = sigma * randn(size(fBlurry)) * max(fBlurry(:));
fNoisy = fBlurry + noise;
figure
imshow(fNoisy)


% compute MumfordShah solution
clear opts;
opts.verbose = 1;
opts.method = 'L2';
opts.isotropic = 2;
opts.groundTruth = img;
opts.muSeq = @(k) (1.05)^k * 1e-6;
%opts.muSeq = @(k) 2^k * 1e-6; % this progession is faster but gives lower quality

% set parameters
j = 0.1;
alpha = 0.15;
gamma = alpha*j^2;

proxHandle = makeProxL2Linop( fNoisy, A);
u = mumfordShah2D(gamma, alpha, proxHandle, opts);

% show solution
subplot(1,2,1)
imshow(fNoisy)
title('Blurred and noisy image')
subplot(1,2,2)
imshow(u)
title('Mumford-Shah result')
