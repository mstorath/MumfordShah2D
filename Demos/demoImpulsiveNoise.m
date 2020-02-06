% Demo for restortation of a color image from impulsive noise 
% using the Mumford-Shah model

img = double(((imread('manWithHat.jpg'))))/255;

% add Salt and Pepper noise
rng(123);
sigma = 0.5;
imgNoisy = imnoise(img, 'salt & pepper', sigma);

% MumShah restoration
clear opts;
s = 0.3 ;
alpha = 40  ;
gamma = s^2 * alpha;
opts.verbose = 1;
opts.method = 'L2';
opts.muSeq = @(k) 2^k * 1e-6;
opts.verbose = true;
prox = makeProxL0w(imgNoisy, ones(size(imgNoisy)));
u = mumfordShah2D(gamma, alpha, prox, opts);

% show solution
figure
subplot(1,2,1)
imshow(imgNoisy)
title('Noisy image')
subplot(1,2,2)
imshow(u)
title('Mumford-Shah result')

