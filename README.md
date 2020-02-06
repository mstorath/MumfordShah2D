# Edge preserving image restoration based on the Mumford-Shah model
This repository contains algorithms for edge preserving image restoration based on the Mumford-Shah model.
This is an implementation of the methods described in the paper

K. Hohm, M. Storath, A. Weinmann,
An algorithmic framework for Mumford–Shah regularization of inverse problems in imaging,
Inverse Problems 31 (11), 115011

## Applications examples
### Edge preserving smoothing of vector-valued images

   - Supports smoothing of vector-valued images (e.g. multispectral images, feature images)
   - Linear complexity in number of color channels
   - No discretization of color space required

![Denoising](/Docs/img_salt_pepper.png)

### Regularization for deconvolution

![Deconvolution](/Docs/img_deconv.png)

### Inpainting

![Deconvolution](/Docs/img_inpainting.png)

## Quickstart:
   - Run the script "setPath.m", it should set all necessary paths
   - For best performance, increase Java heap space in the Matlab preferences (MATLAB - General - Java heap memory)
   - Run a demo from the Demos folder


