# FISH-In-Situ-Hybridization
Fluorescence In Situ Hybridization Registration

The scripts provide 2D registration between large fluorescence in-situ hybridization images.
The registration is piecewise translation/rigid/affine in chunks. After the chunks are registered
between images, they are stitched back to get a smooth representation of the image.

<p align="center">
  <img src="https://github.com/SNIR-NIMH/FISH-In-Situ-Hybridization/blob/main/imgs/original.png" height="300"/>  
 <img src="https://github.com/SNIR-NIMH/FISH-In-Situ-Hybridization/blob/main/imgs/after_reg.png" height="300"/>  
</p>
Zoomed view
<p align="center">
  <img src="https://github.com/SNIR-NIMH/FISH-In-Situ-Hybridization/blob/main/imgs/zoom1.png" height="300"/>  
 <img src="https://github.com/SNIR-NIMH/FISH-In-Situ-Hybridization/blob/main/imgs/zoom2.png" height="300"/>  
</p>

## Installation
The scripts are written in MATLAB. Source codes are provided.
Compiled versions are also provided if MATLAB license isn't
available. To run the compiled code, download MCR for R2022a (v912),
```
https://www.mathworks.com/products/compiler/matlab-runtime.html
```
install some suitable location (no administrator access needed), replace the MCRROOT variable
in the following line `MCRROOT=/usr/local/matlab-compiler/v912` withing the shellscripts
to the locally installed v912 folder.

## Usage
The registration between similar channels (e.g. DAPI or DAPI-like) with overlapping information
can be done using `FISH_registration.sh`

```./FISH_registration.sh FIXED  MOVING  OUTPUTDIR  NAME_1 VALUE_1 ...
 
 FIXED       Source/target/fixed image, or a folder containing multiple 2D tifs
 MOVING      Moving image, to be registered to the target fixed image. It can
             also be a folder containing same number of 2D tifs as FIXED, where
             each of them will be registered to the corresponding tif in FIXED
 OUTPUTDIR   Output directory where all output files will be written
 
 
 Optional name/value parameters
 CHUNKSIZE   (Optional) Number of chunks in height x width, a string separated by x.
             Example: 10x12 means 10 chunks in height and 12 chunks in width are used.
             If not mentioned, it will automatically be computed.
 DSFACTOR    (Optional) Initial downsampling factor to try the coarse
             registration. Too high downsampling factor cause unstable
             registration. Enter a range in Matlab notation in increasing
             order, such as 50:-2:10. Default is 31:-2:3.
 Overlap     Overlap factor between chunks. Default is 0.05 (5% overlap). Enter
             a number between 0 and 1.
 MaxTranslationPixels
             (Optional) Maximum translation in pixels that will be tolerated
             during the chunk-wise registration. If not mentioned, default is 
             2x overlap, i.e., 10% of the chunk size. Increase this tolerance
             if there is significant translation in the image.
 WriteOutput (Optional) Either 1 or 0. Default 1, indicating the registered
             image and transformation matrix will be written to disk as a tif
             file and a mat file, respectively. 0 indicates only their value
             will be returned (useful in interactive mode).
 InitReg     (Optional) Initial registration type to estimate the proper
             downsampling factor. Default is rigid. Options are rigid/translation.
             This is the registration done on the whole image.             
 FinalReg    (Optional) Registration type to use to register chunks. Default is
             similarity (i.e. affine). Options are similarity (affine), rigid, 
             or translation.             
 Modality    (Optional) Either monomodal or multimodal. Default monomodal.
 NumCPU      (Optional) When both FIXED and MOVING are folders, they can be
             processed in parallel using NumCPU processes. Default 12.
 Pyramidlevel  (Optional) Number of pyramid levels to use for registraton.
             Default is 2. Use higher (4) if the images are large. Using small
             (1) may incur in bad registration.           
 Mask        (Optional) A mask image to introduce better registration  when
             there are a lot of background noise
 Quantile    (Optional) Default 0.99. An upper quantile to clamp the outlier
             intensities. It is useful if there are a lot of outliers that can
             cause registration problems.
```

Example command:
```
./FISH_registration.sh myfixedimg_round1_dapi.tif mymovingimg_round2_dapi.tif /home/user/outputfolder/   \
 chunksize 20x40 dsfactor 10:-1:2 finalreg translation initreg translation Modality monomodal overlap 0.2 quantile 0.95
```
A transformation matrix in MATLAB .mat format will be created.
Once a good registration between rounds are obtained, the other channels can be transformed using
`FISH_transformation.sh`,
```
./FISH_transformation.sh FIXED  MOVING  TFORM_MAT   OUTPUT
 
 FIXED      Fixed image used for registration. It is also called "reference"
            image. The output image will have same dimension as this one. It is
            the first argument of the FIST_registration code.
 MOVING     Moving image (input) to be transformed to the fixed image
 TFORM      Transformation matrix obtained from the FISH_registration code
 OUTPUT     Output (.tif) image, transformed version of moving image
```

## Notes:
1. "Initreg" is done on the whole image to estimate the initial linear transformation.
   If the images are very big, more often that not translation or rigid works better than
   similarity (i.e., affine)
2. If there is too much deformation, then choose higher "Overlap" and "Chunksize".
3. For similar looking images (e.g. DAPI), monomodal using phase correlation (`imregcorr`) 
   as a metric works fine. If the images are not very similar, multimodal registration
   (`imregtform`) using mutual information works better


<!-- LICENSE -->
## License

Distributed under the MIT License. See `LICENSE.txt` for more information.

<!-- CONTACT -->
## Contact
Snehashis Roy - email@snehashis.roy@nih.gov
