#!/bin/bash
exe_name=$0
exe_name=`readlink -f ${exe_name}`
exe_dir=`dirname ${exe_name}`

if [ $# -lt "3" ]; then
  echo "
  
 ./FISH_registration.sh FIXED  MOVING  OUTPUTDIR  NAME_1 VALUE_1 ...
 
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
             
  "
  exit 1
fi

export MCR_CACHE_ROOT=/tmp/mcr_${USER}_${RANDOM}
mkdir -p ${MCR_CACHE_ROOT}

MCRROOT=/usr/local/matlab-compiler/v912
LD_LIBRARY_PATH=.:${MCRROOT}/runtime/glnxa64 ;
LD_LIBRARY_PATH=${LD_LIBRARY_PATH}:${MCRROOT}/bin/glnxa64 ;
LD_LIBRARY_PATH=${LD_LIBRARY_PATH}:${MCRROOT}/sys/os/glnxa64;
LD_LIBRARY_PATH=${LD_LIBRARY_PATH}:${MCRROOT}/sys/opengl/lib/glnxa64;
export LD_LIBRARY_PATH;
export MCR_INHIBIT_CTF_LOCK=1


args=
while [ $# -gt 0 ]; do
  token=$1
  args="${args} ${token}" 
  shift
done
echo ${exe_dir}/FISH_registration $args
${exe_dir}/FISH_registration $args
rm -rf ${MCR_CACHE_ROOT}
