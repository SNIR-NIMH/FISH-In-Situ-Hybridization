#!/bin/bash
exe_name=$0
exe_name=`readlink -f ${exe_name}`
exe_dir=`dirname ${exe_name}`

if [ $# -lt "4" ]; then
  echo "
  
 ./FISH_transformation.sh FIXED  MOVING  TFORM_MAT   OUTPUT
 
 FIXED      Fixed image used for registration. It is also called \"reference\"
            image. The output image will have same dimension as this one. It is
            the first argument of the FIST_registration code.
 MOVING     Moving image (input) to be transformed to the fixed image
 TFORM      Transformation matrix obtained from the FISH_registration code
 OUTPUT     Output (.tif) image, transformed version of moving image

  "
  exit 1
fi


MCRROOT=/usr/local/matlab-compiler/v912
LD_LIBRARY_PATH=.:${MCRROOT}/runtime/glnxa64 ;
LD_LIBRARY_PATH=${LD_LIBRARY_PATH}:${MCRROOT}/bin/glnxa64 ;
LD_LIBRARY_PATH=${LD_LIBRARY_PATH}:${MCRROOT}/sys/os/glnxa64;
LD_LIBRARY_PATH=${LD_LIBRARY_PATH}:${MCRROOT}/sys/opengl/lib/glnxa64;
export LD_LIBRARY_PATH;

export MCR_INHIBIT_CTF_LOCK=1
export MCR_CACHE_ROOT=/tmp/mcr_${USER}_${RANDOM}
mkdir -p ${MCR_CACHE_ROOT}

args=
while [ $# -gt 0 ]; do
  token=$1
  args="${args} ${token}" 
  shift
done
echo ${exe_dir}/FISH_transformation $args
${exe_dir}/FISH_transformation $args
rm -rf ${MCR_CACHE_ROOT}
