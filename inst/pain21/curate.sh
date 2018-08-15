#!/bin/bash

# unzip nidm files
for i in $(ls *zip); do unzip -o $i -d $(echo $i | sed "s+.zip++g"); done

# scale SPM contrast maps
for i in $(ls pain_*.nidm/Contrast.nii.gz); do
  echo $i
  out=$(dirname $i | sed "s+.nidm++g");
  out=pain21/contrast_$out.nii.gz
  fslmaths $i -nan -mul 100 $out;
done
# scale and square SPM sd maps
for i in $(ls pain_*.nidm/ContrastStandardError.nii.gz); do
  echo $i
  out=$(dirname $i | sed "s+.nidm++g");
  out=pain21/var_$out.nii.gz
  fslmaths $i -nan -mul 100 -sqr $out;
done

# copy fsl contrast maps
for i in $(ls pain_*.nidm/Contrast_T001.nii.gz); do
  echo $i
  out=$(dirname $i | sed "s+.nidm++g");
  out=pain21/contrast_$out.nii.gz
  cp $i $out;
done

# create mask
# square fsl sd maps
for i in $(ls pain_*.nidm/ContrastStandardError_T001.nii.gz); do
  echo $i
  out=$(dirname $i | sed "s+.nidm++g");
  out=pain21/var_$out.nii.gz
  fslmaths $i -sqr $out;
done

# create mask with over 10 observations
cp pain_01.nidm/Mask.nii.gz pain21/mask.nii.gz
for i in $(ls pain_*/Mask.nii.gz); do fslmaths pain21/mask.nii.gz -add $i pain21/mask.nii.gz; done
fslmaths pain21/mask.nii.gz -thr 10.5 -bin pain21/mask.nii.gz

for i in $(ls pain_*/DesignMatrix.csv); do
  out=$(dirname $i | sed "s+.nidm++g");
  nlines=$(cat $i | wc -l)
  echo $out,$nlines
done > pain21/samplesize.csv
