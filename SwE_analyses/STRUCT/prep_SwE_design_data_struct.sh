#!/bin/bash

cd /data/project/BIPP/dHCP_lucy/FSL_SwE/STRUCT

dos2unix <scans_bothPCs.txt |
while IFS= read s 
do
 echo ${s}
 
cp /data/project/BIPP/dHCP_lucy/longitudinal_preterm_jacobians/${s} bothPCs_IMD_SPS/


done

cd bothPCs_IMD_SPS
echo "merging to 4D"
fslmerge -t bothPCs_IMD_SPS_subs.nii.gz sub*.nii

echo "removing individual scans"
rm sub*.nii

# echo "downsampling merged file to 1.5mm"
# flirt -in bothPCs_IMD_SPS_rad5excl_subs.nii.gz -ref bothPCs_IMD_SPS_rad5excl_subs.nii.gz -out bothPCs_IMD_SPS_rad5excl_subs_15.nii.gz -interp nearestneighbour -applyisoxfm 1.5