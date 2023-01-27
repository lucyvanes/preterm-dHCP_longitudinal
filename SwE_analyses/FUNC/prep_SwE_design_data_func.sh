#!/bin/bash

cd /data/project/BIPP/dHCP_lucy/FSL_SwE/FUNC

dos2unix <scans_bothPCs.txt |
while IFS= read s 
do
 echo ${s}
 
cp /data/project/dHCP/personal/Lucy/DC_weighted_Z_Lucy/${s} bothPCs_IMD_SPS/
 
done

cd bothPCs_IMD_SPS
echo "merging to 4D"
fslmerge -t bothPCs_IMD_SPS_subs.nii.gz sub*.nii

echo "removing individual scans"
rm sub*.nii