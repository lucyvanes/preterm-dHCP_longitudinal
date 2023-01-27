#!/bin/bash

cd /data/project/BIPP/dHCP_lucy/FSL_SwE/STRUCT/bothPCs_IMD_SPS

Text2Vest design.txt design.mat
Text2Vest contrasts.txt design.con
Text2Vest subjects.txt design.sub

/software/system/fsl/6.0.5.2/bin/swe -i bothPCs_IMD_SPS_subs.nii.gz -o bothPCs_IMD_SPS -d design.mat -t design.con -s design.sub -m /data/project/BIPP/dHCP_lucy/FSL_SwE/STRUCT/extended_gm_cort_sub_mask.nii.gz --wb -n 999 --modified --uncorrp --corrp -c 0.001 -T