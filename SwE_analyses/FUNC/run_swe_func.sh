#!/bin/bash

cd /data/project/BIPP/dHCP_lucy/FSL_SwE/FUNC/bothPCs_IMD_SPS

Text2Vest design.txt design.mat
Text2Vest contrasts.txt design.con
Text2Vest subjects.txt design.sub

/software/system/fsl/fsl-6.0.5/bin/swe -i bothPCs_IMD_SPS_subs.nii.gz -o bothPCs_IMD_SPS -d design.mat -t design.con -s design.sub -m /data/project/BIPP/dHCP_lucy/FSL_SwE/FUNC/week40_gm_mask_1_5.nii.gz --wb -n 999 --modified --uncorrp --corrp -c 0.001 -T