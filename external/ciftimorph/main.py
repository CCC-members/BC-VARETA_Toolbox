#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Jul 15 18:15:23 2023

@author: dpaz
"""

from ciftimorph.surface_solver import SurfaceSolver

# Create an instance of the SurfaceSolver class
obj = SurfaceSolver(
    input_dir='/mnt/Cloud/OneDrive/Data/CHBM/Preprocessing/Ciftify',
    subj='sub-CBM00005',
    space='MNINonLinear',
    hemi ='L',
    geom='pial.164k_fs_LR.surf.gii',
    field='atlasroi.164k_fs_LR.shape.gii',
    output_dir='/mnt/Cloud/OneDrive/Data/CHBM/Preprocessing/CiftiMorph',
    sigma = -0.01,
    nefunc = 1000,
    )
obj = obj.LaBel()
obj = obj.Write()


# from ciftimorph.surface_solver import SurfaceSolver

# # Create an instance of the SurfaceSolver class
# obj = SurfaceSolver(
#     input_dir='/mnt/Cloud/OneDrive/Data/CHBM/Preprocessing/Ciftify',
#     subj='sub-CBM00005',
#     space='MNINonLinear',
#     hemi ='L',
#     geom='pial.164k_fs_LR.surf.gii',
#     field='atlasroi.164k_fs_LR.shape.gii',
#     output_dir='/mnt/Cloud/OneDrive/Data/CHBM/Preprocessing/CiftiMorph',
#     sigma = -0.01,
#     nefunc = 1000,
#     )
# obj = obj.LaBel()
# obj = obj.Write()