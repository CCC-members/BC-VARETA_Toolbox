#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jul 17 22:20:23 2023

@author: dpaz
"""

import nibabel  # neuroimaging data formats
import os  # operating system interfaces
from scipy.io import savemat
from scipy.sparse.linalg import LinearOperator, eigsh # scipy sparse linear algebra libraries
from sksparse.cholmod import cholesky # scikit-sparse cholesky decomposition
from nilearn.plotting import plot_surf_stat_map # neuroimaging visualization
from matplotlib import pyplot # matplotlib visualization
from lapy import TriaMesh, Solver # lapy like surface triangulation, and finite elmnt method for Laplace-Beltrami

class SurfaceSolver:
    def __init__(
        self, 
        input_dir: str, 
        subj: str, 
        space: str, 
        hemi: str,
        geom: str,
        field: str,
        output_dir: str,
        nefunc: int,
        sigma: float,
    ):
        # Initialize geometrical atributes
        self.geometry = type('geometry', (), {})()
        ## metadata
        self.geometry.metadata = type('metadata', (), {})()
        self.geometry.metadata.subj = subj
        self.geometry.metadata.hemi = hemi
        self.geometry.metadata.geom = geom
        self.geometry.metadata.path = input_dir + '/' + subj + '/' + space + '/' + subj + '.' + hemi + '.' + geom
        ## data
        self.geometry.data = type('data', (), {})()
        self.geometry.data.img = nibabel.load(self.geometry.metadata.path)
        # Initialize field atributes
        self.field = type('field', (), {})()
        ## metadata
        self.field.metadata = type('metadata', (), {})()
        self.field.metadata.subj = subj
        self.field.metadata.hemi = hemi
        self.field.metadata.geom = geom
        self.field.metadata.field = field
        self.field.metadata.path = input_dir + '/' + subj + '/' + space + '/' + subj + '.' + hemi + '.' + field
        ## data
        self.field.data = type('data', (), {})()
        self.field.data.img = nibabel.load(self.field.metadata.path)
        # Initialize ciftimorph atributes
        self.ciftimorph = type('ciftimorph', (), {})()
        ## metadata
        self.ciftimorph.metadata = type('metadata', (), {})()
        self.ciftimorph.metadata.subj = subj
        self.ciftimorph.metadata.hemi = hemi
        self.ciftimorph.metadata.geom = geom
        self.ciftimorph.metadata.field = field
        self.ciftimorph.metadata.path = output_dir + '/' + subj + '/' + space
        ## data
        self.ciftimorph.data = type('data', (), {})()
        self.ciftimorph.data.nefunc = nefunc
        self.ciftimorph.data.sigma = sigma
        self.ciftimorph.data.vert = self.geometry.data.img.agg_data('pointset')
        self.ciftimorph.data.face = self.geometry.data.img.agg_data('triangle')
        self.ciftimorph.data.tria = TriaMesh(
            self.ciftimorph.data.vert, 
            self.ciftimorph.data.face
            )

    def LaBel(self):      
        ## lapy computation via finite element method (in weak integral form) of the Laplace-Beltrami operator on surface 
        ## returns stifness and mass matrix for eigenvalue problem Stiff*efunc = eval*Mass*efunc
        self.ciftimorph.data.FEM = Solver(
            geometry=self.ciftimorph.data.tria, 
            lump=False, 
            use_cholmod=True
            )
        self.ciftimorph.data.Stiff =  self.ciftimorph.data.FEM.stiffness
        self.ciftimorph.data.Mass =  self.ciftimorph.data.FEM.mass
        print('...Computed stiffness and mass matrices Stiff and Mass via lapy finite element method')     
        ## eigenfunctions for Laplace-Beltrami operator via the shift invert method with shift operator Stiff - sigma*Mass
        ## the new eigenfunction problem turns inv(Stiff - sigma*Mass)*Mass*efunc = eval*efunc, the eigenvalues then turn 
        ## eigenvalues = 1/(eigenvalues - sigma)
        self.ciftimorph.data.chol = cholesky(
            self.ciftimorph.data.Stiff - self.ciftimorph.data.sigma * self.ciftimorph.data.Mass
            )
        self.ciftimorph.data.opinv = LinearOperator(
            matvec=self.ciftimorph.data.chol, 
            shape=self.ciftimorph.data.Stiff.shape, 
            dtype=self.ciftimorph.data.Stiff.dtype
            )
        self.ciftimorph.data.eval, self.ciftimorph.data.efunc = eigsh(
            self.ciftimorph.data.Stiff, 
            self.ciftimorph.data.nefunc, 
            self.ciftimorph.data.Mass, 
            sigma=self.ciftimorph.data.sigma, 
            OPinv=self.ciftimorph.data.opinv, 
            which='LM'
            )
        print('...Computed eigenfunctions and eigenvalues via generalized eigendecomposition for Stiff*efunc = eval*Mass*efunc')
        return self

    def Write(self):
        # path to write files and figures
        if not os.path.exists(self.ciftimorph.metadata.path):
            os.makedirs(self.ciftimorph.metadata.path)
            print('...Created directory for output surface data: ', self.ciftimorph.metadata.path)
        ## write eigenfunctions in gifti
        efunc = nibabel.gifti.GiftiDataArray(
            self.ciftimorph.data.efunc
            )
        efunc.intent = nibabel.nifti1.intent_codes['NIFTI_INTENT_TTEST']
        efunc.datatype = nibabel.nifti1.data_type_codes['NIFTI_TYPE_FLOAT32']
        img = self.geometry.data.img
        img.add_gifti_data_array(efunc)
        nibabel.save(img,
            self.ciftimorph.metadata.path 
            + '/' + self.ciftimorph.metadata.subj 
            + '.' + 'Eigenfunctions' 
            + '.' + self.ciftimorph.metadata.hemi + '.' 
            + self.ciftimorph.metadata.geom)
        del efunc
        del img
        ## write eigenfunctions in mat        
        savemat(
            self.ciftimorph.metadata.path 
            + '/' + self.ciftimorph.metadata.subj 
            + '.' + 'Eigenfunctions' 
            + '.' + self.ciftimorph.metadata.hemi + '.' 
            + self.ciftimorph.metadata.geom + '.mat', 
            {'nefunc': self.ciftimorph.data.nefunc,
             'sigma': self.ciftimorph.data.sigma,
             'vert': self.ciftimorph.data.vert,
             'face': self.ciftimorph.data.face,
             'Stiff': self.ciftimorph.data.Stiff,
             'Mass': self.ciftimorph.data.Mass,
             'efunc': self.ciftimorph.data.efunc})
        # visualization of the Stiffness and Mass matrix
        fig1, axes1 = pyplot.subplots(1, 2, figsize=(12, 6))
        fig1.suptitle('Stiffness and Mass matrix so that Stiff*efunc = eval*Mass*efunc')
        axes1[0].spy(self.ciftimorph.data.Stiff)
        axes1[0].set_title('Stiff')
        axes1[1].spy(self.ciftimorph.data.Mass)
        axes1[1].set_title('Mass')
        ## save figure
        pyplot.savefig(
            self.ciftimorph.metadata.path 
            + '/' + self.ciftimorph.metadata.subj 
            + '.' + 'Stiff&Mass' 
            + '.' + self.ciftimorph.metadata.hemi + '.' 
            + self.ciftimorph.metadata.geom + '.png'
            )
        pyplot.close(fig1)
        # visualization of the eigenfunctions
        if self.ciftimorph.metadata.hemi == 'L':
            hemi = 'left'
        elif self.ciftimorph.metadata.hemi == 'R':
            hemi = 'right'
        fig2 = pyplot.figure(figsize=(24, 12))
        fig2.suptitle('Eigenfunctions of the Laplace Beltrami operator')
        ## iterate for eigenfunctions
        start = 1
        end = 10
        for e in range(start, end + 1):
            ax = pyplot.subplot(2, 5, e, projection='3d')
            title = 'efunc' + str(e) + '(' + str(self.ciftimorph.data.eval[e]) + ')'
            plot_surf_stat_map(
                self.geometry.metadata.path, 
                self.ciftimorph.data.efunc[:,e], 
                hemi=hemi, 
                title=title, 
                colorbar=False,
                threshold=0., 
                bg_map=self.field.metadata.path, 
                bg_on_data=True, 
                engine='matplotlib',
                axes=ax
                )
        ## save figure
        pyplot.savefig(
            self.ciftimorph.metadata.path 
            + '/' + self.ciftimorph.metadata.subj 
            + '.' + 'Eigenfunctions' 
            + '.' + self.ciftimorph.metadata.hemi + '.' 
            + self.ciftimorph.metadata.geom + '.png'
            )
        pyplot.close(fig2)