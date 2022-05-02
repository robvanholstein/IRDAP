'''
This file contains all ADI+PCA auxiliary functions used by IRDAP.
File written by Julien Milli (julien.milli@univ-grenoble-alpes.fr).

IRDAP is a Python package to accurately reduce SPHERE-IRDIS polarimetric data.
Copyright (C) 2019 R.G. van Holstein

Full documentation: https://irdap.readthedocs.io
Feedback, questions, comments: rob.vanholstein@eso.org

When publishing data reduced with IRDAP, please cite van Holstein et al.
(2020): https://ui.adsabs.harvard.edu/abs/2020A%26A...633A..64V/abstract.
For data in pupil-tracking mode please additionally cite van Holstein et al.
(2017): https://ui.adsabs.harvard.edu/abs/2017SPIE10400E..15V.

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <https://www.gnu.org/licenses/>.
'''

import numpy as np

class pca_adi(object):
    """
    Wrapper for the pca class that can handle cubes of images.
    It flattens the images into 1D vectors referred to attributes. The pixels of
    the images are referred to as objects. We therefore end up with a 2D
    maxtrix of Nobj objects (rows) by Katt attributes (columns) that is the starting
    point for the PCA, implemented in another class called pca.
    """

    def __init__(self,datacube,method='cor',verbose=True,radii=None):
        """
        Constructor of the pca_imagecube class.
        Input:
            - datacube: a 3d numpy array
            - method: 'cov' for covariance (default option), 'cor' for correlation
                or 'ssq' for sum of squares
            - verbose: True or False if you want some information printed on the terminal
            - radii: an array containing the radii in pixels of the annuli in
                which the PCA must be calculated. For instance: radii=[10,100,200] means
                the PCA will be computed in 2 annuli defined by 10px-100px and 100px-200px.
                By default, assumes te whole image is used.
        """
        if datacube.ndim != 3:
            raise IndexError('The input datacube must be a 3D  numpy array !')
        self.nframes,self.ny,self.nx = datacube.shape
        if radii is None:
            radii = [0,int(np.round(np.sqrt((self.ny//2)**2+(self.nx//2)**2)))]
        self.method=method
        #For IRDAP, the center of a 1024x1024 frame is in between 4 pixels at (511.5,511.5)
        x_array = np.arange(self.nx)-(self.nx/2.-0.5)
        y_array = np.arange(self.ny)-(self.ny/2.-0.5)
        xx_array,yy_array=np.meshgrid(x_array,y_array)
        distarr = np.abs(xx_array+1j*yy_array)
        self.region_map = np.zeros((self.ny,self.nx),dtype=int)
        self.nb_annuli = len(radii)-1
        self.Nobj_array = np.ndarray(self.nb_annuli)
        self.pca_array = []
        self.x_indices_array = []
        self.y_indices_array = []
        if verbose:
            print('There are {0:d} frames and {1:d} regions.'.format(self.nframes,\
                  self.nb_annuli))
        for i in range(self.nb_annuli):
            y_indices,x_indices = np.where(np.logical_and(distarr>=radii[i],\
                                                          distarr<radii[i+1]))
            self.y_indices_array.append(y_indices)
            self.x_indices_array.append(x_indices)
            self.Nobj_array[i] = len(y_indices)
            self.region_map[y_indices,x_indices] = i+1
            data = datacube[:,y_indices,x_indices].T # Transpose is used to get a shape (Nobj x Katt) where Katt is the number of frames of the datacube
            self.pca_array.append(pca(data,method=method,verbose=verbose))
            if verbose:
                self.pca_array[i].print_explained_inertia(modes=5)

    def compute_residuals(self,truncation=None):
        """
        Reconstructs the datacube, by applying a (truncated) pca and subtract
        the reconstructed data cube from the data themselves to obtain the
        residuals.
        Input:
            - truncation: integer that should be smaller than the number of frames
                            to perform the truncation of the data. If none, use
                            all the frames.
        Output:
            - residuals_datacube:
        """
        residuals_datacube = np.zeros((self.nframes,self.ny,self.nx))*np.nan
        for i in range(self.nb_annuli):
            residuals_data = self.pca_array[i].project_and_subtract(truncation=truncation)
            residuals_datacube[:,self.y_indices_array[i],self.x_indices_array[i]] = \
                residuals_data.T
        return residuals_datacube

class pca(object):
    """
    Object pca.
    The attributes of the object are:
        - matrix: the data matrix (after centring and normalisation if desired)
                Its shape is (Nobj,Katt) i.e. Nobj objects (lines) and Katt attributes
                (columns)
        - Nobj: the number of objects of the PCA
        - Katt: the number of attributes of the PCA
        - method: the type of PCA to perform: 'cov' (covariance, by default),
                'cor' for correlation or 'ssq' (sum of squares)
        - inertia_matrix: the inertia matrix (e.g. covariance matrix, correlation
                        matrix or sum of squares matrix, depending of the PCA
                        method chosen). It is a square matrix of shape is (Katt,Katt)
        - total_inertia: the trace of the inertia matrix (scalar value)
        - eigenval: eigenvalues of the inertia matrix, sorted by descending order.
                    This is a vector of size Katt.
        - eigenvect: eigenvectors of the inertia matrix. The shape of this matrix
                    is (Katt,Katt)
        - pc: the principal components. Its shape is (Nobj,Katt)
    """

    def __init__(self,data,method='cov',verbose=True):
        """
        Constructor of the pca class. It checks the format of the input data matrix
        and builds the covariance/correlation/sum of squares matrix.
        Then it computes the principal components.
        Input:
            - data: a 2d numpy array of shape (Nobj,Katt) i.e with Katt columns and
                    Nobj rows. Katt is the number of attributes and Nobj is the number
                    of objects.
            - method: 'cov' (default), 'cor' or 'ssq'
            - verbose: True or False if you want some information printer on the screen
        """
        if data.ndim != 2:
            raise IndexError('The PCA input matrix is not 2D !')
        if not np.all(np.isfinite(data)):
            raise TypeError('The PCA input matrix must not contain NaN values')
        self.Nobj,self.Katt = data.shape
        self.method=method
        if self.method == 'cov':
            self.mean_att = np.nanmean(data,axis=0) #mean of each attribute, size Katt
            self.matrix = data-self.mean_att
        elif self.method == 'cor':
            self.mean_att = np.nanmean(data,axis=0) #mean of each attribute, size Katt
            self.std_att = np.nanstd(data,axis=0) #std of each attribute, size Katt
            self.matrix = (data-self.mean_att)/self.std_att
        elif self.method == 'ssq':
            self.matrix = data
        else:
            raise TypeError('method not understood. Must be cov (default),\
                             cor or ssq. Got {0:s}'.format(self.method))
        if verbose:
            print('Method chosen: {0:s}'.format(self.method))
            print('There are {0:d} objects (rows) and {1:d} attributes (columns)'.format(self.Nobj,self.Katt))
            print('Computing the {0:d}x{0:d} inertia matrix...'.format(self.Katt))
        self.inertia_matrix = np.dot(self.matrix.T,self.matrix)/self.Nobj
        if verbose:
            print('Done')
        self.total_inertia = np.trace(self.inertia_matrix)
        eigenval, eigenvect = np.linalg.eigh(self.inertia_matrix)    # eigenvalues and eigenvectors
        self.eigenval = eigenval[::-1] # we re-order the eigenvalues by decreasing order
        self.eigenvect = eigenvect[:,::-1] # we re-order the eigenvectors by decreasing order
        if verbose:
            print('Computing the principal components')
        self.pc = self._compute_principal_components()
        if verbose:
            print('Done')

    def get_eigenvect(self,truncation=None):
        """
        Returns the eigenvectors, truncated if desired
        Input:
            - trunaction: None (by default) to return all eigenvectors (matrix of
            shape (Katt,Katt)) or integer smaller or equal that Katt to return a
            truncated matrix of shape (trnucation,Katt).
        """
        if truncation is None:
            return self.eigenvect
        else:
            if truncation<=self.Katt:
                return self.eigenvect[:,0:truncation]
            else:
                print("Can't truncate by more than {0:d}".format(self.Katt))
                return

    def _compute_principal_components(self):
        """
        Compute and return the principal components. The principal components
        is a matrix of shape (Nobj,Katt)
        """
        return np.dot(self.matrix,self.eigenvect)

    def get_principal_components(self,truncation=None):
        """
        Returns the principal components, optionnally after truncating
        a certain number of modes. The principal components is a matrix of shape
        (Nobj,Katt) if truncation is None, or of shape (Nobj,truncation) if truncation
        is an integer.
        Input:
            - truncation: None (by default) to use all vectors or an integer smaller
                        than Katt to truncate the number of modes.
        """
        if truncation is None:
            return self.pc
        else:
            if truncation>self.Katt:
                print("Can't truncate by more than {0:d}".format(self.Katt))
                return
            else:
                return self.pc[:,0:truncation]

    def project_matrix(self,truncation=None):
        """
        Projects a matrix or the data themselves on the (truncated) eigenvectors.
        Input:
            - truncation: None (by default) to use all vectors or an integer smaller
                        than Katt to truncate the number of modes.
        Output:
            - the projected matrix of shape (Ndata,Katt)
        """
        matrix_projected  = np.dot(\
                    self.get_principal_components(truncation=truncation),\
                    self.get_eigenvect(truncation=truncation).T)
        return matrix_projected

    def project_data(self,truncation=None):
        """
        Function that calls the project_matrix function but preliminarily subtracts
        the mean or divide by the standard deviation to normalize the data.
        Input:
        - truncation: None (by default) to use all vectors or an integer smaller
                        than Katt to truncate the number of modes.
        Output:
            - the projected data of shape (Ndata,Katt)
        """
        reconstructed_matrix = self.project_matrix(truncation=truncation)
        if self.method == 'cov':
            return reconstructed_matrix+self.mean_att
        elif self.method == 'cor':
            return reconstructed_matrix*self.std_att+self.mean_att
        elif self.method == 'ssq':
            return reconstructed_matrix

    def project_and_subtract(self,truncation=None):
        """
        Function that calls the project_data function to get the projected
        data points and then subtract this result from the data to obtain the
        residuals.
        Input:
            - truncation: None (by default) to use all vectors or an integer smaller
                        than Katt to truncate the number of modes.
        Output:
            - the residuals of shape (Ndata,Katt)
        """
        if self.method == 'cor':
            return (self.matrix - self.project_matrix(truncation=truncation))*self.std_att # the mean disappears in the difference
        elif self.method == 'cov' or self.method == 'ssq':
            return self.matrix - self.project_matrix(truncation=truncation)

    def print_explained_inertia(self,modes=10):
        """
        Prints the part of inertia explained by each mode, as well as the cumulative
        part of inertia. By default uses 10 modes
        Input:
            - modes: the numer of modes to describe
        """
        if modes>self.Katt:
            modes=self.Katt
        cumulative_explained_inertia = 0.
        print('Total inertia: {0:4.2e} (should be {1:d} if method==cor)'.format(self.total_inertia,self.Katt))
        print('Mode | Eigenvalue | Inertia explained (%) | Cumulative (%)')
        print('----------------------------------------------------------')
        for k in range(modes):
            eigenval = self.eigenval[k]
            explained_inertia = eigenval/self.total_inertia
            cumulative_explained_inertia += explained_inertia
            print('{0:3d}  |  {1:4.2e}  |        {2:6.3f}         |   {3:6.3f}'.format(\
                  k+1,eigenval,explained_inertia*100,cumulative_explained_inertia*100))
        return
