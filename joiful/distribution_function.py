### Tools for binning particles into distribution function from 
### individual particle data.

######################################################################
# Copyright 2019-2020 ANDRÉAS SUNDSTRÖM
#
# This file is part of Joiful.
#
# Joiful is free software: you can redistribute it and/or modify it
# under the terms of the GNU General Public License (version 3) as
# published by the Free Software Foundation.
#
# Joiful is distributed in the hope that it will be useful, but
# WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
# General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with Joiful.  If not, see <http://www.gnu.org/licenses/>.
######################################################################

import numpy as np
import scipy as sp
import warnings


def inRange(x,Min,Max,shift_coordinates=False,delta=0):
    # Helperfunction for if x is in range
    # Min <= x < Max
    with warnings.catch_warnings():
        warnings.filterwarnings('ignore')
        
        if shift_coordinates or delta==0: return np.logical_and(x>=Min, x<Max)
        else: return np.logical_and(x>Min-0.5*delta, x<Max-0.5*delta)



def get_dist1D(x_data,weights, xlim,Nx,
               scale='lin',shift_coordinates=False, **kwargs):
    
    ## Checks the keys
    for key in kwargs.keys():
        if key=='mask':
            mask=kwargs[key]
        else:
            raise KeyError('Key: %s not supported'%key)
    ## Default values
    if 'mask' not in locals():
        mask=np.full(weights.size, True, dtype=bool)
    ##end default values

    # Scale functions
    if scale=='log':
        sf=lambda x: np.log(x)
        sf_inv=lambda y: np.exp(y)
    elif scale=='lin':
        sf=lambda x: x
        sf_inv=lambda y: y

    # Limits and bin sizes
    ulim=sf(xlim[1])
    llim=sf(xlim[0])
    dx=(ulim-llim)/Nx;
    
    mask=np.logical_and(mask, inRange(sf(x_data),llim,ulim,shift_coordinates,delta=dx))
    IND=np.where(mask)[0]
    ## end rangeMask

    
    # Output init:
    f=np.zeros(Nx, dtype=np.float64);
    
    if shift_coordinates:
        x_bins=np.linspace(llim,ulim,num=Nx,endpoint=False) + 0.5*dx
        x_data_ind=np.array(np.floor((sf(x_data)-llim)/dx), dtype=np.int)
    else:
        x_bins=np.linspace(llim,ulim,num=Nx);
        x_data_ind=np.array(np.around((sf(x_data)-llim)/dx), dtype=np.int)
        
    # Loop over every particle, adding its weight to the distribution
    # binning given by the indices given in x_data_ind.
    #
    # If a particle is outside the bounds set by xlim, then it is
    # ignored.  Each particle data is converted to an interger index
    # to be used in the distribution function matrix. The integer
    # conversion is done such that e.g. a position, x_data, will fall
    # under the bin, x_bin, where
    #   (x_bin-dx/2) <= x_data < (x_bin+dx/2).
    for i in IND:
        f[x_data_ind[i]] += weights[i];
    ##end for
    return f, sf_inv(x_bins)
##end def


def get_dist2D(x_data,p_data,weights, xlim,Nx, plim,Np,
               get_density=False,shift_coordinates=True, **kwargs):
    # Calcualtes the distribution function on a uniform x-p bin-grid based on
    # each particle phase-space coordinates and weights.
    # 
    # Input arguments:
    # x_data - [N*1]
    #    containing each particle's configuration space coord.
    # p_data - [N*1]
    #    N element vector containing each particle's momentum space coord.
    # weights - [N*1]
    #    N element vector containing each particle's weight
    # xlim - [2*1]
    #    Then centers of the limiting bnis in configuration space
    # Nx
    #    The number of bins in configuration space
    # plim - [2*1]
    #    Then centers of the limiting bnis in momentum space
    # Np
    #    The number of bins in momentum space
    # 
    # Output arguments:
    # f - [Nx*Np]
    #    The distribution function based on the particles and their weights
    # x_bins - [Nx*1]
    #    The centers of each bin in configuration space
    # p_bins - [Np*1]
    #    The centers of each bin in momentum space
    #
    #
    #
    
    ## Checks the keys
    for key in kwargs.keys():
        if key=='mask':
            mask=kwargs[key]
        else:
            raise KeyError('Key: %s not supported'%key)
    ## Default values
    if 'mask' not in locals(): mask=np.full(weights.size, True, dtype=bool)
    ##end default values

    # The size of the bins in x and p:
    dx=(xlim[1]-xlim[0])/Nx;
    dp=(plim[1]-plim[0])/Np;

    if get_density:
        ## When we want the density, we do not want to exclude the
        ## particles outside the momentum range.
        rangeMask=inRange(x_data,xlim[0],xlim[1],shift_coordinates,dx)
    else:
        rangeMask=np.logical_and(inRange(x_data,xlim[0],xlim[1],shift_coordinates,dx),
                                 inRange(p_data,plim[0],plim[1],shift_coordinates,dp))
    ## end rangeMask
    mask=np.logical_and(mask, rangeMask)
    IND=np.where(mask)[0]
    #N=weights.size #Number of data points

    
    # Output init:
    f=np.zeros((Np,Nx), dtype=np.float64)

    ## Shift coordinates to be used with pcolor or other "pixel"
    ## drawing function, where the coordinates are the corners of the
    ## pixels.
    if shift_coordinates:
        x_bins=np.linspace(xlim[0],xlim[1],num=Nx+1)
        p_bins=np.linspace(plim[0],plim[1],num=Np+1)
        x_data_ind=np.array(np.floor((x_data-xlim[0])/dx), dtype=np.int)
        p_data_ind=np.array(np.floor((p_data-plim[0])/dp), dtype=np.int)
    else:
        x_bins=np.linspace(xlim[0],xlim[1],num=Nx)
        p_bins=np.linspace(plim[0],plim[1],num=Np)
        x_data_ind=np.array(np.floor((x_data-xlim[0])/dx), dtype=np.int)
        p_data_ind=np.array(np.floor((p_data-plim[0])/dp), dtype=np.int)


        
    # Loop over every particle, adding its weight to the distribution function
    # matrix element given by the indices given in x_data_ind and p_data_ind.
    # If a particle is outside the bounds set by xlim and plim, then it is
    # ignored.
    # Each particle position and momentum is converted to an interger index to
    # be used in the distribution function matrix.
    # The integer conversion is done such that e.g. a position, x_data, will
    # fall under the bin, x_bin, where 
    #   (x_bin-dx/2) <= x_data < (x_bin+dx/2).
    if get_density:
        n=np.zeros(Nx, dtype=np.float64);
        for i in IND:
            n[x_data_ind[i]]+=weights[i]
            if inRange(p_data_ind[i],0,Np):
                f[p_data_ind[i], x_data_ind[i]] += weights[i];
        ##end for
        return f, n/dx, x_bins, p_bins
    else: # not get_density
        for i in IND:
            f[p_data_ind[i], x_data_ind[i]] += weights[i];
        ##end for
        return f, x_bins, p_bins
##end def


def get_1moment_rel(x_data,y_data,weights, xlim, Nx,
                    fill_zeros=True,xscale=1, **kwargs):
    # 
    # 
    #

    ## Checks the keys
    for key in kwargs.keys():
        if key=='mask':
            mask=kwargs[key]
        else:
            raise KeyError('Key: %s not supported'%key)
    ## Default values

    if 'mask' not in locals():
        mask=np.full(weights.size, True, dtype=bool)
    ##end default values
    if len(y_data.shape)>1:
        Y=np.zeros((Nx,y_data.shape[1]))
    else: Y=np.zeros(Nx)
    n=np.zeros(Nx)
    dx=(xlim[1]-xlim[0])/Nx;
    
    rangeMask=inRange(x_data,xlim[0],xlim[1],shift_coordinates=False,delta=dx)
    mask=np.logical_and(mask, rangeMask)
    IND=np.where(mask)[0]
    
    for i in IND:
        # Each particle's position is converted to an interger index to
        # be used in the energy and density arrays.
        x_data_ind=int(np.around((x_data[i]-xlim[0])/dx))
        n[x_data_ind] += weights[i]
        Y[x_data_ind] += weights[i] * y_data[i]
    ##end for
    n=np.ma.masked_where(n<=0, n/(dx*xscale))
    if len(y_data.shape)>1:
        nY=n.copy()
        nY=np.tile(nY,(Y.shape[1],1)).transpose()
        Y=Y/(dx*xscale * nY)
    else:
        Y=Y/(dx*xscale * n)
    if fill_zeros:
        return Y.filled(fill_value=0.0), n.filled(fill_value=0.0)
    else:
        return Y, n
##end def

def get_E_moment_rel(x_data,p_data,weights, xlim, Nx, **kwargs):
    # 
    # 
    #
    # E=mc^2 * sqrt(1+(p/mc)^2)

    ## Checks the keys
    for key in kwargs.keys():
        if key=='xscale':
            xscale=kwargs[key]
        elif key=='Escale':
            Escale=kwargs[key]
        elif key=='mask':
            mask=kwargs[key]
        else:
            raise KeyError('Key: %s not supported'%key)
    ## Default values
    if 'xscale' not in locals():
        xscale=1
    if 'Escale' not in locals():
        Escale=1
    if 'mask' not in locals():
        mask=np.full(weights.size, True, dtype=bool)
    ##end default values

    E=np.zeros(Nx)
    n=np.zeros(Nx)
    dx=(xlim[1]-xlim[0])/Nx;
    
    rangeMask=inRange(x_data,xlim[0],xlim[1],shift_coordinates=False,delta=dx)
    mask=np.logical_and(mask, rangeMask)
    IND=np.where(mask)[0]
    

    for i in IND:
        # Each particle's position is converted to an interger index to
        # be used in the energy and density arrays.
        x_data_ind=int(np.around((x_data[i]-xlim[0])/dx))      
        n[x_data_ind] += weights[i]
        E[x_data_ind] += weights[i] * (np.sqrt(1+p_data[i]**2)-1)
        #E[x_data_ind] += weights[i] * p_data[i]**2
    ##end for
    n=n/(dx*xscale)
    E=np.ma.masked_where(n==0, E*(Escale/(dx*xscale)))/n
    E=E.filled(fill_value=0.0)
    return E, n
##end def



