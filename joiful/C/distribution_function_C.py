import numpy as np
from binning import bin2Dw

#bin2Dw(x_vals, p_vals, weights, xmin, xmax, pmin, pmax, nxbins, npxbins)

def get_dist2D(x_data,p_data,weights, xlim,Nx, plim,Np,
               get_density=False,shift_coordinates=True, **kwargs):
    
    # The size of the bins in x and p:
    dx=(xlim[1]-xlim[0])/Nx;
    dp=(plim[1]-plim[0])/Np;
    
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

    if get_density:
        raise NotImplementedError(""" 
You have selected `get_density`, however, this option 
has not been implemeted yet.""")
    
    f=bin2Dw(x_data, p_data, weights, xlim[0], xlim[1], plim[0], plim[1], Nx, Np)
    return f, x_bins, p_bins
