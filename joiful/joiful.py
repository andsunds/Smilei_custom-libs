### Joiful -- Just anOther Independent FUnction Library
### My custom package for processing Smilei data extracted using Happi
import numpy as np
import scipy as sp

def inRange(x,Min,Max):
    # Helperfunction for if x is in range
    # Min <= x < Max
    return x>=Min and x<Max


def get_dist2D(x_data,p_data,weights, xlim,Nx, plim,Np):
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
    
    #print("joiful test") ##DEBUG
    
    N=x_data.size #Number of data points

    # The size of the bins in x and p:
    dx=(xlim[1]-xlim[0])/Nx;
    dp=(plim[1]-plim[0])/Np;
    
    # Output init:
    f=np.zeros((Np,Nx), dtype=np.float64);
    x_bins=np.linspace(xlim[0],xlim[1],num=Nx);
    p_bins=np.linspace(plim[0],plim[1],num=Np);
    
    # Loop over every particle, adding its weight to the distribution function
    # matrix element given by the indices given in x_data_ind and p_data_ind.
    # If a particle is outside the bounds set by xlim and plim, then it is
    # ignored.
    for i in range(N):
        # Each particle position and momentum is converted to an interger index to
        # be used in the distribution function matrix.
        # The integer conversion is done such that e.g. a position, x_data, will
        # fall under the bin, x_bin, where 
        #   (x_bin-dx/2) <= x_data < (x_bin+dx/2).
        x_data_ind=int((x_data[i]-xlim[0])/dx); 
        p_data_ind=int((p_data[i]-plim[0])/dp);
        if inRange(x_data_ind,0,Nx) and inRange(p_data_ind,0,Np):
            f[p_data_ind, x_data_ind] += weights[i];
        #end if
    #end for
    return f, x_bins, p_bins
#end def


def get_density(x_data,weights, xlim):
    # Calcualtes the density of particles within the range given in xlim.
    #
    # Input arguments:
    # x_data - [N*1]
    #    containing each particle's configuration space coord.
    # weights - [N*1]
    #    N element vector containing each particle's weight
    # xlim - [2*1]
    #    Then centers of the limiting bnis in configuration space
    #
    
    N=x_data.size #Number of data points
    n=0
    
    for i in range(N):
        if inRange(x_data[i],xlim[0],xlim[1]):
            n+=weights[i]
    # end for
    return n/(xlim[1]-xlim[0])

def get_p1_moment_nonRel(x_data,p_data,weights, xlim):
    # Calculates the first momentum moment of the particles within
    # the rage given in xlim.
    #
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
    #
    
    N=x_data.size #Number of data points
    n=0
    P=0
    for i in range(N):
        if inRange(x_data[i],xlim[0],xlim[1]):
            n+=weights[i]
            P+=weights[i]*p_data[i]
    # end for
    dx=(xlim[1]-xlim[0])
    n=n/dx
    P=(P/dx)/n
    return P, n

def get_p2_moment_nonRel(x_data,p_data,weights, xlim):
    # Calculates the second momentum moment of the particles within
    # the rage given in xlim.
    #
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
    #
    
    N=x_data.size #Number of data points
    P, n = get_p1_moment_nonRel(x_data,p_data,weights, xlim)
    T=0
    for i in range(N):
        if inRange(x_data[i],xlim[0],xlim[1]):
            T+=weights[i]*(p_data[i]-P)**2
    #end for
    dx=(xlim[1]-xlim[0])
    T=(T/dx)/n
    return T, P, n
