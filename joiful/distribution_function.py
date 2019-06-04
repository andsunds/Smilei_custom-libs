import numpy as np
import scipy as sp

def inRange(x,Min,Max):
    # Helperfunction for if x is in range
    # Min <= x < Max
    return np.logical_and(x>=Min, x<Max)

def get_dist1D(x_data,weights, xlim,Nx, **kwargs):
    
    ## Checks the keys
    for key in kwargs.keys():
        if key=='mask':
            mask=kwargs[key]
        elif key=='scale':
            scale=kwargs[key]
        else:
            raise KeyError('Key: %s not supported'%key)
    ## Default values
    if 'mask' not in locals():
        mask=np.full(weights.size, True, dtype=bool)
    if 'scale' not in locals():
        scale='lin'
    ##end default values
    
    mask=np.logical_and(mask, inRange(x_data,xlim[0],xlim[1]))
    IND=np.where(mask)[0]
    
    if scale=='log':
        sf=lambda x: np.log(x)
        sf_inv=lambda y: np.exp(y)
    elif scale=='lin':
        sf=lambda x: x
        sf_inv=lambda y: y
    ## end rangeMask

    ulim=sf(xlim[1])
    llim=sf(xlim[0])
    
    # The size of the bins in x and p:
    dx=(ulim-llim)/Nx;
    
    # Output init:
    f=np.zeros(Nx, dtype=np.float64);
    x_bins=np.linspace(llim,ulim,num=Nx);

    x_data_ind=np.array((sf(x_data)-llim)/dx, dtype=np.int)
        
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


def get_dist2D(x_data,p_data,weights, xlim,Nx, plim,Np, **kwargs):
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
        if key=='get_density':
            get_density=kwargs[key]
        elif key=='mask':
            mask=kwargs[key]
        else:
            raise KeyError('Key: %s not supported'%key)
    ## Default values
    if 'get_density' not in locals():
        get_density=False
    if 'mask' not in locals():
        mask=np.full(weights.size, True, dtype=bool)
        #IND=np.where(mask)[0]
    #else:
        #IND=range(weights.size)
    ##end default values

    if get_density:
        ## When we want the density, we do not want to exclude the
        ## particles outside the momentum range.
        rangeMask=inRange(x_data,xlim[0],xlim[1])
    else:
        rangeMask=np.logical_and(inRange(x_data,xlim[0],xlim[1]),
                                 inRange(p_data,plim[0],plim[1]))
    ## end rangeMask
    mask=np.logical_and(mask, rangeMask)
    IND=np.where(mask)[0]
    #N=weights.size #Number of data points
    

    # The size of the bins in x and p:
    dx=(xlim[1]-xlim[0])/Nx;
    dp=(plim[1]-plim[0])/Np;
    
    # Output init:
    f=np.zeros((Np,Nx), dtype=np.float64);
    x_bins=np.linspace(xlim[0],xlim[1],num=Nx);
    p_bins=np.linspace(plim[0],plim[1],num=Np);

    x_data_ind=np.array((x_data-xlim[0])/dx, dtype=np.int)
    p_data_ind=np.array((p_data-plim[0])/dp, dtype=np.int)
        
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


def get_p_moment_rel(x_data,p_data,weights, xlim, Nx, **kwargs):
    return get_1moment_rel(x_data,p_data,weights, xlim, Nx, **kwargs)

def get_1moment_rel(x_data,y_data,weights, xlim, Nx, **kwargs):
    # 
    # 
    #

    ## Checks the keys
    for key in kwargs.keys():
        if key=='xscale':
            xscale=kwargs[key]
        elif key=='mask':
            mask=kwargs[key]
        elif key=='fill_zeros':
            fill_zeros=kwargs[key]
        else:
            raise KeyError('Key: %s not supported'%key)
    ## Default values
    if 'xscale' not in locals():
        xscale=1
    if 'mask' not in locals():
        mask=np.full(weights.size, True, dtype=bool)
    if 'fill_zeros' not in locals():
        fill_zeros=True 
    ##end default values
    rangeMask=inRange(x_data,xlim[0],xlim[1])
    mask=np.logical_and(mask, rangeMask)
    IND=np.where(mask)[0]
    
    
    Y=np.zeros(Nx)
    n=np.zeros(Nx)
    dx=(xlim[1]-xlim[0])/Nx;
    for i in IND:
        # Each particle's position is converted to an interger index to
        # be used in the energy and density arrays.
        x_data_ind=int((x_data[i]-xlim[0])/dx); 
        n[x_data_ind] += weights[i]
        Y[x_data_ind] += weights[i] * y_data[i]
    ##end for
    n=np.ma.masked_where(n<=0, n/(dx*xscale))
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
    rangeMask=inRange(x_data,xlim[0],xlim[1])
    mask=np.logical_and(mask, rangeMask)
    IND=np.where(mask)[0]
    
    
    E=np.zeros(Nx)
    n=np.zeros(Nx)
    dx=(xlim[1]-xlim[0])/Nx;
    for i in IND:
        # Each particle's position is converted to an interger index to
        # be used in the energy and density arrays.
        x_data_ind=int((x_data[i]-xlim[0])/dx);       
        n[x_data_ind] += weights[i]
        E[x_data_ind] += weights[i] * (np.sqrt(1+p_data[i]**2)-1)
        #E[x_data_ind] += weights[i] * p_data[i]**2
    ##end for
    n=n/(dx*xscale)
    E=np.ma.masked_where(n==0, E*(Escale/(dx*xscale)))/n
    E=E.filled(fill_value=0.0)
    return E, n
##end def



######################################################################
################ Rather simple and outdated functions ################
######################################################################

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


