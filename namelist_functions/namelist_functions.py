### Package with some useful functions for smilei simulations

import numpy as np
import numpy.random as rnd

def random_select(prob=0.01):
    if prob<=0.0:
        raise Exception("You have set %f pobability of selection <=0" % prob)
    elif prob>=1.0:
        raise Exception("You have set %f probability of selection >=1" %prob)
    else:
        def _random_select(particles):
            return rnd.choice(a=[True, False], size=particles.weight.shape, p=[prob, 1-prob])
        return _random_select
## end random_select


def spatial_filter(xmin=None,xmax=None,ymin=None,ymax=None,zmin=None,zmax=None):
    return lambda p: _in_range(p.x,xmin,xmax) \
        * _in_range(p.y,ymin,ymax)\
        * _in_range(p.z,zmin,zmax)

def _in_range(s, minimum,maximum):
    if (minimum is not None) and (maximum is not None):
        return np.logical_and(minimum<=s, s<maximum)
    elif (minimum is not None) and (maximum is None):
        return minimum <= s
    elif (minimum is None) and (maximum is not None):
        return s < maximum
    else:
        return True # numpy.full(s.shape, True)
