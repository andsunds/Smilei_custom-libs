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
