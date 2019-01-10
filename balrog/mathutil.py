import numpy as np

#-------------------------------------------------------------------------------
# Unit conversions and misc functions

def deg2arcmin(val):
    return 60.0 * val

def arcmin2deg(val):
    return val / 60.0

def sample_uniform_ra(r1, r2, N=None, boundary_cross=False):
    '''
    Sample N random RA values from r1 to r2, where r1<r2. If boundary_cross
    is set to True, will adjust inputs to allow sampling over 360/0 boundary.
    '''

    # Adjust r1, r2 if they cross the 360/0 ra boundary
    if boundary_cross is True:
        # NOTE: ramin/ramax are saved the same way as the geometry file catalog,
        # which is to say ramin/ramax are always 'left'/'right' boundary, even if
        # ramax < ramin across boundary.
        # This means that r1, while 'min', is larger and needs to be adjusted down.
        r1_shift = r1 - 360.0
        shifted_dist = sample_uniform(r1_shift, r2, N)
        shifted_dist[shifted_dist < 0.0] += 360.0
        return shifted_dist

    else:
        # NOTE: Equivalent to below DEC procedure for ra = P(r2-r1)+r1 for P in (0,1)
        return sample_uniform(r1, r2, N)

def sample_uniform_dec(d1, d2, N=None, unit='deg'):
    '''
    Sample N random DEC values from d1 to d2, accounting for curvature of sky.
    '''

    if N is None: N=1

    # Convert dec to radians if not already
    if unit is 'deg' or unit is 'degree':
        d1, d2 = np.deg2rad(d1), np.deg2rad(d2)

    # Uniform sampling from 0 to 1
    P = np.random.rand(N)

    # Can't use `sample_uniform()` as dec needs angular weighting
    delta = np.arcsin(P * (np.sin(d2) - np.sin(d1)) +np.sin(d1))

    # Return correct unit
    if unit is 'deg' or unit is 'degree':
        return np.rad2deg(delta)
    elif unit is 'rad' or unit is 'radian':
        return delta
    else:
        raise TypeError('Only deg/degree or rad/radian are allowed as inputs!')

def sample_uniform_indx(n1, n2, N=None):
    'Samples N random indices from n1 to n2'

    return np.random.randint(n1, high=n2, size=N)

def sample_uniform(v1, v2, N=None):
    'Samples N random values from v1 to v2'

    if N is None: N=1

    return np.random.uniform(low=v1, high=v2, size=N)
