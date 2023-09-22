"""
This is where you set the run initial condition using get_ic()
based on an experiment name passed by the calling code.

Thre are also some utility functions useful for making different
common release patterns.

"""

import numpy as np
    
def get_ic(TR):
    # routines to set particle initial locations, all numpy arrays
    
    # NOTE: "pcs" refers to fractional depth, and goes linearly from -1 to 0
    # between the local bottom and free surface.  It is how we keep track of
    # vertical position, only converting to z-position when needed.
    
    exp_name = TR['exp_name']
    gridname = TR['gridname']
    fn00 = TR['fn00']
        
    if exp_name == 'jdf0': # Mid-Juan de Fuca
        lonvec = np.linspace(-123.85, -123.6, 20)
        latvec = np.linspace(48.2, 48.4, 20)
        pcs_vec = np.array([0])
        plon00, plat00, pcs00 = ic_from_meshgrid(lonvec, latvec, pcs_vec)

    elif exp_name == 'hc_dolph': # EAB: For MURI project detecting dolphins
        lon0 = -122.729779; lat0 = 47.742773 # center of the circle= dolphin pen
        lon_pad = 0.0001*0.5 # desired width of rectangle divided by 2 to center at lon0, lat0
        lat_pad = 0.000025*0.5
        lonvec = np.linspace(lon0-lon_pad,lon0+lon_pad,100)
        latvec = np.linspace(lat0-lat_pad,lat0+lat_pad,100)
        pcs_vec = np.array([0])
        plon00, plat00, pcs00 = ic_from_meshgrid(lonvec, latvec, pcs_vec)
        #radius_km = 0.005 # radius of the circle km, for 10m pen, radius 5m?
        #N = 10000 # number of particles - TEST, next time increase to 10k
        # make random scattering of points in a circle
        #plon00, plat00 = ic_random_in_circle(lon0, lat0, radius_km, N)
        #pcs00 = 0.0* np.ones(N) #surface release
        
    return plon00, plat00, pcs00
    
def ic_from_meshgrid(lonvec, latvec, pcs_vec):
    # First create three vectors of initial locations (as done in some cases above).
    # plat00 and plon00 should be the same length, and the length of pcs00 is
    # as many vertical positions you have at each lat, lon
    # (expressed as fraction of depth -1 < pcs < 0).
    # Then we create full output vectors (each has one value per point).
    # This code takes each lat, lon location and then assigns it to NSP points
    # corresponding to the vector of pcs values.
    lonmat, latmat = np.meshgrid(lonvec, latvec)
    plon_vec = lonmat.flatten()
    plat_vec = latmat.flatten()
    if len(plon_vec) != len(plat_vec):
        print('WARNING: Problem with length of initial lat, lon vectors')
    NSP = len(pcs_vec)
    NXYP = len(plon_vec)
    plon_arr = plon_vec.reshape(NXYP,1) * np.ones((NXYP,NSP))
    plat_arr = plat_vec.reshape(NXYP,1) * np.ones((NXYP,NSP))
    pcs_arr = np.ones((NXYP,NSP)) * pcs_vec.reshape(1,NSP)
    plon00 = plon_arr.flatten()
    plat00 = plat_arr.flatten()
    pcs00 = pcs_arr.flatten()
    return plon00, plat00, pcs00
    
def ic_from_list(lonvec, latvec, pcs_vec):
    # Like ic_from_meshgrid() but treats the lon, lat lists like lists of mooring locations.
    plon_vec = lonvec
    plat_vec = latvec
    if len(plon_vec) != len(plat_vec):
        print('WARNING: Problem with length of initial lat, lon lists')
    NSP = len(pcs_vec)
    NXYP = len(plon_vec)
    plon_arr = plon_vec.reshape(NXYP,1) * np.ones((NXYP,NSP))
    plat_arr = plat_vec.reshape(NXYP,1) * np.ones((NXYP,NSP))
    pcs_arr = np.ones((NXYP,NSP)) * pcs_vec.reshape(1,NSP)
    plon00 = plon_arr.flatten()
    plat00 = plat_arr.flatten()
    pcs00 = pcs_arr.flatten()
    return plon00, plat00, pcs00
    
def ic_random_in_circle(lon0, lat0, radius_km, npoints):
    # Makes lon and lat of npoints scattered randomly in a circle.
    # I think the np.sqrt() used in calculating the radius makes these
    # evenly distributed over the whole circle.
    earth_r = 6371 # average earth radius [km]
    # radius of the circle km
    circle_r = radius_km
    # center of the circle (x, y)
    circle_x = lon0
    circle_y = lat0
    N = npoints # number of particles
    # random angle
    alpha = 2 * np.pi * np.random.rand(N)
    # random radius
    r = (circle_r/earth_r) * (180/np.pi) * np.sqrt(np.random.rand(N))
    # calculating coordinates
    plon00 = r * np.cos(alpha) / np.cos(circle_y*np.pi/180) + circle_x
    plat00 = r * np.sin(alpha) + circle_y
    # we leave it to the user to make pcs00
    return plon00, plat00
    
def ic_from_TEFsegs(fn00, gridname, seg_list, DZ, NPmax=10000):
    import pickle
    import sys
    # select the indir
    from lo_tools import Lfun, zrfun
    Ldir = Lfun.Lstart()
    indir = Ldir['LOo'] / 'tef' / ('volumes_' + gridname)
    # load data
    j_dict = pickle.load(open(indir / 'j_dict.p', 'rb'))
    i_dict = pickle.load(open(indir / 'i_dict.p', 'rb'))
    G = zrfun.get_basic_info(fn00, only_G=True)
    h = G['h']
    xp = G['lon_rho']
    yp = G['lat_rho']
    plon_vec = np.array([])
    plat_vec = np.array([])
    hh_vec = np.array([])
    for seg_name in seg_list:
        jjj = j_dict[seg_name]
        iii = i_dict[seg_name]
        # untested 2021.10.05
        hh_vec = np.append(hh_vec, h[jjj,iii])
        plon_vec = np.append(plon_vec, xp[jjj,iii])
        plat_vec = np.append(plat_vec, yp[jjj,iii])
        # ji_seg = ji_dict[seg_name]
        # for ji in ji_seg:
        #     plon_vec = np.append(plon_vec, xp[ji])
        #     plat_vec = np.append(plat_vec, yp[ji])
        #     hh_vec = np.append(hh_vec, h[ji])
    plon00 = np.array([]); plat00 = np.array([]); pcs00 = np.array([])
    for ii in range(len(plon_vec)):
        x = plon_vec[ii]
        y = plat_vec[ii]
        hdz = DZ*np.floor(hh_vec[ii]/DZ) # depth to closest DZ m (above the bottom)
        if hdz >= DZ:
            zvec = np.arange(-hdz,DZ,DZ) # a vector that goes from -hdz to 0 in steps of DZ m
            svec = zvec/hh_vec[ii]
            ns = len(svec)
            if ns > 0:
                plon00 = np.append(plon00, x*np.ones(ns))
                plat00 = np.append(plat00, y*np.ones(ns))
                pcs00 = np.append(pcs00, svec)
    # subsample the I.C. vectors to around max length around NPmax
    NP = len(plon00)
    print(len(plon00))
    nstep = max(1,int(NP/NPmax))
    plon00 = plon00[::nstep]
    plat00 = plat00[::nstep]
    pcs00 = pcs00[::nstep]
    print(len(plon00))
    sys.stdout.flush()
    return plon00, plat00, pcs00
    
def ic_sect(fn00, lons, lats, NPmax=10000):
    """
    This distributes NPmax particles evenly on a section defined by endpoints
    (lon0, lat0) - (lon1, lat1).
    
    For simplicity we force the section to be NS or EW, we put particles
    only on rho points.
    """
    from lo_tools import Lfun, zfun, zrfun

    Ldir = Lfun.Lstart()

    G = zrfun.get_basic_info(fn00, only_G=True)
    h = G['h']
    m = G['mask_rho']
    xr = G['lon_rho']
    yr = G['lat_rho']
    X = xr[0,:]
    Y = yr[:,0]

    lon0 = lons[0]; lon1 = lons[1]
    lat0 = lats[0]; lat1 = lats[1]

    ix0 = zfun.find_nearest_ind(X, lon0)
    ix1 = zfun.find_nearest_ind(X, lon1)
    iy0 = zfun.find_nearest_ind(Y, lat0)
    iy1 = zfun.find_nearest_ind(Y, lat1)

    # adjust indices to make it perfectly zonal or meridional
    dix = np.abs(ix1 - ix0)
    diy = np.abs(iy1 - iy0)
    if dix > diy: # EW section
        iy1 = iy0
    elif diy > dix: # NS section
        ix1 = ix0
    
    hvec = h[iy0:iy1+1, ix0:ix1+1].squeeze()
    mvec = m[iy0:iy1+1, ix0:ix1+1].squeeze()
    xvec = xr[iy0:iy1+1, ix0:ix1+1].squeeze()
    yvec = yr[iy0:iy1+1, ix0:ix1+1].squeeze()

    # add up total depth of water
    hnet = 0
    for ii in range(len(hvec)):
        if mvec[ii] == 1:
            hnet += hvec[ii]
    p_per_meter = NPmax/hnet
        
    # initialize result arrays
    plon00 = np.array([]); plat00 = np.array([]); pcs00 = np.array([])
    for ii in range(len(hvec)):
        if mvec[ii] == 1:
            this_h = hvec[ii]
            this_np = int(np.floor(p_per_meter * this_h))
            plon00 = np.concatenate((plon00,xvec[ii]*np.ones(this_np)))
            plat00 = np.concatenate((plat00,yvec[ii]*np.ones(this_np)))
            pcs00 = np.concatenate((pcs00,np.linspace(-1,0,this_np)))
    return plon00, plat00, pcs00
    
