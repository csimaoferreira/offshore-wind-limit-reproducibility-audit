from wind_stats import get_gwc_data
import pandas as pd
import xarray as xr
import numpy as np
import re
import sys
from scipy.special import gamma, gammainc
from scipy.optimize import fsolve
import windkit as wk
from pyproj import Proj, CRS
from matplotlib.patches import Polygon as Polygonplt
import alphashape
from shapely.geometry import Polygon
from shapely.geometry.polygon import orient
from shapely.geometry import LineString, Point
import matplotlib.pyplot as plt
import matplotlib.style
import matplotlib as mpl
mpl.style.use('classic')
# Latex font
plt.rcParams['font.family'] = 'STIXGeneral'
plt.title(r'ABC123 vs $\mathrm{ABC123}^{123}$')
plt.rcParams["legend.scatterpoints"] = 1
plt.rcParams["legend.numpoints"] = 1
  
      
# From PyWake v2.6.18 
# Added:
# Nfree = np.minimum(Nfree, Nturb)
# get_phi solve function
class MinimalisticPredictionModel():
    """Sørensen, J.N.; Larsen, G.C.
    A Minimalistic Prediction Model to Determine Energy Production and Costs of Offshore Wind Farms.
    Energies 2021, 14, 448. https://doi.org/10.3390/en14020448"""

    def __init__(self, correction_factor, latitude, CP, Uin, Uout, rho):
        """
        Parameters
        ----------
        correction_factor : int, float or function
            Finite-size wind farm corrrection which multiplied with sqrt(Nturb) gives
            the number of wind turbines exposed to the free wind
        latitude : int or float
            latitude [deg] used to calculate the coriolis parameter
        CP : float, optional
            Wind turbine power coefficient
        Uin : int or float, optional
            Wind turbine cut-in wind speed
        Uout : int or float, optional
            Wind turbine cut-out wind speed
        rho : float, optional
            Air density
        """

        self.CP = CP
        self.Uin = Uin
        self.Uout = Uout
        omega = 2 * np.pi / (24 * 60 * 60)  # earth rotation speed
        self.f = 2 * omega * np.sin(np.deg2rad(latitude))
        self.correction_factor = correction_factor
        self.rho = rho
        
    def predict(self, Pg, CT, D, H, z0, Aw, kw, Nturb, Area):
        """
        Inputs:
            Pg    - [W] Nameplate capacity (generator power)
            CT    - [-] Thrust coefficient
            D     - [m] Rotor diameter
            H     - [m] Tower height
            z0    - [m] roughness length
            Aw    - [m/s] Weibull scale parameter
            kw    - [-] Weibull shape parameter
            Nturb - [-] Number of turbines
            Area  - [m2] Area of wind farm

        Outputs:
            power - [Wh] Annual energy production of the wind farm
            ws_eff - [m/s] Effective mean wind speed including wakes
        """

        kappa = 0.4  # [-] Von Karman constant
        Uin, Uout = self.Uin, self.Uout

        # factor defined by Frandsen, should be used instead of f in eq 13 and 19 (typos in paper)
        fm = self.f * np.exp(4)
        delta = np.log(H / z0)  # eq 19

        # Mean spacing between wt in diameters, eq 8
        S = np.sqrt(Area) / (D * (np.sqrt(Nturb) - 1))

        # Rated wind speed [m/s], eq 4
        Ur = (8 * Pg / (self.rho * np.pi * D**2 * self.CP))**(1 / 3)

        # Power modeled as P = alpha * U^3 + beta, eq 1
        alpha = Pg / (Ur**3 - Uin**3)  # [(m/s)^-3] eq 2
        beta = -Pg * Uin**3 / (Ur**3 - Uin**3)  # [-], eq 2

        Uh0 = Aw * gamma(1 + 1 / kw)  # [m/s] Mean velocity at hub height
        Ctau = np.pi * CT / (8 * S * S)  # [-] Wake parameter, rotor ct smeared on WT area
        nu = np.sqrt(0.5 * Ctau) * D / (kappa**2 * H) * delta  # [-] wake eddy viscosity

        # Finite-size wind farm corrrection, section 2.5
        correction_factor = self.correction_factor
        if hasattr(correction_factor, '__call__'):
            correction_factor = correction_factor(Uh0, S, Nturb)
        Nfree = correction_factor * np.sqrt(Nturb)  # Number of wt exposed to the free wind
        Nfree = np.minimum(Nfree, Nturb)  # To make sure Nfree/Nturb is not larger than one, not yet implemented in PyWake

        # Geostrophic wind speed
        G_last = Uh0
        for n in range(10):
            G = Uh0 * (1 + np.log(G_last / (fm * H)) / delta)
            dG = abs(G - G_last)
            if dG < 1e-5:
                break
            G_last = G

        gam = np.log(G / (fm * H))  # eq 19

        # Mean velocity at hub height without wake effects from geostrophic wind
        Uh0 = G / (1 + gam / kappa * np.sqrt((kappa / delta)**2))  # eq 13, ct=0

        # Power without wake effects, eq 16 modified by
        # - add gamma(1 + 3 / kw) to cancel out normalization in scipy's gammainc
        # - gammainc terms swapped (typo in paper) 
        def get_Py(Aw, Aw_out):  # Yearly power
            return alpha * Aw**3 * gamma(1 + 3 / kw) * (gammainc(1 + 3 / kw, (Ur / Aw)**kw) - gammainc(1 + 3 / kw, (Uin / Aw)**kw)) +\
                beta * (np.exp(-(Uin / Aw)**kw) - np.exp(-(Ur / Aw)**kw)) + \
                Pg * (np.exp(-(Ur / Aw)**kw) - np.exp(-(Uout / Aw_out)**kw))

        P_y = get_Py(Aw, Aw)
 
        # Without cutin and cutout
        # def get_Py(Aw):  # Yearly power
        #     x = (Ur / Aw) ** kw
        #     ks = 1 + 3 / kw
        #     return Pg * (x ** (1.0 - ks) * gamma(ks) * gammainc(ks, x) + np.exp(-x))
        # def get_Py(Aw):  # Yearly power
        #     x = (Ur / Aw)
        #     ks = 1 + 3 / kw
        #     return Pg * (x ** (-3) * gamma(ks) * gammainc(ks, x ** kw) + np.exp(-x ** kw))                      
        # P_y = get_Py(Aw)

        # Mean velocity at hub height with wake effects
        z0_lo = z0  # / (1 - D / (2 * H))**(nu / (1 + nu))  # ???
        Uh = G / (1 + gam * np.sqrt(Ctau + (kappa / np.log(H / z0_lo))**2) / kappa)

        # eq 18. The paper states 3/2 instead of 3.2 which is either a typo or an initial guess
        # eps2 corresponds to eps(Uout) in paper and eps2(Ur)=eps1
        eps1 = (1 + gam / delta) / (1 + gam / kappa * np.sqrt(Ctau + (kappa / delta)**2))
        eps2 = (1 + gam / delta) / (1 + gam / kappa * np.sqrt(Ctau * (Ur / Uh)**3.2 + (kappa / delta)**2))

        # print('e', eps1, correction_factor)
        # Power production with wake effects
        P_WFy = get_Py(eps1 * Aw, eps2 * Aw)

        power = ((Nturb - Nfree) * P_WFy + Nfree * P_y)
        ws_eff = ((Nturb - Nfree) * Uh + Nfree * Uh0) / Nturb
        
        # Phi without cutin and cutout 
        def get_phi(x):  # Yearly power
            # x = U_r / (Uh0 eps)
            ks = 1 + 3 / kw      
            Gm = gamma(1 + 1 / kw)
            return Pg * ((x * Gm) ** (-3) * gamma(ks) * gammainc(ks, (x * Gm) ** kw) + np.exp(-(x * Gm) ** kw)) - power / Nturb
            
        root = fsolve(get_phi, x0=1)  # initial guess        
        phi = root[0]

        return power, ws_eff, P_y, P_WFy, Uh, G, phi
   
       
def calculate_minimalistic_model(inputdata, lat, CP, CT, z0, kw, Uin, Uout, rho, Href, Nwt_row_wf_neighbors):
    nWF = len(inputdata['Name'])
    wfs = inputdata['Name']
    Pgs = inputdata['PG'] * 1e6
    Ds = inputdata['D']
    Hs = inputdata['Ht']
    # It seems that Simão Ferreira et al. (2026) already corrected for this
    # Awref = inputdata['lambda']    
    # Aws = Awref * np.log(Hs / z0) / np.log(Href / z0)
    Aws = inputdata['lambda']
    Nturbs = inputdata['Nt']
    Areas = inputdata['A'] * 1e6  
    WFPrated = inputdata['MW install'] * 1e6           
    # Finite wind farm correction factor ref 15
    NFree = inputdata['Nrowscale parameter'] * inputdata['Nturb_frontal_area']
    correction_factor = NFree / np.sqrt(Nturbs)    
    correction_factor_redo = 2.5 * np.asarray(Nwt_row_wf_neighbors) / np.sqrt(Nturbs)
    Cf_Free = np.zeros((nWF))
    Cf_Inf = np.zeros((nWF))
    Cf_Model = np.zeros((nWF))    
    U_Free = Aws * gamma(1 + 1 / kw)      
    U_Inf = np.zeros((nWF))    
    U_Model = np.zeros((nWF))  
    G_Model = np.zeros((nWF))     
    phi_Model = np.zeros((nWF))    
    phi_Model_redo = np.zeros((nWF))     
    phi_Model_a5p3 = np.zeros((nWF))     
    Cf_Model_redo = np.zeros((nWF))  
    Cf_Model_a5p3 = np.zeros((nWF))        
    for i in range(nWF):  
        wfm = MinimalisticPredictionModel(correction_factor[i], lat, CP, Uin, Uout, rho)
        power, wseff, P_y, P_WFy, Uh, G, phi = wfm.predict(Pgs[i], CT, Ds[i], Hs[i], z0, Aws[i], kw, Nturbs[i], Areas[i])
        Cf_Free[i] = P_y / Pgs[i]
        Cf_Inf[i] = P_WFy / Pgs[i]
        Cf_Model[i] = power / WFPrated[i]        
        U_Inf[i] = Uh
        U_Model[i] = wseff            
        G_Model[i] = G       
        phi_Model[i] = phi
        
        # With my freestream wts method
        wfm = MinimalisticPredictionModel(correction_factor_redo[i], lat, CP, Uin, Uout, rho)
        power, wseff, P_y, P_WFy, Uh, G, phi = wfm.predict(Pgs[i], CT, Ds[i], Hs[i], z0, Aws[i], kw, Nturbs[i], Areas[i])      
        phi_Model_redo[i] = phi  
        Cf_Model_redo[i] = power / WFPrated[i]
         
        # Using fixed correction factor, a=5.3            
        wfm = MinimalisticPredictionModel(5.3, lat, CP, Uin, Uout, rho)
        power, wseff, P_y, P_WFy, Uh, G, phi = wfm.predict(Pgs[i], CT, Ds[i], Hs[i], z0, Aws[i], kw, Nturbs[i], Areas[i])      
        phi_Model_a5p3[i] = phi  
        Cf_Model_a5p3[i] = power / WFPrated[i]          
    return Cf_Free, Cf_Inf, Cf_Model, U_Free, U_Inf, U_Model, G_Model, phi_Model, phi_Model_redo, phi_Model_a5p3, Cf_Model_redo, Cf_Model_a5p3

      
# Mean wind speed from Weibull distribution
def weibull_mean(A, k):
    return A * gamma(1 + 1 / k)
    
    
def get_cf_an(Ur_A, k):
    # Capacity factor of a single turbine, analytic solution, without cutin and cutout ws
    x = Ur_A ** k
    ks = 1 + 3 / k
    return x ** (1 - ks) * gamma(ks) *  gammainc(ks, x) + np.exp(-x)

            
def dms_string_to_decimal(dms_string):
    pattern = r"(\d+)°(\d+)′(\d+)″([NSEW])"
    matches = re.findall(pattern, dms_string)
    
    coords = []
    for deg, minute, sec, direction in matches:
        decimal = int(deg) + int(minute)/60 + int(sec)/3600
        if direction in ['S', 'W']:
            decimal = -decimal
        coords.append(decimal)
    
    return coords[0], coords[1]


def get_latlon_decimal(coords):
    nWF = len(coords)
    latlon = np.zeros((nWF, 2)) 
    for i in range(nWF):
        latlon[i, 0] , latlon[i, 1]  = dms_string_to_decimal(inputdata['Coordinates'][i])   
    return latlon


def get_gwa_windresource(latlon, wfnames, rerun=False, filename='data/GWA.nc'):
    if rerun:
        nWF = len(wfnames)
        # nWF = 1
        WFindex = np.linspace(0, nWF - 1, nWF, dtype=int)
        gwc_data = []
        for i in range(nWF): 
            gwc_data.append(get_gwc_data(latlon[i, 0], latlon[i, 1]))
   
        # Store in netCDF file 
        ds = xr.concat(gwc_data, dim='wf').assign_coords(wf=WFindex)
        ds['wfname'] = xr.DataArray(wfnames[0:nWF], [('wf', WFindex)], attrs={'long_name': 'Wind farm name'})
        ds['lat'] = xr.DataArray(latlon[:, 0], [('wf', WFindex)], attrs={'long_name': 'Latitude'})
        ds['lon'] = xr.DataArray(latlon[:, 1], [('wf', WFindex)], attrs={'long_name': 'Longitude'})
        ds.attrs.pop('coordinates', None)
        ds.to_netcdf(path=filename)
    else:
        ds = xr.open_dataset(filename)
    return ds


def mean_wind_direction(directions_deg, freq):
    '''
    Mean wind direction based on wind rose
    '''
    freq = freq / np.sum(freq)
    directions_rad = np.deg2rad(directions_deg)
    
    u = np.cos(directions_rad)
    v = np.sin(directions_rad)
    
    u_mean = np.mean(u * freq)
    v_mean = np.mean(v * freq)
    
    mean_dir = np.arctan2(v_mean, u_mean)
    mean_dir_deg = np.rad2deg(mean_dir)
    
    return (mean_dir_deg + 360) % 360
    

def calc_windresource(ds, z0, zRef, zH):
    nWF =  len(ds.wfname)
    ARef = np.zeros((nWF))
    AH = np.zeros((nWF))
    kRef = np.zeros((nWF))
    windrose = np.zeros((nWF, 12))
    for i in range(nWF):         
        # Calculate wind rose frequency weighted A and k
        ds_sub = ds.sel(wf=i, height=zRef).interp(roughness=0.0)
        ds_sub = ds_sub.rename({'frequency': 'wdfreq'})             
        ds_sub = ds_sub.assign_coords(west_east=ds_sub['lon'], 
                                      south_north=ds_sub['lat'], 
                                      sector=wk.create_sector_coords(12))
        ds_sub = wk.spatial._point._from_scalar(ds_sub)        
        A_combined, k_combined = wk.weibull_combined(ds_sub)
        ARef[i] = A_combined
        kRef[i] = k_combined              
        windrose[i, :] = ds.sel(wf=i, height=zRef).interp(roughness=0.0)['frequency'] / ds.sel(wf=i, height=zRef).interp(roughness=0.0)['frequency'].sum()
    # Log interpolate A at hub height
    AH = ARef * np.log(zH / z0) / np.log(zRef / z0) 
    return ARef, kRef, AH, windrose


def get_utm_proj(lat, lon):
    zone = int((lon + 180) // 6) + 1
    south = lat < 0
    return Proj(proj="utm", zone=zone, datum="WGS84", south=south)

   
def get_layouts(inputdata, wtdata):

    # We need some renaming
    wtdata["wind_farm"] = wtdata["wind_farm"].replace("Arkona-Becken Südost", "Arkona")
    wtdata["wind_farm"] = wtdata["wind_farm"].replace('EnBW Windpark Baltic 1', 'Baltic 1')
    wtdata["wind_farm"] = wtdata["wind_farm"].replace('EnBW Windpark Baltic 2', 'Baltic 2')
    wtdata["wind_farm"] = wtdata["wind_farm"].replace('Nysted', 'Rodsand 1')
    wtdata["wind_farm"] = wtdata["wind_farm"].replace('Rodsand II', 'Rodsand 2')
    wtdata["wind_farm"] = wtdata["wind_farm"].replace('Anholt', 'Anholt 1')
    wtdata["wind_farm"] = wtdata["wind_farm"].replace('Horns Rev I', 'Horns Rev 1')
    wtdata["wind_farm"] = wtdata["wind_farm"].replace('Horns Rev II', 'Horns Rev 2')
    wtdata["wind_farm"] = wtdata["wind_farm"].replace('Horns Rev III', 'Horns Rev 3')
    wtdata["wind_farm"] = wtdata["wind_farm"].replace('Dan Tysk', 'DanTysk')
    wtdata["wind_farm"] = wtdata["wind_farm"].replace('Meerwind Sued/Ost', 'Meerwind Sud/Ost')
    wtdata["wind_farm"] = wtdata["wind_farm"].replace('EnBW Hohe See', 'Hohe See')
    wtdata["wind_farm"] = wtdata["wind_farm"].replace('Bard Offshore 1', 'BARD')
    wtdata["wind_farm"] = wtdata["wind_farm"].replace('Gode Wind 1 and 2', 'Gode 1 and 2')
    wtdata["wind_farm"] = wtdata["wind_farm"].replace('Borkum Riffgrund 1', 'Borkum Riffgrund I')
    wtdata["wind_farm"] = wtdata["wind_farm"].replace('Borkum Riffgrund 2', 'Borkum Riffgrund II')
    wtdata["wind_farm"] = wtdata["wind_farm"].replace('Trianel Windpark Borkum 1', 'Trianel I and II')
    wtdata["wind_farm"] = wtdata["wind_farm"].replace('Trianel Windpark Borkum 2', 'Trianel I and II')
    wtdata["wind_farm"] = wtdata["wind_farm"].replace('Merkur Offshore', 'Merkur')
    wtdata["wind_farm"] = wtdata["wind_farm"].replace('OWF Prinses Amalia', 'Princess Amalia')
    wtdata["wind_farm"] = wtdata["wind_farm"].replace('OWF Luchterduinen', 'Luchterduinen' )
    wtdata["wind_farm"] = wtdata["wind_farm"].replace('Borssele Kavel I and II', 'Borssele I-II')
    wtdata["wind_farm"] = wtdata["wind_farm"].replace('Thornton Bank - phase I', 'Thortonbank')
    wtdata["wind_farm"] = wtdata["wind_farm"].replace('Thornton Bank - phase II and III', 'Thortonbank')
    wtdata["wind_farm"] = wtdata["wind_farm"].replace('Belwind phase 1', 'Belwind Phase 1')
    wtdata["wind_farm"] = wtdata["wind_farm"].replace('Belwind phase 2', 'Belwind Phase 1')
    wtdata["wind_farm"] = wtdata["wind_farm"].replace('Mermaid', 'Seamade Mermaid')
    wtdata["wind_farm"] = wtdata["wind_farm"].replace('SeaStar', 'Seamade Seastar')
    wtdata["wind_farm"] = wtdata["wind_farm"].replace('Gunfleet Sands', 'Gunfleet Sand')
    wtdata["wind_farm"] = wtdata["wind_farm"].replace('Lynn', 'Lynn and Inner Dowsing')
    wtdata["wind_farm"] = wtdata["wind_farm"].replace('Inner Dowsing', 'Lynn and Inner Dowsing')
    wtdata["wind_farm"] = wtdata["wind_farm"].replace('Hornsea Project 1', 'Hornsea 1')
    wtdata["wind_farm"] = wtdata["wind_farm"].replace('Hornsea Project 2', 'Hornsea 2')
    wtdata["wind_farm"] = wtdata["wind_farm"].replace('Beatrice Offshore Wind Farm', 'Beatrice extension')
    wtdata["wind_farm"] = wtdata["wind_farm"].replace('Walney Extension 3', 'Walney Extension')
    wtdata["wind_farm"] = wtdata["wind_farm"].replace('Walney Extension 4', 'Walney Extension')
    wtdata["wind_farm"] = wtdata["wind_farm"].replace('Gwynt y Mor', 'Gwynt y Môr')
    wtdata["wind_farm"] = wtdata["wind_farm"].replace('Robin Rigg West', 'Robin Rigg')

    # Not available in open turbine database
    # Here it is added to the open data based        
    skip = ['Fryslan']
    
    #print(wtdata)
    #df = wtdata.where(wind_farm='Anholt')
    nWF = len(inputdata['Name'])    
    wf_lats = []
    wf_lons = []
    wt_lats = []
    wt_lons = []    
    wt_xs = []
    wt_ys = []
    ps = []
    for i in range(nWF):
        if inputdata['Name'][i] not in skip:
            if (wtdata["wind_farm"] == inputdata['Name'][i]).any():          
                wf = wtdata.loc[wtdata["wind_farm"] ==  inputdata['Name'][i], :]                
            else:
                print('ERROR: wind farm= %s does not exist in open turbine database' % inputdata['Name'][i])
                sys.exit()           
            wt_lat = wf['latitude']
            wt_lon = wf['longitude']
        elif inputdata['Name'][i] == 'Fryslan':            
            layout = np.genfromtxt('data/fryslan_layout.dat')
            wt_lat = layout[:, 0]
            wt_lon = layout[:, 1]            
            if not len(wt_lat) == inputdata['Nt'][i]:
                print('ERROR: number of turbines for wind farm=%s does not match between input and open turbine database' % inputdata['Name'][i], len(wt_lat), inputdata['Nt'][i])
                sys.exit()              

        # Convert to UTM
        wf_lat = wt_lat.mean()
        wf_lon = wt_lon.mean()
        p = get_utm_proj(wf_lat, wf_lon)
        ps.append(p)
        wt_x, wt_y = p(wt_lon, wt_lat)

        wf_lats.append(wf_lat)
        wf_lons.append(wf_lon)
        wt_lats.append(wt_lat)
        wt_lons.append(wt_lon)         
        wt_xs.append(wt_x)
        wt_ys.append(wt_y) 
           
    return wf_lats, wf_lons, wt_xs, wt_ys, wt_lats, wt_lons, ps


def unique_points_tol(points, tol=1e-3):
    unique = []
    for p in points:
        if not any(np.allclose(p, q, atol=tol) for q in unique):
            unique.append(p)
            
    return np.array(unique)


def calculate_wf_area_and_edge(wt_x, wt_y, ConcaveHull_alpha=None, plot=False, outfile='wflayout_polygon.pdf'):
    """
    Estimate wind farm area, and edge turbines using ConcaveHull. ConcaveHull requires
    alphashape package to be installed
    """
    points = np.transpose(np.vstack((np.asarray(wt_x), np.asarray(wt_y))))       
    if not alphashape:
        print('ConcaveHull requires alphashape package: pip install alphashape')
        sys.exit()
    if not Polygon:
        print('ConcaveHull requires shapely package: pip install shapely')
        sys.exit()
    if ConcaveHull_alpha is None:
        ConcaveHull_alpha = 0.95 * alphashape.optimizealpha(points)
    hull = alphashape.alphashape(points, ConcaveHull_alpha)
    if hull.is_empty:
        print('ConcaveHull failed, try different ConcaveHull_alpha or set ConcaveHull_alpha=None')
        sys.exit()
    pnts = [x for x in hull.boundary.coords]
    p = Polygon(pnts)
    wf_area = p.area
    boundary_points = np.asarray([list(t) for t in pnts])
    if plot:
        fig = plt.figure()
        plt.plot(points[:, 0], points[:, 1], 'o')     
        plt.gca().add_patch(Polygonplt(pnts, fill=False, color='k'))
        plt.xlabel('x')
        plt.ylabel('y')
        plt.axis('equal')
        fig.savefig(outfile)        
    return wf_area, boundary_points, hull


def norm_vec(vec):
    norm = np.linalg.norm(vec)
    return vec / norm if norm > 0 else vec


def calculate_front_row_wts(boundary_wts, wd, hull):
    # Calculate normals for each edge
    hull = orient(hull, sign=1.0)
   
    coords = np.array(hull.exterior.coords[:-1])  # remove duplicate closing point

    normals = []

    for i in range(len(coords)):
        p1 = coords[i]
        p2 = coords[(i + 1) % len(coords)]

        edge = p2 - p1
        normal = np.array([edge[1], -edge[0]])
        normal = norm_vec(normal)
        normals.append(normal)

    normals = np.array(normals)

    # Vertex normals (mean of adjacent edge normals)
    normals_mean = np.zeros_like(normals)

    for i in range(len(coords)):
        n_sum = normals[i] + normals[i - 1]   # previous + current edge
        normals_mean[i] = norm_vec(n_sum)    
           
    # Angle with wd
    deg = (270.0 - wd) / 180.0 * np.pi
    wd_vec = np.array([np.cos(deg), np.sin(deg)]) 

    dot = normals_mean @ wd_vec   # vectorized dot products
    inlet  = dot < 0
    outlet = dot > 0

    # Remove inlet points that are in the shadow of the wind farm polygon using ray casting
    # length large enough to escape hull
    L = 1e6  
    shadow = []
    for p in coords:
        start = Point(p)
        end = Point(p + -wd_vec * L)
        ray = LineString([start, end])
        # count intersections with hull boundary
        intersection = ray.intersection(hull.boundary)
        # If more than one intersection → blocked
        if intersection.geom_type == "MultiPoint":
            shadow.append(True)
        else:
            shadow.append(False)
    shadow = np.array(shadow)
    inlet = inlet & ~shadow

    pnts = [x for x in hull.boundary.coords]  
    boundary_points = np.asarray([list(t) for t in pnts])
    return inlet, boundary_points[:-1], normals_mean


def get_wf_neighbor(wf_lats, wf_lons, wfnames, lat_dist=0.5, lon_dist=0.5):
    # Find wind farm neighbors within lat_dist and lon_dist
    nWF = len(wfnames)
    wf_neighbors = []
    wf_neighbors_flag = np.asarray([False] * nWF)
    for i in range(nWF):
        wf_neighbor = wfnames[(np.abs(wf_lats[i] - wf_lats) < lon_dist) & (np.abs(wf_lons[i] - wf_lons) < lon_dist)]
        if len(wf_neighbor) > 1:
            wf_neighbors_flag[i] = True
        wf_neighbors.append(wf_neighbor.drop([i]))   
    return wf_neighbors, wf_neighbors_flag


def remove_inlet_turbines_wf_neighbors(i, wfname, inlet, boundary_points, hulls, inputdata, wd, wf_neighbors, wf_neighbors_flag, ps, wf_norms, Ldist):
    wf_shadow = False
    wf_neighbor_coords = None
    if wf_neighbors_flag:
        wf_neighbor_coords = []  
        shadow = np.array([False] * len(boundary_points))
        for j in range(len(wf_neighbors)):
            # Get hull of neigbooring farm
            k = wf_neighbors.index[j]
            hull = orient(hulls[k], sign=1.0)             
            
            # Convert boundary coordinates to farm of interest
            coords = np.array(hull.exterior.coords)              
            wt_lon, wt_lat = ps[k](coords[:, 0] + wf_norms[k][0], coords[:, 1] + wf_norms[k][1], inverse=True)            
            wt_x, wt_y = ps[i](wt_lon, wt_lat)
            coords2 = np.zeros(coords.shape)
            coords2[:, 0] = wt_x - wf_norms[i][0]
            coords2[:, 1] = wt_y - wf_norms[i][1]
            wf_neighbor_coords.append(coords2)            
            hull = Polygon(coords2)            
            hull = orient(hull, sign=1.0)       
            
            # Remove inlet points that are in the shadow of a neighboring wind farm polygon using ray casting
            # length large enough to escape hull
            L = 1e6  
            # and are within Ldist
            # Also remove if the inlet point is inside a neighboring wind farm polygon
            # Angle with wd
            deg = (270.0 - wd) / 180.0 * np.pi
            wd_vec = np.array([np.cos(deg), np.sin(deg)])
            #print(boundary_points)           
            for p in range(len(boundary_points)):
                start = Point(boundary_points[p])
                end = Point(boundary_points[p] + -wd_vec * L)
                ray = LineString([start, end])
                # count intersections with hull boundary
                intersection = ray.intersection(hull.boundary)
                 
                length = 1e6          
                if intersection.geom_type == "MultiPoint":
                    dists = [start.distance(pt) for pt in intersection.geoms]
                    length = min(dists)   # first hit along ray
             
                # If more than one intersection → blocked or the point is inside the neighboring wind farm
                if (intersection.geom_type == "MultiPoint" and length < Ldist) or start.within(hull):
                    shadow[p] = True
        if True in shadow:           
            wf_shadow = True                 
        inlet = inlet & ~shadow 
    return inlet, wf_shadow, wf_neighbor_coords


def calculate_layout_metrics(wt_xs, wt_ys, alphas, inputdata, wf_neighbors, wf_neighbors_flag, wt_lats, wt_lons, ps, windrose, plot=False, plot_farm_example=''):
    nWF = len(inputdata['Name'])    
    #nWF = 3
    Areas = []
    Nwt_rows = []
    Nwt_row_wf_neighbors = []
    inlet_all = []
    boundary_points_all = []
    normals_mean_all = []
    hulls = []
    wf_norms = []
    wd_windrose = np.array([0, 30, 60, 90, 120, 150, 180, 210, 240, 270, 300, 330])     
    for i in range(nWF):
        wt_x = np.asarray(wt_xs[i])
        wt_y = np.asarray(wt_ys[i])

        wf_norm = [wt_x.min(), wt_y.min()]
        wf_norms.append(wf_norm)
        wt_x = wt_x - wf_norm[0]
        wt_y = wt_y - wf_norm[1]      
       
        Dref = inputdata['D'][i]
        
        # Wind farm area and edge turbines
        area, boundary_points, hull = calculate_wf_area_and_edge(wt_x, wt_y, ConcaveHull_alpha=alphas[i], plot=False, outfile='wflayout_polygon.pdf')
        Areas.append(area)
                
        # inlet edge turbine and turbine normals for plotting
        inlet_wd = []
        for l in range(12):        
            inlet, boundary_points, normals_mean = calculate_front_row_wts(boundary_points, wd_windrose[l], hull)            
            inlet_wd.append(inlet)
        inlet_all.append(inlet_wd)
        boundary_points_all.append(boundary_points)
        normals_mean_all.append(normals_mean)    
        hulls.append(hull)
        
    for i in range(nWF):
        # Remove inlet turbines due to neighbooring wind farms with a length scale x times the turbine spacing
        S = np.sqrt(inputdata['A'][i] * 1e6) / (inputdata['D'][i] * (np.sqrt(inputdata['Nt'][i])-1))
        L = 10 * S * inputdata['D'][i]
        # print('L', L)
        wf_shadows = []
        inlet_with_wf_neighborss = []
        for l in range(12):
            inlet_with_wf_neighbors, wf_shadow, wf_neighbor_coords = remove_inlet_turbines_wf_neighbors(i, inputdata['Name'][i], inlet_all[i][l], boundary_points_all[i], hulls, inputdata, wd_windrose[l], wf_neighbors[i], wf_neighbors_flag[i], ps, wf_norms, L)
            wf_shadows.append(wf_shadow)
            inlet_with_wf_neighborss.append(inlet_with_wf_neighbors)
   
        Nwt_row = 0
        for l in range(12): 
            Nwt_row = Nwt_row + len(boundary_points_all[i][:, 0][inlet_all[i][l]]) * windrose[i, l]
        Nwt_rows.append(Nwt_row)

        Nwt_row_wf_neighbor = 0        
        for l in range(12): 
            Nwt_row_wf_neighbor = Nwt_row_wf_neighbor + len(boundary_points_all[i][:, 0][inlet_with_wf_neighborss[l]]) * windrose[i, l]         
        Nwt_row_wf_neighbors.append(Nwt_row_wf_neighbor)         
        Dref = inputdata['D'][i]        
        if plot:
            fig, ax = plt.subplots(12, 1, figsize=(10.0, 5*12.0))        
            for l in range(12):        
                ax[l].scatter((wt_xs[i] - wf_norms[i][0]) / Dref, (wt_ys[i] - wf_norms[i][1]) / Dref)
                ax[l].set_aspect('equal')
                deg = (270.0 - wd_windrose[l]) / 180.0 * np.pi

                # Edge turbines    
                ax[l].add_patch(Polygonplt(boundary_points_all[i] / Dref, fill=False, color='k'))

                # Boundary turbines
                ax[l].scatter(boundary_points_all[i][:, 0][inlet_all[i][l]] / Dref, boundary_points_all[i][:, 1][inlet_all[i][l]] / Dref, marker='s', color='r', facecolors='none', s=50)#
                            
                # wd
                farmsize = area ** 0.5
                l1 = 0.1 * farmsize / Dref
                ax[l].arrow(-2*l1, -2*l1, np.cos(deg) * l1, np.sin(deg) * l1,  head_width=0.5 * l1, color='m')

                # Normals
                for j in range(len(normals_mean_all[i])):
                    ax[l].arrow(boundary_points_all[i][j, 0] / Dref, boundary_points_all[i][j, 1] / Dref, normals_mean_all[i][j, 0] * l1, normals_mean_all[i][j, 1] * l1, head_width=0.5 * l1, color='k')    

                ax[l].set_title(r'Nr front wts: %s, $N_{tot}=%s$, $\theta=%3.1f^\circ$, $f=%3.2f$, $D=%s$' % (len(boundary_points_all[i][:, 0][inlet_all[i][l]]), inputdata['Nt'][i], wd_windrose[l], windrose[i][l], inputdata['D'][i]))    
                ax[l].set_xlabel(r'$x/D$')
                ax[l].set_ylabel(r'$y/D$')
            filename = 'plots_frontrowwts/layout_%s.png' % (inputdata['Name'][i].replace('/', '').replace(' ', '_'))
            fig.savefig(filename)
            plt.close()
                                     
            # WF neighboor WTs
            if wf_neighbors_flag[i]:                                                     
                fig, ax = plt.subplots(12, 1, figsize=(10.0, 5*12.0))               
                for l in range(12):
                    ax[l].scatter((wt_xs[i] - wf_norms[i][0]) / Dref, (wt_ys[i] - wf_norms[i][1]) / Dref)
                    ax[l].set_aspect('equal')
                    deg = (270.0 - wd_windrose[l]) / 180.0 * np.pi

                    hull_buffered = hulls[i].buffer(L)             
                    ax[l].add_patch(Polygonplt(hull_buffered.exterior.coords / Dref, fill=False, color='k', alpha=0.5))   

                    # Edge turbines    
                    ax[l].add_patch(Polygonplt(boundary_points_all[i] / Dref, fill=False, color='k'))       
 
                    for j in range(len(wf_neighbor_coords)):
                        # Edge turbines    
                        ax[l].add_patch(Polygonplt(wf_neighbor_coords[j] / Dref, fill=True, color='gray'))

                    ax[l].scatter(boundary_points_all[i][:, 0][inlet_with_wf_neighborss[l]] / Dref, boundary_points_all[i][:, 1][inlet_with_wf_neighborss[l]] / Dref, marker='s', color='r', facecolors='none', s=50)               
                          
                    # wd
                    farmsize = area ** 0.5
                    l1 = 0.5 * farmsize / Dref
                    ax[l].arrow(-2*l1, -2*l1, np.cos(deg) * l1, np.sin(deg) * l1,  head_width=0.5 * l1, color='m')

                    # Normals
                    for j in range(len(normals_mean_all[i])):
                        ax[l].arrow(boundary_points_all[i][j, 0] / Dref, boundary_points_all[i][j, 1] / Dref, normals_mean_all[i][j, 0] * l1, normals_mean_all[i][j, 1] * l1, head_width=0.5 * l1, color='k')    

                    ax[l].set_title(r'Nr front wts: %g, Nr front wts w. wf neighbors: %s, $N_{tot}=%s$, $\theta=%3.1f^\circ,  f=%3.2f, L/D=%g, D=%s$' % (len(boundary_points_all[i][:, 0][inlet_with_wf_neighborss[l]]), len(boundary_points_all[i][:, 0][inlet_with_wf_neighborss[l]]), inputdata['Nt'][i], wd_windrose[l], windrose[i][l], L/Dref, inputdata['D'][i]))    
                    ax[l].set_xlabel(r'$x/D$')
                    ax[l].set_ylabel(r'$y/D$')
                filename = 'plots_frontrowwts/layout_%s_wf_neighbor.png' % inputdata['Name'][i].replace('/', '').replace(' ', '_')
                fig.savefig(filename)
                plt.close()
                    

        if inputdata["Name"][i] == plot_farm_example:
            l = 8
            # Make an example plot                   
            fig = plt.figure(figsize=(10.0, 5))
            ax1 = fig.add_subplot(2,2,3)
            ax3 = fig.add_subplot(1,2,2)
            ax = [ax1, ax3]

            ax[0].scatter((wt_xs[i] - wf_norms[i][0]) / Dref, (wt_ys[i] - wf_norms[i][1]) / Dref)
            ax[0].set_aspect('equal')
            deg = (270.0 - wd_windrose[l]) / 180.0 * np.pi

            # Edge turbines    
            ax[0].add_patch(Polygonplt(boundary_points_all[i] / Dref, fill=False, color='k', linestyle='--'))

            # Boundary turbines
            ax[0].scatter(boundary_points_all[i][:, 0][inlet_all[i][l]] / Dref, boundary_points_all[i][:, 1][inlet_all[i][l]] / Dref, marker='s', color='r', facecolors='none', s=50)
                            
            # wd
            farmsize = area ** 0.5
            l1 = 0.1 * farmsize / Dref
            ax[0].arrow(-2*l1, -2*l1, np.cos(deg) * l1, np.sin(deg) * l1,  head_width=0.5 * l1, color='m')

            # Normals
            for j in range(len(normals_mean_all[i])):
                ax[0].arrow(boundary_points_all[i][j, 0] / Dref, boundary_points_all[i][j, 1] / Dref, normals_mean_all[i][j, 0] * l1, normals_mean_all[i][j, 1] * l1, head_width=0.5 * l1, color='k')    

            ax[0].set_title(r'$M_{\rm turbines} = %g$' % (len(boundary_points_all[i][:, 0][inlet_all[i][l]])))    
            ax[0].set_xlabel(r'$x/D$')
            ax[0].set_ylabel(r'$y/D$')
            
            if wf_neighbors_flag[i]:
                ax[1].scatter((wt_xs[i] - wf_norms[i][0]) / Dref, (wt_ys[i] - wf_norms[i][1]) / Dref, label='Amrumbank West layout')
                ax[1].set_aspect('equal')
                deg = (270.0 - wd_windrose[l]) / 180.0 * np.pi

                hull_buffered = hulls[i].buffer(L)             
                ax[1].add_patch(Polygonplt(hull_buffered.exterior.coords / Dref, fill=False, color='g', alpha=0.5, label=r'Inter wind farm distance $L=10S$'))   

                # Edge turbines    
                ax[1].add_patch(Polygonplt(boundary_points_all[i] / Dref, fill=False, color='k', label='Concave polygon', linestyle='--'))       
 
                for j in range(len(wf_neighbor_coords)):
                    # Edge turbines
                    if j == 0:
                        ax[1].add_patch(Polygonplt(wf_neighbor_coords[j] / Dref, fill=True, color='gray', label='Neighbor wind farms'))
                    else:
                        ax[1].add_patch(Polygonplt(wf_neighbor_coords[j] / Dref, fill=True, color='gray'))                        

                ax[1].scatter(boundary_points_all[i][:, 0][inlet_with_wf_neighborss[l]] / Dref, boundary_points_all[i][:, 1][inlet_with_wf_neighborss[l]] / Dref, marker='s', color='r', facecolors='none', s=50, label='Inlet turbines')              
                          
                # wd
                farmsize = area ** 0.5
                l1 = 0.2 * farmsize / Dref
                ax[1].arrow(-2*l1, -2*l1, np.cos(deg) * l1, np.sin(deg) * l1,  head_width=0.5 * l1, color='m')

                # Normals
                for j in range(len(normals_mean_all[i])):
                    ax[1].arrow(boundary_points_all[i][j, 0] / Dref, boundary_points_all[i][j, 1] / Dref, normals_mean_all[i][j, 0] * l1, normals_mean_all[i][j, 1] * l1, head_width=0.5 * l1, color='k') 
                ax[1].set_title(r'$M_{\rm turbines} = %g$' % (len(boundary_points_all[i][:, 0][inlet_with_wf_neighborss[l]])))    
                ax[1].set_xlabel(r'$x/D$')
                ax[1].set_ylabel(r'$y/D$') 
                ax[1].set_ylim(-150, 100)                            
            # Get handles and labels
            handles, labels = plt.gca().get_legend_handles_labels()
            order = [0, 4, 2, 1, 3]
            ax[0].legend([handles[ii] for ii in order], [labels[ii] for ii in order], bbox_to_anchor=(1.0, 2.5))          
            ax[0].text(-16, 32, '(a)', fontsize=14)
            ax[1].text(-90, 80, '(b)', fontsize=14)            
            filename = 'fig1_Nwt_edge_example_%s.pdf' % inputdata['Name'][i].replace('/', '').replace(' ', '_')
            fig.savefig(filename)
            plt.close()                
           
    return Areas, Nwt_rows, Nwt_row_wf_neighbors

        
def eps_inf_simple(s, c1, c2):
    # A very simple model for epsilon = U_infinity / U_free, for constant CT=0.75, delta and gamma
    return (1.0 + c1) / (1.0 + np.sqrt(c2 / s ** 2 + c1 ** 2))

            
def phi_simple(s, Ur_Aw, Nfree_Ntot, c1, c2, c3):
    # A very simple model for phi
    f_epsinf =  eps_inf_simple(s, c1, c2)
    return c3 * Ur_Aw / (Nfree_Ntot * (1 - f_epsinf) + f_epsinf)


def fit_line(x, y):
    # Fit linear model y = m x + b
    m, b = np.polyfit(x, y, 1)
    # Predicted values
    y_pred = m * x + b
    # R^2 calculation
    ss_res = np.sum((y - y_pred)**2)
    ss_tot = np.sum((y - np.mean(y))**2)
    r2 = 1 - ss_res / ss_tot
    return m, b, r2


if __name__ == '__main__': 

    # Main input parameters
    z0 = 10 ** (-4)  # Roughness length [m]
    zRef = 100.0  # Reference height at which GWA data is taken [m]
    kw = 2.4  # Weibull shape parameter [-]
    rho = 1.225  # Air density [kg/m^3]
    lat = 55  # Latitude [degree]
    CP = 0.46  # Turbine power coefficient [-]
    CT = 0.75  # Turbine thrust coefficient [-]
    Uin = 0.0  # Cut-in wind speed [m/s]
    Uout = 1e6  # Cut-out wind speed [m/s]    
    # Uin = 3.0
    # Uout = 25.0

    # Run flags        
    # Re create output
    recreate_output = False
    # recreate_output = True
    # Rerun Global windatlas api    
    rerun_gwa = False
    # plot all layouts per wd sector including result of the inflow edge turbines
    # Note that this generates many plots!
    plot_layouts = False
    # plot_layouts = True
    if recreate_output:
        # Wind farms input data send by email from Carlos Simao Ferreira
        infile = 'data/2026-02-11_input_vars_send_by_CarlosSimaoFerreira.csv'
        inputdata = pd.read_csv(infile)    

        # Get list of lat lon in decimal format
        latlon = get_latlon_decimal(inputdata['Coordinates'])

        # Get Global Wind Atlas wind resource for each wind farm
        # get_gwa_windresource(latlon, inputdata['Name'], rerun=True)    
        ds = get_gwa_windresource(latlon, inputdata['Name'], rerun=rerun_gwa)    
        
        # Calculate A, k and mean wd
        ARef, kRef, AH, windrose = calc_windresource(ds, z0, zRef, np.asarray(inputdata['Ht']))

        # Load wind farm layouts from open turbine database https://zenodo.org/records/17311571
        infile = 'data/TurbineDatabase/20251218_eww_opendatabase.csv'
        wtdata = pd.read_csv(infile)
    
        # Remove two turbines from Robin Rigg which have been decommissioned
        # https://www.sciencedirect.com/science/article/pii/S0267726123000489
        # 5026                   54.7715     -3.7006
        # 5028                   54.7692     -3.7093    
        wtdata = wtdata.drop([5026, 5028])     
        # print(wtdata)       
        wf_lats, wf_lons, wt_xs, wt_ys, wt_lats, wt_lons, ps = get_layouts(inputdata, wtdata)    
    
        # Find wind farm neighbors
        wf_neighbors, wf_neighbors_flag = get_wf_neighbor(wf_lats, wf_lons, inputdata['Name'], lat_dist=0.5, lon_dist=0.5)
        # print('Isolated wind farms:', inputdata["Name"][~wf_neighbors_flag])
        
        # Calculate wind farm area, edges and "inlet" turbines
        # A Concave hull method is used to determine the farm edge, which is not unique and depends on alpha.
        # A larger alpha means more concave. A too large alpha can result in strange results.
        alphas = [0.0001] * len(inputdata['Nt'])
        alphas[inputdata.loc[inputdata["Name"] == 'Amrumbank West', :].index[0]] = 0.0005        
        alphas[inputdata.loc[inputdata["Name"] == 'Kriegers Flak', :].index[0]] = 0.00025
        alphas[inputdata.loc[inputdata["Name"] == 'Beatrice extension', :].index[0]] = 0.0005    
        alphas[inputdata.loc[inputdata["Name"] == 'Belwind Phase 1', :].index[0]] = 0.0005        
        alphas[inputdata.loc[inputdata["Name"] == 'Borssele I-II', :].index[0]] = 0.0005       
        alphas[inputdata.loc[inputdata["Name"] == 'DanTysk', :].index[0]] = 0.0005    
        alphas[inputdata.loc[inputdata["Name"] == 'Dudgeon', :].index[0]] = 0.0005
        alphas[inputdata.loc[inputdata["Name"] == 'Galloper', :].index[0]] = 0.00015  
        alphas[inputdata.loc[inputdata["Name"] == 'Gode 1 and 2', :].index[0]] = 0.0005
        alphas[inputdata.loc[inputdata["Name"] == 'Greater Gabbard', :].index[0]] = 0.00025        
        alphas[inputdata.loc[inputdata["Name"] == 'Gunfleet Sand', :].index[0]] = 0.0005     
        alphas[inputdata.loc[inputdata["Name"] == 'Gwynt y Môr', :].index[0]] = 0.001
        alphas[inputdata.loc[inputdata["Name"] == 'Hornsea 1', :].index[0]] = 0.00025
        alphas[inputdata.loc[inputdata["Name"] == 'Hornsea 2', :].index[0]] = 0.0002
        alphas[inputdata.loc[inputdata["Name"] == 'Horns Rev 3', :].index[0]] = 0.0005    
        alphas[inputdata.loc[inputdata["Name"] == 'Kaskasi', :].index[0]] = 0.0005
        alphas[inputdata.loc[inputdata["Name"] == 'Lillgrund', :].index[0]] = 0.0005
        alphas[inputdata.loc[inputdata["Name"] == 'Lincs', :].index[0]] = 0.001
        alphas[inputdata.loc[inputdata["Name"] == 'Meerwind Sud/Ost', :].index[0]] = 0.0005    
        alphas[inputdata.loc[inputdata["Name"] == 'Merkur', :].index[0]] = 0.0005 
        alphas[inputdata.loc[inputdata["Name"] == 'Moray East', :].index[0]] = 0.0005
        alphas[inputdata.loc[inputdata["Name"] == 'Nobelwind', :].index[0]] = 0.0005
        alphas[inputdata.loc[inputdata["Name"] == 'Nordsee One', :].index[0]] = 0.0005                
        alphas[inputdata.loc[inputdata["Name"] == 'Norther', :].index[0]] = 0.0005
        alphas[inputdata.loc[inputdata["Name"] == 'Northwester 2', :].index[0]] = 0.0005    
        alphas[inputdata.loc[inputdata["Name"] == 'Princess Amalia', :].index[0]] = 0.001  
        alphas[inputdata.loc[inputdata["Name"] == 'Race Bank', :].index[0]] = 0.001 
        alphas[inputdata.loc[inputdata["Name"] == 'Rampion', :].index[0]] = 0.00025 
        alphas[inputdata.loc[inputdata["Name"] == 'Veja Mate', :].index[0]] = 0.00025
        alphas[inputdata.loc[inputdata["Name"] == 'Walney 2', :].index[0]] = 0.00025
        alphas[inputdata.loc[inputdata["Name"] == 'Walney Extension', :].index[0]] = 0.00025
        Areas, Nwt_rows, Nwt_row_wf_neighbors = calculate_layout_metrics(wt_xs, wt_ys, alphas, inputdata, wf_neighbors, wf_neighbors_flag, wt_lats, wt_lons, ps, windrose, plot=plot_layouts,         plot_farm_example='Amrumbank West')        
    
        # Calculate capacity factor and wind speed reduction from model of Jens and Gunner                                                                
        Cf_Free, Cf_Inf, Cf_Model, U_Free, U_Inf, U_Model, G_Model, Phi_Model, Phi_Model_redo, phi_Model_a5p3, Cf_Model_redo, Cf_Model_a5p3 = calculate_minimalistic_model(inputdata, lat, CP, CT, z0, kw, Uin, Uout, rho, zRef, Nwt_row_wf_neighbors) 
    
        # Rated wind speed
        Ur = (inputdata['PG'] * 1e6 / (0.125 * rho * inputdata['D'] ** 2 * np.pi * CP)) ** (1.0 / 3.0)
        # Save to output to csv
        df = inputdata.copy()
        df["k_ref"] = kRef
        df["A_ref"] = ARef
        df["A_h"] = AH
        df["lat"] = latlon[:, 0]
        df["lon"] = latlon[:, 1]
        df["lat2"] = wf_lats
        df["lon2"] = wf_lons
        df["Area"] = np.asarray(Areas) * 1e-6
        df["Nwt_rows"] = Nwt_rows
        df["Nwt_row_wf_neighbors"] = Nwt_row_wf_neighbors        
        df["Cf_Free"] = Cf_Free
        df["Cf_Inf"] = Cf_Inf
        df["Cf_Model"] = Cf_Model
        df["U_Free"] = U_Free
        df["U_Inf"] = U_Inf    
        df["U_Model"] = U_Model
        df["G_Model"] = G_Model        
        df["U_rated"] = Ur
        df["Phi_Model"] = Phi_Model             
        df["Phi_Model_redo"] = Phi_Model_redo  
        df["Phi_Model_a5p3"] = phi_Model_a5p3                      
        df["Cf_Model_redo"] = Cf_Model_redo                 
        df["Cf_Model_a5p3"] = Cf_Model_a5p3         
        df.to_csv("output.csv", index=False)
 
    ##############
    # Make plots
    ##############
    
    # Load wind farm input data from Carlos Simao Ferreira and added output from this script
    infile = 'output.csv'
    outputdata = pd.read_csv(infile)
    
    # Output data from reference 15 + column with phi from table s3
    infile = 'data/ref15/Analysis_result.csv'         
    outputdata_ref15 = pd.read_csv(infile)     
    # Drop last 4 rows
    outputdata_ref15 = outputdata_ref15.drop([72, 73, 74, 75])  
    
    # Digitized phi from table s1
    infile = 'data/digitized_phi_table_s1.csv'
    outputdata_ref15_phi = pd.read_csv(infile)   
    outputdata_ref15["phi_table_s1"] = outputdata_ref15_phi["phi_table_s1"]    
        
    # Capacity factor
    f_loss_factor = 0.9
    fig, ax = plt.subplots(1, 1, figsize=(10.0, 5.0))
    Ur_A = np.arange(0.5, 2.0 + 0.01, 0.01)    
    Cf_an = get_cf_an(Ur_A, kw)
    x =  Ur_A /  gamma(1 + 1 / kw)         
    ax.plot(x, Cf_an, '-b', label=r'Normalized analytic gross AEP', lw=2)
    ax.plot(x, f_loss_factor * Cf_an, '--b', label=r'$f_{\rm loss}=0.9$', lw=2)       
    # Drop wind farm for which we do not have individual capacity factors   
    outputdata_drop = outputdata.drop([4, 5, 35, 36])      
    outputdata_ref15_drop = outputdata_ref15.drop([4, 5, 35, 36])        

    ax.scatter(outputdata_ref15_drop['phi_table_s1'], outputdata_ref15_drop['Capacity Factor real'] * 0.01, color='r', label=r'Ferreira et al. (2026), Table S1', zorder=2)      
    phi0 = outputdata_drop['U_rated'] / outputdata_drop['U_Free']
    phi1 = outputdata_drop['U_rated'] / outputdata_drop['U_Inf']
    phic = 0.5 * (phi0 + phi1)
    xerr = np.zeros((2, len(phi0)))
    xerr[0, :] = phic - phi0
    xerr[1, :] = phi1 - phic
    ax.errorbar(phic, outputdata_ref15_drop['Capacity Factor real'] * 0.01, xerr=xerr, color='gray', fmt='none', elinewidth=3, capsize=3, capthick=True, zorder=-1, label=r'Model range: $N_{\rm free}=0, N_{\rm free}=N_{\rm tot}$', alpha=0.5)    
    ax.scatter(outputdata_drop['Phi_Model'], outputdata_ref15_drop['Capacity Factor real'] * 0.01, color='g', label=r'PyWake, $N_{\rm free}$ from Ferreira et al. (2026)', marker='s', facecolor='None')    
    ax.scatter(outputdata_drop['Phi_Model_redo'], outputdata_ref15_drop['Capacity Factor real'] * 0.01, color='m', label=r'PyWake, $N_{\rm free}$ from automated script', marker='^')
    ax.scatter(outputdata_drop['Phi_Model_a5p3'], outputdata_ref15_drop['Capacity Factor real'] * 0.01, color='k', label=r'PyWake, $N_{\rm free}=min(5.3\sqrt{N_{\rm tot}}, N_{\rm tot})$', marker='v') 
    # ax.scatter(outputdata_drop['Phi_Model_a5p3'], outputdata_ref15_drop['Capacity Factor real'] * 0.01, color='k', label=r'PyWake, $N_{\rm free}=min(a\sqrt{N_{\rm tot}}, N_{\rm tot})$,' + '\n' +  r'$a=5.3$, from Sørensen et al. (2024)', marker='v')     
    ax.grid(True)
    ax.set_xlim(1.0, 2.0)
    ax.set_ylim(0.3, 0.6)
    ax.set_ylabel(r'Measured capacity factor')
    ax.set_xlabel(r'$\phi=U_r / (U_0 \varepsilon)$')
    handles, labels = plt.gca().get_legend_handles_labels()
    order = [0, 1, 6, 2, 3, 4, 5]
    ax.legend([handles[ii] for ii in order], [labels[ii] for ii in order], bbox_to_anchor=(1.1, 1.1))       
    filename = 'fig5_CapacityFactor.pdf'
    fig.savefig(filename)    

    # Compare Nwt_rows for all wind farms
    fig, ax = plt.subplots(2, 1, sharey=True, figsize=(10.0, 10.0))  
    xbin = np.linspace(0, 36 - 1, 36)   
    Nfree_Ntot_a5p3 = np.minimum(5.3 / np.sqrt(outputdata['Nt']), 1.0)
    ax[0].bar(xbin, np.asarray(2.5 * outputdata['Nturb_frontal_area'] / outputdata['Nt'])[0:36], 0.2, color='r', edgecolor='r', label='Manual method from Ferreira et al. (2026)')        
    ax[0].bar(xbin+0.2, np.asarray(2.5 * outputdata['Nwt_rows'] / outputdata['Nt'])[0:36], 0.2, color='g', edgecolor='g', label='Automated method without upstream wind farms', alpha=0.3)        
    ax[0].bar(xbin+0.2, np.asarray(2.5 * outputdata['Nwt_row_wf_neighbors'] / outputdata['Nt'])[0:36], 0.2, color='g', edgecolor='g', label='Automated method with upstream wind farms') 
            
    ax[0].bar(xbin+0.4, np.asarray(Nfree_Ntot_a5p3)[0:36], 0.2, color='k', edgecolor='k', alpha=0.5, label=r'$N_{\rm free}=min(5.3\sqrt{N_{\rm tot}}, N_{\rm tot})$')
    # ax[0].bar(xbin+0.4, np.asarray(Nfree_Ntot_a5p3)[0:36], 0.2, color='k', edgecolor='k', alpha=0.5, label=r'$N_{\rm free}=min(a\sqrt{N_{\rm tot}}, N_{\rm tot})$,' + '\n' +  r'$a=5.3$, from Sørensen et al. (2024)')         
    ax[0].set_xticks(xbin, np.asarray(outputdata['Name'])[0:36], rotation=90)    
    xbin = np.linspace(0, len(outputdata) - 1 - 36, len(outputdata) -36) 
    ax[1].bar(xbin, np.asarray(2.5 * outputdata['Nturb_frontal_area'] / outputdata['Nt'])[36:], 0.2, color='r', edgecolor='r', label='Manual method from Ferreira et al. (2026)')        
    ax[1].bar(xbin+0.2, np.asarray(2.5 * outputdata['Nwt_rows'] / outputdata['Nt'])[36:], 0.2, color='g', edgecolor='g', label='Automated method without upstream wind farms', alpha=0.3)        
    ax[1].bar(xbin+0.2, np.asarray(2.5 * outputdata['Nwt_row_wf_neighbors'] / outputdata['Nt'])[36:], 0.2, color='g', edgecolor='g', label='Automated method with upstream wind farms')  
    ax[1].bar(xbin+0.4, np.asarray(Nfree_Ntot_a5p3)[36:], 0.2, color='k', edgecolor='k', alpha=0.5, label=r'$N_{\rm free}=min(5.3\sqrt{N_{\rm tot}}, N_{\rm tot})$')
    # ax[1].bar(xbin+0.4, np.asarray(Nfree_Ntot_a5p3)[36:], 0.2, color='k', edgecolor='k', alpha=0.5, label=r'$N_{\rm free}=min(a\sqrt{N_{\rm tot}}, N_{\rm tot})$,' + '\n' +  r'$a=5.3$, from Sørensen et al. (2024)')                                                  
    ax[1].set_xticks(xbin, np.asarray(outputdata['Name'])[36:], rotation=90)  
    for i in range(2):
        ax[i].set_ylabel(r'$N_{\rm free}/N_{\rm tot}$')    
        ax[i].grid(True)        
        ax[i].set_xlim(-1, 36)        
    ax[0].set_ylim(0.0, 1.0)
    fig.tight_layout(rect=[0.0, 0.0, 1.0, 0.85])    
    ax[0].legend(bbox_to_anchor=(0.75, 1.75))
    filename = 'fig2_Nwtrows_all.pdf'
    fig.savefig(filename)
       
    # A very simple model for U_inf
    fig, ax = plt.subplots(1, 1, figsize=(8.0, 4.0))      
    S = np.sqrt(outputdata['A'] * 1e6) / (outputdata['D'] * (np.sqrt(outputdata['Nt'])-1))   
    ax.scatter(S, outputdata['U_Inf'] / outputdata['U_Free'], color='k', label=r'Model results') 
    c1 = 0.218  # gamma / delta
    c2 = 16.5  # pi CT / 8 / (gamma/kappa) ** 2
    x = np.arange(0.01, 14 + 0.01, 0.01)
    eps_inf_simple1 = eps_inf_simple(x, c1, c2)
    ax.plot(x, eps_inf_simple1, '-r', label=r'$f(s)=\frac{1 + c_1}{1 + \sqrt{c_2 / s^2 + c_1^2}}, c_1=0.218, c_2=16.5$')             
    ax.grid(True)
    ax.set_xlim(0.0, 14.0)
    ax.set_ylim(0.0, 1.0)
    ax.set_ylabel(r'$\varepsilon_\infty$')
    ax.set_xlabel(r'$s$')
    ax.legend(loc=0)
    fig.tight_layout()    
    filename = 'fig3_SuperSimpleModel_epsinf.pdf'
    fig.savefig(filename)    

    # A very simple model for phi
    fig, ax = plt.subplots(1, 1, figsize=(8.0, 4.0))      
    phi = outputdata['Phi_Model']
    c3 = 1.128  # 1 / gamma(1+1/kw)    
    phi_simple1 = phi_simple(S, outputdata['U_rated'] / outputdata['lambda'], 2.5 * outputdata['Nturb_frontal_area'] / outputdata['Nt'], c1, c2, c3)
    # error = (phi_simple1 - phi) / phi
    # print('error', error.max(), error.min())  
    ax.scatter(phi, phi_simple1, color='k', label=r'Model results')
    ax.plot([1.0, 2.0], [1.0, 2.0], '-r', label=r'$y=x$')    
    ax.grid(True)
    ax.set_xlim(1.0, 2.0)
    ax.set_ylim(1.0, 2.0)
    ax.set_ylabel(r'$\phi$')
    ax.set_xlabel(r'$c_0 \frac{U_r}{A_w}/(\frac{N_{free}}{N_{tot}} (1-\varepsilon_{\infty}(s)) + \varepsilon_{\infty}(s)), c_0=1.128$')
    ax.legend(loc=0)
    fig.tight_layout()
    filename = 'fig4_SuperSimpleModel_phi.pdf'
    fig.savefig(filename)          

    # Measured vs modelled capacity factor
    # Fit linear model y = m x + b
    m1, b1, r21 = fit_line(outputdata_ref15_drop['Capacity Factor model'] * 0.01, outputdata_ref15_drop['Capacity Factor real'] * 0.01)
    m2, b2, r22 = fit_line(outputdata_drop['Cf_Model'], outputdata_ref15_drop['Capacity Factor real'] * 0.01)
    m3, b3, r23 = fit_line(outputdata_drop['Cf_Model_redo'], outputdata_ref15_drop['Capacity Factor real'] * 0.01)    
    m4, b4, r24 = fit_line(outputdata_drop['Cf_Model_a5p3'], outputdata_ref15_drop['Capacity Factor real'] * 0.01)      

    fig, ax = plt.subplots(2, 2, sharex=True, sharey=True, figsize=(8.5, 8.5))   
    ax[0, 0].scatter(outputdata_ref15_drop['Capacity Factor model'] * 0.01, outputdata_ref15_drop['Capacity Factor real'] * 0.01, color='r', label='68 wind farms')
    ax[0, 1].scatter(outputdata_drop['Cf_Model'], outputdata_ref15_drop['Capacity Factor real'] * 0.01, color='g', marker='s', facecolor='None', label='68 wind farms')
    ax[1, 0].scatter(outputdata_drop['Cf_Model_redo'], outputdata_ref15_drop['Capacity Factor real'] * 0.01, color='m', marker='^', label='68 wind farms')    
    ax[1, 1].scatter(outputdata_drop['Cf_Model_a5p3'], outputdata_ref15_drop['Capacity Factor real'] * 0.01, color='k', marker='v', label='68 wind farms')    
    ax[0, 0].plot([0.2, 0.6], np.array([0.2, 0.6]) * m1 + b1, '--r', label=r'$y=%2.2f x + %2.2f, r^2=%2.2f$' % (m1, b1, r21))
    ax[0, 1].plot([0.2, 0.6], np.array([0.2, 0.6]) * m2 + b2, '--g', label=r'$y=%2.2f x + %2.2f, r^2=%2.2f$' % (m2, b2, r22))       
    ax[1, 0].plot([0.2, 0.6], np.array([0.2, 0.6]) * m3 + b3, '--m', label=r'$y=%2.2f x + %2.2f, r^2=%2.2f$' % (m3, b3, r23))
    ax[1, 1].plot([0.2, 0.6], np.array([0.2, 0.6]) * m4 + b4, '--k', label=r'$y=%2.2f x + %2.2f, r^2=%2.2f$' % (m4, b4, r24))        
    ax[0, 0].set_xlim(0.2, 0.6)
    ax[0, 0].set_ylim(0.2, 0.6)
    labels = [[r'Ferreira et al. (2026), Table S1',  r'PyWake, $N_{\rm free}$ from Ferreira et al. (2026)',],
    #          [r'PyWake, $N_{\rm free}$ from automated script', r'PyWake, $N_{\rm free}=min(a\sqrt{N_{\rm tot}}, N_{\rm tot})$,' + '\n' + '$a=5.3$ from Sørensen et al. (2024)']]
              [r'PyWake, $N_{\rm free}$ from automated script', r'PyWake, $N_{\rm free}=min(5.3\sqrt{N_{\rm tot}}, N_{\rm tot})$']]              
    labels2 = [['(a)', '(b)'], ['(c)', '(d)']]
    for i in range(2):       
        for j in range(2):          
            ax[i, j].text(0.205, 0.57, labels2[i][j], fontsize=14)
            ax[i, j].set_title(labels[i][j])
            ax[i, j].legend(loc=0)  
            ax[i, j].plot([0.2, 0.6], [0.2, 0.6], '-k', label=r'$y=x$')   
            ax[i, j].grid(True)                    
    fig.text(0.025, 0.45, 'Measured capacity factor', rotation=90, fontsize=14)
    fig.text(0.45, 0.025, 'Modeled capacity factor', fontsize=14)
    fig.tight_layout(rect=[0.05, 0.05, 1.0, 1.0])
    filename = 'fig6_Capacity_factor_model.pdf'
    fig.savefig(filename)    
