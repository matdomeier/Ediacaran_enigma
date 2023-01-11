import pygplates as pygp
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
from cartopy.geodesic import Geodesic
from shapely.geometry import Polygon as Pygon

def sph_car(v):
    """
    converts a unit vector [lat, lon] in spherical coordinates to cartesian [x,y,z] coordinates
    """
    
    x = np.cos(np.radians(v[0])) * np.cos(np.radians(v[1]))  
    y = np.cos(np.radians(v[0])) * np.sin(np.radians(v[1]))  
    z = np.sin(np.radians(v[0]))
    
    return [x,y,z]

def get_gcd(pt1, pt2):
    """
    calculates great circle distance between two points [lat, lon] in spherical coordinates
    """
    
    gcd = np.arccos(np.dot(sph_car(pt1), sph_car(pt2)))
    
    return np.degrees(gcd)


def get_rand_ranges(min_time, max_time, min_events, max_events, min_duration, max_duration):
    
    n_ranges = np.random.randint(min_events, max_events)

    starts, durations = [], []
    for i in range(n_ranges):
        t_dur = np.random.randint(min_duration, max_duration)
        t_sta = np.random.randint(min_time, max_time-t_dur)
        while any((t_sta >= starts[j] and t_sta <= starts[j]+durations[j]) 
                  or (t_sta+t_dur >= starts[j] and t_sta+t_dur <= starts[j]+durations[j])
                  for j in range(len(starts))):
            t_dur = np.random.randint(min_duration, max_duration)
            t_sta = np.random.randint(min_time, max_time-t_dur)
        durations.append(t_dur)
        starts.append(t_sta)

    ranges = [[starts[i], starts[i]+durations[i]] for i in range(len(starts))]
    sorted_ranges = sorted(ranges)
    
    return (sorted_ranges)


def generate_epoles(changepoint_ages, min_velocity, max_velocity, min_epole_distance):
    """
    Generates sequence of Euler poles; returned as list of pole locations and list of omegas.
    Input:  list of changepoint ages (sorted and ascending from 0), min and max plate velocities (mm/yr)
            and minimum angular distance between Euler poles and plate centroid
    """
    
    euler_poles, omegas = [], []
    for i in range(len(changepoint_ages)):
        if i == 0: 
            euler_poles.append([0, 0])
            omegas.append(0)
            
        else:
            # Euler poles generated from random points on sphere
            lon = np.degrees(2*np.pi * np.random.rand())
            cutoff = np.sin(np.radians(90-min_epole_distance))
            lat = np.degrees(np.arcsin(np.random.uniform(-cutoff, cutoff)))
            euler_poles.append([lat,lon])

            # get duration of each rotation
            time = changepoint_ages[i]-changepoint_ages[i-1]
            
            # compute omegas from randomly generated velocities
            velocity = np.random.randint(min_velocity, max_velocity)
            omega = velocity/(6371 * np.sin(np.radians(90-abs(lat)))) * time
            if np.random.rand(1) < 0.5: omegas.append(omega)
            else: omegas.append(omega*-1)
            
    return euler_poles, omegas


def assemble_rotations(centroid, euler_poles, omegas):
    """
    Assemble Euler poles into a rotation history
    Input: centroid, list of Euler poles, list of corresponding omegas
    """
    
    rotation_sequence = []
    init_rot = pygp.FiniteRotation([0, centroid[1]-90], -np.radians(90-centroid[0]))
    total_rotation = pygp.FiniteRotation([0, 0], 0)

    for i in range(len(euler_poles)):
        point = pygp.PointOnSphere(euler_poles[i])
        init_point = init_rot * point # rotate Euler pole to reference of plate at t=0
        
        # now apply total reconstruction sequence preceeding this rotation
        rotated_point = total_rotation * init_point
        stage_rotation = pygp.FiniteRotation(rotated_point, omegas[i])
        updated_total_rotation = stage_rotation * total_rotation ### NOTE: pygplates describes this operation in their tutorials the wrong way (the first rot is second)!!
        rotation_sequence.append(updated_total_rotation.get_lat_lon_euler_pole_and_angle_degrees())
        total_rotation = updated_total_rotation
        
    return rotation_sequence


def get_motion_path(time_range, time_step, plate_id, location, rotation_model):
    """
    Generates motion path of plate location from time t=0 back to some previous time specified by time range
    Input: time range (to reconstruct back to) and time-step size, plate id, location of point and pygplates rotation model
    """
    
    # assign starting point to pygplates point
    tracer = pygp.PointOnSphere(location)
    
    # generate motion path by reconstructing backward in time
    motion_path = []
    for i in range(0, time_range, time_step):
        # reconstruct tracer at each time step
        eulerpole = rotation_model.get_rotation(i, plate_id, 0)
        rot_tracer = eulerpole * tracer
        motion_path.append(rot_tracer.to_lat_lon())
    
    return motion_path


def get_apwp(time_range, time_step, plate_id, rotation_model):
    """
    Generates synthetic APWP
    Input: start and stop time (going forward toward t=0) and time-step size, plate id and pygplates rotation model
    """
    
    # generate APWP means as tracers at some specified time step
    apwp = []
    for i in range(0, time_range, time_step):
        # generate point at pole
        paleopole = pygp.PointOnSphere([90, 0])    

        # get stage rotation and rotate point to time=t0
        eulerpole = rotation_model.get_rotation(0, plate_id, i)
        rot_paleopole = eulerpole * paleopole
        apwp.append(rot_paleopole.to_lat_lon())
    
    return apwp


def norm_pdf(mu, sigma, x):
    
    term1 = 1/(sigma*np.sqrt(2*np.pi))
    term2 = np.exp(-1*(x-mu)**2/(2*sigma)**2)
    
    return term1*term2


def generate_tpw(tpw_onset, tpw_cease, tpw_lon, tpw_max, init_ang, pltid, refid):
    """
    Generates true polar wander episode (assuming IITPW)
    """ 
    
    stepsize = 0.1
    rel_rate = []
    xrange = np.arange(tpw_cease, tpw_onset, stepsize)
    mu = tpw_cease + (tpw_onset-tpw_cease)/2
    sigma = 2.8
    for x in xrange:
        rel_rate.append(norm_pdf(mu, sigma, x))
    
    # normalize rates and add to rotation file
    tpw_rate, tpw_history = [], []
    tpw_ang = init_ang
    for i in range(len(rel_rate)): 
        tpw_rate = (rel_rate[i]-min(rel_rate))/(max(rel_rate)-min(rel_rate)) * tpw_max * stepsize
        tpw_ang += tpw_rate
        tpw_history.append([pltid, xrange[i], 0, tpw_lon, tpw_ang, refid])
    
    return tpw_history


def update_rot_model(rotation_features, pltid, time, reference_position, desired_position):
    
    # original recon model
    original_rotation_model = pygp.RotationModel(rotation_features.get_features())
    
    # reconstruction time
    recon_time = pygp.GeoTimeInstant(time)
    
    # cast present_position to point feature
    pt_feature = pygp.Feature()
    pt_feature.set_reconstruction_plate_id(pltid)
    pt_feature.set_geometry(reference_position)

    # reconstruct point feature
    recon_feature = []
    pygp.reconstruct(pt_feature, original_rotation_model, recon_feature, recon_time)
    recon_pt = recon_feature[0].get_reconstructed_geometry()

    if recon_pt != desired_position:
        adjustment = pygp.FiniteRotation(recon_pt, desired_position)
        
    # now update plate model
    for rotation_feature in rotation_features.get_features():
        
        total_reconstruction_pole = rotation_feature.get_total_reconstruction_pole()
        if not total_reconstruction_pole: continue
            
        fixed_plate_id, moving_plate_id, rotation_sequence = total_reconstruction_pole
        if moving_plate_id != pltid: continue
        
        enabled_rotation_samples = rotation_sequence.get_enabled_time_samples()
        if not enabled_rotation_samples: continue

        if not (enabled_rotation_samples[0].get_time() <= recon_time and
                enabled_rotation_samples[-1].get_time() >= recon_time): continue

        rotation_property_value = rotation_sequence.get_value(recon_time)
        if not rotation_property_value: continue
        
        rotation = rotation_property_value.get_finite_rotation()
        
        fixed_plate_frame = original_rotation_model.get_rotation(recon_time, fixed_plate_id)
        fixed_plate_frame_adjustment = fixed_plate_frame.get_inverse() * adjustment * fixed_plate_frame
        adjusted_rotation = fixed_plate_frame_adjustment * rotation

        rotation_sequence.set_value(pygp.GpmlFiniteRotation(adjusted_rotation), recon_time)

    return rotation_features


def interpolate_rot_model(rotation_features, pltid, time, refid):
    
    # get rotation model
    original_rotation_model = pygp.RotationModel(rotation_features.get_features())
    
    # reconstruction time
    recon_time = pygp.GeoTimeInstant(time)
    
    for rotation_feature in rotation_features.get_features():

        total_reconstruction_pole = rotation_feature.get_total_reconstruction_pole()
        if not total_reconstruction_pole: continue

        fixed_plate_id, moving_plate_id, rotation_sequence = total_reconstruction_pole
        if moving_plate_id != pltid: continue

        enabled_rotation_samples = rotation_sequence.get_enabled_time_samples()
        if not enabled_rotation_samples: continue

        if not (enabled_rotation_samples[0].get_time() <= recon_time and
                enabled_rotation_samples[-1].get_time() >= recon_time): continue

        rotation_property_value = rotation_sequence.get_value(recon_time)
        if not rotation_property_value: continue

        rotation = rotation_property_value.get_finite_rotation()
        rotation_sequence.set_value(pygp.GpmlFiniteRotation(rotation), recon_time)

    return rotation_features
    

def plot_motion_paths(central_lon, central_lat, polygons, paths, colors, title):

    fig = plt.figure(figsize=(6, 6))
    proj = ccrs.Orthographic(central_longitude=central_lon, central_latitude=central_lat)
    proj._threshold = proj._threshold/50.  # this setting is necessary to keep sparsesly sampled great-circle segments smooth
    ax = plt.axes(projection = proj)
    ax.set_global()
    gl = ax.gridlines()
    gl.xlocator = mpl.ticker.FixedLocator([-180, -150, -120, -90, -60, -30, 0, 30, 60, 90, 120, 150, 180])
    gl.ylocator = mpl.ticker.FixedLocator([-90, -60, -30, 0, 30, 60, 90])

    # plot polygons if desired
    if polygons:
        for polygon in polygons:
            for geometry in polygon.get_geometries():
                vertices = geometry.to_lat_lon_array()
                plt.plot(vertices[:, 1], vertices[:, 0], transform=ccrs.Geodetic(), c='k', linewidth=0.5)

    # add motion paths
    for i, path in enumerate(paths):
        pathlats = [x[0] for x in path]
        pathlons = [x[1] for x in path]
        plt.plot(pathlons, pathlats, transform=ccrs.Geodetic(), zorder=1, linewidth=1.5, c=colors[i])
        plt.scatter(pathlons, pathlats, s=0.1, transform=ccrs.PlateCarree(), zorder=2, c='k')
        
    plt.title(title)
    
        
def plot_apwps(central_lon, central_lat, polygons, paths, A95s, colors, title):
    
    # setup plot for APWPs
    fig = plt.figure(figsize=(6, 6))
    proj = ccrs.Orthographic(central_longitude=central_lon, central_latitude=central_lat)
    proj._threshold = proj._threshold/50.  # this setting is necessary to keep sparsesly sampled great-circle segments smooth
    ax = plt.axes(projection = proj)
    ax.set_global()
    gl = ax.gridlines()
    gl.xlocator = mpl.ticker.FixedLocator([-180, -150, -120, -90, -60, -30, 0, 30, 60, 90, 120, 150, 180])
    gl.ylocator = mpl.ticker.FixedLocator([-90, -60, -30, 0, 30, 60, 90])

    # plot polygons if desired
    if polygons:
        for polygon in polygons:
            for geometry in polygon.get_geometries():
                vertices = geometry.to_lat_lon_array()
                plt.plot(vertices[:, 1], vertices[:, 0], transform=ccrs.Geodetic(), c='k', linewidth=0.5)

    # add APWPs
    for i, path in enumerate(paths):
        pathlats = [x[0] for x in path]
        pathlons = [x[1] for x in path]
        plt.plot(pathlons, pathlats, transform=ccrs.Geodetic(), zorder=1, linewidth=1.5, c=colors[i])
        plt.scatter(pathlons, pathlats, s=0.1, transform=ccrs.PlateCarree(), zorder=2, c='k')

    for i, poles in enumerate(A95s):
        for pole in poles:
            A95 = Pygon(Geodesic().circle(lon=pole[1], lat=pole[0], radius=3*111000))
            ax.add_geometries((A95,), crs=ccrs.PlateCarree().as_geodetic(), facecolor='none', linewidth=0.75, edgecolor=colors[i], alpha=0.5)
            
    plt.title(title)