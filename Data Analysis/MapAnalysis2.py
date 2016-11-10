'''
!!!!
This script now only opens the new file structure. All analysis is performed by MapAnalysis.pyplot
!!!!

Change_log

6-9-16 Added gridfit interpolation option to test if it is appropriate for data
8-9-16 Calculate interpolation boundary from initial global centre rather than transformed centre- adjusts for distortion
5-10-16 All functions and data processing methods now use MapAnalysis- this script is for loading
data collected after the missing EMG bug fix
'''

import numpy as np
from scipy import interpolate
from scipy import optimize
import matplotlib.pyplot as plt
from matplotlib import mlab
from mpl_toolkits.mplot3d import Axes3D
import amcmorl_py_tools.vecgeom.coords as sphere
import amcmorl_py_tools.vecgeom.rotations as rotations
from amcmorl_py_tools.vecgeom.measure import angle_between
import pickle
import matplotlib.colors as colors
import csv
import matlab.engine, matlab

pretrigger_sample = int((0.15-0.005)*2000)
time = np.arange(-0.15, 1.1, 1.0/2000)

class RMS(object):
    pass


class MuscleMap(object):
    
    def  __init__(self, EMG_fname, pos_fname, target, action, target_RMS, control_RMS):
        self.target = target
        self.RMS = RMS()
        EMG_data = np.genfromtxt(EMG_fname, delimiter=',')
        self.pos_data = np.genfromtxt(pos_fname, delimiter=',')[:,1:]
        self.AD_EMG = EMG_data[np.where(EMG_data[:,2]==1)[0], 3:]*1000
        self.FDI_EMG = EMG_data[np.where(EMG_data[:,2]==0)[0], 3:]*1000        
        if target == 'AD':
            self.RMS.AD_low = target_RMS*0.05
            self.RMS.AD_high = target_RMS*0.15
            self.RMS.FDI_low = control_RMS*0.05
            self.RMS.FDI_high = control_RMS*0.15

        elif target == 'FDI':
            self.RMS.FDI_low = target_RMS*0.05
            self.RMS.FDI_high = target_RMS*0.15
            self.RMS.AD_low = control_RMS*0.05
            self.RMS.AD_high = control_RMS*0.15
        else:
            raise ValueError('Enter appropriate target') 
        self.task = action
        print self.FDI_EMG.shape[0],  self.pos_data.shape[0]
        if np.any(np.isnan(self.pos_data)):
            raise ValueError('nan in position data')

'''
def preprocess_EMG(MuscleMap):
    global time, pretrigger_sample
    if MuscleMap.target == 'AD':
        data = MuscleMap.AD_EMG
        target_RMS_low = MuscleMap.RMS.AD_low
        target_RMS_high = MuscleMap.RMS.AD_high
        alternative_data = MuscleMap.FDI_EMG
        alternative_RMS_low = MuscleMap.RMS.FDI_low
        alternative_RMS_high = MuscleMap.RMS.FDI_high
    elif MuscleMap.target == 'FDI':
        data = MuscleMap.FDI_EMG
        target_RMS_low = MuscleMap.RMS.FDI_low
        target_RMS_high = MuscleMap.RMS.FDI_high
        alternative_data = MuscleMap.AD_EMG
        alternative_RMS_low = MuscleMap.RMS.AD_low
        alternative_RMS_high = MuscleMap.RMS.AD_high
    
    EMG_mask = np.ones(data.shape[0]) + 1.
    
    for i, volt in enumerate(data):
        plt.ion()
        plt.plot(time, volt)
        plt.draw()
        while EMG_mask[i] == 2.:
            print i
            plt.pause(1)
            trial_good = raw_input('Is this trial good? [y]/[n]').lower()
            if trial_good == 'y':
                EMG_mask[i] = True
            elif trial_good =='n':
                EMG_mask[i] = False
            else:
                print 'Invalid: Try Again'
        plt.clf()
    
    EMG_mask[np.where(np.logical_or(np.sqrt(np.mean(np.square(data[:, :pretrigger_sample]), axis=1)) > (target_RMS_high),
                       np.sqrt(np.mean(np.square(data[:, :pretrigger_sample]), axis=1)) < (target_RMS_low)))[0]] = False
    if MuscleMap.task == 'both':
        EMG_mask[np.where(np.logical_or(np.sqrt(np.mean(np.square(alternative_data[:, :pretrigger_sample]), axis=1)) > (alternative_RMS_high),
                    np.sqrt(np.mean(np.square(alternative_data[:, :pretrigger_sample]), axis=1)) < (alternative_RMS_low)))[0]] = False
    else:
        EMG_mask[np.where(np.sqrt(np.mean(np.square(alternative_data[:, :pretrigger_sample]), axis=1)) > (0.01))[0]] = False
    MuscleMap.Amplitude, MuscleMap.MEParea = get_amplitude(data)
    MuscleMap.mask = EMG_mask.astype(bool)
    print MuscleMap.mask

def get_amplitude(data):
    windowstart = int((0.15+0.01)*2000)
    windowend = int((0.15+0.045)*2000)
    P2P = np.zeros(data.shape[0])
    for i, trial in enumerate(data):
        P2P[i] = trial[windowstart:windowend].max() - trial[windowstart:windowend].min()
    P2P[np.where(P2P < 0.2)] = 0.
    MEParea = np.trapz(np.abs(data[:, windowstart:windowend]), dx=1, axis=1) - np.trapz(np.abs(data[:, (0.15-0.04)*2000:(0.15-0.005)*2000]), dx=1, axis=1)
    return P2P, MEParea    

def fit_head(data, AP, ML):
    #create sphere bound
    u = np.linspace(0, 2*np.pi, 100)
    v = np.linspace(0, np.pi, 100)

    p0 = [0, 0, 0, 1] #initial parameters

    errfunc = lambda p, x: _fitfunc(p, x) - p[3]

    p1, flag = optimize.leastsq(errfunc, p0, args=(data,))
    
    R= p1[3] #radius of sphere
    
    X = R * np.outer(np.cos(u), np.sin(v)) 
    Y = R * np.outer(np.sin(u), np.sin(v)) 
    Z = R * np.outer(np.ones(np.size(u)), np.cos(v))  #convert sphere bounds to cartesian

    
    shifted_data = data - p1[:-1] #shift data to be at origin 0,0,0
    Cm = shifted_data.mean(axis=0) #find centre of data
    Cm_sphere = (R/(np.sqrt(sum(Cm**2))))*Cm #find projection of centre of data on sphere
    
    theta, phi  = sphere.cart2pol_v2(Cm_sphere/R) #convert CM xyz to polar
    map_AP, map_ML = AP/R, ML/R
    map_boundary = np.array([[theta-map_ML/2, theta+map_ML/2], 
                            [phi-map_AP/2, phi+map_AP/2]])
    print theta, phi,  map_AP, map_ML, map_boundary
    #find Cartesian boundary based on calculated polar map boundaries 
    mapboundaryinitial = sphere.pol2cart(np.array([[map_boundary[0,0],phi],
                                [map_boundary[0,1],phi],
                                [theta,map_boundary[1,0]],
                                [theta,map_boundary[1,1]],
                                ]))*R
    rot = np.array([theta- np.pi/2., phi]) #create rotation array
    #new_pos = rotations.rotate_by_angles(shifted_data.T, rot[0], rot[1], reverse_order=True).T #find new pos after rotation
    new_Cm = rotations.rotate_by_angles(Cm_sphere, rot[0], rot[1], reverse_order=True) #shift CM to new rotation
    polars =  sphere.cart2pol_v2(shifted_data/R)- rot #convert data points xyz to polar
    polars[:,0] = _zerotopi(polars[:,0])
    polars[:,1] = _pitopi(polars[:,1])
    new_pos = sphere.pol2cart(polars)*R
    map_boundary[0] -= rot[0]
    map_boundary[1] -= rot[1]
    mapboundarynew = sphere.pol2cart(np.array([[map_boundary[0,0],phi-rot[1]],
                                [map_boundary[0,1], phi-rot[1]],
                                [theta-rot[0], map_boundary[1,0]],
                                [theta-rot[0], map_boundary[1,1]],
                                ]))*R
    
    fig = plt.figure()
    ax = fig.gca(projection='3d')
    ax.plot_surface(X, Y, Z, rstride=1, cstride=1, alpha=0.2)
    ax.scatter(shifted_data[:,0], shifted_data[:,1], shifted_data[:,2], c='r', s=50)
    ax.scatter(new_Cm[0], new_Cm[1], new_Cm[2], c='b', s=50)
    ax.scatter(Cm_sphere[0], Cm_sphere[1], Cm_sphere[2], c='y', s=50)
    ax.scatter(new_pos[:, 0], new_pos[:, 1], new_pos[:, 2], c='g', s=50)
    ax.scatter(mapboundaryinitial[:,0], mapboundaryinitial[:,1], mapboundaryinitial[:,2],  c='y', s=50)
    ax.scatter(mapboundarynew[:,0], mapboundarynew[:,1], mapboundarynew[:,2],  c='b', s=50)
    plt.xlabel('X')
    plt.ylabel('Y')
    ax.set_zlabel('Z')
    plt.show()
    print rot
    return p1, rot, R, sphere.cart2pol_v2(new_Cm/R), map_boundary
 
def project_points(MuscleMap, p, t):
    u = np.linspace(0, 2*np.pi, 100)
    v = np.linspace(0, np.pi, 100)
    R= p[3]
    
    X = R * np.outer(np.cos(u), np.sin(v)) 
    Y = R * np.outer(np.sin(u), np.sin(v)) 
    Z = R * np.outer(np.ones(np.size(u)), np.cos(v)) 

    data = MuscleMap.pos_data
    shifted_data = data - p[:-1]
    MuscleMap.R = R
    
    Cm = shifted_data.mean(axis=0) #find centre of data
    Cm_sphere = (R/(np.sqrt(sum(Cm**2))))*Cm #find projection of centre of data on sphere
    projected_data = (R/np.sqrt(np.sum(shifted_data**2, axis=1)))[:,None]*shifted_data    

    #new_pos = rotations.rotate_by_angles(shifted_data.T, t[0], t[1], reverse_order=True).T
    new_Cm = rotations.rotate_by_angles(Cm_sphere, t[0], t[1], reverse_order=True)
    polars =  sphere.cart2pol_v2(projected_data/R) - t #convert data points xyz to polar(polar, azimuth)
    polars[:,0] = _zerotopi(polars[:,0])
    polars[:,1] = _pitopi(polars[:,1]) #corrections for near 2pi values of phi
    new_pos = sphere.pol2cart(polars)*R
    #polars = sphere.cart2pol_v2(projected_data/R)
   
    fig = plt.figure()
    ax = fig.gca(projection='3d')
    ax.plot_surface(X, Y, Z, rstride=1, cstride=1, alpha=0.2)
    ax.scatter(new_pos[:,0], new_pos[:,1], new_pos[:,2], c='r', s=50)
    ax.scatter(Cm_sphere[0], Cm_sphere[1], Cm_sphere[2], c='y', s=50)
    ax.scatter(new_Cm[0], new_Cm[1], new_Cm[2], c='b', s=50)
    plt.xlabel('X')
    plt.ylabel('Y')
    ax.set_zlabel('Z')
    ax.axis('equal')
    ax.axis('tight')
    plt.show()
    MuscleMap.polars = polars

def interpolate_data(MuscleMap, global_centre, mapboundary):
    p2p = MuscleMap.Amplitude
    position =  MuscleMap.polars
    MuscleMap.mask[np.isnan(p2p)] = False
    print np.sum(MuscleMap.mask)
    amp = p2p[np.where(MuscleMap.mask)]
    theta_values = position[:, 0][np.where(MuscleMap.mask)]
    phi_values = position[:, 1][np.where(MuscleMap.mask)]
    
    xi, theta_step = np.linspace(mapboundary[0,0], mapboundary[0,1], 750, endpoint=True, retstep=True)
    yi, phi_step = np.linspace(mapboundary[1,0], mapboundary[1,1], 750, endpoint=True, retstep=True)
    #X is theta which is the ML axis, Y is phi which is the AP axis!
    f = interpolate.Rbf(theta_values, phi_values, amp, function='gaussian', smooth=-0.05)
    X, Y = np.meshgrid(xi,yi)
    Z = f(X,Y)
    MuscleMap.area_multiplier = np.zeros((Z.shape))
    for i, line in enumerate(X):
        for j, val in enumerate(line):
            MuscleMap.area_multiplier[i,j] = (MuscleMap.R**2)*np.sin(val*theta_step*phi_step)
            #polar area = r^2 * sin(Theta * dtheta *dphi)
            #X is theta which is the ML axis, Y is phi which is the AP axis!
            
    fig, ax = plt.subplots(2, 1)
    plot = ax[0].pcolormesh(Y,X,Z, linewidths=0.5, cmap=plt.get_cmap('gnuplot2'))
    scatter = ax[1].scatter(phi_values,theta_values, amp)
    fig.colorbar(plot, ax=ax[0], extend='both')
    fig.show()
    plt.show()
    #plt.pcolormesh(np.where(Z >= amp.max()*0.1, Z, np.zeros_like(Z)))
    #plt.show()
    MuscleMap.theta = X
    MuscleMap.phi = Y
    MuscleMap.Z = Z

def interpolate_data_gridfit(MuscleMap, global_centre, mapboundary):
    p2p = MuscleMap.Amplitude
    position =  MuscleMap.polars
    MuscleMap.mask[np.isnan(p2p)] = False
    print np.sum(MuscleMap.mask)
    amp = p2p[np.where(MuscleMap.mask)]
    theta_values = position[:, 0][np.where(MuscleMap.mask)]
    phi_values = position[:, 1][np.where(MuscleMap.mask)]
    
    eng = matlab.engine.start_matlab()
    
    xi, theta_step = np.linspace(mapboundary[0,0], mapboundary[0,1], 750, endpoint=True, retstep=True)
    yi, phi_step = np.linspace(mapboundary[1,0], mapboundary[1,1], 750, endpoint=True, retstep=True)
    #X is theta which is the ML axis, Y is phi which is the AP axis!
    f = interpolate.Rbf(theta_values, phi_values, amp, function='gaussian', smooth=-0.05)
    X, Y = np.meshgrid(xi,yi)
    
    Z = np.array(eng.gridfit(matlab.double(theta_values.tolist()), matlab.double(phi_values.tolist()),\
    matlab.double(amp.tolist()), matlab.double(xi.tolist()),matlab.double(yi.tolist()), 'smoothness', matlab.double([2,2]))).reshape((750,750))

    MuscleMap.area_multiplier = np.zeros((Z.shape))
    for i, line in enumerate(X):
        for j, val in enumerate(line):
            MuscleMap.area_multiplier[i,j] = (MuscleMap.R**2)*np.sin(val*theta_step*phi_step)
            #polar area = r^2 * sin(Theta * dtheta *dphi)
            #X is theta which is the ML axis, Y is phi which is the AP axis!
            
    fig, ax = plt.subplots(2, 1)
    plot = ax[0].pcolormesh(Y,X,Z, linewidths=0.5, cmap=plt.get_cmap('gnuplot2'))
    scatter = ax[1].scatter(phi_values,theta_values, amp)
    fig.colorbar(plot, ax=ax[0], extend='both')
    fig.show()
    plt.show()
    #plt.pcolormesh(np.where(Z >= amp.max()*0.1, Z, np.zeros_like(Z)))
    #plt.show()
    MuscleMap.theta = X
    MuscleMap.phi = Y
    MuscleMap.Z = Z
    
def view_data(MuscleMap):
    p2p = MuscleMap.Amplitude
    amp = p2p[np.where(MuscleMap.mask)]
    phi_values = MuscleMap.polars[:,1][np.where(MuscleMap.mask)]
    theta_values = MuscleMap.polars[:,0][np.where(MuscleMap.mask)]

    fig = plt.figure()
    ax = fig.add_subplot(2,1,2, projection='3d')
    ax.plot_surface(MuscleMap.phi, MuscleMap.theta, MuscleMap.Z,cmap=plt.get_cmap('gnuplot2'))   
    plot2 = ax.scatter3D(phi_values, theta_values, amp,  c=amp, s=50, cmap=plt.get_cmap('RdPu'))
    fig.colorbar(plot2, ax=ax, extend='both')
    ax = fig.add_subplot(2, 1, 1)
    plot1 = ax.pcolormesh(MuscleMap.phi, MuscleMap.theta, MuscleMap.Z, linewidths=0.5, cmap=plt.get_cmap('gnuplot2'))
    fig.colorbar(plot1, ax=ax, extend='both')
    plt.show()

def get_map_parameters(AD_map, FDI_map, R):
    p2p_AD = AD_map.Amplitude[np.where(AD_map.mask)]
    p2p_FDI = FDI_map.Amplitude[np.where(FDI_map.mask)]
    AD_max = np.nanmax(p2p_AD)
    FDI_max = np.nanmax(p2p_FDI)
    AD_map.COG = _find_COG(AD_map.theta, AD_map.phi, AD_map.Z, AD_max, R)
    FDI_map.COG = _find_COG(FDI_map.theta, FDI_map.phi, FDI_map.Z, FDI_max, R)
    #X is theta which is the ML axis, Y is phi which is the AP axis!
    Distance =  angle_between(AD_map.COG,  FDI_map.COG)*R
    
    AD_map.theta_COG, AD_map.phi_COG = sphere.cart2pol_v2(AD_map.COG/R)
    FDI_map.theta_COG, FDI_map.phi_COG = sphere.cart2pol_v2(FDI_map.COG/R)
    
    Overlap_Volume = np.nansum(np.where(np.logical_and((AD_map.Z >= 0.1*AD_max),(FDI_map.Z >= 0.1*FDI_max)),
                        AD_map.Z, np.zeros_like(AD_map.Z)) \
        + np.where(np.logical_and((AD_map.Z >= 0.1*AD_max), (FDI_map.Z >= 0.1*FDI_max)),
                        FDI_map.Z, np.zeros_like(FDI_map.Z)))

    Total_Volume = np.nansum(np.where(AD_map.Z >= 0.1*AD_max, AD_map.Z, np.zeros_like(AD_map.Z)) \
        + np.where(FDI_map.Z >= 0.1*FDI_max, FDI_map.Z, np.zeros_like(FDI_map.Z)))

    Volume_ratio = Overlap_Volume/Total_Volume
    
    area_T = np.nansum(np.where(np.logical_or((AD_map.Z >= 0.1*AD_max), (FDI_map.Z >= 0.1*FDI_max)),
            AD_map.area_multiplier, np.zeros_like(AD_map.area_multiplier)))
            
    area_O = np.nansum(np.where(np.logical_and((AD_map.Z >= 0.1*AD_max), (FDI_map.Z >= 0.1*FDI_max)), \
            AD_map.area_multiplier, np.zeros_like(AD_map.area_multiplier)))
        
    area_R = area_O / area_T
    
    Overlay = np.where(AD_map.Z >= 0.1*AD_max, AD_map.Z, np.zeros_like(AD_map.Z)) \
        + np.where(FDI_map.Z >= 0.1*FDI_max, FDI_map.Z, np.zeros_like(FDI_map.Z))
    
    plt.pcolormesh(FDI_map.phi, FDI_map.theta, np.where(np.logical_and((AD_map.Z >= 0.1*AD_max), (FDI_map.Z > 0.1*FDI_max)), Overlay, np.zeros_like(Overlay)), cmap=plt.get_cmap('gnuplot2'))
    plt.colorbar(extend='both')
    plt.show()

    parameter = [AD_map.task, np.nanmax(p2p_AD), np.nanmax(p2p_FDI), AD_map.theta_COG, AD_map.phi_COG,  FDI_map.theta_COG, FDI_map.phi_COG,\
    Distance, area_T, area_O, area_R, Total_Volume, Overlap_Volume, Volume_ratio]
    print np.all(AD_map.phi == FDI_map.phi)
    print 'COG Distance = ', Distance, 
    print 'Total Volume = ', Total_Volume, 'Overlap_Volume = ', Overlap_Volume, 'Volume Ratio = ', Volume_ratio 
    print 'Area Total = ', area_T, 'Area Overlap = ', area_O, 'Area Ratio = ', area_R
    print  'AD Peak = ', np.nanmax(p2p_AD), 'FDI Peak = ', np.nanmax(p2p_FDI)
    return parameter
     
def save_data(MuscleMap, Folder):
    File = open(Folder + '\\' + '_'.join((MuscleMap.target,MuscleMap.task)) + '.pickle', 'wb')
    pickle.dump(MuscleMap, File, pickle.HIGHEST_PROTOCOL)
    File.close() 
    
def _fitfunc(p, coords):
    x0, y0, z0, R = p
    x, y, z = coords.T
    return np.sqrt((x-x0)**2 + (y-y0)**2 + (z-z0)**2)  #equation of sphere

def _find_COG(X, Y, Z, amp, R):
    theta_area = np.where(Z >=0.1*amp, X, np.zeros_like(X))
    phi_area = np.where(Z >= 0.1*amp, Y, np.zeros_like(Y))
    p = sphere.pol2cart(np.dstack((theta_area, phi_area)))*R
    print p.shape
    Z_area = np.where(Z >=0.1*amp, Z, np.zeros_like(Z))
    return np.array([np.divide(np.nansum(p[:,:, 0]*Z_area), np.nansum(Z_area)), np.divide(np.nansum(p[:,:,1]*Z_area),np.nansum(Z_area)), np.divide(np.nansum(p[:,:,2]*Z_area),np.nansum(Z_area))])
    
def _pitopi(radians):
    return radians - 2*np.pi*np.floor((radians+np.pi)/(2*np.pi))

def _zerotopi(radians):
    return radians - np.pi*np.floor(radians/(np.pi))
        
'''