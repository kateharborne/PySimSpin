#!/usr/bin/env python
import h5py
import math as math
import numpy as np

class PtypeNotValidException(Exception):
    pass

class M2lNotValidException(Exception):
    pass

class ApShapeNotValidException(Exception):
    pass

class Simulation:

    def __init__(self, filename, ptype="All", m2l="Solar", ssp=None):

        if ssp != None: # if ssp is supplied, assign the contained properties
            sspf = h5py.File(ssp, mode="r")
            self.ssp_wave = np.array(sspf[("/PartType4/Wavelength")])
            self.ssp_lum = np.array(sspf[("/PartType4/Luminosity")])
            sspf.close()

        f = h5py.File(filename, mode="r")
        ppart = list(f.keys())  # list of present particle types in file
        if ptype=="All":
            ptype = ppart   # if all particles are requested, set to ppart
        else:
            rptype = [] # placeholder for requested ptypes
            for i in ptype:
                rptype += ["PartType{}".format(i)]  # change requested ptypes to comparable format
            ptest = np.in1d(rptype, ppart)  # do all requested ptypes exist in ppart?
            if not all(ptest):  # if not, raise exception
                nppart = np.array([i for i, val in enumerate(ptest) if not val])
                raise PtypeNotValidException("The requested ptype(s) {} are not present in the simulation provided.".format(str([rptype[i] for i in nppart])))
            else:
                ptype = rptype  # if all requested types are present, set ptype to requested list
        
        if m2l == "Solar":  
                m2l = [1]*len(ptype)    # if m2l is default, make all m2l 1 for ppart
                if ssp!=None:
                    m2l.insert(ptype.index("PartType4"), np.nan)    # if SSP present, make m2l for stars NA

        pids = [pp[-1:] for pp in ptype] # getting the present particle types (number alone)
        if len(pids) != len(m2l):
            if ssp != None:
                m2l.insert(pids.index("4"), np.nan) # fixing the case in which ssp is provided and no place holder m2l value is given from ptype = 4
            else:
                raise M2lNotValidException("The requested particle types do not have matching mass to light ratios.")
        
        ds_x, ds_y, ds_z, ds_vx, ds_vy, ds_vz, ds_mass = ([] for i in range(7))
        part_type, x, y, z, vx, vy, vz, mass, luminosity = ([] for i in range(9))

        for i in ptype: # generating lists of the contained dataset names within HDF5 file
            ds_x += ["/{}/x".format(i)]
            ds_y += ["/{}/y".format(i)]
            ds_z += ["/{}/z".format(i)]
            ds_vx += ["/{}/vx".format(i)]
            ds_vy += ["/{}/vy".format(i)]
            ds_vz += ["/{}/vz".format(i)]
            ds_mass += ["/{}/Mass".format(i)]
        
        for idx in range(len(ptype)):   # extracting the data from each of those named datasets
            part_type.append([pids[idx]]*len(np.array(f[(ds_x[idx])])))
            x.append(np.array(f[(ds_x[idx])]))
            y.append(np.array(f[(ds_y[idx])]))
            z.append(np.array(f[(ds_z[idx])]))
            vx.append(np.array(f[(ds_vx[idx])]))
            vy.append(np.array(f[(ds_vy[idx])]))
            vz.append(np.array(f[(ds_vz[idx])]))
            mass.append(np.array(f[(ds_mass[idx])]))
            luminosity.append(np.array(f[(ds_mass[idx])])*(1e10 / m2l[idx]))

        self.ptype = np.concatenate(part_type, axis=None)   # assign all particle data to simulation properties
        self.x = np.concatenate(x, axis=None)   
        self.y = np.concatenate(y, axis=None)
        self.z = np.concatenate(z, axis=None)
        self.vx = np.concatenate(vx, axis=None)
        self.vy = np.concatenate(vy, axis=None)
        self.vz = np.concatenate(vz, axis=None)
        self.Mass = np.concatenate(mass, axis=None)
        self.Lum = np.concatenate(luminosity, axis=None)

        f.close()

        # open hdf5 file, read which particle types it contains, read data from requested particle types,
        # error if requested types don't exist

class Telescope:
    
    def __init__(self, fov, ap_shape, central_wvl, lsf_fwhm, pixel_sscale, pixel_vscale):
        self.fov = fov  # field of view of the IFU, diameter in arcsec
        self.ap_shape = ap_shape    # shape of the fov, options "Circular", "Hexagonal" or "Square"
        self.central_wvl = central_wvl  # central filter wavelength used, in angstroms
        self.lsf_fwhm = lsf_fwhm    # line spread function FWHM, in angstroms
        self.pixel_sscale = pixel_sscale    # spatial pixel scale, in arcsec
        self.pixel_vscale = pixel_vscale    # velocity puxel scale, in angstroms

class Observation:

    def __init__(self, telescope, z, inc_deg, cosmo):
        self.ang_size = 1 / cosmo.arcsec_per_kpc_proper(z)  # angular size for a given z, in kpc/arcsec
        self.lum_dist = cosmo.luminosity_distance(z)    # luminosity distance, in Mpc
        self.inc_rad = inc_deg * (math.pi / 180)    # galaxy inclination, in radians
        self.ap_size = self.ang_size * telescope.fov    # diameter of aperture, in kpc
        self.sbin = math.floor(telescope.fov / telescope.pixel_sscale)  # bin size in spatial plane, in pixels
        self.sbinsize = self.ap_size / self.sbin    # kpc per pixel
        self.vbinsize = (telescope.pixel_vscale / telescope.central_wvl) * (3e8 / 1e3)  # km/s per pixel
        self.lsf_size = ((telescope.lsf_fwhm / telescope.central_wvl) * (3e8 / 1e3)) / (2 * math.sqrt(2*math.log(2)))   # velocity uncertainty
        self.ap_region = np.zeros((self.sbin, self.sbin))   # empty matrix for aperture mask

        xcentre = self.sbin/2 + 0.5
        ycentre = self.sbin/2 + 0.5

        if telescope.ap_shape == "Circular":
            x = np.tile(np.arange(1, self.sbin+1), (self.sbin, 1))
            y = x.T
            xx = x - xcentre
            yy = y - ycentre
            rr = np.sqrt(np.square(xx) + np.square(yy))
            self.ap_region[rr<=self.sbin/2] = 1.    # circular aperture mask
        elif telescope.ap_shape == "Square":
            self.ap_region = np.ones((self.sbin, self.sbin))    # square aperture mask
        elif telescope.ap_shape == "Hexagonal":
            for x in list(range(1,self.sbin+1)):
                for y in list(range(1,self.sbin+1)):
                    xx = x - xcentre
                    yy = y - ycentre
                    rr = (2 * (self.sbin / 4) * (self.sbin * math.sqrt(3) / 4)) - ((self.sbin / 4) ) * abs(yy) - ((self.sbin * math.sqrt(3) / 4)) * abs(xx)
                    if rr >= 0 and abs(xx) < self.sbin/2 and abs(yy) < (self.sbin  * math.sqrt(3) / 4):
                        self.ap_region[x,y] = 1.    # hexagonal aperture mask
        else:
            raise ApShapeNotValidException("The ap_shape {0} is not valid. Please specify ap_shape = 'Circular', 'Square' or 'Hexagonal'".format(telescope.ap_shape))
            # error raised when invalid aperture shape is requested
            
class ObsGalaxy(Simulation, Observation):

    def __init__(self, filename, telescope, z, inc_deg, cosmo, ptype="All", m2l="Solar", ssp=None, centre=True):
        Simulation.__init__(self, filename, ptype="All", m2l="Solar", ssp=None)
        Observation.__init__(self, telescope, z, inc_deg, cosmo)

        if centre:
            xcen = np.median(self.x)
            ycen = np.median(self.y)
            zcen = np.median(self.z)
            vxcen = np.median(self.vx)
            vycen = np.median(self.vy)
            vzcen = np.median(self.vz)
            self.x = self.x - xcen
            self.y = self.y - ycen
            self.z = self.z - zcen
            self.vx = self.vx - vxcen
            self.vy = self.vy - vycen
            self.vz = self.vz - vzcen

        self.r, self.z_obs, self.vy_obs, self.r_obs = [],[],[],[]
    
        for i in range(len(self.x)):
            self.r.append(math.sqrt((self.x[i] * self.x[i]) + (self.y[i] * self.y[i]) + (self.z[i] * self.z[i])))
            self.z_obs.append((math.sin(self.inc_rad) * self.z[i]) + (math.cos(self.inc_rad) * self.y[i])) 
            self.vy_obs.append((math.cos(self.inc_rad) * self.vz[i]) - (math.sin(self.inc_rad) * self.vy[i]))
            self.r_obs.append(math.sqrt((self.x[i] * self.x[i]) + (self.z_obs[i] * self.z_obs[i])))
        
