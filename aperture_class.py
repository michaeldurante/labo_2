import numpy as np
import matplotlib.pyplot as plt
import pickle 
from astropy.io import fits
from astropy import coordinates as coord, units as u
from astropy.time import Time, TimeDelta


import sys

if not sys.warnoptions:
    import warnings
    warnings.simplefilter("ignore")

class AperturePhotometry:
    def __init__(self):
        self.data_path = './group09_HAT-P-12_20200422/'
        
        self.readout_noise = 7.1  
        self.gain = 1.91 
        
        self.median_bias = pickle.load(open(self.data_path + 'median_bias.p', 'rb'))
        self.median_bias_errors = pickle.load(open(self.data_path + 'median_bias_errors.p', 'rb'))
        
        self.median_normalized_flat = pickle.load(open(self.data_path + 'median_normalized_flat.p', 'rb'))
        self.median_normalized_flat_errors = pickle.load(open(self.data_path + 'median_normalized_flat_errors.p', 'rb'))
        
        self.science_path = self.data_path + 'science/'
        self.science_list = np.genfromtxt(self.science_path + 'science.list', dtype=str)
        self.science_size = len(self.science_list)
        
        ylen, xlen  = np.shape(self.median_bias)
        self.X_axis = np.arange(0, xlen, 1)
        self.Y_axis = np.arange(0, ylen, 1)
        self.X, self.Y = np.meshgrid(self.X_axis, self.Y_axis)

        lon, lat, height = self.get_observatory_coordinates(self)
        
        self.observatory = coord.EarthLocation(
            lon=lon * u.deg,
            lat=lat * u.deg,
            height=height * u.m
        )
        
        with fits.open(self.science_path + self.science_list[0]) as first_fits:
            jd_observation = first_fits[0].header['JD']

        target_J2000 = coord.SkyCoord(
            '13 57 33.47', '+43 29 36.60',
            unit = (u.hourangle, u.deg), frame='icrs', obstime=Time('J2000'),
            pm_ra_cosdec=-134.793 * u.mas/u.yr, pm_dec=-44.229*u.mas/u.yr,
            distance=(7.0448*u.mas).to(u.pc, equivalencies=u.parallax()),
            radial_velocity=-40.806702*u.km/u.s
            )
        
        reference_epoch = Time(jd_observation, format='jd', scale='utc', location=self.observatory)
        self.target = target_J2000.apply_space_motion(new_obstime=reference_epoch)

    def get_observatory_coordinates(self):
        science_fits = fits.open(self.science_path + self.science_list[0])
        longitude = science_fits[0].header['SITELONG']
        latitude = science_fits[0].header['SITELAT']
        elevation = science_fits[0].header['SITEELEV']

        return longitude, latitude, elevation
    
    def provide_aperture_parameters(self, sky_inner_radius, sky_outer_radius, aperture_radius, x_initial, y_initial):
            self.sky_inner_radius = sky_inner_radius 
            self.sky_outer_radius = sky_outer_radius 
            self.aperture_radius = aperture_radius 
            self.x_initial = x_initial
            self.y_initial = y_initial

    def correct_science_frame(self, science_frame):
        
        science_debiased = science_frame - self.median_bias
        science_corrected = science_debiased / self.median_normalized_flat
        
        science_debiased_errors = np.sqrt(
            self.readout_noise**2 + self.median_bias_errors**2 + np.abs(science_debiased)
        )
         
        science_corrected_errors = science_corrected * np.sqrt(
            (science_debiased_errors / science_debiased)**2 + 
            (self.median_normalized_flat_errors / self.median_normalized_flat)**2
        )
        
        return science_corrected, science_corrected_errors

    def compute_centroid(self, science_frame, x_target_initial, y_target_initial, maximum_number_of_iterations=20):

        x_target_refined = x_target_initial
        y_target_refined = y_target_initial
        
        for i_iter in range(maximum_number_of_iterations):
            if i_iter == 0:
                x_target_previous = x_target_initial
                y_target_previous = y_target_initial
            else:
                x_target_previous = x_target_refined
                y_target_previous = y_target_refined
                
            target_distance = np.sqrt((self.X - x_target_previous)**2 + (self.Y - y_target_previous)**2)
            annulus_sel = (target_distance < self.sky_inner_radius)
            
            weighted_X = np.sum(science_frame[annulus_sel] * self.X[annulus_sel])
            weighted_Y = np.sum(science_frame[annulus_sel] * self.Y[annulus_sel])
            total_flux = np.sum(science_frame[annulus_sel])
            
            x_target_refined = weighted_X / total_flux
            y_target_refined = weighted_Y / total_flux
            
            difference_x = np.abs(x_target_refined - x_target_previous)
            difference_y = np.abs(y_target_refined - y_target_previous)
            
            if np.abs(difference_x) < 0.1 and np.abs(difference_y) < 0.1:
                break
                
        return x_target_refined, y_target_refined
        
    
    def compute_sky_background(self, science_frame, science_frame_error, x_pos, y_pos):
        target_distance = np.sqrt((self.X - x_pos)**2 + (self.Y - y_pos)**2)
        annulus_selection = (target_distance > self.sky_inner_radius) & (target_distance <= self.sky_outer_radius)
        sky_flux_average = np.sum(science_frame[annulus_selection]) / np.sum(annulus_selection)
        sky_flux_median = np.median(science_frame[annulus_selection])
        sky_flux_error = np.sqrt(np.sum(science_frame_error[annulus_selection]**2)/(np.sum(annulus_selection)**2))
        #sky_flux_error = np.std(science_frame[annulus_selection]) / np.sqrt(len(science_frame[annulus_selection]))
        
        return sky_flux_median, sky_flux_error

    def aperture_photometry(self):
            self.airmass = np.empty(self.science_size)
            self.exptime = np.empty(self.science_size)
            self.julian_date = np.empty(self.science_size)
            
            self.aperture = np.empty(self.science_size)
            self.aperture_error = np.empty(self.science_size)
            self.x_refined = np.empty(self.science_size)
            self.y_refined = np.empty(self.science_size)
            self.sky_background  = np.empty(self.science_size)
            
            #Optional
            # self.x_fwhm = np.empty(self.science_size)
            # self.y_fwhm = np.empty(self.science_size)
            
            for ii_science, science_name in enumerate(self.science_list):
                science_fits = fits.open(self.science_path + science_name)
                
                self.airmass[ii_science] = science_fits[0].header['AIRMASS']
                self.exptime[ii_science] = science_fits[0].header['EXPTIME']
                self.julian_date[ii_science] = science_fits[0].header['JD']
                
                science_data = science_fits[0].data * self.gain
                science_fits.close()
                
                science_corrected, science_corrected_errors = self.correct_science_frame(science_data)
                x_refined, y_refined = self.compute_centroid(science_corrected, self.x_initial, self.y_initial, maximum_number_of_iterations=10)
                sky_background_value, sky_background_value_error = self.compute_sky_background(science_corrected,
                                                            science_corrected_errors, x_refined, y_refined)
                
                science_sky_corrected = science_corrected - sky_background_value
                science_sky_corrected_errors = np.sqrt(science_corrected_errors**2 + sky_background_value_error**2)
                x_refined, y_refined = self.compute_centroid(science_sky_corrected,  x_refined, y_refined)
                
                target_distance = np.sqrt((self.X-x_refined)**2 + (self.Y-y_refined)**2)
                aperture_selection = (target_distance < self.aperture_radius)
                self.aperture[ii_science] =  np.sum(science_sky_corrected[aperture_selection])
                self.aperture_error[ii_science] = np.sqrt(np.sum(science_sky_corrected_errors[aperture_selection]**2))
                self.x_refined[ii_science] = x_refined
                self.y_refined[ii_science] = y_refined
                self.sky_background[ii_science] = sky_background_value
    
                # try:
                #     self.x_fwhm[ii_science], self.y_fwhm[ii_science] = self.compute_fwhm(science_sky_corrected, x_refined, y_refined)
                # except ValueError:
                #     self.x_fwhm[ii_science] = -1.0
                #     self.y_fwhm[ii_science] = -1.0
    
            #jd = self.julian_date + self.exptime/86400./2.
            tm_start = Time(self.julian_date, format='jd', scale='utc', location=self.observatory)
            tm = tm_start + TimeDelta(self.exptime / 2, format='sec') #Light travel time correction to the barycentre or heliocentre.
            ltt_bary = tm.light_travel_time(self.target, kind='barycentric') 
            self.bjd_tdb = tm.tdb + ltt_bary

        # def _determine_FWHM_axis(self, reference_axis, normalized_cumulative_distribution):
        #     ...
        #     return FWHM

        # def compute_fwhm(self, science_frame, x_pos, y_pos):
        #  #flux with the annulus
        #  #2D array with distanceof each pixel from target

        #  ...

        #     return FWHM_x, FWHM_y