'''this file .py needs to be submitted with the report for the exam'''

import numpy as np
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
        self.data_path = './group00_GJ3470_20220328/'
        
        self.readout_noise = 7.4  # [e] photoelectrons
        self.gain = 1.91 # [e/ADU]
        
        self.median_bias = pickle.load(open(self.data_path + 'median_bias.p', 'rb'))
        self.median_bias_errors = pickle.load(open(self.data_path + 'median_bias_error.p', 'rb'))
        # put the value here if you decided to use a constant value
        self.median_normalized_flat = pickle.load(open(self.data_path + 'median_normalized_flat.p', 'rb'))
        self.median_normalized_flat_errors = pickle.load(open(self.data_path + 'median_normalized_flat_errors.p', 'rb'))
        
        self.science_path = self.data_path + 'science/'
        self.science_list = np.genfromtxt(self.science_path + 'science.list', dtype=str)
        self.science_size = len(self.science_list)
        
        ylen, xlen  = np.shape(self.median_bias)
        self.X_axis = np.arange(0, xlen, 1)
        self.Y_axis = np.arange(0, ylen, 1)
        self.X, self.Y = np.meshgrid(self.X_axis, self.Y_axis)

        ## Here: read the first science frame and extract the information about the observatory
        # telescope_longitude, telescope_latitude, telescope_elevation, and  jd_observation (JD of the observation)
        self.observatory = coord.EarthLocation.from_geodetic(
            lon=telescope_longitude * u.deg,
            lat=telescope_latitude * u.deg,
            height=telescope_elevation * u.m)

        # SIMBAD / Gaia EDR3, ICRS coordinates at the epoch J2000
        target_j2000 = coord.SkyCoord( ...)
        reference_epoch = Time(jd_observation, format='jd', scale='utc', location=observatory)
        self.target = target_j2000.apply_space_motion(new_obstime=reference_epoch)

        
    def provide_aperture_parameters(self, sky_inner_radius, sky_outer_radius, aperture_radius, x_initial, y_initial):
        #You can decided to include a method to obtain the the new aperture photometry parameters
        #Or include them directly ad initialization
        self.sky_inner_radius = sky_inner_radius 
        self.sky_outer_radius = sky_outer_radius 
        self.aperture_radius = aperture_radius 
        self.x_initial = x_initial
        self.y_initial = y_initial


    
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
        self.x_fwhm = np.empty(self.science_size)
        self.y_fwhm = np.empty(self.science_size)
        
        for ii_science, science_name in enumerate(self.science_list):
            science_fits = fits.open(self.science_path + science_name)
            
            #You must read the info from the header before closing the file
            self.airmass[ii_science] = science_fits[0].header['AIRMASS']
            self.exptime[ii_science] = science_fits[0].header['EXPTIME']
            self.julian_date[ii_science] = science_fits[0].header['JD']
            
            science_data = science_fits[0].data * self.gain # save the data from the first HDU 
            science_fits.close()
            
            science_corrected, science_corrected_errors = self.correct_science_frame(science_data)
            x_refined, y_refined = self.compute_centroid(science_corrected, self.x_initial, self.y_initial, maximum_number_of_iterations=10)
            sky_background_value, sky_background_value_error = self.compute_sky_background(science_corrected,
                                                        science_corrected_errors, x_refined, y_refined)
            
            science_sky_corrected = ...
            science_sky_corrected_errors = ...
            x_refined, y_refined = self.compute_centroid(science_sky_corrected,  x_refined, y_refined)
            
            target_distance = np.sqrt((self.X-x_refined)**2 + (self.Y-y_refined)**2)
            aperture_selection = (target_distance < self.aperture_radius)
            self.aperture[ii_science] =  np.sum(science_sky_corrected[aperture_selection])
            self.aperture_error[ii_science] =  ...
            self.x_refined[ii_science] = x_refined
            self.y_refined[ii_science] = y_refined
            self.sky_background[ii_science] = sky_background_value

            try:
                self.x_fwhm[ii_science], self.y_fwhm[ii_science] = self.compute_fwhm(science_sky_corrected, x_refined, y_refined)
            except ValueError:
                self.x_fwhm[ii_science] = -1.0
                self.y_fwhm[ii_science] = -1.0

        tm = Time(self.julian_date, format='jd', scale='utc', location=('45.8472d', '11.569d')) 

        jd = self.julian_date + array_exptime/86400./2.
        tm_start = Time(array_jd, format='jd', scale='utc', location=self.observatory)
        tm = tm_start + TimeDelta(array_exptime / 2, format='sec')
        #Light travel time correction to the barycentre or heliocentre.
        ltt_bary = tm.light_travel_time(target, kind='barycentric')
        self.bjd_tdb = tm.tdb + ltt_bary
                
    
    def correct_science_frame(self, science_frame):
        science_debiased = science_frame - self.median_bias
        science_corrected = science_debiased / self.median_normalized_flat

        ## Error associated to the science corrected frame
        science_debiased_errors = ...
        science_corrected_errors = ...

        return science_corrected, science_corrected_errors
    
    
    def compute_centroid(self, science_frame, x_target_initial, y_target_initial, maximum_number_of_iterations=20):

        for i_iter in range(0, maximum_number_of_iterations):

            if i_iter == 0:
                # first iteration
                x_target_previous = x_target_initial
                y_target_previous = y_target_initial
            else:
                # using the previous result as starting point
                x_target_previous = x_target_refined
                y_target_previous = y_target_refined

            # 2D array with the distance of each pixel from the target star 
            target_distance = np.sqrt((self.X-x_target_previous)**2 + (self.Y-y_target_previous)**2)

            # Selection of the pixels within the inner radius
            annulus_sel = (target_distance < self.sky_inner_radius)

            # Weighted sum of coordinates
            weighted_X = np.sum(science_frame[annulus_sel]*self.X[annulus_sel])
            weighted_Y = np.sum(science_frame[annulus_sel]*self.Y[annulus_sel])

            # Sum of the weights
            total_flux = np.sum(science_frame[annulus_sel])

            # Refined determination of coordinates
            x_target_refined = weighted_X/total_flux
            y_target_refined = weighted_Y/total_flux

            variance_x = (x_target_refined-x_target_previous)
            variance_y = (y_target_refined-y_target_previous) 
            # exit condition: coordinates changed less than 0.1 pixels 
            if np.abs(percent_variance_x)<0.1 and  np.abs(percent_variance_y)<0.1:
                  break

        return x_target_refined, y_target_refined
        
    def _determine_FWHM_axis(self, reference_axis, normalized_cumulative_distribution):

        ...

        return FWHM

        
    def compute_fwhm(self, science_frame, x_refined, y_refined):
        # flux within the annulus 
        # 2D array with the distance of each pixel from the target star 
        
        ...

        return FWHM_x, FWHM_y
    
    def compute_sky_background(self, science_frame, science_frame_error, x_pos, y_pos):
        target_distance = np.sqrt((self.X-x_pos)**2 + (self.Y-y_pos)**2)

        annulus_selection = (target_distance > self.sky_inner_radius) & (target_distance<=self.sky_outer_radius)

        sky_flux_average = np.sum(science_frame[annulus_selection]) / np.sum(annulus_selection)
        sky_flux_median = np.median(science_frame[annulus_selection])
        sky_flux_error = ...
        
        return sky_flux_median, sky_flux_error