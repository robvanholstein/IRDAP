from astropy.coordinates import EarthLocation
import numpy as np
from astropy.time import Time

def convert_ESOHeaderRADec_to_Degrees_PynPoint(tmp_ra, tmp_dec):
    """Convert the R.A. and Decl. values from ESO header `ESO INS4 DROT2 RA` and `ESO INS4 DROT2 DEC` to degrees. 
    # Code copied from https://github.com/PynPoint/PynPoint/blob/main/pynpoint/processing/psfpreparation.py#L542-L558
    Input:
        (1) ra_header, float. Should be direct reading from `ESO INS4 DROT2 RA` whose format is likely `(h)hmmss.ssss`.
        (2) dec_header, float. Should be direct reading from `ESO INS4 DROT2 DEC` whose format is likely `ddmmss.ssss`.
    Output:
        (1) ra_deg, float.
        (2) dec_deg, float.
        
    """
    # get sign of declination
    tmp_dec_sign = np.sign(tmp_dec)
    tmp_dec = np.abs(tmp_dec)

    # parse RA
    tmp_ra_s = tmp_ra % 100
    tmp_ra_m = ((tmp_ra - tmp_ra_s) / 1e2) % 100
    tmp_ra_h = ((tmp_ra - tmp_ra_s - tmp_ra_m * 1e2) / 1e4)

    # parse DEC
    tmp_dec_s = tmp_dec % 100
    tmp_dec_m = ((tmp_dec - tmp_dec_s) / 1e2) % 100
    tmp_dec_d = ((tmp_dec - tmp_dec_s - tmp_dec_m * 1e2) / 1e4)

    # get RA and DEC in degree
    ra = (tmp_ra_h + tmp_ra_m / 60. + tmp_ra_s / 3600.) * 15.
    dec = tmp_dec_sign * (tmp_dec_d + tmp_dec_m / 60. + tmp_dec_s / 3600.)

    return ra, dec
    
    
def parallactic_angles_IRDIS_correct(header_input):
    """Correct the START and END parallactic angles for IRDIS observations considering
        (1) DATE-OBS header
        (2) Actual latitude, longitude, elevation (actually not needed) from SPHERE
        (3) Overhead times from IRDIS
    Input: Original header
    Output: Updated header with new 'ESO TEL PARANG START' and 'ESO TEL PARANG END' values
    """
    # Adopted from https://github.com/PynPoint/PynPoint/blob/893ca12c414789b20b5a05dea0e51bbb6c907106/pynpoint/processing/psfpreparation.py#L448
    # overheads in cube mode (several NDITS) in hours
    param_O_START = 0.3 / 3600.        # Overhead for start, see SPHERE manual page 98 (Issue p107.2)
    param_ROT = 0.838 / 3600.          # Readout delay for IRDIS, see SPHERE manual page 98 (Issue p107.2)
    param_DIT_DELAY = 0.1 / 3600.      # Additional readout delay for IRDIS, see SPHERE manual page 98 (Issue p107.2)

    # rotator offset in degrees
    param_rot_offset = 0.              # no offset here
            
    # Adjusted from https://github.com/PynPoint/PynPoint/blob/893ca12c414789b20b5a05dea0e51bbb6c907106/pynpoint/processing/psfpreparation.py#L582
    
    # Load cube sizes
    steps = header_input['NAXIS3']      # should be the same as the next line
    ndit = header_input['ESO DET NDIT'] # should be the same as the above line
    if steps != ndit:
        printandlog('Warning: Likely corrupted files (NDIT is not equal to the number of slices in NAXIS3)! Continuing anyway.\n')

    # Load exposure time [hours]
    exptime = header_input['EXPTIME']/3600.

    # Load start times of exposures
    obs_date = header_input['DATE-OBS'] # actually both date and time are loaded and used later
    
    # Load telescope location
    tel_lat = header_input['ESO TEL GEOLAT']
    tel_lon = header_input['ESO TEL GEOLON']
    tel_elev = header_input['ESO TEL GEOELEV']
    location_SPHERE = EarthLocation(lat=tel_lat, lon=tel_lon, height = tel_elev)
    
    # Load temporary target position
    tmp_ra = header_input['ESO INS4 DROT2 RA']
    tmp_dec = header_input['ESO INS4 DROT2 DEC']
    
    ra, dec = convert_ESOHeaderRADec_to_Degrees_PynPoint(tmp_ra, tmp_dec) #do the RA Dec conversion

    t = Time(obs_date, location=location_SPHERE)

    sid_time = t.sidereal_time('apparent').value

    #start and end time arrays for IRDAP to read
    sid_time_arr = (sid_time + param_O_START) +  np.array([0 , (exptime + param_DIT_DELAY + param_ROT) * steps])
                                                           #start,  # end
    
    # Convert to degrees
    sid_time_arr_deg = sid_time_arr * 15.

    # Calculate hour angle in degrees
    hour_angle = sid_time_arr_deg - ra

    # Conversion to radians:
    hour_angle_rad = np.deg2rad(hour_angle)
    dec_rad = np.deg2rad(dec)
    lat_rad = np.deg2rad(tel_lat)

    p_angle = np.rad2deg(np.arctan2(np.sin(hour_angle_rad),
                         np.cos(dec_rad) * np.tan(lat_rad) - np.sin(dec_rad) * np.cos(hour_angle_rad)
                        ))
    header_input['ESO TEL PARANG START'] = p_angle[0] #new start parang.
    header_input['ESO TEL PARANG END'] = p_angle[1] #new end parang
    
    return header_input