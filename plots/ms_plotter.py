"""Creates various plots for Balrog validation testing.


To run: $python ms_plotter.py base_path_to_catalogs output_directory realization tile
Example: $python ms_plotter.py /data/des71.a/data/kuropat/des2247-4414_sof/y3v02/ /Users/mspletts/BalValPlots/ all DES2247-4414

Relies on ms_matcher. User may need to replace `/data/des71.a/data/mspletts/balrog_validation_tests/scripts/BalVal/ms_matcher` with the correct path to ms_matcher.
FOF analysis relies on fof_matcher. User may need to replace `...` with the correct path to fof_matcher.

Plot attributes are specified with constants (many of them booleans) at the top of this script.
Constants that the user may wish to change are indicated by: # !!!!! {description and/or warnings} #. For example, user may wish to set `PRINTOUTS = False` or comment out `notice`.

# Comments are ABOVE the code they correspond to (with the exception of FIXMEs and TODOs). #

Megan Splettstoesser mspletts@fnal.gov"""


# astropy is needed only if analyzing a coadd catalog #
'''
from astropy.io import fits
from astropy.table import Column
from astropy.table import Table
'''
import matplotlib.pyplot as plt
import matplotlib
import numpy as np
import os
import pandas as pd
import subprocess
import sys




### Command line args ###
BASEPATH, OUTDIR, realizations, tiles = sys.argv[1], sys.argv[2], sys.argv[3].split(','), sys.argv[4].split(',')
# Catch error from inadequate number of command line args #
if len(sys.argv) != 5:
        sys.exit("Args: basepath (location of catalogs), output directory, realizations (can be 'all'), tiles (can be 'all') \n")


BALROG_RUN = BASEPATH[BASEPATH.replace('/', ';', 4).find('/')+1:BASEPATH.replace('/', ';', 5).find('/')]

ALL_FILTERS = [ 'g', 'r', 'i', 'z' ]

# !!!!! Number of realizations depends on the tile # 
ALL_REALIZATIONS = [ '0', '1', '2', '3', '4', '5', '6', '7', '8', '9' ]
#ALL_REALIZATIONS = [ '0', '1', '2' ]

# !!!!! Available tiles depends on run #
ALL_TILES = [ 'DES0347-5540', 'DES2329-5622', 'DES2357-6456' ]

if tiles != 'all':
	ALL_TILES = tiles
if realizations != 'all':
	ALL_REALIZATIONS = realizations





################################################################### Specify plot and catalog attributes ###################################################################

### Colorbar ###
# !!!!! Add colorbar according to one of the following (only one can be True at a time). If all are False a scatter-plot is made. Colorbars cannot be used if NORMALIZE is True. #
HEXBIN = False
CM_T_S2N_COLORBAR = False
CM_T_ERR_COLORBAR = False
CM_T_COLORBAR = False
BIN_CM_T_S2N = False
# Normalizes plot to 1-sigma magnitude error. If NORMALIZE is True, PLOT_1SIG must be True else errors will not be computed and normalization cannot be performed #
NORMALIZE = False 

# Use quality cuts introduced by Eric Huff? Link: https://github.com/sweverett/Balrog-GalSim/blob/master/plots/balrog_catalog_tests.py. Can only be performed if catalog has all the necessary headers: cm_s2n_r, cm_T, cm_T_err, and psfrec_T. #
EH_CUTS = False

# !!!!! Plot 1-sigma magnitude error curve? Must be True if NORMALIZE is True. If NORMALIZE is True will also plot the 68th percentile of the data in each error bin. #
PLOT_1SIG = True

# !!!!! What to do with the plot? #
SAVE_PLOT = False 
SHOW_PLOT = True

# !!!!! Limits for the vertical axis. 'None' is an allowed value and will result in default scaling #
YLOW, YHIGH = None, None

# Swap horizontal axis? Default is magnitude1. Matching script ms_matcher determines which catalog is 1 and which is 2. Generally SWAP_HAX does not need to be changed unless the truth catalog values are not on the horizontal axis. #
SWAP_HAX = False


### Catalog attributes ###
# !!!!! Allowed values: sof, mof, star_truth, gal_truth, coadd. Both can be 'sof' and both can be 'mof' if INJ1 and INJ2 are different. Note that truth catalogs always have INJ=True. #
MATCH_CAT1, MATCH_CAT2 = 'gal_truth', 'sof'
# !!!!! Booleans. Examine injected catalogs? #
INJ1, INJ2 = True, True
# !!!!! Must be used with realization=all at command line #
STACK_REALIZATIONS = False

# !!!!! Make 2x2 subplots of each griz filter? Or make individual plots? #
SUBPLOT = True

# !!!!! If directories do no exist, make them or force sys.exit() to edit dirs within script? NO_DIR_EXIT is invoked first, so sys.exit() will exit if NO_DIR_EXIT = True and NO_DIR_MAKE = True #
NO_DIR_EXIT = False
NO_DIR_MAKE = True




### For FOF analysis. To ignore FOF analysis set `RUN_TYPE=None` ###
# !!!!! Allowed values: 'ok' 'rerun' None. 'ok' refers to FOF groups unchanged after Balrog-injection. 'rerun' refers to groups changed after Balrog-injection. #
RUN_TYPE = None

# !!!!! Only used if RUN_TYPE is not None #
MOF = False
SOF = True

if RUN_TYPE is not None:
	print 'Doing FOF analysis ... \n '
        # Overwrite MATCH_CATs #
        if MOF:
                MATCH_CAT1, MATCH_CAT2 = 'mof', 'mof'
        if SOF:
                MATCH_CAT1, MATCH_CAT2 = 'sof', 'sof'
        INJ1, INJ2 = True, False


### !!!!! Make region files? #
MAKE_REG = False 

### Miscellaneous ###
# Print progress? #
PRINTOUTS = True
# Only refers to printouts within get_floats_from_string(), get_matrix_diagonal_element(), and bin_and_cut_measured_magnitude_error() # 
PRINTOUTS_MINOR = False

# Not currently in use or under constructrion #
LOG_FLAGS = False
PLOT_FLAGGED_OBJS = True
SHOW_FLAG_TYPE = False




# Catch errors from plot attributes #
if ((CM_T_S2N_COLORBAR and CM_T_ERR_COLORBAR) or (CM_T_S2N_COLORBAR and HEXBIN) or (CM_T_ERR_COLORBAR and HEXBIN) or (CM_T_COLORBAR and CM_T_ERR_COLORBAR) or (CM_T_COLORBAR and CM_T_S2N_COLORBAR) or (CM_T_COLORBAR and HEXBIN)) or (NORMALIZE and PLOT_1SIG is False) or (YLOW is not None and YHIGH is not None and YHIGH == YLOW) or (NORMALIZE and (CM_T_S2N_COLORBAR or CM_T_ERR_COLORBAR or CM_T_COLORBAR)) or (STACK_REALIZATIONS and realization != 'all') or (MOF and SOF):
	sys.exit('ERROR: Only of of the following may be True at a time: CM_T, CM_T_S2N_COLORBAR, CM_T_ERR_COLORBAR, HEXBIN. Otherwise, colorbar will be overwritten. \nERROR: If NORMALIZE is True so must be PLOT_1SIG. \nERROR: YHIGH cannot be equal to YLOW. \n NORMALIZE must be False if any of: CM_T_S2N_COLORBAR CM_T_ERR_COLORBAR CM_T_S2N_COLORBAR are True.\nERROR: STACK_REALIZATIONS is True must be used with realization = all. \nERROR: MOF and SOF cannot be simultaneously True. \n')


# !!!!! Check that plots will not be overwritten, etc #
NOTICE = raw_input(' \n !! CHECK BEFORE RUNNING !! \n Save plot(s) -- ' + str(SAVE_PLOT) + '\n Showing plot(s) -- ' + str(SHOW_PLOT) + '\n Normalize plot(s) -- ' + str(NORMALIZE) + '\n Hexbin -- ' + str(HEXBIN) + '\n cm_T colorbar -- ' + str(CM_T_COLORBAR) + '\n cm_T_err colorbar -- ' + str(CM_T_ERR_COLORBAR) + '\n cm_T_s2n -- ' + str(CM_T_S2N_COLORBAR) + '\n Plot limits -- ' + str(YLOW) + ', ' + str(YHIGH) + '\n Plotting 1-sigma curve -- ' + str(PLOT_1SIG) +'\n Plotting flagged objects -- ' + str(PLOT_FLAGGED_OBJS) + '\n Print flags and flag types -- ' + str(SHOW_FLAG_TYPE) + '\n Logging flags -- ' + str(LOG_FLAGS) + '\n --> Press enter to proceed, control+c to stop...\n')









################################################################### Store catalog information ###################################################################

class CoaddCat():
	"""Declare headers for coadd catalogs .../coadd/{tile}_{filter}_cat.fits. There is a separate catalog for each filter."""

        # Once matched, headers will have form 'hdr_1' or 'hdr_2' with a suffix (suf) #
	def __init__(self, inj, suf):
		"""Declare headers.

		Args:
			inj (bool) -- Balrog injected catalog?
			suf (int) -- Refers to order in which catalog was matched in ms_matcher. Allowed values: '_1' '_2'
		"""

		# For plot title #
		if inj:
			self.title_piece = 'Inj Coadd Cat'
		if inj is False:
			self.title_piece = 'Coadd Cat'
		self.axlabel = 'meas'
		# Magnitude, is one number #
		self.mag_hdr = 'MAG_AUTO' + str(suf)
		self.mag_axlabel = 'MAG_AUTO_meas'
		# For error calculation #
		self.mag_err_hdr = 'MAGERR_AUTO' + str(suf)
		self.cm_flux_hdr = None
		self.cm_flux_cov_hdr = None
		# Size #
		self.cm_t_hdr = None
		self.cm_t_err_hdr = None
		self.cm_t_s2n_axlabel = None
		# Flags #
		self.flags_hdr = 'FLAGS' + str(suf)
		self.obj_flags_hdr = None
		self.psf_flags_hdr = None
		self.cm_flags_hdr = None
		self.cm_max_flags_hdr = None
		self.cm_mof_flags_hdr = None
		self.cm_flags_r_hdr = None
		# For region file #
		self.ra_hdr = 'ALPHAWIN_J2000' + str(suf)
		self.dec_hdr = 'DELTAWIN_J2000' + str(suf)
		# Units: pixels #
		self.a_hdr = 'A_IMAGE' + str(suf)
		self.b_hdr = 'B_IMAGE' + str(suf)
		# For Eric Huff (EH) quality cuts #
                self.cm_s2n_r_hdr = None
                self.psfrec_t_hdr = None









class SOFGalTruthCat():
	"""Declare headers and axes labels for galaxy truth catalogs in the sof directory /data/des71.a/data/kuropat/des2247-4414_sof/."""

	# Once matched, headers will have form 'hdr_1' or 'hdr_2' with a suffix (suf) #
	def __init__(self, inj, suf):
		"""Declare constants.

                Args:
                        inj (bool) -- Balrog injected catalog?
			suf (int) -- Refers to order in which catalog was matched in ms_matcher. Allowed values: '_1' '_2'
                """

		if inj:
			self.title_piece = 'Inj Gal Truth Cat'
		if inj is False:
			self.title_piece = 'Gal Truth Cat'
		self.axlabel = 'true'
		# Headers are the same as MOFCat class. Reproduced below for clarity in ms_plotter.py #
		# Magnitude, is string of form '(mag_g, mag_r, mag_i, mag_z)' #
	    	self.mag_hdr = 'cm_mag' + str(suf)
		self.mag_axlabel = 'cm_mag_true'
                # For error calculation #
                self.mag_err_hdr = None
                self.cm_flux_hdr = 'cm_flux' + str(suf)
                self.cm_flux_cov_hdr = 'cm_flux_cov' + str(suf)
                # Size. cm_T units: arcseconds squared. #
                self.cm_t_hdr = 'cm_T'  + str(suf)
                self.cm_t_err_hdr = 'cm_T_err'  + str(suf)
		self.cm_t_s2n_axlabel = 'cm_T_s2n_true'
                # Flags #
                self.flags_hdr = 'flags' + str(suf)
                self.obj_flags_hdr = 'obj_flags' + str(suf)
                self.psf_flags_hdr = 'psf_flags' + str(suf)
                self.cm_flags_hdr = 'cm_flags' + str(suf)
                self.cm_max_flags_hdr = 'cm_max_flags' + str(suf)
                self.cm_flags_r_hdr = 'cm_flags_r' + str(suf)
		self.cm_mof_flags_hdr = 'cm_mof_flags' + str(suf)
                # For region file #
                self.ra_hdr = 'ra' + str(suf)
                self.dec_hdr = 'dec' + str(suf)
                self.a_hdr, self.b_hdr = None, None
                # For Eric Huff (EH) quality cuts #
		self.cm_s2n_r_hdr = 'cm_s2n_r' + str(suf)
                self.psfrec_t_hdr = 'psfrec_T' + str(suf)



	





class SOFCat():
        """Declare headers and axis labels for sof catalog in /data/des71.a/data/kuropat/des2247-4414_sof/y3v02/balrog_images/{realization}/{tile}/sof/{tile}_sof.fits and /data/des71.a/data/kuropat/sof_stars/y3v02/balrog_images/{realization}/{tile}/mof/{tile}_mof.fits"""

	# Once matched, headers will have form 'hdr_1' or 'hdr_2' with a suffix (suf) #
        def __init__(self, inj, suf):
		"""Declare constants.

                Args:
                        inj (bool) -- Balrog injected catalog?
			suf (int) -- Refers to order in which catalog was matched in ms_matcher. Allowed values: '_1' '_2'
                """

                if inj:
                        self.title_piece = 'Inj SOF Cat'
                if inj is False:
                        self.title_piece = 'SOF Cat'
		self.axlabel = 'meas'
                # Headers are the same as MOFCat class with the exception of cm_mof_flags_hdr. Reproduced below for clarity #
		# Magnitude, is string of form '(mag_g, mag_r, mag_i, mag_z)' #
		self.mag_hdr = 'cm_mag' + str(suf)
		self.mag_axlabel = 'cm_mag_meas'
                # For error calculation #
                self.mag_err_hdr = None
                self.cm_flux_hdr = 'cm_flux' + str(suf)
                self.cm_flux_cov_hdr = 'cm_flux_cov' + str(suf)
                # Size #
                self.cm_t_hdr = 'cm_T'  + str(suf)
                self.cm_t_err_hdr = 'cm_T_err'  + str(suf)
                self.cm_t_s2n_axlabel = 'cm_T_s2n_meas'
		# Flags #
                self.flags_hdr = 'flags' + str(suf)
                self.obj_flags_hdr = 'obj_flags' + str(suf)
                self.psf_flags_hdr = 'psf_flags' + str(suf)
                self.cm_flags_hdr = 'cm_flags' + str(suf)
                self.cm_max_flags_hdr = 'cm_max_flags' + str(suf)
                self.cm_flags_r_hdr = 'cm_flags_r' + str(suf)
                self.cm_mof_flags_hdr = None
                # For region file #
                self.ra_hdr = 'ra' + str(suf)
                self.dec_hdr = 'dec' + str(suf)
                self.a_hdr, self.b_hdr = None, None
                # For Eric Huff (EH) quality cuts #
		self.cm_s2n_r_hdr = 'cm_s2n_r' + str(suf)
                self.psfrec_t_hdr = 'psfrec_T' + str(suf)









class MOFCat():
	"""Declare headers and axes labels for MOF catalog. Currently, the galaxy truth catalogs are created using MOF and have the same headers. Works (mostly) with /data/des71.a/data/kuropat/sof_stars/y3v02/balrog_images/{realization}/{tile}/{tile}_{realization}_balrog_truth_cat_gals.fits, /data/des71.a/data/kuropat/sof_stars/y3v02/balrog_images/{realization}/{tile}/mof/{tile}_mof.fits, ..."""

	# Once matched, headers will have form 'hdr_1' or 'hdr_2' with a suffix (suf) #
        def __init__(self, inj, suf):
		"""Declare constants.

		Args:
			inj (bool) -- Balrog injected catalog?
			suf (int) -- Refers to order in which catalog was matched in ms_matcher. Allowed values: '_1' '_2'
		"""

		# For plot title #
		if inj:
			self.title_piece = 'Inj MOF Cat'
		if inj is False:
			self.title_piece = 'MOF Cat'
		self.axlabel = 'meas'
		# Magnitude, is string of form (mag_g, mag_r, mag_i, mag_z)  #
		self.mag_hdr = 'cm_mag' + str(suf)
		self.mag_axlabel = 'cm_mag_meas'
		# For error calculation #
		self.mag_err_hdr = None 
		self.cm_flux_hdr = 'cm_flux' + str(suf)
		self.cm_flux_cov_hdr = 'cm_flux_cov' + str(suf)
		# Size #
		self.cm_t_hdr = 'cm_T'  + str(suf)
		self.cm_t_err_hdr = 'cm_T_err'  + str(suf)
		self.cm_t_s2n_axlabel = 'cm_T_s2n_meas'
		# Flags #
		self.flags_hdr = 'flags' + str(suf)
                self.obj_flags_hdr = 'obj_flags' + str(suf)
                self.psf_flags_hdr = 'psf_flags' + str(suf)
                self.cm_flags_hdr = 'cm_flags' + str(suf)
                self.cm_max_flags_hdr = 'cm_max_flags' + str(suf)
                self.cm_flags_r_hdr = 'cm_flags_r' + str(suf)
		self.cm_mof_flags_hdr = 'cm_mof_flags' + str(suf)
		# For region file #
		self.ra_hdr = 'ra' + str(suf)
		self.dec_hdr = 'dec' + str(suf)
		self.a_hdr, self.b_hdr = None, None
		# For Eric Huff (EH) quality cuts #
                self.cm_s2n_r_hdr = 'cm_s2n_r' + str(suf) 
                self.psfrec_t_hdr = 'psfrec_T' + str(suf) 









class StarTruthCat(): #are there sep headers for MOFStarTruthCat and SOFStarTruthCat?
        """Declare headers and axes labels for star truth catalogs in /data/des71.a/data/kuropat/sof_stars/."""

	# Once matched, headers will have form 'hdr_1' or 'hdr_2' with a suffix (suf) #
	def __init__(self, inj, suf):
		"""Declare constants.
	
		Args:
			inj (bool) -- Balrog injected catalog?
			suf (int) -- Refers to order in which catalog was matched in ms_matcher. Allowed values: '_1' '_2'
		"""
		
		if inj:
			self.title_piece = 'Inj Star Truth Cat'
		if inj is False:
			self.title_piece = 'Star Truth Cat'
		self.axlabel = 'true'
		# Magnitude #
		self.mag_hdr = 'g_Corr' + str(suf) 
		self.mag_axlabel = 'mag_true' 
		# For error calculation #
                # Is of form 'PSF_MAG_ERR_{filter}' + str(suf) #
		self.mag_err_hdr = 'PSF_MAG_ERR' + str(suf)
                self.cm_flux_hdr = None
                self.cm_flux_cov_hdr = None
		# Size #
		self.cm_t_hdr = None
		self.cm_t_err_hdr = None
		self.cm_t_s2n_axlabel = None
		# Flags #
		self.flags_hdr = None
		self.obj_flags_hdr = None
		self.psf_flags_hdr = None
		self.cm_flags_hdr = None
		self.cm_max_flags_hdr = None
		self.cm_mof_flags_hdr = None
		self.cm_flags_r_hdr = None
		# For region file #
		self.ra_hdr = 'RA_new' + str(suf)
		self.dec_hdr = 'DEC_new' + str(suf)
		self.a_hdr, self.b_hdr = None, None # Use cm_T
		# For Eric Huff (EH) quality cuts #
		self.cm_s2n_r_hdr = None
		self.psfrec_t_hdr = None









def get_class(cat_type, inj, suf):
        """Get the appropriate class for the catalog type.

        Args:
                cat_type -- Catalog type. Allowed values: 'gal_truth', 'mof', 'star_truth', 'sof', 'coadd'.
                inj (bool) -- Balrog injected catalog?
		suf (int) -- Refers to order in which catalog was matched in ms_matcher. Allowed values: '_1' '_2'
        Returns:
                cat_type_class -- Points to the appropriate class which contains constants.
        """

        if cat_type == 'gal_truth':
                cat_type_class = SOFGalTruthCat(inj=inj, suf=suf)

        if cat_type == 'mof':
                cat_type_class = MOFCat(inj=inj, suf=suf)

        if cat_type == 'star_truth':
                cat_type_class = StarTruthCat(inj=inj, suf=suf)

        if cat_type == 'sof':
                cat_type_class = SOFCat(inj=inj, suf=suf)

        if cat_type == 'coadd':
                cat_type_class = CoaddCat(inj=inj, suf=suf)

        return cat_type_class









def get_match_type(title_piece1, title_piece2):
        """Transform plot title of form 'Inj MOF Cat & Truth Cat' to 'inj_mof_cat_truth_cat'.

        Args:
                title_piece1, title_piece2 (str) -- Ex: Injected MOF
        Return:
                match_type (str) -- Ex: injected_mof_truth_catalog
        """
        title_piece1, title_piece2 = title_piece1.lower(), title_piece2.lower()
        title_piece1, title_piece2 = title_piece1.replace(' ', '_'), title_piece2.replace(' ', '_')
        match_type = str(title_piece1)+'_'+str(title_piece2)

        return match_type









def get_fd_names():
        """Generate names for the following log files: flags, magnitude bins, number of objects plotted, number of objects within one sigma.
        Relies on directory structure: outdir/log_files/`BALROG_RUN`/`MATCH_TYPE`/

        Args:
                outdir (str) -- Output directory. Files will be saved here.
        Returns:
                fn1, fn2, fn3, fn4 (str) -- Filenames for flag log file, magnitude log, number of objects plotted log, number of objects within 1-sigma, respectively.
        """

        # !!!!! User may wish to edit directory structure #
        ### Check for directory existence ###
        if RUN_TYPE is None:
                log_dir = os.path.join(OUTDIR, 'log_files', BALROG_RUN, MATCH_TYPE)
        if RUN_TYPE is not None:
                log_dir = os.path.join(OUTDIR, 'log_files', BALROG_RUN, MATCH_TYPE, 'fof_analysis')

        if os.path.isdir(log_dir) is False:
                if NO_DIR_EXIT:
                        sys.exit('Directory ' + str(log_dir) + ' does not exist. \n Change directory structure in ms_plotter.get_fd_names() or set `NO_DIR_MAKE=True`')
                if NO_DIR_MAKE:
                        print 'Making directory ', log_dir, '...\n'
                        os.makedirs(log_dir)

        fn1 = os.path.join(log_dir, 'flag_log_'+str(MATCH_TYPE)+'.csv')
        fn2 = os.path.join(log_dir, 'magnitude_bins_'+str(MATCH_TYPE)+'.txt')
        fn3 = os.path.join(log_dir, 'num_objs_plotted_'+str(MATCH_TYPE)+'.txt')
        fn4 = os.path.join(log_dir, 'one_sigma_objects_'+str(MATCH_TYPE)+'.txt')

        if RUN_TYPE is not None:
		fn1 = fn1[:-4] + '_' + str(RUN_TYPE) + fn1[-4:]; fn2 = fn2[:-4] + '_' + str(RUN_TYPE) + fn2[-4:]
		fn3 = fn3[:-4] + '_' + str(RUN_TYPE) + fn3[-4:]; fn4 = fn4[:-4] + '_' + str(RUN_TYPE) + fn4[-4:]

        print '-----> Saving log file for flags as: ', fn1, '\n'
        print '-----> Saving log file for magnitude and error bins as: ', fn2, '\n'
        print '-----> Saving log file for number of objects plotted as: ', fn3, '\n'
        print '-----> Saving log file for number of objects within 1-sigma as: ', fn4, '\n'

	return fn1, fn2, fn3, fn4










def get_reg_names(tile_name, realization_number):
	"""Generate names for region files of different join types in STILTS script ms_matcher or fof_matcher.

        Args:
                outdir (str) -- Output directory. Files will be saved here.
        Returns:
                fn1, fn2, fn3, (str) -- Filenames for join=1and2, join=1not2, join=2not1, respectively. 
        """

        # !!!!! User may wish to edit directory structure #
        ### Check for directory existence ###
        if RUN_TYPE is None:
                reg_dir = os.path.join(OUTDIR, 'region_files', BALROG_RUN, MATCH_TYPE, tile_name, realization_number)
        if RUN_TYPE is not None:
                reg_dir = os.path.join(OUTDIR, 'region_files', BALROG_RUN, MATCH_TYPE, tile_name, realization_number, 'fof_analysis')


        if os.path.isdir(reg_dir) is False:
                if NO_DIR_EXIT:
                        sys.exit('Directory ' + str(reg_dir) + ' does not exist. \n Change directory structure in ms_plotter.get_fd_names() or set `NO_DIR_MAKE=True`')
                if NO_DIR_MAKE:
                        print 'Making directory ', reg_dir, '...\n'
                        os.makedirs(reg_dir)


        fn1 = os.path.join(reg_dir, str(tile_name) + '_' + str(realization_number) + '_' + str(MATCH_TYPE)+'_match1and2.reg')
        fn2 = os.path.join(reg_dir, str(tile_name) + '_' + str(realization_number) + '_' + str(MATCH_TYPE)+'_match1not2.reg')
        fn3 = os.path.join(reg_dir, str(tile_name) + '_' + str(realization_number) + '_' + str(MATCH_TYPE)+'_match2not1.reg')


	if RUN_TYPE is not None:
		fn1 = fn1[:-15] + '_' + str(RUN_TYPE) + fn1[-15:]; fn2 = fn2[:-15] + '_' + str(RUN_TYPE) + fn2[-15:]; fn3 = fn3[:-15] + '_' + str(RUN_TYPE) + fn3[-15:]


	return fn1, fn2, fn3







################################################################### Declare necessary constants ###################################################################

### For data analysis ###
# CLASS1 refers to in1 in ms_matcher. in1 appends _1 to all the headers, hence suf=1. fof_matcher is done such that injected catalogs have no suffix #
if RUN_TYPE is not None:
	CLASS1 = get_class(cat_type=MATCH_CAT1, inj=INJ1, suf='')
if RUN_TYPE is None:
	CLASS1 = get_class(cat_type=MATCH_CAT1, inj=INJ1, suf='_1')
CLASS2 = get_class(cat_type=MATCH_CAT2, inj=INJ2, suf='_2')

# Get arguments to pass to ms_matcher. Need to transform header of form 'ra_1' to 'ra', hence [:-2] #
RA_HDR1, RA_HDR2 = CLASS1.ra_hdr[:-2], CLASS2.ra_hdr[:-2]
DEC_HDR1, DEC_HDR2 = CLASS1.dec_hdr[:-2], CLASS2.dec_hdr[:-2]
# For plot labels #
AXLABEL1, AXLABEL2 = CLASS1.axlabel, CLASS2.axlabel

# Magnitudes #
M_HDR1, M_HDR2 = CLASS1.mag_hdr, CLASS2.mag_hdr
M_AXLABEL1, M_AXLABEL2 = CLASS1.mag_axlabel, CLASS2.mag_axlabel

# For magnitude error calculation #
M_ERR_HDR1, M_ERR_HDR2 = CLASS1.mag_err_hdr, CLASS2.mag_err_hdr
CM_FLUX_HDR1, CM_FLUX_HDR2 = CLASS1.cm_flux_hdr, CLASS2.cm_flux_hdr
CM_FLUX_COV_HDR1, CM_FLUX_COV_HDR2 = CLASS1.cm_flux_cov_hdr, CLASS2.cm_flux_cov_hdr

# For signal to noise calculation #
CM_T_HDR1, CM_T_HDR2 = CLASS1.cm_t_hdr, CLASS2.cm_t_hdr
CM_T_ERR_HDR1, CM_T_ERR_HDR2 = CLASS1.cm_t_err_hdr, CLASS2.cm_t_err_hdr
CM_T_S2N_AXLABEL1, CM_T_S2N_AXLABEL2 = CLASS1.cm_t_s2n_axlabel, CLASS2.cm_t_s2n_axlabel

# Flags #
FLAGS_HDR1, FLAGS_HDR2 = CLASS1.flags_hdr, CLASS2.flags_hdr
OBJ_FLAGS_HDR1, OBJ_FLAGS_HDR2 = CLASS1.obj_flags_hdr, CLASS2.obj_flags_hdr
# psf_flags is a string of form '(0,0,0,0)'; must pass through get_floats_from_string() if used. #
PSF_FLAGS_HDR1, PSF_FLAGS_HDR2 = CLASS1.psf_flags_hdr, CLASS2.psf_flags_hdr
CM_FLAGS_HDR1, CM_FLAGS_HDR2 = CLASS1.cm_flags_hdr, CLASS2.cm_flags_hdr
CM_MAX_FLAGS_HDR1, CM_MAX_FLAGS_HDR2 = CLASS1.cm_max_flags_hdr, CLASS2.cm_max_flags_hdr
CM_FLAGS_R_HDR1, CM_FLAGS_R_HDR2 = CLASS1.cm_flags_r_hdr, CLASS2.cm_flags_r_hdr
CM_MOF_FLAGS_HDR1, CM_MOF_FLAGS_HDR2 = CLASS1.cm_mof_flags_hdr, CLASS2.cm_mof_flags_hdr

# For quality cuts introduced by Eric Huff #
CM_S2N_R_HDR1, CM_S2N_R_HDR2 = CLASS1.cm_s2n_r_hdr, CLASS2.cm_s2n_r_hdr
PSFREC_T_HDR1, PSFREC_T_HDR2 = CLASS1.psfrec_t_hdr, CLASS2.psfrec_t_hdr

# For region file #
MAJOR_AX_HDR1, MAJOR_AX_HDR2 = CLASS1.a_hdr, CLASS2.a_hdr 
MINOR_AX_HDR1, MINOR_AX_HDR2 = CLASS1.b_hdr, CLASS2.b_hdr

FLAG_HDR_LIST = [ FLAGS_HDR1, FLAGS_HDR2, CM_FLAGS_HDR1, CM_FLAGS_HDR2, CM_MOF_FLAGS_HDR1, CM_MOF_FLAGS_HDR2, OBJ_FLAGS_HDR1, OBJ_FLAGS_HDR2, PSF_FLAGS_HDR1, PSF_FLAGS_HDR2, CM_MAX_FLAGS_HDR1, CM_MAX_FLAGS_HDR2, CM_FLAGS_R_HDR1, CM_FLAGS_R_HDR2 ]

# Used if LOG_FLAGS is True #
flag_idx = []


### For plot names, plot titles, log file names ###
TITLE_PIECE1, TITLE_PIECE2 = CLASS1.title_piece, CLASS2.title_piece
MATCH_TYPE = get_match_type(title_piece1=TITLE_PIECE1, title_piece2=TITLE_PIECE2)

### Names for file directors (fd) ###
FD_FLAG_NAME, FD_MAG_BINS_NAME, FD_NOP_NAME, FD_1SIG_NAME = get_fd_names()

# Create log file for number of objects plotted (nop) #
FD_NOP = open(FD_NOP_NAME, 'w')
# Create log file for number of objects within 1-sigma #
FD_1SIG = open(FD_1SIG_NAME, 'w')
# Create log file for magnitude bins #
FD_MAG_BINS = open(FD_MAG_BINS_NAME, 'w')
FD_MAG_BINS.write('NUM_OBJS_IN_BIN, BIN_LHS, BIN_RHS, MEDIAN_HAXIS_MAG, MEDIAN_ERROR \n') 
# Create log file for flags #
FD_FLAG = open(FD_FLAG_NAME, 'w')
FD_FLAG.write('TILE, FILTER, TYPE, REALIZATION, FLAG1_HEADER, FLAG2_HEADER, FLAG1_VALUE, FLAG2_VALUE, MAG1, MAG2 \n')
if LOG_FLAGS is False:
        FD_FLAG.write('Flags not logged because LOG_FLAGS is False.')









################################################################### Analysis ###################################################################
def get_floats_from_string(df, filter_name, hdr):
	"""Transform a list of strings of form '[ (1, 2, 3, 4), (1, 2, 3, 4), ... ]' to a list of floats of form '[1,1,...]' (if filter_name="g"), '[2,2,...]' ("r"), '[3,3,...]' ("i"), or '[4,4,...]' ("z").

	Args:
            df (pandas DataFrame)
            filter_name (str) -- Allowed values: 'g' 'r' 'i' 'z'.
            hdr (str) -- Header refers to a column name in the matched catalog. Must refer to a list of strings where each element is of form '(1,2,3,4)'.
        Returns:
            list_a (list of floats) -- Collection of the numbers corresponding to a particular index in a list of form '[ (1, 2, 3, 4), (1, 2, 3, 4), ... ]. 
	"""

	strings = df[hdr]; list_a = []

	# Each element (elmt) is of form '(1, 2, 3, 4)' #
	for elmt in strings:

		if filter_name == 'g':
			i = 1
			idx1 = elmt.find('(') + i
			idx2 = elmt.find(',')

		if filter_name == 'r':
			i = 2
			idx1 = elmt.replace(',', ';', 0).find(',') + i
			idx2 = elmt.replace(',', ';', 1).find(',')

		if filter_name == 'i':
			i = 2
			idx1 = elmt.replace(',', ';', 1,).find(',') + i
			idx2 = elmt.replace(',', ';', 2).find(',')

		if filter_name == 'z':
			i = 2
			idx1 = elmt.replace(',', ';', 2).find(',') + i
			idx2 = elmt.find(')')

		list_a.append(float(elmt[idx1:idx2]))

	if PRINTOUTS_MINOR:
		print 'Got ', hdr, ' for filter ', filter_name, '...'
		print ' Check: ', strings[0], ' & ', list_a[0], '\n'


	return list_a









def get_matrix_diagonal_element(df, filter_name, hdr):
	"""Transforms a list of 4x4 matrices where each element is a string of form '((11,12,13,14), (21,22,23,24), (31,32,33,34), (41,42,43,44))' into a list of either the 11 (if filter_name is "g"), 22 ("r"), 33 ("i"), or 44 ("z") matrix elements.

	Args:
            df (pandas DataFrame)
            filter_name (str) -- Allowed values: 'g' 'r' 'i' 'z'.
            hdr (str) -- Header refers to a column name in the matched catalog. Must refer to a list of strings where each element is of form '((11,12,13,14), (21,22,23,24), (31,32,33,34), (41,42,43,44))'.
        Returns:
            list_aa (list of floats) -- Collection of the numbers corresponding to a particular diagonal element in a list of 4-by-4 matrices.
	"""

	matrices = df[hdr]; list_aa = []

	# Each element in `matrices` is a matrix of form '((11,12,13,14), (21,22,23,24), (31,32,33,34), (41,42,43,44))' #
	for matrix in matrices:

		if filter_name == 'g':
 			i, j = 2, 0
			idx1 = 0
			idx2 = matrix.find(',')

		if filter_name == 'r':
			i, j = 2, 0
			idx1 = matrix.replace(',', ';', 4).find(',')
			idx2 = matrix.replace(',', ';', 5).find(',')

		if filter_name == 'i':
			i, j = 2, 0
			idx1 = matrix.replace(',', ';', 9).find(',')
			idx2 = matrix.replace(',', ';', 10).find(',')

		if filter_name == 'z':
			i, j = 2, -1
			idx1 = matrix.replace(',', ';', 14).find(',')
			idx2 = matrix.replace(',', ';', 15).find(',')

		list_aa.append(float(matrix[idx1+i:idx2+j]))

	if PRINTOUTS_MINOR:
		print 'Got ', hdr, ' for filter ', filter_name
		print ' Check: ', matrices[0], ' & ', list_aa[0], '\n'


	return list_aa









def get_good_index_using_primary_flags(df, full_magnitude1, full_magnitude2, cm_flag_hdr1, cm_flag_hdr2, flag_hdr1, flag_hdr2):
	"""Get indices of objects without flags as indicated by the headers 'flags' and 'cm_flags'. Also get indices of objects with  magnitudes not equal to +/- 99, +/- 9999, and 37.5. Store the bad indices as well (if PLOT_FLAGGED_OBJS is True).

	Args:
		df (pandas DataFrame)
		full_magnitude1, full_magnitude2 (list of floats) -- Uncleaned lists containing magnitudes. 
	Returns:
		idx_good (list of ints)
		idx_bad (list of ints) -- Is empty if PLOT_FLAGGED_OBJS is False.
	"""

	if cm_flag_hdr2 is None and cm_flag_hdr1 is None and flag_hdr1 is None and flag_hdr2 is None:
                sys.exit('No headers to clean flags with...')

	# If one catalog does not have the appropriate header, check it twice in the catalog that does have it so code still runs #
	if cm_flag_hdr2 is None:
                cm_flag_hdr2 = cm_flag_hdr1

        if cm_flag_hdr1 is None:
                cm_flag_hdr1 = cm_flag_hdr2

        if flag_hdr1 is None:
                flag_hdr1 = flag_hdr2

        if flag_hdr2 is None:
                flag_hdr2 = flag_hdr1

        ### Get flags ###
        flag1, flag2 = df[flag_hdr1], df[flag_hdr2]
        cm_flag1, cm_flag2 = df[cm_flag_hdr1], df[cm_flag_hdr2]


	# Make arrays to take absolute value in next step #
	full_magnitude1, full_magnitude2 = np.array(full_magnitude1), np.array(full_magnitude2)

	# Get rid of these objects 37.5 corresponds to a negative flux #
	idx_good= np.where( (abs(full_magnitude1) != 9999.0) & (abs(full_magnitude1) != 99.0) & (abs(full_magnitude1) != 37.5) & (abs(full_magnitude2) != 9999.0) & (abs(full_magnitude2) != 99.0) & (abs(full_magnitude2) != 9999.0) & (abs(full_magnitude2) != 99.0) & (abs(full_magnitude2) != 37.5) & (flag1 == 0) & (flag2 == 0) & (cm_flag1 == 0) & (cm_flag2 == 0) )[0]


	if PLOT_FLAGGED_OBJS:
		idx_bad = np.where( (abs(full_magnitude1) != 9999.0) & (abs(full_magnitude1) != 99.0) & (abs(full_magnitude1) != 37.5) & (abs(full_magnitude2) != 9999.0) & (abs(full_magnitude2) != 99.0) & (abs(full_magnitude2) != 9999.0) & (abs(full_magnitude2) != 99.0) & ((flag2 != 0) | (flag1 != 0) | (cm_flag1 != 0) | (cm_flag2 != 0)) )[0]

	if PLOT_FLAGGED_OBJS is False:
		idx_bad = None

		
        if PRINTOUTS:
                print 'Eliminated ', len(full_magnitude1) - len(idx_good), ' objects with magnitudes equal to +/- 9999, +/- 99, and 37.5 and objects with nonzero flags for: ', flag_hdr1, ', ', flag_hdr2, ', ', cm_flag_hdr1, ', ', cm_flag_hdr2, ' ... \n'


	return idx_good, idx_bad	









def get_good_index_using_quality_cuts(df, full_magnitude1, full_magnitude2):
	"""Get indices of objects that satisfy quality cuts introduced by Eric Huff. Also get indices of objects without flags as indicated by the headers 'flags' and 'cm_flags'. Also get indices of objects with  magnitudes not equal to +/- 99, +/- 9999, and 37.5. Store the bad indices as well (if PLOT_FLAGGED_OBJS is True).

	Args:
		df (pandas DataFrame)
	        *_hdr (str) -- Headers refer to column names in the matched catalog.
		full_magnitude1, full_magnitude2 (list of floats) -- Values read directly from pandas DataFrame or passed through `get_floats_from_string()`; no flags removed. 		
        Returns:
		idx_good (list of ints) -- Indices of objects without flags and objects which met criteria for quality cuts.
		idx_bad (list of ints) -- Is empty if PLOT_FLAGGED_OBJS is False.
	"""

	if 'true' in AXLABEL1 and 'true' in AXLABEL2:
		sys.exit('ERROR. Cuts should be performed on measured catalog, not truth catalog.')

	# Ignore truth catalog if one is present. Preserve ability to check both catalogs. #
	if 'true' in AXLABEL1 and 'meas' in AXLABEL2:
		CM_T_HDR1, CM_T_ERR_HDR1, CM_S2N_R_HDR1, PSFREC_T_HDR1 = CM_T_HDR2, CM_T_ERR_HDR2, CM_S2N_R_HDR2, PSFREC_T_HDR2

	if 'true' in AXLABEL2 and 'meas' in AXLABEL1:
		CM_T_HDR2, CM_T_ERR_HDR2, CM_S2N_R_HDR2, PSFREC_T_HDR2 = CM_T_HDR1, CM_T_ERR_HDR1, CM_S2N_R_HDR1, PSFREC_T_HDR1

	idx_good, idx_bad = [], []

	# If one catalog does not have the appropriate header, check it twice in the catalog that does have it so code still runs #
	if cm_flag_hdr2 is None:
		cm_flag_hdr2 = cm_flag_hdr1

	if cm_flag_hdr1 is None:
		cm_flag_hdr1 = cm_flag_hdr2


	### Define flags ###
	flag1, flag2 = df[flag_hdr1], df[flag_hdr2]
	cm_flag1, cm_flag2 = df[cm_flag_hdr1], df[cm_flag_hdr2]

	### Define parameters needed for quality cuts ###
	# Size squared of object #
	cm_t1, cm_t2 = df[CM_T_HDR1], df[CM_T_HDR2]
	# Size error #
	cm_t_err1, cm_t_err2 = df[CM_T_ERR_HDR1], df[CM_T_ERR_HDR2]
	# Signal to noise #
	cm_s2n_r1, cm_s2n_r2 = df[CM_S2N_R_HDR1], df[CM_S2N_R_HDR2]
	# PSF size #
	psfrec_t1, psfrec_t2 = df[PSFREC_T_HDR1], df[PSFREC_T_HDR2]

	# Cast into array to take absolute value in next step #
	full_magnitude1 = np.array(full_magnitude1)
	full_magnitude2 = np.array(full_magnitude2)

	idx_good = np.where( (abs(full_magnitude1) != 9999.0) & (abs(full_magnitude1) != 99.0) & (abs(full_magnitude1) != 37.5) & (abs(full_magnitude2) != 9999.0) & (abs(full_magnitude2) != 99.0) & (abs(full_magnitude2) != 9999.0) & (abs(full_magnitude2) != 99.0) & (abs(full_magnitude2) != 37.5) & (flag1 == 0) & (flag2 == 0) & (cm_flag1 == 0) & (cm_flag2 == 0) & (cm_s2n_r1 > 10) & (cm_s2n_r2 > 10) & (cm_t1/cm_t_err1 > 0.5) & (cm_t2/cm_t_err2 > 0.5) & (cm_t1/psfrec_t1 > 0.5) & (cm_t1/psfrec_t2 > 0.5) )[0]


	idx_bad = []


	if PLOT_FLAGGED_OBJS:
		counter_bad_mag = 0
		for i in np.arange(0, len(full_magnitude1)):

			# Get rid of these objects #
			if abs(full_magnitude1[i]) != 9999.0 and abs(full_magnitude2[i]) != 9999.0 and full_magnitude1[i] != 37.5 and full_magnitude2[i] != 37.5 and full_magnitude1[i] != 99.0 and full_magnitude2[i] != 99:
				counter_bad_mag += 1
				if flag1[i] != 0 or flag2[i] != 0 or cm_flag1[i] != 0 or cm_flag2[i] != 0 or cm_s2n_r1[i] < 10 or cm_s2n_r2[i] < 10 or cm_t1[i]/cm_t_err1[i] < 0.5 or cm_t2[i]/cm_t_err2[i] < 0.5 or cm_t1[i]/psfrec_t1[i] < 0.5 or cm_t2[i]/psfrec_t2[i] < 0.5:
					idx_bad.append(i)

	if PRINTOUTS:
		print 'Eliminated objects with magnitudes equal to +/- 9999, +/- 99, and 37.5 and objects with nonzero flags for: ', flag_hdr1, ', ', flag_hdr2, ', ', CM_FLAGS_HDR1, ', ', cm_flag_hdr2, ' ...'
		print ' Eliminated objects with signal-to-noise ratio < 10 ...'
		print ' Eliminated objects with cm_T/cm_T_err < 0.5 ...'
		print ' Eliminated objects with cm_T/psfrec_T < 0.5 ...'

	if PLOT_FLAGGED_OBJS is False:
		idx_bad = None


	return idx_good, idx_bad









def handle_flags(df, flag_hdr1, flag_hdr2, filter_name, full_magnitude1, full_magnitude2, realization_number, tile_name):
	"""Examine a particular flag and write to log file. Can also be used to check all flags in a list of flags.

	Args:
		df (pandas DataFrame)
		filter_name (str) -- Allowed values: 'g' 'r' 'i' 'z'.
		full_magnitude1, full_magnitude2 (numpy.ndarray if directly from `df[hdr]` OR list of floats if from `get_floats_from_string()`) -- Values read directly from pandas DataFrame via `df[hdr]`; no objects removed using nonzero flag values and no quality cuts performed.
		realization_number (int) -- Allowed values: 0 1 2 None. Refers to Balrog injection and None refers to a one-realization run.
		tile_name (str) -- 
        Returns:
		idx_good (list of ints) -- Indices of objects with flags values of zero.
		idx_bad (list of ints) -- Indices of objects with nonzero flag values.
	"""

	idx_good, idx_bad = [], []; counter_idx_bad = 0

	# If one catalog does not have the appropriate header, check it twice in the catalog that does have it so code still runs #
	if flag_hdr1 is None:
		flag_hdr1 = flag_hdr2

	if flag_hdr2 is None:
		flag_hdr2 = flag_hdr1


	### psf_flags are strings of form '(0,0,0,0)' ###
	if flag_hdr1 is not None and flag_hdr2 is not None and 'psf' not in flag_hdr1 and 'psf' not in flag_hdr2:
		flag1 = df[flag_hdr1]
		flag2 = df[flag_hdr2]

	if 'psf' in flag_hdr1 and 'psf' in flag_hdr2:
		flag1 = get_floats_from_string(df=df, hdr=flag_hdr1, filter_name=filter_name)
		flag2 = get_floats_from_string(df=df, hdr=flag_hdr2, filter_name=filter_name)


	### Check for flags ###
	for i in np.arange(0, len(full_magnitude1)):

		if abs(full_magnitude1[i]) != 9999.0 and abs(full_magnitude2[i]) != 9999.0 and full_magnitude1[i] != 37.5 and full_magnitude2[i] != 37.5 and full_magnitude1[i] != 99.0 and full_magnitude2[i] != 99:
			
			if flag1[i] == 0 and flag2[i] == 0:
				idx_good.append(i)
	
			if flag1[i] != 0 or flag2[i] != 0:
				idx_bad.append(i)
                                counter_idx_bad += 1

				### Write flags to file with headers TILE, FILTER, TYPE, REALIZATION, FLAG1_HEADER, FLAG2_HEADER, FLAG1_VALUE, FLAG2_VALUE, MAG1, MAG2 ###
				FD_FLAG.write(str(tile_name) + '\t' + str(filter_name) + '\t' + str(RUN_TYPE) + '\t' + str(realization_number) + '\t' + str(flag_hdr1) + '\t' + str(flag_hdr2) + '\t' + str(flag1[i]) + '\t' + str(flag2[i]) + '\t' + str(full_magnitude1[i]) + '\t' + str(full_magnitude2[i]) +'\n')



	if PRINTOUTS:
		print 'For tile: ', str(tile_name), ' and filter: ', str(filter_name), ', checked flags: ', flag_hdr1, ' & ', flag_hdr2, '...'

	### Check if flags were found ###
	if counter_idx_bad > 0 and PRINTOUTS:
		print ' Number of flags for magnitudes values 9999, 99, 37.5 and flags ', str(flag_hdr1), ' and ', str(flag_hdr2), ' : ', counter_idx_bad, '\n'


	return idx_good, idx_bad









def calculate_total_fractional_magnitude_error(cov_hdr, df, filter_name, flux_hdr, idx_good):
	"""Calculate the magnitude error via 1.08 * (flux_cov[i][i])^0.5 / flux[i] and ignore flagged objects in error calculation.

	Args:
            cov_hdr (str) -- Header for flux covariance matrix in the matched catalog.
            df (pandas DataFrame)
            filter_name (str) -- Allowed values: 'g' 'r' 'i' 'z'.
            flux_hdr (str) -- Headers refer to column names in the matched catalog.
            idx_good (list of ints) -- Indices with flag values equal to zero.
        Returns:
            error (list of floats) -- The magnitude error corresponding to each object.
	"""

	# Uncleaned lists for flux and flux covariance #
	full_flux = get_floats_from_string(df=df, hdr=flux_hdr, filter_name=filter_name)
	full_flux_cov = get_matrix_diagonal_element(df=df, hdr=cov_hdr, filter_name=filter_name)

	error, flux, fluxcov = [], [], []; counter_neg = 0

	# 'Safe' indices #
	for i in idx_good:
		flux.append(full_flux[i])
		fluxcov.append(full_flux_cov[i])

	# Calculations #
	for i in np.arange(0, len(flux)):

		# Throw out negative fluxcov (error calculation involves taking the square root of a negative) #
		if fluxcov[i] < 0:
			error.append(0)
			counter_neg += 1

		if fluxcov[i] == 0:
			print 'cm_flux_cov is 0'

		if fluxcov[i] > 0:
			err = 1.08 * fluxcov[i]**0.5 / flux[i] # Pogsons number = 1.08
			error.append(err)

	if PRINTOUTS:
		print 'Calculated the magnitude error for filter: ', filter_name
		print ' Number of negative cm_flux_cov: ', counter_neg, ' / ', len(flux), '\n'


	return error









def calculate_and_bin_cm_T_signal_to_noise(cm_t_hdr, cm_t_err_hdr, df, idx_good, clean_magnitude1, clean_magnitude2):
	"""Calculate measured signal-to-noise ratio via cm_T/cm_T_err (cuts performed on truth catalogs).

	Args:
		cm_t_hdr (str) -- Header for the size squared of object. Headers refer to column names in the matched catalog.
		cm_t_err_hdr (str) -- Header for error on the size squared of object. Headers refer to column names in the matched catalog.
		df (pandas DataFrame)
		idx_good (list of ints) -- Indices where no flags exist. 
		clean_magnitude1, clean_magnitude2 (list of floats) -- Magnitudes with flags removed.
	Returns:
		s2n (list of floats) -- cm_T signal-to-noise at each 'safe' index.
	"""

	cm_t = get_good_data(df=df, hdr=cm_t_hdr, idx_good=idx_good, magnitude=False, filter_name=None)
	cm_t_err = get_good_data(df=df, hdr=cm_t_err_hdr, idx_good=idx_good, magnitude=False, filter_name=None)

	# cm_T signal to noise (s2n) #
	cm_t_s2n = []
	for i in np.arange(0, len(cm_t)):
		cm_t_s2n.append(abs(cm_t[i]/cm_t_err[i]))

	# Bin signal to noise #
	# Bins suggested by Spencer Everett #
	bins = [0, 1, 9, 20, max(cm_t_s2n)]
	if PRINTOUTS:
		print 'Binning cm_T_s1n with bins: ', bins, ' and headers/axlabels:', cm_t_hdr, ', ', cm_t_err_hdr, '...'
		print ' Min and max absolute value of cm_T signal-to-noise: ', min(cm_t_s2n), ' and ', max(cm_t_s2n), '...'

	# idx_list is a list of lists to preserve bin structure #	
	binned_s2n, binned_hax_mag, binned_vax_mag, idx_list = [], [], [], []

	for j in np.arange(0, len(bins)-1):
		idx_temp = []
		for i in np.arange(0, len(cm_t_s2n)):
			if cm_t_s2n[i] > bins[j] and cm_t_s2n[i] < bins[j+1]:
				idx_temp.append(i)	
		if PRINTOUTS:
			print ' For cm_T_s2n, number of objects in bin ', bins[j], '-', bins[j+1], ': ', len(idx_temp)
		idx_list.append(idx_temp)
		#idx_temp = np.where(s2n > bins[j] & (s2n < bins[j+1]))

	if PRINTOUTS:
		print ' '


	return idx_list, bins, cm_t_s2n









def get_68percentile_from_normalized_data(norm_dm_list, bins, hax_mag_list):
	"""Calculate the point on the normalized vertical axis corresponding to the 68th percentile of the data for each bin used in the error calculation.

	Args:
		norm_dm_list (list of list of floats) -- Normalized delta magnitudes. Bin structure preserved.
		bins (list of floats) -- Bins used in error calculation.
		hax_mag_list (list of list of floats) -- Magnitudes on the horizontal axis. Bin structure preserved.
	Returns:
		vax_68percentile (list of floats) -- Point on the vertical axis (vax) corresponding to 68 percentile. Each element in the list corresponds to a different bin.
		bins (list of floats) -- Bins used in error calculation.
	"""

	vax_68percentile, neg_vax_34percentile, pos_vax_34percentile = [], [], []

	PLOT_HIST = False

	# Loop through bins (b) #
	for b in np.arange(0, len(norm_dm_list)):

		if norm_dm_list[b] is None:
			vax_68percentile.append(None)
			neg_vax_34percentile.append(None)
			pos_vax_34percentile.append(None)

		if norm_dm_list[b] is not None:

			### Find 68th percentile about zero ### 	
			# Values in current bin (icb) #
			vax_mag_list_icb = norm_dm_list[b]
			# Take absolute value of each point in bin #
			abs_vax_mag_list_icb = [abs(elmt) for elmt in vax_mag_list_icb]
			# Percentile sorts the data #
			vax_68percentile.append(np.percentile(abs_vax_mag_list_icb, 68, interpolation='lower'))	

			if PRINTOUTS_MINOR:
				# Check the percentile because interpolation='lower' was used #
				num = 0
				for j in np.arange(0, len(norm_dm_list[b])):
					if abs(norm_dm_list[b][j]) <= np.percentile(abs_vax_mag_list_icb, 68, interpolation='lower'):
						num += 1
				print 'Number of objects within 68 percentile via np.percentile(interpolation=lower): ', float(num)/len(norm_dm_list[b]), '...\n'


			### Find 34th percentile of positive and negative values separately ###
			neg_vax, pos_vax = [], []
			counter_neg, counter_pos = 0, 0
			for j in np.arange(0, len(vax_mag_list_icb)):
				if vax_mag_list_icb[j] < 0:
					neg_vax.append(vax_mag_list_icb[j])
					counter_neg += 1
				if vax_mag_list_icb[j] > 0:
					pos_vax.append(vax_mag_list_icb[j])
					counter_pos += 1

			# Check if lists are populated #
			if counter_neg > 0: 
				neg_vax_34percentile.append(np.percentile(neg_vax, 34, interpolation='lower'))
			if counter_pos > 0:
				pos_vax_34percentile.append(np.percentile(pos_vax, 34, interpolation='lower'))
			if counter_neg  == 0:
				neg_vax_34percentile.append(None)
			if counter_pos  == 0:
				pos_vax_34percentile.append(None)


			# Plot histogram to see distrubtion of data (data is not normally distributed) #
                        if PLOT_HIST:
                                plt.figure()
                                norm_dm = [abs(elmt) for elmt in norm_dm_list[b]]
                                plt.hist(norm_dm)
                                plt.title('Bin LHS: ' + str(bins[b]))
                                plt.xlabel(r'$\Delta M$')
                                plt.axvline(x=0, color='black', linestyle=':', linewidth=0.5)
                                plt.show()


	return vax_68percentile, bins, neg_vax_34percentile, pos_vax_34percentile









def bin_and_cut_measured_magnitude_error(clean_magnitude1, clean_magnitude2, error1, error2, filter_name):
        """Remove error values corresponding to objects where |Delta-M| > 3. Do not consider error corresponding to empty bins nor bins with a small number of objects.

        Args:
                clean_magnitude1, clean_magnitude2 (list of floats) -- Objects with flag values of zero and/or quality cuts performed.
                error1, error2 (list of floats) -- 1 and 2 refer to the matched catalogs. 
        Returns:
                binned_hax_mag_median (list of floats) -- List of medians of the horizontal axis magnitude in each bin.
                binned_vax_mag_median (list of floats) -- List of medians of the vertical axis magnitude in each bin. Vertical axis is computed via clean_magnitude1 - clean_magnitude2.
                binned_err_median (list of floats) -- Median of the error in each bin.
                bins (list of floats) -- Bins used. Binned according to horizontal axis.
		binned_hax_mag_list, binned_vax_mag_list, binned_err_list (list of lists of floats) -- Stores values in each bin (horizontal axis magnitude, vertical axis magnitude, error, respectively).
        """

	### !!!!! Comment this block out if errors are to be computed using both catalogs regardless of origin (measured catalog or truth catalog) ###
	if 'meas' in AXLABEL1 and 'meas' not in AXLABEL2:
		error2 = np.zeros(len(error1))
		if PRINTOUTS:
			print 'Using measured catalog (catalog1) for error calculation ... '

	if 'meas' in AXLABEL2 and 'meas' not in AXLABEL1:
		error1 = np.zeros(len(error2))
		if PRINTOUTS:
			print 'Using measured catalog (catalog2) for error calculation ... '

	if 'meas' in AXLABEL1 and 'meas' in AXLABEL2:
		if PRINTOUTS:
			print 'Using measured catalog (catalog1 AND catalog2) for error calculation ... '

	if 'true' in AXLABEL1 and 'true' in AXLABEL2:
		sys.exit('Errors are to be computed using the measured catalog(s), not the truth catalog(s).')


	### Define bins ###
        step = 0.5
        # Find the absolute min and max of the magnitudes in the matched catalog #
        limlow1, limlow2 = min(clean_magnitude1), min(clean_magnitude2)
        limhigh1, limhigh2 = max(clean_magnitude1), max(clean_magnitude2)
        limlow, limhigh = min([limlow1, limlow2]), max([limhigh1, limhigh2])

	# Define bins limits by ints #
	limlow, limhigh = int(limlow), int(limhigh)
	# Introduce magnitude cutoff to tame errors #
	if 'gal' in MATCH_CAT1 or 'gal' in MATCH_CAT2:
		limhigh = 26
	if 'star' in MATCH_CAT1 or 'star' in MATCH_CAT2:
		limhigh = 24

        if PRINTOUTS:
                print 'Forcing magnitudes to be binned with max ', limhigh, '...'

        bins = np.arange(limlow, limhigh, step)

	# Stores median of values in each bin #
        binned_hax_mag_median, binned_vax_mag_median, binned_err_median = [], [], []
	# List of lists. Stores all values in each bin #
	binned_hax_mag_list, binned_vax_mag_list, binned_err_list = [], [], []
        counter_empty_bin = 0


        # Bin magnitude errors according to the magnitude on the horizontal axis #
        if SWAP_HAX:
                hax_mag = clean_magnitude2
        if SWAP_HAX is False:
                hax_mag = clean_magnitude1

	# Magnitude on the vertical axis (vax) #	
	vax_mag = np.array(clean_magnitude1) - np.array(clean_magnitude2)


	### Write filter header to log file ###
	FD_MAG_BINS.write('Filter: ' + str(filter_name) + '\n')


	### Populate each bin ###
	for j in np.arange(limlow, limhigh, step):

                binned_hax_mag_temp, binned_vax_mag_temp, binned_err_temp, counter_err = [], [], [], 0

                for i in np.arange(0, len(clean_magnitude1)):

                        # Do not calculate errors using outlier magnitudes (chosen to be |Delta-M| > 3). Bin magnitude errors according to the magnitude on the horizontal axis of the plot #
                        if hax_mag[i] >= j and hax_mag[i] < j+step and abs(vax_mag[i]) < 3:
                                binned_err_temp.append((error1[i]**2 + error2[i]**2)**0.5)
                                binned_hax_mag_temp.append(hax_mag[i])
                                binned_vax_mag_temp.append(vax_mag[i])
                                counter_err += 1

		# Written in log file, hence 'minor' #
                if PRINTOUTS_MINOR:
                        print ' For magnitude, number of objects in bin ', round(j, 2), '-', round(j+step, 2), ': ', counter_err, '...'


		### Write to log file ###
		if counter_err == 0:
			write_median, write_err = None, None
		if counter_err > 0:
			write_median, write_err = np.median(binned_hax_mag_temp), np.median(binned_err_temp)
		FD_MAG_BINS.write(str(counter_err) + '\t' + str(round(j, 2)) + '\t' + str(round(j+step, 2)) + '\t' + str(write_median)+ '\t' + str(write_err) + '\n')


                ### Tame error calculation and normalization by adding zeros to empty bins and bins with a small number of points ###
		# Define 'small' #
		if STACK_REALIZATIONS:
			CONST = 30
		if STACK_REALIZATIONS is False:
			CONST = 10

		if counter_err <= CONST:
                        counter_empty_bin += 1
                        binned_err_median.append(None)
                        binned_hax_mag_median.append(None)
                        binned_vax_mag_median.append(None)

			# Add to list of lists to keep bin structure #
			binned_err_list.append(None)
                        binned_hax_mag_list.append(None)
                        binned_vax_mag_list.append(None)		

                if counter_err > CONST:
                        binned_err_median.append(np.median(binned_err_temp))
                        binned_hax_mag_median.append(np.median(binned_hax_mag_temp))
                        binned_vax_mag_median.append(np.median(binned_vax_mag_temp))
	
			# Add to list of lists to keep bin structure #	
			binned_err_list.append(binned_err_temp)
                        binned_hax_mag_list.append(binned_hax_mag_temp)
                        binned_vax_mag_list.append(binned_vax_mag_temp)


	if PRINTOUTS:
                if SWAP_HAX:
                        print 'Binned clean_magnitude2 with step size: ', step, ', and minimum: ', limlow, ', and maximum: ', limhigh, '...'
		if SWAP_HAX is False:
                        print 'Binned clean_magnitude1 with step size: ', step, ', and minimum: ', limlow, ', and maximum: ', limhigh, '...'
		print ' Calculated errors using objects where |DeltaM| < 3 ... '
		print ' Excluded bins with less than ', CONST, ' objects ... \n'


        return binned_hax_mag_median, binned_vax_mag_median, binned_err_median, bins, binned_hax_mag_list, binned_vax_mag_list, binned_err_list









def normalize_plot_maintain_bin_structure(clean_magnitude1, clean_magnitude2, error1, error2, filter_name):
	"""Normalize the vertical axis using error and uphold the bin structure. 

	Args:
		clean_magnitude1, clean_magnitude2 (list of floats) --
		error1, error2 (list of floats) --
	Returns:
		norm_dm_list (list of list of floats) --
		bins (list of floats) -- Bins used in error calculation. 
	"""

	# List of lists. Stores all values in each bin #
	norm_dm_list, hax_mag_list = [], []

	# binned_err_median: stores median of vales in bin. *_list: stores all values in each bin #
	binned_err_median, bins, binned_hax_mag_list, binned_vax_mag_list = bin_and_cut_measured_magnitude_error(clean_magnitude1=clean_magnitude1, clean_magnitude2=clean_magnitude2, error1=error1, error2=error2, filter_name=filter_name)[2:-1]


	# Loop through bins (b) #
	for b in np.arange(0, len(binned_vax_mag_list)):

		# Normalized Delta-Magnitudes (dm) in current bin (icb) #
		norm_dm_icb, hax_mag_icb = [], []	


		# 0 is a placeholder for empty bins and bins with few objects #
		if binned_err_median[b] is None:
			norm_dm_list.append(None)	
			hax_mag_list.append(None)


		#if vax_mag_icb != 0:
		if binned_err_median[b] is not None:
			
			vax_mag_icb = binned_vax_mag_list[b]

			for i in np.arange(0, len(vax_mag_icb)):
				norm_dm_icb.append(vax_mag_icb[i]/binned_err_median[b])
				hax_mag_icb.append(binned_hax_mag_list[b][i])

			# List of lists to keep bin structure #
			hax_mag_list.append(hax_mag_icb)
			norm_dm_list.append(norm_dm_icb)

	return norm_dm_list, bins, hax_mag_list









def normalize_plot(norm_delta_mag_list, bins, hax_mag_list):
	"""Normalize plot to 1-sigma curve using tame magnitude errors only (use bin_and_cut_measured_magnitude_error()).

	Args:
		norm_dm_list (list of list of floats) -- Normalized delta magnitudes in each bin. 
                bins (list of floats) -- Bins used in error calculation.
                hax_mag_list (list of list of floats) -- Magnitudes on the horizontal axis. Bin structure preserved.
        Returns:
		norm_dm (list of floats) -- Delta-Magnitude normalized by error. Delta-Magnitude computed via magnitude1 - magnitude2. 
                hax_mag (list of floats) -- Magnitude to be plotted on the horizontal axis.
	"""

	### Remove zeros so that lists can be flattened. Zeros (int) were placeholders for missing lists due to empty or small bin. ###
	norm_delta_mag_list[:] = [temp for temp in norm_delta_mag_list if temp is not None]
	hax_mag_list[:] = [temp for temp in hax_mag_list if temp is not None]

	### Flatten lists ###
	hax_mag = [item for sublist in hax_mag_list for item in sublist]
	norm_dm = [item for sublist in norm_delta_mag_list for item in sublist]


	### Check ###
	idx = []
	for b in np.arange(0, len(bins)-1):
		for j in np.arange(0, len(hax_mag)):
			if hax_mag[j] >= bins[b] and hax_mag[j] <= bins[b+1]:
				idx.append(j)
	norm_dm, hax_mag = np.array(norm_dm), np.array(hax_mag)
	norm_dm, hax_mag = norm_dm[idx], hax_mag[idx]

	return norm_dm, hax_mag, bins









def one_sigma_counter(norm_delta_mag, clean_magnitude1, bins, hax_mag):
	"""Find the number of objects within 1-sigma, where 1-sigma is calculated according to the error. This function is only called if NORMALIZE is True.

	Args:
		norm_delta_mag (list of floats) -- Normalized Delta-Magnitude. 
        Returns:
            counter_1sig (int) -- Number of objects within 1-sigma curve.
	"""

	counter_1sig = 0


	# Cutoffs were introduced in error calculation. Consider only points not cutoff #
	maglow, maghigh = min(bins), max(bins)
	hax_mag, norm_delta_mag = np.array(hax_mag), np.array(norm_delta_mag)
	norm_delta_mag = norm_delta_mag[(hax_mag >= maglow) & (hax_mag <= maghigh)]

	for k in norm_delta_mag:
		if abs(k) < 1.0:
			counter_1sig += 1

	if PRINTOUTS:
		print 'Fraction of objects within 1-sigma: ', counter_1sig, ' / ', len(norm_delta_mag), ' = ', str(float(counter_1sig) / len(norm_delta_mag))
		print ' Fraction of objects considered (objects plotted on normalized plot / objects plotted on scatter plot): ', str(float(len(norm_delta_mag)) / len(clean_magnitude1)), '\n'


	return counter_1sig









def get_flag_type(df, k):
	"""Print the flag type() once.

	Args:
            df (pandas DataFrame)
            k (int) -- Counter implemented so printout is not repeated.
        Returns:
            0
	"""

	if k == 0:
		for flag_hdr in FLAG_HDR_LIST:
			print 'HEADER:', str(flag_hdr), ' -- EXAMPLE:', df[flag_hdr][0], ' -- TYPE:', type(df[flag_hdr][0])
			k += 1
        return 0









def get_color(filter_name):
	"""Color code plot such that each griz band is a different color.

	Args:
		filter_name (str) -- Allowed values: 'g' 'r' 'i' 'z'
	Returns:
		color (str)
		cmap (str) -- Colormap used for Delta-Magnitude colorbars.
	"""

	if filter_name == 'g':
		color, cmap = 'green', 'Greens'

	if filter_name == 'r':
		color, cmap = 'orange', 'Oranges'

	if filter_name == 'i':
		#color, cmap = 'purple', 'Purples'
		color, cmap = 'darkgrey', 'Greys'

	if filter_name == 'z':
		#color, cmap = 'blue', 'Blues'
		color, cmap = 'navy', 'Blues'

	return color, cmap









def logger(delta_mag, tile_name, filter_name, realization_number, clean_magnitude1, full_magnitude1, bins, hax_mag):
	"""Write to log files to record number of objects plotted and number of objects within 1sigma.

	Args:
            filter_name (str) -- Allowed values: 'g' 'r' 'i' 'z'.
            clean_magnitude1 (list of floats) -- Objects with nonzero flags and/or quality cuts removed.
            full_magnitude (list of floats) -- Values read directly from pandas DataFrame via `df[hdr]`; no objects removed using nonzero flag values and no quality cuts performed.
            realization_number (int) -- Allowed values: 0 1 2 None. Refers to Balrog injection and None refers to a one-realization run.
        Returns:
            0
	"""

	if NORMALIZE:
		num_1sig = one_sigma_counter(norm_delta_mag=delta_mag, clean_magnitude1=clean_magnitude1, bins=bins, hax_mag=hax_mag)

		# Record number of objects plotted within 1sigma #
		FD_1SIG.write('Number of objects within 1sigma for tile ' + str(tile_name) + ', filter ' + str(filter_name) + ', type ' + str(RUN_TYPE) + ', realization ' + str(realization_number) + ' : ' + str(num_1sig) + ' / ' + str(len(clean_magnitude1)) + ' = ' + str(float(num_1sig) / len(clean_magnitude1)) + '\n')


	# Record number of objects plotted (nop) #
	FD_NOP.write('Number of objects plotted after flag cuts for tile ' + str(tile_name) + ', filter ' + str(filter_name) + ', type ' + str(RUN_TYPE) + ', realization ' + str(realization_number) + ' : ' + str(len(clean_magnitude1)) + ' / ' + str(len(full_magnitude1)) + ' = ' + str(float(len(clean_magnitude1)) / len(full_magnitude1)) + '\n')


        return 0









def get_colorbar_value(df, cm_t_hdr, cm_t_err_hdr, idx_good, clean_magnitude1, clean_magnitude2, axlabel, inj):
	"""Get data that will be used for the colorbar of plot.

	Args:
		df (pandas DataFrame)
		*_hdr (str) -- Headers refer to columns in the matched catalog.
		inj (bool)  
	Returns:
		cbar_val -- Values used to make colorbar.
		cbar_idx_list -- Can be None
		cbar_bins -- Can be None
		cbar_axlabel (str) -- Label for the colorbar.
	"""

	if 'true' in axlabel:
                sys.exit('ERROR. Colorbars should describe measured catalog values, not truth catalog values.')


        if CM_T_S2N_COLORBAR:
                cbar_idx_list, cbar_bins, cbar_val = calculate_and_bin_cm_T_signal_to_noise(cm_t_hdr=cm_t_hdr, cm_t_err_hdr=cm_t_err_hdr, df=df, idx_good=idx_good, clean_magnitude1=clean_magnitude1, clean_magnitude2=clean_magnitude2)
		cbar_axlabel = 'cm_T_s2n_'+str(AXLABEL)

        if CM_T_ERR_COLORBAR:
                # For measured catalog, cuts performed on truth catalogs #
                cbar_val = get_good_data(df=df, hdr=cm_t_err_hdr, idx_good=idx_good, magnitude=False, filter_name=None)
		cbar_axlabel = str(cm_t_err_hdr[:-2]) + '_' + str(AXLABEL)
		cbar_idx_list, cbar_bins = None, None

        if CM_T_COLORBAR:
                cbar_val = get_good_data(df=df, hdr=cm_t_hdr, idx_good=idx_good, magnitude=False, filter_name=None)
		cbar_axlabel = str(cm_t_hdr[:-2]) + '_' + str(AXLABEL)
                cbar_idx_list, cbar_bins = None, None

	if CM_T_S2N_COLORBAR is False and CM_T_ERR_COLORBAR is False and CM_T_COLORBAR is False:
		cbar_val, cbar_idx_list, cbar_bins, cbar_axlabel = None, None, None, None

	if inj and cbar_axlabel is not None:
		cbar_axlabel = 'inj_' + cbar_axlabel
	

	return cbar_val, cbar_idx_list, cbar_bins, cbar_axlabel









def get_errors(mag_err_hdr1, mag_err_hdr2, df, filter_name, idx_good):
	"""Get errors for plot.

	Args:
		*_hdr (str ) -- Headers refer to columns in the matched catalog. Can be None.
		df (pandas DataFrame)
	Returns:
		err1, err2 (list of floats) -- Will be None if PLOT_1SIG is False.
	"""

        if PLOT_1SIG:
                if mag_err_hdr1 is None:
                        err1 = calculate_total_fractional_magnitude_error(df=df, flux_hdr=CM_FLUX_HDR1, cov_hdr=CM_FLUX_COV_HDR1, filter_name=filter_name, idx_good=idx_good)
                if mag_err_hdr1 is not None:
			if MATCH_CAT1 == 'coadd':
				err1 = get_floats_from_string(df=df, hdr=mag_err_hdr1, filter_name=filter_name)
			if MATCH_CAT1 == 'star_truth':
                                err1 = df[str(mag_err_hdr1[:-2]) + '_' + filter_name.upper() + str(mag_err_hdr1[-2:])]

                if mag_err_hdr2 is None:
                        err2 = calculate_total_fractional_magnitude_error(df=df, flux_hdr=CM_FLUX_HDR2, cov_hdr=CM_FLUX_COV_HDR2, filter_name=filter_name, idx_good=idx_good)
                if mag_err_hdr2 is not None:
			if MATCH_CAT2 == 'coadd':
				err2 = get_floats_from_string(df=df, hdr=mag_err_hdr2, filter_name=filter_name)
			if MATCH_CAT2 == 'star_truth':
				err2 = df[str(mag_err_hdr2[:-2]) + '_' + filter_name.upper() + str(mag_err_hdr2[-2:])] 

        if PLOT_1SIG is False:
                print 'WARNING: Not plotting 1-sigma curve so log file will FALSELY report that ZERO objects are within 1sigma ...\n'
                err1, err2 = None, None

	return err1, err2









def get_good_data(df, hdr, idx_good, magnitude, filter_name):
	"""Get the data corresponding to good indices (no flags or post quality cuts).

	Args:
		df (pandas DataFrame)
		hdr (str) -- Header for the DataFrame.
		idx_good (list of floats) -- Safe indices
		magnitude (bool) -- Get data for magnitudes?
		filter_name (str) -- Only used if magnitude is True.
	"""

	if magnitude:
		full_data = get_floats_from_string(df=df, hdr=hdr, filter_name=filter_name)
	if magnitude is False:
                full_data = df[hdr]

	clean_data = np.array(full_data)[idx_good]

	return clean_data









def get_plot_variables(filter_name, df, mag_hdr1, mag_hdr2, mag_err_hdr1, mag_err_hdr2, realization_number, tile_name, mag_axlabel1, mag_axlabel2):
	"""Get quantities needed for plotter() and subplotter().

	Args:
		df (pandas DataFrame)
		
		*_hdr (str) -- Can be None.
	Returns:
	"""

	 # Rewrite mag_axlabels. Transform, for example, cm_mag_true to cm_mag_{filter}_true or psf_mag_meas to psf_mag_{filter}_meas #
        if RUN_TYPE is None:
                mag_axlabel1 = str(mag_hdr1[:-2]) + '_' + str(AXLABEL1)
        if RUN_TYPE is not None:
                mag_axlabel1 = str(mag_hdr1) + '_' + str(AXLABEL1)
        mag_axlabel2 = str(mag_hdr2[:-2]) + '_' + str(AXLABEL2)

        # Transform, for example, cm_mag_true to cm_mag_{filter}_true, or psf_mag_meas to psf_mag_{filter}_meas #
        mag_axlabel1 = mag_axlabel1[:-4] + str(filter_name) + '_' + mag_axlabel1[-4:]
        mag_axlabel2 = mag_axlabel2[:-4] + str(filter_name) + '_' + mag_axlabel2[-4:]

        if INJ1:
                mag_axlabel1 = 'inj_' + mag_axlabel1
        if INJ2:
                mag_axlabel2 = 'inj_' + mag_axlabel2
	
	# Coadd catalogs. Combined to get '(m_g, m_r, m_i, m_z)' then matched.  #
	if MATCH_CAT1 == 'coadd':
		mag_axlabel1 = 'MAG_AUTO_'+ str(AXLABEL1)
	if MATCH_CAT2 == 'coadd':
		mag_axlabel2 = 'MAG_AUTO_'+ str(AXLABEL2)

        if CM_T_HDR1 is not None:
                cm_t_axlabel1 = str(CM_T_HDR1[:-2]) + '_' + str(AXLABEL1)
        if CM_T_HDR2 is not None:
                cm_t_axlabel2 = str(CM_T_HDR2[:-2]) + '_' + str(AXLABEL2)

        if CM_T_ERR_HDR1 is not None:
                cm_t_err_axlabel1 = str(CM_T_ERR_HDR1[:-2]) + '_' + str(AXLABEL1)
        if CM_T_ERR_HDR2 is not None:
                cm_t_err_axlabel2 = str(CM_T_ERR_HDR2[:-2]) + '_' + str(AXLABEL2)



	### Define variables ###

	# Get magnitude1 #
	fullmag1 = get_floats_from_string(df=df, hdr=mag_hdr1, filter_name=filter_name)

	# Get magnitude2 #
	fullmag2 = get_floats_from_string(df=df, hdr=mag_hdr2, filter_name=filter_name)



	### Clean the data: removed flags and/or perform quality cuts ###
	if EH_CUTS:
		idx_good = get_good_index_using_quality_cuts(df, full_magnitude1=fullmag1, full_magnitude2=fullmag2, cm_flag_hdr1=CM_FLAGS_HDR1, cm_flag_hdr2=CM_FLAGS_HDR2, flag_hdr1=FLAGS_HDR1, flag_hdr2=FLAGS_HDR2)[0]
	
	if EH_CUTS is False:
                idx_good = get_good_index_using_primary_flags(df=df, full_magnitude1=fullmag1, full_magnitude2=fullmag2, cm_flag_hdr1=CM_FLAGS_HDR1, cm_flag_hdr2=CM_FLAGS_HDR2, flag_hdr1=FLAGS_HDR1, flag_hdr2=FLAGS_HDR2)[0]

	cleanmag1 = get_good_data(df=df, hdr=mag_hdr1, idx_good=idx_good, magnitude=True, filter_name=filter_name)
	cleanmag2 = get_good_data(df=df, hdr=mag_hdr2, idx_good=idx_good, magnitude=True, filter_name=filter_name)


	# Some variables set to None because must pass to plotter() #
	cbar_val, cbar_idx_list, cbar_bins, cbar_axlabel = get_colorbar_value(df=df, cm_t_hdr=CM_T_HDR2, cm_t_err_hdr=CM_T_ERR_HDR2, idx_good=idx_good, clean_magnitude1=cleanmag1, clean_magnitude2=cleanmag2, axlabel=AXLABEL2, inj=INJ2)


	### Define errors ###
	err1, err2 = get_errors(mag_err_hdr1=mag_err_hdr1, mag_err_hdr2=mag_err_hdr2, df=df, filter_name=filter_name, idx_good=idx_good)


	### Write flags to file ###
	if LOG_FLAGS:
		for i in np.arange(0, len(FLAG_HDR_LIST), 2):
			# Bad index #
			temp_idx = handle_flags(df=df, filter_name=f, realization_number=realization_number, flag_hdr1=FLAG_HDR_LIST[i], flag_hdr2=FLAG_HDR_LIST[i+1], full_magnitude1=fullmag1, full_magnitude2=fullmag2, tile_name=tile_name)[1]
		#flag_idx.append(temp_idx)
		FLAG_HDR_LIST.extend(temp_idx)
		if SHOW_PLOT is False and SAVE_PLOT is False:
		# Reset to avoid errors #
			counter_subplot = 1


	### Print out the type() for each flag ###
	if SHOW_FLAG_TYPE:
		get_flag_type(df=df, k=counter_flag_type_printout)
		counter_flag_type_printout += 1


	return cbar_val, cbar_idx_list, cbar_bins, err1, err2, cleanmag1, cleanmag2, idx_good, cbar_axlabel, fullmag1, mag_axlabel1, mag_axlabel2	









def plotter(mag_hdr1, mag_hdr2, cbar_val, error1, error2, filter_name, clean_magnitude1, full_magnitude1, mag_axlabel1, clean_magnitude2, mag_axlabel2, plot_title, realization_number, tile_name, idx_list, bins, cbar_axlabel, plot_name):
	"""Plot a single magnitude versus delta-magnitude plot.

	Args:
            full_magnitude1, full_magnitude2 (numpy.ndarray if directly from `df` OR list of floats if from `get_floats_from_string()`) -- Values read directly from pandas DataFrame via `df[hdr]`; no objects removed using nonzero flag values and no quality cuts performed.
            realization_number (str) -- Allowed values: 0 1 2 None. Refers to Balrog injection and None refers to a one-realization run.
        Returns:
		0
	"""

	### Get labels for vertical and horizontal axes. ###
	'''
   	# Rewrite mag_axlabels. Transform, for example, cm_mag_true to cm_mag_{filter}_true or psf_mag_meas to psf_mag_{filter}_meas #
        if RUN_TYPE is None:
                mag_axlabel1 = str(mag_hdr1[:-2]) + '_' + str(AXLABEL1)
        if RUN_TYPE is not None:
                mag_axlabel1 = str(mag_hdr1) + '_' + str(AXLABEL1)
        mag_axlabel2 = str(mag_hdr2[:-2]) + '_' + str(AXLABEL2)

	# Transform, for example, cm_mag_true to cm_mag_{filter}_true, or psf_mag_meas to psf_mag_{filter}_meas #
        mag_axlabel1 = mag_axlabel1[:-4] + str(filter_name) + '_' + mag_axlabel1[-4:]
        mag_axlabel2 = mag_axlabel2[:-4] + str(filter_name) + '_' + mag_axlabel2[-4:]
	if INJ1:
		mag_axlabel1 = 'inj_' + mag_axlabel1
	if INJ2:
		mag_axlabel2 = 'inj_' + mag_axlabel2

	'''


	### Values to plot for normalized plot ###
	if NORMALIZE:
		# Args needed to call normalize_plot() #
		norm_dm_list, bins, hax_mag_list = normalize_plot_maintain_bin_structure(clean_magnitude1=clean_magnitude1, clean_magnitude2=clean_magnitude2, error1=error1, error2=error2, filter_name=filter_name) 


		PLOT_68P, PLOT_34P_SPLIT = True, True

		if PLOT_1SIG:

			### Plot 1-sigma curve according to error calculation ###
                        plt.axhline(y=1.0, color='red', linestyle='--', linewidth=0.7, label='$1 \sigma_{mag\_meas}$')
                        plt.axhline(y=-1.0, color='red', linestyle='--', linewidth=0.7)


			# Line width for top and sides of 68th percentile bins #
			lwt = 1.1; lws = 0.7
			
			if PLOT_34P_SPLIT:
				### Plot the 68th percentile calculated from np.percentile() ###
				vax_68percentile_list, bins, neg_vax_34percentile, pos_vax_34percentile = get_68percentile_from_normalized_data(norm_dm_list=norm_dm_list, bins=bins, hax_mag_list=hax_mag_list)
				counter_legend1 = 0; color1 = 'cyan'
	
				for b in np.arange(0, len(neg_vax_34percentile)-1):
					# Horizontal bar bounds #
					x_hbound = np.array([bins[b], bins[b+1]])
					x_vbound1 = np.array([bins[b], bins[b]])
					x_vbound2 = np.array([bins[b+1], bins[b+1]])

					if neg_vax_34percentile[b] is not None:
						# Horizontal bar bounds #
						neg_y_hbound = np.array([neg_vax_34percentile[b], neg_vax_34percentile[b]])
						# Vertical bar bounds #
						y_vbound = np.array([neg_vax_34percentile[b], 0])
						# Plot #
						# Plot legend once #
						if counter_legend1 == 0:
							plt.plot(x_hbound, neg_y_hbound, color=color1, linewidth=lwt, label='$\pm P_{34}$')
							counter_legend1 = 1
						if counter_legend1 == 1:
							plt.plot(x_hbound, neg_y_hbound, color=color1)
							plt.plot(x_vbound1, y_vbound, color=color1, linewidth=lws, linestyle=':')
							plt.plot(x_vbound2, y_vbound, color=color1, linewidth=lws, linestyle=':')

					if pos_vax_34percentile[b] is not None:
						# Horizontal bar bounds #
						pos_y_hbound = np.array([pos_vax_34percentile[b], pos_vax_34percentile[b]])
						# Vertical bar bounds #
						y_vbound = np.array([0, pos_vax_34percentile[b]])
						# Plot #
						plt.plot(x_hbound, pos_y_hbound, color=color1, linewidth=lwt)
						plt.plot(x_vbound1, y_vbound, color=color1, linewidth=lws, linestyle=':')
						plt.plot(x_vbound2, y_vbound, color=color1, linewidth=lws, linestyle=':')
				

	
			if PLOT_68P:

				counter_legend2 = 0; color2 = 'fuchsia'

				for b in np.arange(0, len(vax_68percentile_list)-1):

					if vax_68percentile_list[b] is not None:

						# Horizontal bar bounds #
						x_hbound = np.array([bins[b], bins[b+1]])
						y_hbound = np.array([vax_68percentile_list[b], vax_68percentile_list[b]])
						# Vertical bar bounds #
						x_vbound1, x_vbound2 = np.array([bins[b], bins[b]]), np.array([bins[b+1], bins[b+1]])
						y_vbound = np.array([-1*vax_68percentile_list[b], vax_68percentile_list[b]])

						# Plot legend once #
						if counter_legend2 == 0:
							plt.plot(x_hbound, y_hbound, color=color2, label='$P_{68}$', linewidth=lwt)	
							counter_legend2 += 1

						if counter_legend2 == 1:
							# Horizontal bar #
							plt.plot(x_hbound, y_hbound, color=color2, linewidth=lwt)
							plt.plot(x_hbound, -1.0*y_hbound, color=color2, linewidth=lwt)
							# Vertical bar #
							plt.plot(x_vbound1, y_vbound, color=color2, linewidth=lws, linestyle=':')
							plt.plot(x_vbound2, y_vbound, color=color2, linewidth=lws, linestyle=':')
			

		### Values to plot ###
		deltam, hax_mag, bins = normalize_plot(norm_delta_mag_list=norm_dm_list, bins=bins, hax_mag_list=hax_mag_list)
		# Labels and appearance #
                plt.ylabel('('+str(mag_axlabel1) + ' - ' + str(mag_axlabel2)+') / '+ '$\sigma$', fontsize=8)



	### For scatter plot ###
	if NORMALIZE is False:
		# Values to plot #
		deltam = np.array(clean_magnitude1) - np.array(clean_magnitude2)
		if SWAP_HAX:
			hax_mag = clean_magnitude2
		if SWAP_HAX is False:
                        hax_mag = clean_magnitude1
		# Labels and appearance #
		plt.ylabel(str(mag_axlabel1) + ' - ' + str(mag_axlabel2), fontsize=9)
		### 1-sigma curve ###
		if PLOT_1SIG:
			hax, vax, err, bins = bin_and_cut_measured_magnitude_error(error1=error1, error2=error2, clean_magnitude1=clean_magnitude1, clean_magnitude2=clean_magnitude2, filter_name=filter_name)[:4]

			### Remove zeros from x, y, and err (zeros were placeholders for instances in which there were no objects in a particular magnitude bin) ###
			err[:] = [temp for temp in err if temp is not None]
			hax[:] = [temp for temp in hax if temp is not None]
			vax[:] = [temp for temp in vax if temp is not None]
			### Plot 1-sigma curve ###
			plt.plot(hax, np.array(vax) + np.array(err), color='red', linestyle='-', linewidth=0.7, label='$1 \sigma_{mag\_meas}$')
			plt.plot(hax, np.array(vax) - np.array(err), color='red', linestyle='-', linewidth=0.7)


	### Write to log files to record the number of objects plotted and the number of objects within 1sigma ###
	logger(delta_mag=deltam, filter_name=filter_name, clean_magnitude1=clean_magnitude1, full_magnitude1=full_magnitude1, realization_number=realization_number, tile_name=tile_name, bins=bins, hax_mag=hax_mag)


	if PRINTOUTS:
                print 'Plotting ', len(clean_magnitude1), ' objects ... \n'

	### Plot ###
	# One colorbar at a time. This error is caught at beginning of script #
	if CM_T_S2N_COLORBAR is False and BIN_CM_T_S2N is False and CM_T_COLORBAR is False and CM_T_ERR_COLORBAR is False:
		plt.scatter(hax_mag, deltam, color=get_color(filter_name=filter_name)[0], alpha=0.25, s=0.25)

	
	if CM_T_S2N_COLORBAR or CM_T_ERR_COLORBAR or CM_T_COLORBAR:
		'''To plot only the worst (smallest) s2n ratio:
		plt.scatter(np.array(hax_mag)[idx_list[0]], np.array(deltam)[idx_list[0]], color='purple', alpha=1, s=1, label='%1.f'%bins[0]+'<cm_T_s2n<%1.f'%bins[1])
		'''
		plt.scatter(hax_mag, deltam, c=cbar_val, alpha=0.25, s=0.25, norm=matplotlib.colors.LogNorm(), cmap='gist_rainbow')
		plt.colorbar(label=cbar_axlabel)


	if BIN_CM_T_S2N:
		colors = ['green', 'purple', 'cyan', 'orange', 'pink', 'yellow', 'black', 'blue']
		for i in np.arange(0, len(idx_list)):
			plt.scatter(np.array(hax_mag)[idx_list[i]], np.array(deltam)[idx_list[i]], color=colors[i], alpha=0.25, s=0.25, label='%1.f'%bins[i]+'<cm_T_s2n<%1.f'%bins[i+1])


	if HEXBIN:
		if NORMALIZE:
			grid = (100, 1000)
			if PRINTOUTS:
				print ' Normalized hexbin has a large number of grid cells. Will take a moment to plot ... \n'
		if NORMALIZE is False:
			grid = 100
		plt.hexbin(hax_mag, deltam, gridsize=grid, cmap=get_color(filter_name=filter_name)[1], bins='log')
		plt.colorbar(label='log(counts)')


	# Labels and appearance #
	if SWAP_HAX:
		plt.xlabel(str(mag_axlabel2), fontsize=9)
	if SWAP_HAX is False:
                plt.xlabel(str(mag_axlabel1), fontsize=9)
	plt.axhline(y=0.0, color='k', linestyle=':', linewidth=0.5)
        if YLOW is not None and YHIGH is not None:
            plt.ylim([YLOW, YHIGH])


	### Plot legend ###
	if PLOT_1SIG and BIN_CM_T_S2N is False:
		plt.legend(fontsize=8).draggable()
	if BIN_CM_T_S2N:
		# Increase marker size and opacity in legend #
		lgnd = plt.legend(markerscale=4, fontsize=8)
		for l in lgnd.legendHandles:
			l.set_alpha(1)


	if SUBPLOT is False:

		plot_name = plot_name.replace('griz', filter_name)

		### Save plot ###
		if SAVE_PLOT:
			print '-----> Saving plot as: ', plot_name
			plt.savefig(plot_name)

		if SHOW_PLOT:
			plt.title(plot_title)
			plt.show()


        return 0









def subplotter(df, flag_idx, mag_hdr1, mag_hdr2, mag_err_hdr1, mag_err_hdr2, plot_name, plot_title, realization_number, tile_name):
	"""Combine four subplots into a single plot with four panels (2-by-2). Declare variables needed for plotting.

	Args:
            *_hdr (str) -- Headers refer to columns in the matched catalog.
            df (pandas DataFrame)
            plot_name (str) -- Path and name for the plot. Used when save_plot is True and normalize is False.
            realization_number (int) -- Allowed values: 0 1 2 None. Refers to Balrog injection and None refers to a one-realization run.
        Returns:
            flag_idx (list of ints) -- If log_flags is True, will check for all nonzero flag values in `FLAG_HDR_LIST` and `flag_idx` will contain indices that have nonzero flag values. Will be empty if LOG_FLAGS is False.
	"""



	# Counter for flag type() printout #
	counter_flag_type_printout = 0


        ### Create 4-by-4 subplot ###
	counter_subplot = 1
	# Figure size units: inches #
	plt.figure(figsize=(10, 8))


        ### Create one subplot for each griz filter ###
	for f in ALL_FILTERS:

		### Define variables ###
		cbar_val, cbar_idx_list, cbar_bins, err1, err2, cleanmag1, cleanmag2, index_good, cbar_axlabel, fullmag1, mag_axlabel1, mag_axlabel2 = get_plot_variables(filter_name=f, df=df, mag_hdr1=mag_hdr1, mag_hdr2=mag_hdr2, mag_err_hdr1=mag_err_hdr1, mag_err_hdr2=mag_err_hdr2, realization_number=realization_number, tile_name=tile_name, mag_axlabel1=M_AXLABEL1, mag_axlabel2=M_AXLABEL2)



		### Subplot ###
		if SUBPLOT:
			plt.subplot(2, 2, counter_subplot)

		plotter(mag_hdr1=mag_hdr1, mag_hdr2=mag_hdr2, cbar_val=cbar_val, plot_title=plot_title, error1=err1, error2=err2, filter_name=f, full_magnitude1=fullmag1, clean_magnitude1=cleanmag1, clean_magnitude2=cleanmag2, mag_axlabel1=mag_axlabel1, mag_axlabel2=mag_axlabel2, realization_number=realization_number, tile_name=tile_name, idx_list=cbar_idx_list, bins=cbar_bins, cbar_axlabel=cbar_axlabel, plot_name=plot_name)

		counter_subplot += 1


	if SUBPLOT:

		### Show or save the plot once all four subplots have been filled ###
		plt.subplots_adjust(hspace=0.4)
		plt.subplots_adjust(wspace=0.3)
		plt.tight_layout(pad=3, h_pad=2.5)


		### Title ###
		plt.suptitle(plot_title)


		### Save plot ###
		if SAVE_PLOT:
			print '-----> Saving plot as: ', plot_name
			plt.savefig(plot_name)

		### Show plot ###
		if SHOW_PLOT:
			plt.show()

	
	return flag_idx









def get_plot_suptitle(realization_number, tile_name):
	"""Generate plot title.

	Args:
		match_type (str) -- Ex: inj_mof_vs_truth_cat 
		realization_number (str) -- Allowed values: '0' '1' '2' ... 'stacked'.
		tile_name (str)
	Returns:
		title (str) -- Ex: 'Inj MOF Cat & Truth Cat' 
	"""

	title = str(TITLE_PIECE1) + ' & ' + str(TITLE_PIECE2) +'. Tile: ' + str(tile_name) + '. Realization: ' + str(realization_number) + '.'

	if RUN_TYPE == 'ok': 
		title = title + ' Unchanged FOF groups.'
	if RUN_TYPE == 'rerun':
		title = title + ' Changed FOF groups.'

	if NORMALIZE:
		title = 'Normalized. ' + title

	return title









def get_plot_save_name(realization_number, tile_name):
        """Generate name of the plot that will be used in plt.savefig().
	Relies on directory structure: outdir/plots/`BALROG_RUN`/`MATCH_TYPE`/{tile}/{plot_type}/{realization}/ where allowed values for plot_type are: 'normalized' 'scatter'. 

        Args:
                outdir (str) -- Output directory
                realization_number (str) -- Allowed values: '0' '1' '2' 'stacked'
                tile_name (str)
        Returns:
                fn (str) -- The complete filename which includes path.
        """

	### Get filename ###
	if YLOW is None and YHIGH is None:
		# Default scale for the vertical axis (vax) is used #
		ylim = 'defaultvax'
	if YLOW is not None and YHIGH is not None:
		ylim = str(YLOW)+'y'+str(YHIGH)

	if RUN_TYPE is None:	
		endname = str(tile_name) + '_' + str(realization_number) + '_griz_' + str(MATCH_TYPE) + '_' + str(ylim) + '.png'
	if RUN_TYPE is not None:
		endname = str(tile_name) + '_' + str(realization_number) + '_griz_' + str(MATCH_TYPE) + '_' + str(RUN_TYPE) + '_' + str(ylim) + '.png'

	# dm = delta magnitude #	
        if CM_T_S2N_COLORBAR:
		outname = 'm_vs_dm_cm_t_s2n_' + endname
        if CM_T_COLORBAR:
		outname = 'm_vs_dm_cm_t_' + endname
	if CM_T_ERR_COLORBAR:
		outname = 'm_vs_dm_cm_t_err_' + endname
	if HEXBIN:
		outname = 'm_vs_dm_hexbin_' + endname
	if CM_T_S2N_COLORBAR is False and CM_T_COLORBAR is False and CM_T_ERR_COLORBAR is False and HEXBIN is False:
		outname = 'm_vs_dm_' + endname


	# !!!!! User may wish to edit directory structure #
	### Check for directory existence ###
	if RUN_TYPE is not None:
		plot_dir_pref = os.path.join(OUTDIR, 'plots', BALROG_RUN, MATCH_TYPE, tile_name, realization_number, 'fof_analysis') 
	if RUN_TYPE is None:
		plot_dir_pref = os.path.join(OUTDIR, 'plots', BALROG_RUN, MATCH_TYPE, tile_name, realization_number)

	if NORMALIZE:
		plot_dir = os.path.join(plot_dir_pref, 'normalized')
	if NORMALIZE is False:
		plot_dir = os.path.join(plot_dir_pref, 'scatter')
	
	if os.path.isdir(plot_dir) is False:
		if NO_DIR_EXIT:
			sys.exit('Directory ' + str(plot_dir) + ' does not exist. \n Change directory structure in ms_plotter.get_plot_save_name() or set `NO_DIR_MAKE=True`')
		if NO_DIR_MAKE:
			print 'Making directory ', plot_dir, '...\n'
			os.makedirs(plot_dir)


	### Get filename and path ###
        if NORMALIZE:
		fn = os.path.join(plot_dir, 'norm_' + str(outname))
        if NORMALIZE is False:
		fn = os.path.join(plot_dir, outname)


        return fn









def get_coadd_mag_and_mag_err(fn_g, fn_r, fn_i, fn_z, mag_hdr, err_hdr):
	"""Solely for use with coadd catalogs. Creates a list of magnitudes of form '(mag_g, mag_r, mag_i, mag_z)' from four catalogs.

	Args:
		fn -- Filenames. Must be FITS files.
		hdr (str) -- Header for the magnitude. Headers refer to columns in the matched catalog.
	Returns:
		m_griz (list of str) -- Stores magnitude of each filter in form '(mag_g, mag_r, mag_i, mag_z)'
		m_err_griz (list of str) -- Stores error in magnitude of each filter in form '(mag_g, mag_r, mag_i, mag_z)'
	"""

	# Files have not yet been matched, and do not have hdr_1 #
	mag_hdr = mag_hdr[:-2]
	err_hdr = err_hdr[:-2]

	# Open FITS files #
	hdu_g = fits.open(fn_g); hdu_r = fits.open(fn_r); hdu_i = fits.open(fn_i); hdu_z = fits.open(fn_z)
	
	# Read data #
	data_g = hdu_g[1].data; data_r = hdu_r[1].data; data_i = hdu_i[1].data; data_z = hdu_z[1].data

	# Get magnitudes #
	m_g = data_g[mag_hdr]; m_r = data_r[mag_hdr]; m_i = data_i[mag_hdr]; m_z = data_z[mag_hdr]

	# Get magnitude errors #
	err_g = data_g[err_hdr]; err_r = data_r[err_hdr]; err_i = data_i[err_hdr]; err_z = data_z[err_hdr]

	m_griz, m_err_griz = [], []

        for i in np.arange(0, len(m_g)):
                m_griz.append("'("+ str(m_g[i]) + ', ' + str(m_r[i]) + ', ' + str(m_i[i]) + ', ' + str(m_z[i]) + ")'")
		m_err_griz.append("'("+ str(err_g[i])+ ', ' + str(err_r[i])+ ', ' + str(err_i[i]) + ', ' + str(err_z[i]) + ")'")

        return m_griz, m_err_griz









def get_star_mag(df):
	"""Solely for use with star truth catalogs. Computes and creates a list of magnitudes of form '(mag_g, mag_r, mag_i, mag_z)'.

	Args:
		df (pandas DataFram)
	Returns:
		m_griz (list of str) -- Stores magnitudes of each filter in form '(mag_g, mag_r, mag_i, mag_z)'.
	"""

	m_g = df['g_Corr_1']

	m_r = df['g_Corr_1'] - df['gr_Corr_1']

	m_i = df['g_Corr_1'] - df['gr_Corr_1'] - df['ri_Corr_1']

	m_z = df['g_Corr_1'] - df['gr_Corr_1'] - df['ri_Corr_1'] - df['iz_Corr_1']

	m_griz = []

	for i in np.arange(0, len(m_g)):
		m_griz.append("'("+ str(m_g[i]) + ', ' + str(m_r[i]) + ', ' + str(m_i[i]) + ', ' + str(m_z[i]) + ")'")	
        
	return m_griz









def get_catalog(cat_type, inj, realization_number, tile_name, filter_name):
        """Get catalog to analyze.
	
        Args:
                cat_type -- Catalog type. Allowed values: 'gal_truth', 'mof', 'star_truth', 'sof', 'coadd'.
                inj (bool)
                realization_number (str) -- Allowed values: '0' '1' '2' ...
                tile_name -- Different allowed values depending on catalog.
		filter_name (str) -- Only used with coadd catalogs.
        Returns:
                fn -- Filename
        """

        if cat_type == 'gal_truth' and inj:
                fn = os.path.join(BASEPATH, 'balrog_images', realization_number, tile_name, tile_name+'_'+realization_number+'_balrog_truth_cat_gals.fits')
        if cat_type == 'gal_truth' and inj is False:
                sys.exit('No non-injected truth catalog exists.')

        if cat_type == 'star_truth' and inj:
                fn = os.path.join(BASEPATH, 'balrog_images', realization_number, tile_name, tile_name+'_'+realization_number+'_balrog_truth_cat_stars.fits')
        if cat_type == 'star_truth' and inj is False:
		sys.exit('No non-injected truth catalog exists.')

        if cat_type == 'sof' and inj:
                fn = os.path.join(BASEPATH, 'balrog_images', realization_number, tile_name, 'sof', tile_name+'_sof.fits')
        if cat_type == 'sof' and inj is False:
                fn = os.path.join(BASEPATH, tile_name, 'sof', tile_name+'_sof.fits')

        if cat_type == 'mof' and inj:
                fn = os.path.join(BASEPATH, 'balrog_images', realization_number, tile_name, 'mof', tile_name+'_mof.fits')
        if cat_type == 'mof' and inj is False:
                fn = os.path.join(BASEPATH, tile_name, 'mof', tile_name+'_mof.fits')

        if cat_type == 'coadd' and inj:
                fn = os.path.join(BASEPATH, 'balrog_images', realization_number, tile_name, 'coadd', tile_name+'_'+filter_name+'_cat.fits')
        if cat_type == 'coadd' and inj is False:
                fn = os.path.join(BASEPATH, tile_name, 'coadd', tile_name+'_'+filter_name+'_cat.fits')

        return fn









def matcher(realization_number, tile_name, filter_name):
        """Match two catalogs on RA and DEC with a tolerance of 1 arcsecond via STILTS.

        Args:
                outdir (str) -- Path to where matched catalogs are saved.
                realization_number (str or int) -- Currently allowed values: 0 1 2 3 4 5 6 7 8 9 depending on the basepath.
                tile_name (str) -- Currently allowed values: DES0347-5540  DES2329-5622  DES2357-6456 DES2247-4414 depending on the basepath.
        Returns:
                outname (str) -- Name of matched catalog. Headers will have _1 appended for truth catalog and _2 appended for mof catalog.
        """

        ### Get arguments to pass to ms_matcher ###

        # Input catalogs for STILTS #
	if MATCH_CAT1 is not 'coadd':
		in1 = get_catalog(cat_type=MATCH_CAT1, inj=INJ1, realization_number=realization_number, tile_name=tile_name, filter_name=filter_name)
	if MATCH_CAT1 == 'coadd':
		in1 =  get_coadd_matcher_catalog(cat_type=MATCH_CAT1, inj=INJ1, realization_number=realization_number, tile_name=tile_name, mag_hdr=M_HDR1, err_hdr=mag_err_hdr1)

	if MATCH_CAT2 is not 'coadd':
		in2 = get_catalog(cat_type=MATCH_CAT2, inj=INJ2, realization_number=realization_number, tile_name=tile_name, filter_name=filter_name)
	if MATCH_CAT2 == 'coadd':
		in2 =  get_coadd_matcher_catalog(cat_type=MATCH_CAT2, inj=INJ2, realization_number=realization_number, tile_name=tile_name, mag_hdr=M_HDR2, err_hdr=mag_err_hdr2)

        # !!!!! User may wish to edit directory structure. Output catalog name for STILTS #
	### Check for directory existence ###
	match_dir = os.path.join(OUTDIR, 'catalog_compare', BALROG_RUN, MATCH_TYPE)
	if os.path.isdir(match_dir) is False:
		if NO_DIR_EXIT:
			sys.exit('Directory ' + str(match_dir) + ' does not exist. \n Change directory structure in ms_plotter.matcher() or set `NO_DIR_MAKE=True`')
		if NO_DIR_MAKE:
			print 'Making directory ', match_dir, '...\n'
			os.makedirs(match_dir)


        outname_match = os.path.join(match_dir, tile_name+'_'+realization_number+'_'+str(MATCH_TYPE)+'_match1and2.csv')
	outname_1not2 = os.path.join(match_dir, tile_name+'_'+realization_number+'_'+str(MATCH_TYPE)+'_match1not2.csv')
	outname_2not1 = os.path.join(match_dir, tile_name+'_'+realization_number+'_'+str(MATCH_TYPE)+'_match2not1.csv')

	# Overwrite matched catalogs if one already exists? # 
        overwrite = False

	# Check `outname` existence #
	if os.path.isfile(outname_2not1) is False or (os.path.isfile(outname_2not1) and overwrite):

		print '\nMatching ', in1, in2, '...\n'

		### Matching done in ms_matcher. Args: in1, in2, out, RA_HDR1, DEC_HDR1, RA_HDR2, DEC_HDR2, overwrite ###
		# !!!!! Ensure that path to ms_matcher is correct #
		subprocess.call(['/data/des71.a/data/mspletts/balrog_validation_tests/scripts/BalVal/ms_matcher', in1, in2, outname_match, outname_1not2, outname_2not1, RA_HDR1, DEC_HDR1, RA_HDR2, DEC_HDR2])

        return outname_match, outname_1not2, outname_2not1









def fof_matcher(realization_number, tile_name):
        """Get catalogs to analyze. Return FOF-analysed catalogs.

        Args:
		outdir (str) 
                realization_number (str) -- Allowed values: '0' '1' '2' ...
                tile_name -- Different allowed values depending on catalog.
        Returns:
                fn_* -- Filenames
        """


        ### Filenames for input catalogs used in fof_matcher ###
	# FOF #
        fof = os.path.join(BASEPATH, tile_name, 'mof', tile_name+'_fofslist.fits')
        inj_fof = os.path.join(BASEPATH, 'balrog_images', realization_number, tile_name, 'mof', tile_name+'_fofslist.fits')
        # MOF or SOF #
        if MOF:
                mof = os.path.join(BASEPATH, tile_name, 'mof', tile_name+'_mof.fits')
                inj_mof = os.path.join(BASEPATH, 'balrog_images', realization_number, tile_name, 'mof', tile_name+'_mof.fits')
                fof = os.path.join(BASEPATH, tile_name, 'mof', tile_name+'_fofslist.fits')
                inj_fof = os.path.join(BASEPATH, 'balrog_images', realization_number, tile_name, 'mof', tile_name+'_fofslist.fits')
        if SOF:
                mof = os.path.join(BASEPATH, tile_name, 'sof', tile_name+'_sof.fits')
                inj_mof = os.path.join(BASEPATH, 'balrog_images', realization_number, tile_name, 'sof', tile_name+'_sof.fits')
                fof = os.path.join(BASEPATH, tile_name, 'sof', tile_name+'_fofslist.fits')
                inj_fof = os.path.join(BASEPATH, 'balrog_images', realization_number, tile_name, 'sof', tile_name+'_fofslist.fits')
        # Coadds. Using i-band #
        coadd = os.path.join(BASEPATH, tile_name, 'coadd', tile_name+'_i_cat.fits')
        inj_coadd = os.path.join(BASEPATH, 'balrog_images', realization_number, tile_name, 'coadd', tile_name+'_i_cat.fits')


        ### Filenames for outputs of fof_matcher ###
        #TODO relies on certain directory structure. make note of this in README.md #
        outdir = os.path.join(OUTDIR, 'fof_analysis_catalog_compare', BALROG_RUN) 
        # Repeated #
        inj_outdir = os.path.join(outdir, tile_name, realization_number)
        inj_outname = tile_name + '_' + realization_number


        ### Check directory existence. Make directories if not present. ###
        if os.path.isdir(inj_outdir) is False:
                os.makedirs(inj_outdir)
        if os.path.isdir(os.path.join(outdir, tile_name)) is False:
                os.makedirs(os.path.join(outdir, tile_name))


        fofcoadd = os.path.join(outdir, tile_name, tile_name+ '_num-match_fof_coadd.csv')
        fofgroups = os.path.join(outdir, tile_name, tile_name+ 'fofgroups.csv')
	inj_fofcoadd = os.path.join(inj_outdir, inj_outname + '_num-match_inj_fof_inj_coadd.csv')
        inj_fofgroups = os.path.join(inj_outdir, inj_outname + '_inj_fofgroups.csv')
        origfof_injfof = os.path.join(inj_outdir, inj_outname + '_inj_fofgroup_fofgroup_match1and2.csv')
        ok = os.path.join(inj_outdir, inj_outname + '.ok')
        rerun = os.path.join(inj_outdir, inj_outname + '.rerun')
        ok_inj_mof = os.path.join(inj_outdir, inj_outname + '_ok_inj_mof.csv')
        rerun_inj_mof = os.path.join(inj_outdir, inj_outname + '_rerun_inj_mof.csv')
        ok_mof = os.path.join(inj_outdir, inj_outname + '_ok_mof.csv')
        rerun_mof = os.path.join(inj_outdir, inj_outname + '_rerun_mof.csv')

        ok_match = os.path.join(inj_outdir, inj_outname + '_ok_inj_mof_ok_mof_match1and2.csv')
	ok_1not2 = os.path.join(inj_outdir, inj_outname + '_ok_inj_mof_ok_mof_match1not2.csv')
        ok_2not1 = os.path.join(inj_outdir, inj_outname + '_ok_inj_mof_ok_mof_match2not1.csv')

        rerun_match = os.path.join(inj_outdir, inj_outname + '_rerun_inj_mof_rerun_mof_match1and2.csv')
	rerun_1not2 = os.path.join(inj_outdir, inj_outname + '_rerun_inj_mof_rerun_mof_match1not2.csv')
	rerun_2not1 = os.path.join(inj_outdir, inj_outname + '_rerun_inj_mof_rerun_mof_match2not1.csv')

        # Output directory for files made in par.py #
        parpy_outdir = os.path.join(inj_outdir, inj_outname)


        # May need to overwrite if matching was interupted #
        overwrite = False 

        ### Check file existence of last file made in fof_matcher ###
        if os.path.isfile(rerun_match) is False or (os.path.isfile(rerun_match) and overwrite):

                ### Run fof_matcher ###
                subprocess.call(['/data/des71.a/data/mspletts/balrog_validation_tests/scripts/BalVal/fof_matcher', fof, inj_fof, mof, inj_mof, coadd, inj_coadd, parpy_outdir, fofcoadd, fofgroups, inj_fofcoadd, inj_fofgroups, origfof_injfof, ok, rerun, ok_inj_mof, rerun_inj_mof, ok_mof, rerun_mof, ok_match, rerun_match, ok_1not2, rerun_1not2, ok_2not1, ok_2not1])


	if RUN_TYPE == 'ok':
		return ok_match, ok_1not2, ok_2not1
	if RUN_TYPE == 'rerun':
		return rerun_match, rerun_1not2, ok_2not1









def make_plots(mag_hdr1, mag_hdr2, mag_err_hdr1, mag_err_hdr2):
	"""Makes plots.

	Args:
		mag_hdr1, mag_hdr2 (str) -- Headers for magnitude. May be altered, hence passed as parameters.
		mag_err_hdr1, mag_err_hdr2 (str) -- Headers for magnitude error. May be altered, hence passed as parameters.
	Returns:
		0
	"""

	for t in ALL_TILES:

		### For plotting all realizations at once, stacked ###
		if STACK_REALIZATIONS:

			# Filename #
			fn_stack = os.path.join(OUTDIR, t+'_stacked_'+str(MATCH_TYPE)+'_match1and2.csv')
			
			# Check if stacked realization file already exists #
			if os.path.isfile(fn_stack):
				print 'Stacked realization catalog exists. Not overwriting ... \n'
				df1and2 = pd.read_csv(fn_stack)
			
			# Combine all realizations for one tile into a single catalog. Catalogs combined AFTER matching. #
			if os.path.isfile(fn_stack) is False:
				all_fn = []
				for r in ALL_REALIZATIONS:
					if RUN_TYPE is None:
						fn = matcher(realization_number=r, tile_name=t, filter_name=None)[0]
					if RUN_TYPE is not None:
						fn = fof_matcher(realization_number=r, tile_name=t)[0]
					all_fn.append(fn)

				print 'Stacking realizations. ', len(all_fn), 'files ...'
				df1and2 = pd.concat((pd.read_csv(fn) for fn in all_fn))
				print 'Stacking complete ... \n'

				# Save stacked catalog as DataFrame #
				df1and2.to_csv(fn_stack, sep=',')
				print '-----> Saving stacked realization catalog as ', fn_stack

			# Name for plt.savefig() #
			fn = get_plot_save_name(realization_number='stacked_realizations', tile_name=t) 

			# Title for plot #
			title = get_plot_suptitle(realization_number='stacked '+str(len(ALL_REALIZATIONS)), tile_name=t) 


			### Handle star truth catalogs ###
			if MATCH_CAT1 == 'star_truth' or MATCH_CAT2 == 'star_truth':
				print 'Adding new column to matched csv ...'
				star_mag = get_star_mag(df=df1and2)
				# 'mag_a' short for mag_all. New header must be of the form {base}_x where x is a single character because of the way m_axlabel is created from m_hdr #
				df1and2.insert(len(df1and2.columns), 'mag_a', star_mag)

			if MATCH_CAT1 == 'star_truth':
				mag_hdr1 = 'mag_a'
			if MATCH_CAT2 == 'star_truth':
				mag_hdr2 = 'mag_a'

			 ### Handle coadd catalogs. New column has been added with name 'mag_c'. Catalog combined then matched so has suffix (unlike star) #
                        if MATCH_CAT1 == 'coadd':
                                mag_hdr1 = 'mag_c_1'
                                mag_err_hdr1 = 'mag_err_c_1'
                        if MATCH_CAT2 == 'coadd':
                                mag_hdr2 = 'mag_c_2'
                                mag_err_hdr2 = 'mag_err_c_2'


			subplotter(df=df1and2, flag_idx=flag_idx, mag_hdr1=mag_hdr1, mag_hdr2=mag_hdr2, mag_err_hdr1=mag_err_hdr1, mag_err_hdr2=mag_err_hdr2, plot_name=fn, plot_title=title, realization_number='stacked', tile_name=t) 




	### For plotting one realization at a time ###
        if STACK_REALIZATIONS is False:

		print 'Not stacking realizations...'

                for r in ALL_REALIZATIONS:
			
			# Filename #
			if RUN_TYPE is None:
				fn_match, fn_1not2, fn_2not1 = matcher(realization_number=r, tile_name=t, filter_name=None)
			if RUN_TYPE is not None:
				fn_match, fn_1not2, fn_2not1 = fof_matcher(realization_number=r, tile_name=t)
			# DataFrame #
                        df1and2 = pd.read_csv(fn_match)
			df1not2 = pd.read_csv(fn_1not2)
			df2not1 = pd.read_csv(fn_2not1)

			### ####
			if MAKE_REG:
				make_region_file(df_match=df1and2, df_1not2=df1not2, df_2not1=df2not1, realization_number=r, tile_name=t)


                        # Name for plt.savefig() #
                        fn = get_plot_save_name(realization_number=r, tile_name=t)

                        # Title for plot #
                        title = get_plot_suptitle(realization_number=r, tile_name=t) 


			### Handle star truth catalogs ###
                        if MATCH_CAT1 == 'star_truth' or MATCH_CAT2 == 'star_truth':
                                print 'Adding new column to matched csv ... \n'
                                star_mag = get_star_mag(df=df1and2)
                                df1and2.insert(len(df1and2.columns), 'mag_a', star_mag)

			# Star truth catalogs matched then combined # #FIXME rename to mag_s for mag_star 
                        if MATCH_CAT1 == 'star_truth':
                                mag_hdr1 = 'mag_a'
                        if MATCH_CAT2 == 'star_truth':
                                mag_hdr2 = 'mag_a'

			
			### Handle coadd catalogs. New column has been added with name 'mag_c'. Catalog combined then matched so has suffix (unlike star) #
			if MATCH_CAT1 == 'coadd':
				mag_hdr1 = 'mag_c_1'
				mag_err_hdr1 = 'mag_err_c_1'
			if MATCH_CAT2 == 'coadd':
				mag_hdr2 = 'mag_c_2'
				mag_err_hdr2 = 'mag_err_c_2'


			subplotter(df=df1and2, flag_idx=flag_idx, mag_hdr1=mag_hdr1, mag_hdr2=mag_hdr2, mag_err_hdr1=mag_err_hdr1, mag_err_hdr2=mag_err_hdr2, plot_name=fn, plot_title=title, realization_number=r, tile_name=t) 


	return 0







def get_coadd_matcher_catalog(cat_type, inj, realization_number, mag_hdr, err_hdr, tile_name):
	"""Make FITS file that includes a column of form '(m_g, m_r, m_i, m_z)' where m is magnitude. Column will be added to '..._i_cat.fits'. This will be used in matcher().

	Args:
		cat_type (str) -- Catalog type. Allowed values: mof, sof, coadd, gal_truth, star_truth.
		inj (bool) -- Is the catalog (`cat_type`) injected?
		realization_number (str)
		mag_hdr (str) -- Header for magnitude. Headers refer to columns in the matched catalog.
		err_hdr (str) -- Header for error. Headers refer to columns in the matched catalog.
		tile_name (str)
	Returns:
		fn (str) -- Filename. Is a FITS file.
	"""

	fn_new = os.path.join(OUTDIR, str(tile_name) + '_i_cat.fits')

	# Check if new coadd catalog has already been created #
	if os.path.isfile(fn_new):
		print 'New coadd catalog already exists ...\n'


	if os.path.isfile(fn_new) is False:	
		print 'Adding a column to i-band coadd catalog. Will take a moment ...\n'

		# Get list of filenames #
		fn_griz = []
		for f in ALL_FILTERS:
			fn_griz.append(get_catalog(cat_type=cat_type, inj=inj, realization_number=realization, tile_name=tile_name, filter_name=f))
		fn_g, fn_r, fn_i, fn_z = fn_griz

		# Get coadd magnitude (mag_c) and magnitude error to be of form '(m_g, m_r, m_i, m_z)'. Recall that this is a string #
		mag_c, mag_err_c = get_coadd_mag_and_mag_err(fn_g=fn_g, fn_r=fn_r, fn_i=fn_i, fn_z=fn_z, mag_hdr=mag_hdr, err_hdr=err_hdr)
 
	       # Create new table #
		mag_c = Column(mag_c, name='mag_c')
		mag_err_c = Column(mag_err_c, name='mag_err_c')

		# Add new table to i-band coadd catalog #
		table = Table.read(fn_i)
		table.add_column(mag_c, index=0)
       		table.add_column(mag_err_c, index=1)
 
		# Save new table as FITS #
		table.write(fn_new)

	return fn_new 









#TODO accept idx_good as input param? Will only be used for df_match. 
def make_region_file(df_match, df_1not2, df_2not1, realization_number, tile_name):
	"""Make DS9 region files for catalogs matched via join=1and2, join=1not2, and 2not1.

	Args:
		df_* (pandas DataFrame)
	Returns:
		fn_* -- Filenames
	"""

	### Get filenames and open files ###
	fn_match, fn_1not2, fn_2not1 = get_reg_names(tile_name=tile_name, realization_number=realization_number)

	overwrite = False
	
	if os.path.isfile(fn_2not1) and overwrite is False:
		print 'Region files already exist. Not overwriting ...'
	

	if os.path.isfile(fn_2not1) is False or (os.path.isfile(fn_2not1) and overwrite):
		fd_match = open(fn_match, 'w'); fd_1not2 = open(fn_1not2, 'w'); fd_2not1 = open(fn_2not1, 'w')
		# Write coordinate system #
		fd_match.write('J20000 \n'); fd_1not2.write('J20000 \n'), fd_2not1.write('J20000 \n')

		# Handle matched catalog #
		if RUN_TYPE is None:
			ra1 = RA_HDR1 + str('_1'); dec1 = DEC_HDR1 + str('_1')
			ra2 = RA_HDR2 + str('_2'); dec2 = DEC_HDR2 + str('_2')
		if RUN_TYPE is not None:
			# MOF or SOF catalogs #
			ra1 = 'ra'; dec1 = 'dec'
			ra2 = 'ra_2'; dec2 = 'dec_2' 

		### Get position. Arbitrarily using MATCH_CAT1 for RA and DEC ###
		ra_match, dec_match = df_match[ra2], df_match[dec2] 
		ra_1not2, dec_1not2 = df_1not2[ra1], df_1not2[dec1]
		ra_2not1, dec_2not1 = df_2not1[ra2], df_2not1[dec2]

		### Write to region file for matched catalog. Units are arcseconds. ###
		# Coadds allow for elliptical regions #
		if MATCH_CAT1 == 'coadd' or MATCH_CAT2 == 'coadd':
			### Get semimajor and semiminor axes (a and b, respectively). Coadds have these values. ###
			a_match, b_match = df_match[MAJOR_AX_HDR1], df_match[MINOR_AX_HDR1]
			a_1not2, b_1not2 = df_1not2[MAJOR_AX_HDR1], df_1not2[MINOR_AX_HDR1]
			a_2not1, b_2not2 = df_2not1[MAJOR_AX_HDR2], df_2not1[MINOR_AX_HDR2]

			for i in np.arange(0, len(ra_match)):
				fd_match.write('ellipse ' + str(ra_match[i]) + ' ' + str(dec_match[i]) + ' ' + str(a_match[i]) + '" ' + str(b_match[i]) + '" #color=green width=3 \n')
			for i in np.arange(0, len(ra_1not2)):
				fd_1not2.write('ellipse ' + str(ra_1not2[i]) + ' ' + str(dec_1not2[i]) + ' ' + str(a_1not2[i]) + '" ' + str(b_1not2[i]) + '" #color=yellow width=3 \n')
			for i in np.arange(0, len(ra_2not1)):
				fd_2not1.write('ellipse ' + str(ra_2not1[i]) + ' ' + str(dec_2not1[i]) + ' ' + str(a_2not1[i]) + '" ' + str(b_2not1[i]) + '" #color=blue width=3 \n')


		# Non-coadd catalogs allow for circular regions #
		if MATCH_CAT1 != 'coadd' and MATCH_CAT2 != 'coadd':
			size_sq_match = df_match[CM_T_HDR1]
			size_sq_1not2 = df_1not2[CM_T_HDR1]
			size_sq_2not1 = df_2not1[CM_T_HDR2]

			# Use a typical radius of 2 arcsec? #

			for i in np.arange(0, len(ra_match)):
				if size_sq_match[i] > 0:# and np.isnan(size_sq_match[i]) is False:
					fd_match.write('circle ' + str(ra_match[i]) + ' ' + str(dec_match[i]) + ' ' + str(size_sq_match[i]**0.5) + '" #color=green width=3 \n')
			for i in np.arange(0, len(ra_1not2)):
				if size_sq_1not2[i] > 0: # and np.isnan(size_sq_1not2[i]) is False:
					fd_1not2.write('circle ' + str(ra_1not2[i]) + ' ' + str(dec_1not2[i]) + ' ' + str(size_sq_1not2[i]**0.5) + '" #color=yellow width=3 \n')
			for i in np.arange(0, len(ra_2not1)):
				if size_sq_2not1[i] > 0: # and np.isnan(size_sq_2not1[i]) is False:
					fd_2not1.write('circle ' + str(ra_2not1[i]) + ' ' + str(dec_2not1[i]) + ' ' + str(size_sq_2not1[i]**0.5) + '" #color=blue width=3 \n')

		# Close files #
		fd_match.close(); fd_1not2.close(); fd_2not1.close()

	print '-----> Saving region files as: ', fn_match
	print ' -----> ', fn_1not2
	print ' ----->', fn_2not1

	return 0 








################################################################### Run script. 0 returned when complete. ###################################################################

### !!!!! Run once. Log files are closed once 0 is returned. ###
print make_plots(mag_hdr1=M_HDR1, mag_hdr2=M_HDR2, mag_err_hdr1=M_ERR_HDR1, mag_err_hdr2=M_ERR_HDR2)



### Loop over vertical axis limits. Suggestions: for normalized plot with star truth catalog use y=[3, 10], for normalized plot with galaxy truth catalog use y=[3, 20]. For non-normalized plot with star truth catalog or galaxy truth catalog use y=[0.5, None]. ###
# !!!!! Loop over vertical axis limits? #
YLOOP = False 

if YLOOP:
	for y in [0.5, None]: 
		if y is None:
			YLOW, YHIGH = None, None
		if y is not None:
			YLOW, YHIGH = -1.0*y, y

		# Must pass constants as parameters here because used for plot name #
		print make_plots(mag_hdr1=M_HDR1, mag_hdr2=M_HDR2, mag_err_hdr1=M_ERR_HDR1, mag_err_hdr2=M_ERR_HDR2)



### !!!!! Loop over possible colorbars? ###
CBAR_LOOP = False

if CBAR_LOOP:
	# Possible combinations for colorbars. Only one True allowed at a time and NORMALIZE and HEXBIN must both be False. #
	cbar_bools_list =[[True, False, False, False, False], [False, True, False, False, False], [False, False, True, False, False]]
	for cbar_bools in cbar_bools_list:
		# Reset constants #
		CM_T_COLORBAR, CM_T_S2N_COLORBAR, CM_T_ERR_COLORBAR, NORMALIZE, HEXBIN = cbar_bools
		make_plots(mag_hdr1=M_HDR1, mag_hdr2=M_HDR2, mag_err_hdr1=M_ERR_HDR1, mag_err_hdr2=M_ERR_HDR2)



RUN_TYPE_LOOP = False 

if RUN_TYPE_LOOP:
	run_type_list = [None, 'ok', 'rerun']
	for run_type in run_type_list:
		RUN_TYPE = run_type
		make_plots(mag_hdr1=M_HDR1, mag_hdr2=M_HDR2, mag_err_hdr1=M_ERR_HDR1, mag_err_hdr2=M_ERR_HDR2)


### Close log files ###
FD_1SIG.close()
FD_FLAG.close()
FD_MAG_BINS.close()
FD_NOP.close()

