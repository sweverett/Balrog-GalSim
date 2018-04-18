"""Creates various plots for Balrog validation testing.


To run: $python ms_plotter.py base_path_to_catalogs output_directory realization tile
Example: $python ms_plotter.py /data/des71.a/data/kuropat/des2247-4414_sof/y3v02/ /Users/mspletts/BalValPlots/ all DES2247-4414

Relies on ms_matcher. User may need to replace `/data/des71.a/data/mspletts/balrog_validation_tests/scripts/ms_matcher` with the correct path to ms_matcher.

Plot attributes are specified with constants (many of them booleans) at the top of this script.
Constants that the user may wish to change are indicated by: # !!!!! {description and/or warnings} #. For example, user may wish to set `PRINTOUTS = False` or comment out `notice`.

# Comments are ABOVE the code they correspond to (with the exception of #FIXMEs and #TODOs). #

Megan Splettstoesser mspletts@fnal.gov"""


import matplotlib.pyplot as plt
import matplotlib
import numpy as np
import os
import pandas as pd
import subprocess
import sys




### Command line args ###
basepath, outdir, realization, tilename = sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4]
# Catch error from inadequate number of command line args #
if len(sys.argv) != 5:
        sys.exit("Args: basepath (location of catalogs), output directory, realization (can be 'all'), tile (can be 'all') \n")


ALL_FILTERS = [ 'g', 'r', 'i', 'z' ]

# !!!!! Number of realizations depends on the tile # 
ALL_REALIZATIONS = [ '0', '1', '2', '3', '4', '5', '6', '7', '8', '9' ]
#ALL_REALIZATIONS = [ '0', '1', '2' ]

# !!!!! Available tiles depends on run #
ALL_TILES = [ 'DES0347-5540', 'DES2329-5622', 'DES2357-6456' ]


### Make command line args lists ###
if tilename != 'all':
        ALL_TILES = [tilename]
if realization != 'all':
        ALL_REALIZATIONS = [realization]




################################################################### Specify plot and catalog attributes ###################################################################

### Colorbar ###
# !!!!! Add colorbar according to one of the following (only one can be True at a time). If all are False a scatter-plot is made. Colorbars cannot be used if NORMALIZE is True. #
HEXBIN = False
CM_T_S2N_COLORBAR = False
CM_T_ERR_COLORBAR = False
CM_T_COLORBAR = False
BIN_CM_T_S2N = False
# Normalizes plot to 1-sigma magnitude error. If NORMALIZE is True, PLOT_1SIG must be True else errors will not be computed and normalization cannot be performed #
NORMALIZE = True

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
# !!!!! Allowed values: sof, mof, star_truth, gal_truth, coadd (#TODO under construction). Both can be 'sof' and both can be 'mof' if INJ1 and INJ2 are different. Note that truth catalogs always have INJ=True. #
MATCH_CAT1, MATCH_CAT2 = 'gal_truth', 'sof'
# !!!!! Booleans. Examine injected catalogs? #
INJ1, INJ2 = True, True
# !!!!! Must be used with realization=all at command line #
STACK_REALIZATIONS = False


### Miscellaneous ###
# Print progress? #
PRINTOUTS = True
# Only refers to printouts within get_floats_from_string() and get_matrix_diagonal_element() # 
PRINTOUTS_MINOR = False

# Not currently in use or under constructrion #
LOG_FLAGS = False
PLOT_FLAGGED_OBJS = False
SHOW_FLAG_TYPE = False
SUBPLOT = True #TODO




# Catch errors from plot attributes #
if ((CM_T_S2N_COLORBAR and CM_T_ERR_COLORBAR) or (CM_T_S2N_COLORBAR and HEXBIN) or (CM_T_ERR_COLORBAR and HEXBIN) or (CM_T_COLORBAR and CM_T_ERR_COLORBAR) or (CM_T_COLORBAR and CM_T_S2N_COLORBAR) or (CM_T_COLORBAR and HEXBIN)) or (NORMALIZE and PLOT_1SIG is False) or (YLOW is not None and YHIGH is not None and YHIGH == YLOW) or (NORMALIZE and (CM_T_S2N_COLORBAR or CM_T_ERR_COLORBAR or CM_T_COLORBAR)) or (STACK_REALIZATIONS and realization != 'all'):
	sys.exit('ERROR: Only of of the following may be True at a time: CM_T, CM_T_S2N_COLORBAR, CM_T_ERR_COLORBAR, HEXBIN. Otherwise, colorbar will be overwritten. \nERROR: If NORMALIZE is True so must be PLOT_1SIG. \nERROR: YHIGH cannot be equal to YLOW. \n NORMALIZE must be False if any of: CM_T_S2N_COLORBAR CM_T_ERR_COLORBAR CM_T_S2N_COLORBAR are True.\nERROR: STACK_REALIZATIONS is True must be used with realization = all. \n')


# !!!!! Check that plots will not be overwritten, etc #
NOTICE = raw_input(' \n !! CHECK BEFORE RUNNING !! \n Save plot(s) -- ' + str(SAVE_PLOT) + '\n Showing plot(s) -- ' + str(SHOW_PLOT) + '\n Normalize plot(s) -- ' + str(NORMALIZE) + '\n Hexbin -- ' + str(HEXBIN) + '\n cm_T colorbar -- ' + str(CM_T_COLORBAR) + '\n cm_T_err colorbar -- ' + str(CM_T_ERR_COLORBAR) + '\n cm_T_s2n -- ' + str(CM_T_S2N_COLORBAR) + '\n Plot limits -- ' + str(YLOW) + ', ' + str(YHIGH) + '\n Plotting 1-sigma curve -- ' + str(PLOT_1SIG) +'\n Plotting flagged objects -- ' + str(PLOT_FLAGGED_OBJS) + '\n Print flags and flag types -- ' + str(SHOW_FLAG_TYPE) + '\n Logging flags -- ' + str(LOG_FLAGS) + '\n --> Press enter to proceed, control+c to stop...\n')









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

	strings = df[hdr]

	list_a = []

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

	matrices = df[hdr]

	list_aa = []

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









def get_good_index_using_primary_flags(df, cm_flag_hdr1, cm_flag_hdr2, flag_hdr1, flag_hdr2, full_magnitude1, full_magnitude2):
	"""Get indices of objects without flags as indicated by the headers 'flags' and 'cm_flags'. Also get indices of objects with  magnitudes not equal to +/- 99, +/- 9999, and 37.5. Store the bad indices as well (if PLOT_FLAGGED_OBJS is True).

	Args:
		df (pandas DataFrame)
		cm_flag_hdr1, cm_flag_hdr2, flag_hdr1, flag_hdr2 (str) -- Headers refer to column names in the matched catalog.
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


	idx_good, idx_bad = [], []


	for i in np.arange(0, len(full_magnitude1)):

		# Keep objects with these indices #
		if abs(full_magnitude1[i]) != 9999.0 and abs(full_magnitude2[i]) != 9999.0 and full_magnitude1[i] != 37.5 and full_magnitude2[i] != 37.5 and full_magnitude1[i] != 99.0 and full_magnitude2[i] != 99 and flag1[i] == 0 and flag2[i] == 0 and cm_flag1[i] == 0 and cm_flag2[i] == 0:
			idx_good.append(i)

		if PLOT_FLAGGED_OBJS:
			# Get rid of objects with these indices #
			if abs(full_magnitude1[i]) != 9999.0 and abs(full_magnitude2[i]) != 9999.0 and full_magnitude1[i] != 37.5 and full_magnitude2[i] != 37.5 and full_magnitude1[i] != 99.0 and full_magnitude2[i] != 99 and flag1[i] != 0 or flag2[i] != 0 or cm_flag1[i] != 0 or cm_flag2[i] != 0:
				idx_bad.append(i)


        if PRINTOUTS:
                print 'Eliminated ', len(full_magnitude1) - len(idx_good), ' objects with magnitudes equal to +/- 9999, +/- 99, and 37.5 and objects with nonzero flags for: ', flag_hdr1, ', ', flag_hdr2, ', ', cm_flag_hdr1, ', ', cm_flag_hdr2, ' ... \n'


	return idx_good, idx_bad	









def get_good_index_using_quality_cuts(cm_flag_hdr1, cm_flag_hdr2, df, cm_t_hdr1, cm_t_hdr2, cm_t_err_hdr1, cm_t_err_hdr2, cm_s2n_r_hdr1, cm_s2n_r_hdr2, flag_hdr1, flag_hdr2, full_magnitude1, full_magnitude2, psfrec_t_hdr1, psfrec_t_hdr2, cm_t_err_axlabel1, cm_t_err_axlabel2):
	"""Get indices of objects that satisfy quality cuts introduced by Eric Huff. Also get indices of objects without flags as indicated by the headers 'flags' and 'cm_flags'. Also get indices of objects with  magnitudes not equal to +/- 99, +/- 9999, and 37.5. Store the bad indices as well (if PLOT_FLAGGED_OBJS is True).

	Args:
		df (pandas DataFrame)
	        *_hdr (str) -- Headers refer to column names in the matched catalog.
		cm_t_hdr1, cm_t_hdr2 (str) -- Header for size squared of object.
		cm_t_err_hdr1, cm_t_err_hdr2 (str) -- Header for error on the size squared of object.
		cm_s2n_hdr1, cm_s2n_hdr2 (str) -- Header for signal to noise of object.
		psfrec_t_hdr1, psfrec_t_hdr2 (str) -- Header for PSF size.
		full_magnitude1, full_magnitude2 (list of floats) -- Values read directly from pandas DataFrame or passed through `get_floats_from_string()`; no flags removed. 		
		cm_t_err_axlabel1, cm_t_err_axlabel2 (str) -- Allowed values: 'true' or 'meas'. Value is set in classes and should not be manually set.
        Returns:
		idx_good (list of ints) -- Indices of objects without flags and objects which met criteria for quality cuts.
		idx_bad (list of ints) -- Is empty if PLOT_FLAGGED_OBJS is False.
	"""

	if 'true' in cm_t_axlabel1 and 'true' in cm_t_axlabel2:
		sys.exit('ERROR. Cuts should be performed on measured catalog, not truth catalog.')

	# Ignore truth catalog #
	if 'true' in cm_t_axlabel1 and 'meas' in cm_t_axlabel2:
		cm_t_hdr1, cm_t_err_hdr1, cm_s2n_r_hdr1, psfrec_t_hdr1 = cm_t_hdr2, cm_t_err_hdr2, cm_s2n_r_hdr2, psfrec_t_hdr2

	if 'true' in cm_t_axlabel2 and 'meas' in cm_t_axlabel1:
		cm_t_hdr2, cm_t_err_hdr2, cm_s2n_r_hdr2, psfrec_t_hdr2 = cm_t_hdr1, cm_t_err_hdr1, cm_s2n_r_hdr1, psfrec_t_hdr1

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
	cm_t1, cm_t2 = df[cm_t_hdr1], df[cm_t_hdr2]
	# Size error #
	cm_t_err1, cm_t_err2 = df[cm_t_err_hdr1], df[cm_t_err_hdr2]
	# Signal to noise #
	cm_s2n_r1, cm_s2n2_r = df[cm_s2n_r_hdr1], df[cm_s2n_r_hdr2]
	# PSF size #
	psfrec_t1, psfrec_t2 = df[psfrec_t_hdr1], df[psfrec_t_hdr2]


	### Search for flags, etc ###
	counter1, counter2 = 0, 0


	for i in np.arange(0, len(full_magnitude1)):

	# Keep these objects #
		if abs(full_magnitude1[i]) != 9999.0 and abs(full_magnitude2[i]) != 9999.0 and full_magnitude1[i] != 37.5 and full_magnitude2[i] != 37.5 and full_magnitude1[i] != 99.0 and full_magnitude2[i] != 99 and flag1[i] == 0 and flag2[i] == 0 and cm_flag1[i] == 0 and cm_flag2[i] == 0:

			counter1 += 1

			# From Eric Huff: keep (catalog['cm_s2n_r'] > 10) & (catalog['cm_T']/catalog['cm_T_err'] > .5) & (catalog['cm_T']/catalog['psfrec_T'] > 0.5) via https://github.com/sweverett/Balrog-GalSim/blob/master/plots/balrog_catalog_tests.py #
			if cm_s2n_r1[i] > 10 and cm_s2n2_r[i] > 10 and cm_t1[i]/cm_t_err1[i] > 0.5 and cm_t2[i]/cm_t_err2[i] > 0.5 and cm_t1[i]/psfrec_t1[i] > 0.5 and cm_t2[i]/psfrec_t2[i] > 0.5:
				idx_good.append(i)
				counter2 += 1


		if PLOT_FLAGGED_OBJS:
			# Get rid of these objects #
			if abs(full_magnitude1[i]) != 9999.0 and abs(full_magnitude2[i]) != 9999.0 and full_magnitude1[i] != 37.5 and full_magnitude2[i] != 37.5 and full_magnitude1[i] != 99.0 and full_magnitude2[i] != 99 and flag1[i] != 0 or flag2[i] != 0 or cm_flag1[i] != 0 or cm_flag2[i] != 0:
				if cm_s2n_r1[i] < 10 or cm_s2n2_r[i] < 10 or cm_t1[i]/cm_t_err1[i] < 0.5 or cm_t2[i]/cm_t_err2[i] < 0.5 or cm_t1[i]/psfrec_t1[i] < 0.5 or cm_t2[i]/psfrec_t2[i] < 0.5:
					idx_bad.append(i)

	if PRINTOUTS:
		print 'Eliminated objects with magnitudes equal to +/- 9999, +/- 99, and 37.5 and objects with nonzero flags for: ', flag_hdr1, ', ', flag_hdr2, ', ', cm_flag_hdr1, ', ', cm_flag_hdr2, ' ...'
		print ' Eliminated objects with signal-to-noise ratio < 10 ...'
		print ' Eliminated objects with cm_T/cm_T_err < 0.5 ...'
		print ' Eliminated objects with cm_T/psfrec_T < 0.5 ...'
		print ' Number of objects with primary flags and outlier magnitudes removed: ', counter1
		print " Number of objects when Eric Huff's additional constraints are used: ", counter2, '\n'

	
	return idx_good, idx_bad









def handle_flags(df, fd_flag, filter_name, flag_hdr1, flag_hdr2, full_magnitude1, full_magnitude2, realization_number, run_type, tile_name):
	"""Examine a particular flag and write to log file. Can also be used to check all flags in a list of flags.

	Args:
		df (pandas DataFrame)
		fd_flag (file descriptor) -- File to write flags to with headers TILE, FILTER, TYPE, REALIZATION, FLAG1_HEADER, FLAG2_HEADER, FLAG1_VALUE, FLAG2_VALUE, MAG1, MAG2.
		filter_name (str) -- Allowed values: 'g' 'r' 'i' 'z'.
		flag_hdr1, flag_hdr2 (str) -- Headers refer to column names in the matched catalog.
		full_magnitude1, full_magnitude2 (numpy.ndarray if directly from `df[hdr]` OR list of floats if from `get_floats_from_string()`) -- Values read directly from pandas DataFrame via `df[hdr]`; no objects removed using nonzero flag values and no quality cuts performed.
		realization_number (int) -- Allowed values: 0 1 2 None. Refers to Balrog injection and None refers to a one-realization run.
		run_type (str) -- Allowed values: 'ok' 'rerun' None. 'ok' refers to FOF groups unchanged after Balrog injection. 'rerun' refers to changed FOF groups.
		tile_name (str) -- 
        Returns:
		idx_good (list of ints) -- Indices of objects with flags values of zero.
		idx_bad (list of ints) -- Indices of objects with nonzero flag values.
	"""

	idx_good, idx_bad, counter_idx_bad = [], [], 0

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
				fd_flag.write(str(tile_name) + '\t' + str(filter_name) + '\t' + str(run_type) + '\t' + str(realization_number) + '\t' + str(flag_hdr1) + '\t' + str(flag_hdr2) + '\t' + str(flag1[i]) + '\t' + str(flag2[i]) + '\t' + str(full_magnitude1[i]) + '\t' + str(full_magnitude2[i]) +'\n')



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

	counter_neg, error, flux, fluxcov = 0, [], [], []

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

	'''
	full_cm_t = df[cm_t_hdr]
	full_cm_t_err = df[cm_t_err_hdr]

	cm_t, cm_t_err = [], []

	# Remove objects with flags #
	for i in idx_good:
		cm_t.append(full_cm_t[i])
		cm_t_err.append(full_cm_t_err[i])
	'''

	cm_t = get_good_data(df=df, hdr=cm_t_hdr, idx_good=idx_good, magnitude=False, filter_name=None)
	cm_t_err = get_good_data(df=df, hdr=cm_t_err_hdr, idx_good=idx_good, magnitude=False, filter_name=None)

	# cm_T signal to noise (s2n) #
	s2n = []
	for i in np.arange(0, len(cm_t)):
		s2n_temp = abs(cm_t[i] / cm_t_err[i])
		s2n.append(s2n_temp)

	# Bin signal to noise #
	# Bins suggested by Spencer Everett #
	bins = [0, 1, 9, 20, max(s2n)]
	if PRINTOUTS:
		print 'Binning cm_T_s1n with bins: ', bins, ' and headers/axlabels:', cm_t_hdr, ', ', cm_t_err_hdr, '...'
		print ' Min and max absolute value of cm_T signal-to-noise: ', min(s2n), ' and ', max(s2n), '...'
	
	binned_s2n, binned_hax_mag, binned_vax_mag, idx = [], [], [], []

	for j in np.arange(0, len(bins)-1):
		idx_temp = []
		for i in np.arange(0, len(s2n)):
			if s2n[i] > bins[j] and s2n[i] < bins[j+1]:
				idx_temp.append(i)	
		if PRINTOUTS:
			print ' For cm_T_s2n, number of objects in bin ', bins[j], '-', bins[j+1], ': ', len(idx_temp)
		idx.append(idx_temp)
		#idx_temp = np.where(s2n > bins[j] & (s2n < bins[j+1]))

	if PRINTOUTS:
		print ' '


	return idx, bins, s2n









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

	vax_68percentile = []

	PLOT_HIST = False

	# Loop through bins (b) #
	for b in np.arange(0, len(norm_dm_list)):

		if norm_dm_list[b] == 0:
			vax_68percentile.append(0)

		if norm_dm_list[b] != 0:

			# Plot histogram to see distrubtion of data (data is not normally distributed) #
			if PLOT_HIST:
				plt.figure()
				norm_dm = [abs(elmt) for elmt in norm_dm_list[b]]
				plt.hist(norm_dm)
				plt.title('Bin LHS: ' + str(bins[b]))
				plt.xlabel(r'$\Delta M$')
				plt.axvline(x=0, color='black', linestyle=':', linewidth=0.5)	
				plt.show()
		
			# Values in current bin (icb) #
			vax_mag_list_icb = norm_dm_list[b]
			# Take absolute value of each point in bin #
			vax_mag_list_icb = [abs(elmt) for elmt in vax_mag_list_icb]
			# Percentile sorts the data #
			vax_68p_icb = np.percentile(vax_mag_list_icb, 68, interpolation='lower')
			vax_68percentile.append(vax_68p_icb)
			# Check the percentile because interpolation='lower' was used #
			num = 0
			for j in np.arange(0, len(norm_dm_list[b])):
				if abs(norm_dm_list[b][j]) <= vax_68p_icb:
					num += 1
			if PRINTOUTS:
				print 'Number of objects within 68 percentile via np.percentile(interpolation=lower): ', float(num)/len(norm_dm_list[b]), '...\n'


	return vax_68percentile, bins









def bin_and_cut_measured_magnitude_error_dm3(clean_magnitude1, clean_magnitude2, error1, error2, swap_hax, axlabel1, axlabel2, fd_mag_bins):
        """Remove error values corresponding to objects where |Delta-M| > 3.

        Args:
                clean_magnitude1, clean_magnitude2 (list of floats) -- Objects with flag values of zero and/or quality cuts performed.
                error1, error2 (list of floats) -- 1 and 2 refer to the matched catalogs. 
		swap_hax (bool) -- "swap horizontal axis". If False plots magnitude1 on the x-axis. If True plots magnitude2 on the x-axis.
		fd_mag_bins (file director) -- Log file to record information about the bins used in error calculation.
        Returns:
                b_hax_mag_median (list of floats) -- List of medians of the horizontal axis magnitude in each bin.
                b_vax_mag_median (list of floats) -- List of medians of the vertical axis magnitude in each bin. Vertical axis is computed via clean_magnitude1 - clean_magnitude2.
                b_err_median (list of floats) -- Median of the error in each bin.
                bins (list of floats) -- Bins used. Binned according to horizontal axis.
        """

	### !!!!! Comment this block out if errors are to be computed using both catalogs regardless of origin (measured catalog or truth catalog) ###
	if 'meas' in axlabel1 and 'meas' not in axlabel2:
		error2 = np.zeros(len(error1))
		if PRINTOUTS:
			print 'Using measured catalog (catalog1) for error calculation ... '

	if 'meas' in axlabel2 and 'meas' not in axlabel1:
		error1 =np.zeros(len(error2))
		if PRINTOUTS:
			print 'Using measured catalog (catalog2) for error calculation ... '

	if 'meas' in axlabel1 and 'meas' in axlabel2:
		if PRINTOUTS:
			print 'Using measured catalog (catalog1 AND catalog2) for error calculation ... '

	if 'true' in axlabel1 and 'true' in axlabel2:
		sys.exit('Erros are to be computed using the measured catalog(s), not the truth catalog(s).')


	### Define bins ###
        step = 0.5
        # Find the absolute min and max of the magnitudes in the matched catalog #
        limlow1, limlow2 = min(clean_magnitude1), min(clean_magnitude2)
        limhigh1, limhigh2 = max(clean_magnitude1), max(clean_magnitude2)
        limlow, limhigh = min([limlow1, limlow2]), max([limhigh1, limhigh2])
	# Define bins limits by ints #
	limlow = int(limlow)

	# Introduce magnitude cutoff to tame errors #
	if 'gal' in MATCH_CAT1 or 'gal' in MATCH_CAT2:
		limhigh = 26
	if 'star' in MATCH_CAT1 or 'star' in MATCH_CAT2:
		limhigh = 24

        if PRINTOUTS:
                print 'Forcing magnitudes to be binned with max ', limhigh, '...'

        bins = np.arange(limlow, limhigh, step)


        # b = binned #
        b_hax_mag_median, b_vax_mag_median, b_err_median = [], [], []
	# List of lists. Stores values in each bin #
	b_hax_mag_list, b_vax_mag_list, b_err_list = [], [], []
        counter_empty_bin = 0


        # Bin magnitude errors according to the magnitude on the horizontal axis #
        if swap_hax:
                hax_mag = clean_magnitude2
        if swap_hax is False:
                hax_mag = clean_magnitude1

	# Magnitude on the vertical axis (vax) #	
	vax_mag = np.array(clean_magnitude1) - np.array(clean_magnitude2)


	# Write to log file once outside of loop. FIXME pass filter_name arg to this function #
	fd_mag_bins.write('New filter. First instance: g. Second: r. Third: i. Fourth: z. \n')


	for j in np.arange(limlow, limhigh, step):

                b_hax_mag_temp, b_vax_mag_temp, b_err_temp, counter_err = [], [], [], 0

                for i in np.arange(0, len(clean_magnitude1)):

                        # Do not calculate errors using outlier magnitudes (chosen to be |Delta-M| > 3). Bin magnitude errors according to the magnitude on the horizontal axis of the plot #
                        if hax_mag[i] >= j and hax_mag[i] < j+step and abs(vax_mag[i]) < 3:
                                b_err_temp.append((error1[i]**2 + error2[i]**2)**0.5)
                                b_hax_mag_temp.append(hax_mag[i])
                                b_vax_mag_temp.append(vax_mag[i])
                                counter_err += 1

                if PRINTOUTS:
                        print ' For magnitude, number of objects in bin ', round(j, 2), '-', round(j+step, 2), ': ', counter_err, '...'


		# Write to log file #
		if counter_err == 0:
			write_median, write_err = 'NA', 'NA'
		if counter_err > 0:
			write_median, write_err = np.median(b_hax_mag_temp), np.median(b_err_temp)
		fd_mag_bins.write(str(counter_err) + '\t' + str(round(j, 2)) + '\t' + str(round(j+step, 2)) + '\t' + str(write_median)+ '\t' + str(write_err) + '\n')

                # Add zeros to empty bins and bins with a small number of points #
		# CONST = 0.0002, len(hax_mag) * CONST #
		CONST = 10
		if counter_err <= CONST:
                        counter_empty_bin += 1
                        b_err_median.append(0.0)
                        b_hax_mag_median.append(0.0)
                        b_vax_mag_median.append(0.0)

			b_err_list.append(0)
                        b_hax_mag_list.append(0)
                        b_vax_mag_list.append(0)		

                if counter_err > CONST:
                        b_err_median.append(np.median(b_err_temp))
                        b_hax_mag_median.append(np.median(b_hax_mag_temp))
                        b_vax_mag_median.append(np.median(b_vax_mag_temp))
	
			# Add to list of lists to keep bin structure #	
			b_err_list.append(b_err_temp)
                        b_hax_mag_list.append(b_hax_mag_temp)
                        b_vax_mag_list.append(b_vax_mag_temp)


	#FIXME : Edit `bins` returned by this function... do not consider the bins with <= CONST objs in it

	if PRINTOUTS:
                if swap_hax:
                        print 'Binned clean_magnitude2 with step size: ', step, ', and minimum: ', limlow, ', and maximum: ', limhigh, '...'
		if swap_hax is False:
                        print 'Binned clean_magnitude1 with step size: ', step, ', and minimum: ', limlow, ', and maximum: ', limhigh, '...'
		print ' Calculated errors using objects where |DeltaM| < 3 ... '
		print ' Excluded bins with less than ', CONST, ' objects ... \n'


        return b_hax_mag_median, b_vax_mag_median, b_err_median, bins, b_hax_mag_list, b_vax_mag_list, b_err_list









def normalize_plot_maintain_bin_structure(clean_magnitude1, clean_magnitude2, error1, error2, swap_hax, axlabel1, axlabel2, fd_mag_bins):
	"""Normalize the vertical axis using error and uphold the bin structure. 

	Args:
		clean_magnitude1, clean_magnitude2 (list of floats) --
		error1, error2 (list of floats) --
		swap_hax (bool) -- Swap horizontal axis? Is a global constant.
		fd_mag_bins (file director) --
		axlabel1, axlabel2 -- Allowed values: 
	Returns:
		norm_dm_list (list of list of floats) --
		bins (list of floats) -- Bins used in error calculation. 
	"""

	norm_dm_list, hax_mag_list = [], []

	b_err_median, bins, b_hax_mag_list, b_vax_mag_list = bin_and_cut_measured_magnitude_error_dm3(clean_magnitude1=clean_magnitude1, clean_magnitude2=clean_magnitude2, error1=error1, error2=error2, swap_hax=swap_hax, axlabel1=axlabel1, axlabel2=axlabel2, fd_mag_bins=fd_mag_bins)[2:-1]


	### Identify which magnitude will be plotted on the horizontal axis (hax) ###
        if swap_hax:
                hax_clean_mag = clean_magnitude2
	if swap_hax is False:
                hax_clean_mag = clean_magnitude1

	# Loop through bins (b) #
	for b in np.arange(0, len(b_vax_mag_list)):

		# Normalized Delta-Magnitude (dm) in current bin (icb) #
		norm_dm_icb, hax_mag_icb = [], []	
		vax_mag_icb = b_vax_mag_list[b]

		# 0 is a placeholder for empty bins and bins with few objects #
		if vax_mag_icb == 0:
			norm_dm_list.append(0.0)	
			hax_mag_list.append(0)
	
		#if vax_mag_icb != 0:
		if b_err_median[b] != 0:
			for i in np.arange(0, len(vax_mag_icb)):
				norm = vax_mag_icb[i]/b_err_median[b]
				norm_dm_icb.append(norm)
				hax_mag_icb.append(hax_clean_mag[i])
			# List of lists to keep bin structure #
			hax_mag_list.append(hax_mag_icb)
			norm_dm_list.append(norm_dm_icb)


	return norm_dm_list, bins, hax_mag_list









def normalize_plot(norm_dm_list, bins, hax_mag_list):
	"""Normalize plot to 1-sigma curve using tame magnitude errors only (use bin_and_cut_measured_magnitude_error()).

	Args:
		norm_dm_list (list of list of floats) -- Normalized delta magnitudes. Bin structure preserved.
                bins (list of floats) -- Bins used in error calculation.
                hax_mag_list (list of list of floats) -- Magnitudes on the horizontal axis. Bin structure preserved.
        Returns:
		norm_dm (list of floats) -- Delta-Magnitude normalized by error. Delta-Magnitude computed via magnitude1 - magnitude2. 
                hax_mag (list of floats) -- Magnitude to be plotted on the horizontal axis.
	"""

	### Remove zeros so that lists can be flattened. Zeros (int) were placeholders for missing lists due to empty or small bin. ###
	norm_dm_list[:] = [temp for temp in norm_dm_list if temp != 0]
	hax_mag_list[:] = [temp for temp in hax_mag_list if temp != 0]

	### Flatten lists ###
	hax_mag = [item for sublist in hax_mag_list for item in sublist]
	norm_dm = [item for sublist in norm_dm_list for item in sublist]


	return norm_dm, hax_mag, bins









def one_sigma_counter(norm_delta_mag):
	"""Find the number of objects within 1-sigma, where 1-sigma is calculated according to the error.

	Args:
		norm_delta_mag (list of floats) -- Normalized Delta-Magnitude. This function is only called if NORMALIZE is True.
        Returns:
            counter_1sig (int) -- Number of objects within 1-sigma curve.
	"""

	counter_1sig = 0

	for k in norm_delta_mag:
		if abs(k) < 1.0:
			counter_1sig += 1

	if PRINTOUTS:
		print 'Fraction of objects within 1-sigma: ', counter_1sig, ' / ', len(clean_magnitude1), ' = ', str(float(counter_1sig) / len(clean_magnitude1)), ' \n'

	return counter_1sig









def get_flag_type(df, flag_hdr_list, k):
	"""Print the flag type() once.

	Args:
            df (pandas DataFrame)
            flag_hdr_list (list of str) -- All flag headers to be examined.
            k (int) -- Counter implemented so printout is not repeated.
        Returns:
            0
	"""

	if k == 0:
		for flag_hdr in flag_hdr_list:
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
		color, cmap = 'purple', 'Purples'

	if filter_name == 'z':
		color, cmap = 'blue', 'Blues'

	return color, cmap









def logger(fd_nop, fd_1sig, delta_mag, tile_name, filter_name, run_type, realization_number, clean_magnitude1, full_magnitude1):
	"""Write to log files to record number of objects plotted and number of objects within 1sigma.

	Args:
            fd_nop (file descriptor) -- Write to this file the number of objects plotted (nop).
            fd_1sig (file descriptor) -- Write to this file the number of objects plotted that lie within 1-sigma.            
            filter_name (str) -- Allowed values: 'g' 'r' 'i' 'z'.
            clean_magnitude1 (list of floats) -- Objects with nonzero flags and/or quality cuts removed.
            full_magnitude (list of floats) -- Values read directly from pandas DataFrame via `df[hdr]`; no objects removed using nonzero flag values and no quality cuts performed.
            realization_number (int) -- Allowed values: 0 1 2 None. Refers to Balrog injection and None refers to a one-realization run.
            run_type (str) -- Allowed values: 'ok' 'rerun' None. 'ok' refers to FOF groups unchanged after Balrog injection. 'rerun' refers to changed FOF groups.
        Returns:
            0
	"""

	if NORMALIZE:
		num_1sig = one_sigma_counter(norm_delta_mag=delta_mag)

		# Record number of objects plotted within 1sigma #
		fd_1sig.write('Number of objects within 1sigma for tile ' + str(tile_name) + ', filter ' + str(filter_name) + ', type ' + str(run_type) + ', realization ' + str(realization_number) + ' : ' + str(num_1sig) + ' / ' + str(len(clean_magnitude1)) + ' = ' + str(float(num_1sig) / len(clean_magnitude1)) + '\n')


	# Record number of objects plotted (nop) #
	fd_nop.write('Number of objects plotted after flag cuts for tile ' + str(tile_name) + ', filter ' + str(filter_name) + ', type ' + str(run_type) + ', realization ' + str(realization_number) + ' : ' + str(len(clean_magnitude1)) + ' / ' + str(len(full_magnitude1)) + ' = ' + str(float(len(clean_magnitude1)) / len(full_magnitude1)) + '\n')


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
                # For measured catalog #FIXME in this function add check to see if meas is in axlabel and if not exit and throw error  #
                cbar_idx_list, cbar_bins, cbar_val = calculate_and_bin_cm_T_signal_to_noise(cm_t_hdr=cm_t_hdr, cm_t_err_hdr=cm_t_err_hdr, df=df, idx_good=idx_good, clean_magnitude1=clean_magnitude1, clean_magnitude2=clean_magnitude2)
		cbar_axlabel = 'cm_T_s2n_meas' #FIXME

        if CM_T_ERR_COLORBAR:
                # For measured catalog, cuts performed on truth catalogs #
                cbar_val = get_good_data(df=df, hdr=cm_t_err_hdr, idx_good=idx_good, magnitude=False, filter_name=None)
		cbar_axlabel = str(cm_t_err_hdr[:-2]) + '_' + str(axlabel)
		cbar_idx_list, cbar_bins = None, None

        if CM_T_COLORBAR:
                cbar_val = get_good_data(df=df, hdr=cm_t_hdr, idx_good=idx_good, magnitude=False, filter_name=None)
		cbar_axlabel = str(cm_t_hdr[:-2]) + '_' + str(axlabel)
                cbar_idx_list, cbar_bins = None, None

	if CM_T_S2N_COLORBAR is False and CM_T_ERR_COLORBAR is False and CM_T_COLORBAR is False:
		cbar_val, cbar_idx_list, cbar_bins, cbar_axlabel = None, None, None, None

	if inj and cbar_axlabel is not None:
		cbar_axlabel = 'inj_' + cbar_axlabel
	

	return cbar_val, cbar_idx_list, cbar_bins, cbar_axlabel









def get_errors(mag_err_hdr1, mag_err_hdr2, df, cm_flux_hdr1, cm_flux_hdr2, cm_flux_cov_hdr1, cm_flux_cov_hdr2, filter_name, idx_good):
	"""Get errors for plot.

	Args:
		*_hdr (str ) -- Headers refer to columns in the matched catalog. Can be None.
		df (pandas DataFrame)
	Returns:
		err1, err2 (list of floats) -- Will be None if PLOT_1SIG is False.
	"""

        if PLOT_1SIG:
                if mag_err_hdr1 is None:
                        err1 = calculate_total_fractional_magnitude_error(df=df, flux_hdr=cm_flux_hdr1, cov_hdr=cm_flux_cov_hdr1, filter_name=filter_name, idx_good=idx_good)
                if mag_err_hdr1 is not None:
			if 'star' not in MATCH_CAT1:
                                err1 = df[mag_err_hdr1]
			if 'star' in MATCH_CAT1:
                                err1 = df[str(mag_err_hdr1[:-2]) + '_' + filter_name.upper() + str(mag_err_hdr1[-2:])]

                if mag_err_hdr2 is None:
                        err2 = calculate_total_fractional_magnitude_error(df=df, flux_hdr=cm_flux_hdr2, cov_hdr=cm_flux_cov_hdr2, filter_name=filter_name, idx_good=idx_good)
                if mag_err_hdr2:
                        err2 = df[mag_err_hdr2]

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

	clean_data = []

	for i in idx_good:
		clean_data.append(full_data[i])

	return clean_data









def get_plot_variables(filter_name, df, mag_hdr1, mag_hdr2, flag_hdr_list, cm_s2n_r_hdr1, cm_s2n_r_hdr2, psfrec_t_hdr1, psfrec_t_hdr2, mag_err_hdr1, mag_err_hdr2, cm_flux_hdr1, cm_flux_hdr2, cm_flux_cov_hdr1, cm_flux_cov_hdr2, realization_number, tile_name, run_type, axlabel1, axlabel2, cm_t_hdr1, cm_t_hdr2, cm_t_err_hdr1, cm_t_err_hdr2, flag_list):
	"""Get quantities needed for plotter() and subplotter().

	Args:
		df (pandas DataFrame)
		
		*_hdr (str) -- Can be None.
		run_type (str) -- Can be None.
	Returns:
	"""

	# Rewrite mag_axlabels. Transform, for example, cm_mag_true to cm_mag_{filter}_true or psf_mag_meas to psf_mag_{filter}_meas #
	mag_axlabel1 = str(mag_hdr1[:-2]) + '_' + str(axlabel1)
	mag_axlabel2 = str(mag_hdr2[:-2]) + '_' + str(axlabel2)

        if cm_t_hdr1 is not None:
                cm_t_axlabel1 = str(cm_t_hdr1[:-2]) + '_' + str(axlabel1)
        if cm_t_hdr2 is not None:
                cm_t_axlabel2 = str(cm_t_hdr2[:-2]) + '_' + str(axlabel2)

        if cm_t_err_hdr1 is not None:
                cm_t_err_axlabel1 = str(cm_t_err_hdr1[:-2]) + '_' + str(axlabel1)
        if cm_t_err_hdr2 is not None:
                cm_t_err_axlabel2 = str(cm_t_err_hdr2[:-2]) + '_' + str(axlabel2)



	### Define variables ###

	# Get magnitude1 #
	fullmag1 = get_floats_from_string(df=df, hdr=mag_hdr1, filter_name=filter_name)

	# Get magnitude2 #
	fullmag2 = get_floats_from_string(df=df, hdr=mag_hdr2, filter_name=filter_name)



	### Clean the data: removed flags and/or perform quality cuts ###
	if EH_CUTS:
		idx_good = get_good_index_using_quality_cuts(cm_flag_hdr1=cm_flag_hdr1, cm_flag_hdr2=cm_flag_hdr2, df=df, cm_t_hdr1=cm_t_hdr1, cm_t_hdr2=cm_t_hdr2, cm_t_err_hdr1=cm_t_err_hdr1, cm_t_err_hdr2=cm_t_err_hdr2, cm_s2n_r_hdr1=cm_s2n_r_hdr1, cm_s2n_r_hdr2=cm_s2n_r_hdr2, flag_hdr1=flag_hdr1, flag_hdr2=flag_hdr2, full_magnitude1=fullmag1, full_magnitude2=fullmag2, psfrec_t_hdr1=psfrec_t_hdr1, psfrec_t_hdr2=psfrec_t_hdr2, cm_t_err_axlabel1=cm_t_err_axlabel1, cm_t_err_axlabel2=cm_t_err_axlabel2)[0]
	
	if EH_CUTS is False:
                idx_good = get_good_index_using_primary_flags(df=df, cm_flag_hdr1=flag_list[2], cm_flag_hdr2=flag_list[3], flag_hdr1=flag_list[0], flag_hdr2=flag_list[1], full_magnitude1=fullmag1, full_magnitude2=fullmag2)[0]

	cleanmag1 = get_good_data(df=df, hdr=mag_hdr1, idx_good=idx_good, magnitude=True, filter_name=filter_name)
	cleanmag2 = get_good_data(df=df, hdr=mag_hdr2, idx_good=idx_good, magnitude=True, filter_name=filter_name)


	# Some variables set to None because must pass to plotter() #
	cbar_val, cbar_idx_list, cbar_bins, cbar_axlabel = get_colorbar_value(df=df, cm_t_hdr=cm_t_hdr2, cm_t_err_hdr=cm_t_err_hdr2, idx_good=idx_good, clean_magnitude1=cleanmag1, clean_magnitude2=cleanmag2, axlabel=axlabel2, inj=INJ2)


	### Define errors ###
	err1, err2 = get_errors(mag_err_hdr1=mag_err_hdr1, mag_err_hdr2=mag_err_hdr2, df=df, cm_flux_hdr1=cm_flux_hdr1, cm_flux_hdr2=cm_flux_hdr2, cm_flux_cov_hdr1=cm_flux_cov_hdr1, cm_flux_cov_hdr2=cm_flux_cov_hdr2, filter_name=filter_name, idx_good=idx_good)


	### Write flags to file ###
	if LOG_FLAGS:
		for i in np.arange(0, len(flag_list), 2):
			# Bad index #
			temp_idx = handle_flags(df=df, fd_flag=fd_flag, filter_name=f, realization_number=realization_number, flag_hdr1=flag_list[i], flag_hdr2=flag_list[i+1], full_magnitude1=fullmag1, full_magnitude2=fullmag2, tile_name=tile_name, run_type=run_type)[1]
			#flag_idx.append(temp_idx)
			flag_list.extend(temp_idx)
			if SHOW_PLOT is False and SAVE_PLOT is False:
			# Reset to avoid errors #
				counter_subplot = 1


	### Print out the type() for each flag ###
	if SHOW_FLAG_TYPE:
		get_flag_type(df=df, flag_hdr_list=flag_list, k=counter_flag_type_printout)
		counter_flag_type_printout += 1


	return cbar_val, cbar_idx_list, cbar_bins, err1, err2, cleanmag1, cleanmag2, idx_good, cbar_axlabel, fullmag1, mag_axlabel1, mag_axlabel2	









def plotter(cbar_val, error1, error2, fd_nop, fd_1sig, filter_name, clean_magnitude1, full_magnitude1, mag_axlabel1, clean_magnitude2, mag_axlabel2, plot_title, realization_number, swap_hax, run_type, tile_name, idx_list, bins, cbar_axlabel, fd_mag_bins, ylow, yhigh):
	"""Plot a single magnitude versus delta-magnitude plot.

	Args:
            fd_nop (file descriptor) -- Write to this file the number of objects plotted (nop).
            fd_1sig (file descriptor) -- Write to this file the number of objects plotted that lie within 1-sigma.
            full_magnitude1, full_magnitude2 (numpy.ndarray if directly from `df` OR list of floats if from `get_floats_from_string()`) -- Values read directly from pandas DataFrame via `df[hdr]`; no objects removed using nonzero flag values and no quality cuts performed.
            realization_number (str) -- Allowed values: 0 1 2 None. Refers to Balrog injection and None refers to a one-realization run.
            run_type (str) -- Allowed values: 'ok' 'rerun' None. 'ok' refers to FOF groups unchanged after Balrog injection. 'rerun' refers to changed FOF groups.
            swap_hax (bool) -- "swap horizontal axis". If False plots magnitude1 on the x-axis. If True plots magnitude2 on the x-axis.
        Returns:
		0
	"""

	### Get labels for vertical and horizontal axes. ###
	# Transform, for example, cm_mag_true to cm_mag_{filter}_true, or psf_mag_meas to psf_mag_{filter}_meas #
        mag_axlabel1 = mag_axlabel1[:-4] + str(filter_name) + '_' + mag_axlabel1[-4:]
        mag_axlabel2 = mag_axlabel2[:-4] + str(filter_name) + '_' + mag_axlabel2[-4:]
	if INJ1:
		mag_axlabel1 = 'inj_' + mag_axlabel1
	if INJ2:
		mag_axlabel2 = 'inj_' + mag_axlabel2


	### Values to plot for normalized plot ###
	if NORMALIZE:
		# Args needed to call normalize_plot() #
		norm_dm_list, bins, hax_mag_list = normalize_plot_maintain_bin_structure(clean_magnitude1=clean_magnitude1, clean_magnitude2=clean_magnitude2, error1=error1, error2=error2, swap_hax=swap_hax, axlabel1=axlabel1, axlabel2=axlabel2, fd_mag_bins=fd_mag_bins) 

		deltam, hax_mag, bins = normalize_plot(norm_dm_list=norm_dm_list, bins=bins, hax_mag_list=hax_mag_list)


                # Labels and appearance #
                plt.ylabel('('+str(mag_axlabel1) + ' - ' + str(mag_axlabel2)+') / '+ '$\sigma$', fontsize=8)

		if PLOT_1SIG:
			### Plot the 68th percentile calculated from np.percentile() ###
			vax_68percentile_list, bins = get_68percentile_from_normalized_data(norm_dm_list=norm_dm_list, bins=bins, hax_mag_list=hax_mag_list)
			# Plot legend once #
			counter_legend = 0
			for b in np.arange(0, len(vax_68percentile_list)-1):
				if vax_68percentile_list[b] != 0:
					if counter_legend == 0:	
						plt.plot(np.array([bins[b], bins[b+1]]), np.array([vax_68percentile_list[b], vax_68percentile_list[b]]), color='fuchsia', label='$P_{68}$', linewidth=0.7)
						counter_legend += 1
					if counter_legend == 1:
						# Horizontal bounds #
						plt.plot(np.array([bins[b], bins[b+1]]), np.array([vax_68percentile_list[b], vax_68percentile_list[b]]), color='fuchsia', linewidth=0.7)
						plt.plot(np.array([bins[b], bins[b+1]]), -1.0*np.array([vax_68percentile_list[b], vax_68percentile_list[b]]), color='fuchsia', linewidth=0.7)
						# Vertical bounds #
						plt.plot(np.array([bins[b], bins[b]]), np.array([-1*vax_68percentile_list[b], vax_68percentile_list[b]]), color='fuchsia', linewidth=0.5, linestyle=':')
						plt.plot(np.array([bins[b+1], bins[b+1]]), np.array([-1*vax_68percentile_list[b], vax_68percentile_list[b]]), color='fuchsia', linewidth=0.5, linestyle=':')
			

			### Plot 1-sigma curve according to error calculation ###
			plt.axhline(y=1.0, color='red', linestyle='--', linewidth=0.7, label='$1 \sigma_{mag\_meas}$')
			plt.axhline(y=-1.0, color='red', linestyle='--', linewidth=0.7)



	### For scatter plot ###
	if NORMALIZE is False:
		# Values to plot #
		deltam = np.array(clean_magnitude1) - np.array(clean_magnitude2)
		if swap_hax:
			hax_mag = clean_magnitude2
		if swap_hax is False:
                        hax_mag = clean_magnitude1
		# Labels and appearance #
		plt.ylabel(str(mag_axlabel1) + ' - ' + str(mag_axlabel2), fontsize=9)
		### 1-sigma curve ###
		if PLOT_1SIG:
			hax, vax, err, placehold1, p3, p4, p5 = bin_and_cut_measured_magnitude_error_dm3(error1=error1, error2=error2, clean_magnitude1=clean_magnitude1, clean_magnitude2=clean_magnitude2, swap_hax=swap_hax, axlabel1=mag_axlabel1, axlabel2=mag_axlabel2, fd_mag_bins=fd_mag_bins)
			### Remove zeros from x, y, and err (zeros were placeholders for instances in which there were no objects in a particular magnitude bin) ###
			err[:] = [temp for temp in err if temp != 0]
			hax[:] = [temp for temp in hax if temp != 0]
			vax[:] = [temp for temp in vax if temp != 0]
			### Plot 1-sigma curve ###
			plt.plot(hax, np.array(vax) + np.array(err), color='red', linestyle='-', linewidth=0.7, label='$1 \sigma_{mag\_meas}$')
			plt.plot(hax, np.array(vax) - np.array(err), color='red', linestyle='-', linewidth=0.7)


	if PRINTOUTS:
		print 'Plotting ', len(clean_magnitude1), ' objects ... \n'


	### Write to log files to record the number of objects plotted and the number of objects within 1sigma ###
	logger(delta_mag=deltam, fd_nop=fd_nop, fd_1sig=fd_1sig, filter_name=filter_name, clean_magnitude1=clean_magnitude1, full_magnitude1=full_magnitude1, realization_number=realization_number, run_type=run_type, tile_name=tile_name)


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
		plt.colorbar()#FIXME: label?


	# Labels and appearance #
	if swap_hax:
		plt.xlabel(str(mag_axlabel2), fontsize=9)
	if swap_hax is False:
                plt.xlabel(str(mag_axlabel1), fontsize=9)
	plt.axhline(y=0.0, color='k', linestyle=':', linewidth=0.5)
        if ylow is not None and yhigh is not None:
            plt.ylim([ylow, yhigh])


	### Plot legend ###
	if PLOT_1SIG and BIN_CM_T_S2N is False:
		plt.legend(fontsize=8).draggable()
	if BIN_CM_T_S2N:
		# Increase marker size and opacity in legend #
		lgnd = plt.legend(markerscale=4, fontsize=8)
		for l in lgnd.legendHandles:
			l.set_alpha(1)


	if SUBPLOT is False and SHOW_PLOT:
		plt.title(plot_title)
		plt.show()


        return 0









def subplotter(cm_flux_hdr1, cm_flux_hdr2, cm_flux_cov_hdr1, cm_flux_cov_hdr2, cm_t_hdr1, cm_t_hdr2, cm_t_err_hdr1, cm_t_err_hdr2, cm_s2n_r_hdr1, cm_s2n_r_hdr2, psfrec_t_hdr1, psfrec_t_hdr2, df, fd_1sig, fd_flag, fd_nop, flag_idx, flag_list, mag_hdr1, mag_hdr2, mag_err_hdr1, mag_err_hdr2, plot_name, plot_title, realization_number, run_type, tile_name, swap_hax, stack_realizations, cm_t_s2n_axlabel1, cm_t_s2n_axlabel2, axlabel1, axlabel2, fd_mag_bins, ylow, yhigh):
	"""Combine four subplots into a single plot with four panels (2-by-2). Declare variables needed for plotting.

	Args:
            *_hdr (str) -- Headers refer to columns in the matched catalog.
            df (pandas DataFrame)
            fd_nop (file descriptor) -- Write to this file the number of objects plotted (nop).
            fd_1sig (file descriptor) -- Write to this file the number of objects plotted that lie within 1-sigma.
            fd_flag (file descriptor) -- Write to this file the nonzero flags if LOG_FLAGS is True.
            flag_list (list of str) -- List of all headers corresponding to flags. List must be paired via [ 'flags_1', 'flags_2', 'cm_flags_1', 'cm_flags_2', ... ] but headers without a pair can be set to None. The first two pairs must be 'flags' and 'cm_flags'. 
            plot_name (str) -- Path and name for the plot. Used when save_plot is True and normalize is False.
            realization_number (int) -- Allowed values: 0 1 2 None. Refers to Balrog injection and None refers to a one-realization run.
            run_type (str) -- Allowed values: 'ok' 'rerun' None. 'ok' refers to FOF groups unchanged after Balrog injection. 'rerun' refers to changed FOF groups.
            swap_hax (bool) -- "swap horizontal axis". If False plots magnitude1 on the x-axis. If True plots magnitude2 on the x-axis.
        Returns:
            flag_idx (list of ints) -- If log_flags is True, will check for all nonzero flag values in `flag_list` and `flag_idx` will contain indices that have nonzero flag values. Will be empty if LOG_FLAGS is False.
	"""



	# Counter for flag type() printout #
	counter_flag_type_printout = 0


        ### Create 4-by-4 subplot ###
	counter_subplot = 1
	# Set figure size once. Units: inches #
	plt.figure(figsize=(10, 8))

	# The coadds have mag_err header and a separate _cat.fits for each filter #
	#TODO add filter_dependent_df or combine catalogs


        ### Create one subplot for each griz filter ###
	for f in ALL_FILTERS:

		### Define variables ###
		cbar_val, cbar_idx_list, cbar_bins, err1, err2, cleanmag1, cleanmag2, index_good, cbar_axlabel, fullmag1, mag_axlabel1, mag_axlabel2 = get_plot_variables(filter_name=f, df=df, mag_hdr1=mag_hdr1, mag_hdr2=mag_hdr2, flag_hdr_list=flag_list, cm_s2n_r_hdr1=cm_s2n_r_hdr1, cm_s2n_r_hdr2=cm_s2n_r_hdr2, psfrec_t_hdr1=psfrec_t_hdr1, psfrec_t_hdr2=psfrec_t_hdr2, mag_err_hdr1=mag_err_hdr1, mag_err_hdr2=mag_err_hdr2, cm_flux_hdr1=cm_flux_hdr1, cm_flux_hdr2=cm_flux_hdr2, cm_flux_cov_hdr1=cm_flux_cov_hdr1, cm_flux_cov_hdr2=cm_flux_cov_hdr2, realization_number=realization_number, tile_name=tile_name, run_type=run_type, axlabel1=axlabel1, axlabel2=axlabel2, cm_t_hdr1=cm_t_hdr1, cm_t_hdr2=cm_t_hdr2, cm_t_err_hdr1=cm_t_err_hdr1, cm_t_err_hdr2=cm_t_err_hdr2, flag_list=flag_list)


		### Subplot ###
		plt.subplot(2, 2, counter_subplot)
		plotter(cbar_val=cbar_val, plot_title=plot_title, fd_1sig=fd_1sig, run_type=run_type, error1=err1, error2=err2, fd_nop=fd_nop, filter_name=f, full_magnitude1=fullmag1, clean_magnitude1=cleanmag1, clean_magnitude2=cleanmag2, mag_axlabel1=mag_axlabel1, mag_axlabel2=mag_axlabel2, realization_number=realization_number, tile_name=tile_name, swap_hax=swap_hax, idx_list=cbar_idx_list, bins=cbar_bins, fd_mag_bins=fd_mag_bins, cbar_axlabel=cbar_axlabel, ylow=ylow, yhigh=yhigh)
		counter_subplot += 1


	### Show or save the plot once all four subplots have been filled ###
	plt.subplots_adjust(hspace=0.4)
	plt.subplots_adjust(wspace=0.3)
	plt.tight_layout(pad=3, h_pad=2.5)


	### Title ###
	plt.suptitle(plot_title)


	### Save plot ###
	if SAVE_PLOT:
		filename = plot_name
		print '-----> Saving plot as: ', filename
		plt.savefig(plot_name)

	### Show plot ###
	if SHOW_PLOT:
		plt.show()

	
	return flag_idx









def get_match_type(title_piece1, title_piece2):
        """Transform plot title of form 'Inj MOF Cat & Truth Cat' to 'inj_mof_cat_truth_cat'.

        Args:
                title_piece1, title_piece2 (str) -- Ex: Injected MOF
        Return:
                match_type (str) -- Ex: injected_mof_truth_catalog
        """

        title_piece1, title_piece2  = title_piece1.lower(), title_piece2.lower()
        title_piece1, title_piece2 = title_piece1.replace(' ', '_'), title_piece2.replace(' ', '_')
        match_type = str(title_piece1)+'_'+str(title_piece2)

        return match_type









def get_fd_names(title_piece1, title_piece2, outdir):
        """Generate names for the following log files: flags, magnitude bins, number of objects plotted, number of objects within one sigma.

        Args:
                title_piece1, title_piece2 (str) --
                outdir (str) -- Output directory. Files will be saved here.
        Returns:
                fn1, fn2, fn3, fn4 (str) -- Filenames for flag log file, magnitude log, number of objects plotted log, number of objects within 1-sigma, respectively.
        """

        match_type = get_match_type(title_piece1, title_piece2)

	fn1 = os.path.join(outdir, 'flag_log_'+str(match_type)+'.csv')
	fn2 = os.path.join(outdir, 'magnitude_bins_'+str(match_type)+'.txt')
	fn3 = fn = os.path.join(outdir, 'num_objs_plotted_'+str(match_type)+'.txt')
	fn4 = os.path.join(outdir, 'one_sigma_objects_'+str(match_type)+'.txt')

	return fn1, fn2, fn3, fn4









def get_plot_suptitle(realization_number, tile_name, title_piece1, title_piece2):
	"""Generate plot title.

	Args:
		match_type (str) -- Ex: inj_mof_vs_truth_cat #FIXME should contain numbers for subscript so axes are easier to read. USE TITLE PIECE
		realization_number (str) -- Allowed values: '0' '1' '2' ... 'stacked'.
		tile_name (str)
	Returns:
		title (str) -- Ex: 'Inj MOF Cat & Truth Cat' 
	"""

	title = str(title_piece1) + ' & ' + str(title_piece2) +'. Tile: ' + str(tile_name) + '. Realization: ' + str(realization_number) + '. \n'
	if NORMALIZE:
		title = 'Normalized. ' + title

	return title









def get_plot_save_name(outdir, realization_number, tile_name, title_piece1, title_piece2, ylow, yhigh):
        """Generate name of the plot that will be used in plt.savefig(). Relies on directory structure fn_dir/{tile}/{plot_type}/{realization}/ where allowed values for plot_type are: 'normalized' 'scatter'. 

        Args:
                outdir (str) -- Output directory
                realization_number (str) -- Allowed values: '0' '1' '2' 'stacked'
                tile_name (str)
        Returns:
                fn (str) -- The complete filename which includes path.
        """

        match_type = get_match_type(title_piece1, title_piece2)

	if YLOW is None and YHIGH is None:
		# Default scale for the vertical axis (vax) is used #
		ylim = 'defaultvax'
	if YLOW is not None and YHIGH is not None:
		ylim = str(YLOW)+'y'+str(YHIGH)
	
	endname = str(tile_name) + '_' + str(realization_number) + '_' + str(match_type) + '_'+str(ylim)+'_rm_flags.png'

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
        if NORMALIZE:
                fn = os.path.join(outdir, tile_name, 'normalized', realization_number, 'norm_' + str(outname))
        if NORMALIZE is False:
                fn = os.path.join(outdir, tile_name, 'scatter', realization_number, outname)

        return fn









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









################################################################### Store catalog information ###################################################################

class CoaddCat():
	"""Declare headers for coadd catalogs .../coadd/{tile}_{filter}_cat.fits. There is a separate catalog for each filter."""

        # Once matched, headers will have form 'hdr_1' or 'hdr_2' with a suffix (suf) #
	def __init__(self, inj, suf):
		"""Declare headers.

		Args:
			inj (bool) -- Balrog injected catalog?
			suf (int) -- Refers to order in which catalog was matched in ms_matcher. Allowed values: 1 2
		"""

		# For plot title #
		if inj:
			self.title_piece = 'Inj Coadd Cat (_' + str(suf) +')'
		if inj is False:
			self.title_piece = 'Coadd Cat (_' + str(suf) +')'
		self.axlabel = 'meas'
		# Magnitude, is one number #
		self.mag_hdr = 'MAG_AUTO_' + str(suf)
		self.mag_axlabel = 'MAG_AUTO_meas'
		# For error calculation #
		self.mag_err_hdr = 'MAGERR_AUTO_' + str(suf)
		self.flux_hdr = None
		self.flux_cov_hdr = None
		# Size #
		self.cm_t_hdr = None
		self.cm_t_err_hdr = None
		self.cm_t_s2n_axlabel = None
		# Flags #
		self.flags_hdr = 'FLAGS_' + str(suf)
		self.obj_flags_hdr = None
		self.psf_flags_hdr = None
		self.cm_flags_hdr = None
		self.cm_max_flags_hdr = None
		self.cm_mof_flags_hdr = None
		self.cm_flags_r_hdr = None
		# For region file #
		self.ra_hdr = 'ALPHAWIN_J2000' + str(suf)
		self.dec_hdr = 'DELTAWIN_J2000' + str(suf)
		self.a_hdr = 'A_IMAGE_' + str(suf)
		self.b_hdr = 'B_IMAGE_' + str(suf)









class SOFGalTruthCat():
	"""Declare headers and axes labels for galaxy truth catalogs in the sof directory /data/des71.a/data/kuropat/des2247-4414_sof/."""

	# Once matched, headers will have form 'hdr_1' or 'hdr_2' with a suffix (suf) #
	def __init__(self, inj, suf):
		"""Declare constants.

                Args:
                        inj (bool) -- Balrog injected catalog?
			suf (int) -- Refers to order in which catalog was matched in ms_matcher. Allowed values: 1 2
                """

		if inj:
			self.title_piece = 'Inj Gal Truth Cat'
		if inj is False:
			self.title_piece = 'Gal Truth Cat'
		self.axlabel = 'true'
		# Headers are the same as MOFCat class. Reproduced below for clarity in ms_plotter.py #
		# Magnitude, is string of form '(mag_g, mag_r, mag_i, mag_z)' #
	    	self.mag_hdr = 'cm_mag_' + str(suf)
		self.mag_axlabel = 'cm_mag_true'
                # For error calculation #
                self.mag_err_hdr = None
                self.cm_flux_hdr = 'cm_flux_' + str(suf)
                self.cm_flux_cov_hdr = 'cm_flux_cov_' + str(suf)
                # Size #
                self.cm_t_hdr = 'cm_T_'  + str(suf)
                self.cm_t_err_hdr = 'cm_T_err_'  + str(suf)
		self.cm_t_s2n_axlabel = 'cm_T_s2n_true'
                # Flags #
                self.flags_hdr = 'flags_' + str(suf)
                self.obj_flags_hdr = 'obj_flags_' + str(suf)
                self.psf_flags_hdr = 'psf_flags_' + str(suf)
                self.cm_flags_hdr = 'cm_flags_' + str(suf)
                self.cm_max_flags_hdr = 'cm_max_flags_' + str(suf)
                self.cm_flags_r_hdr = 'cm_flags_r_' + str(suf)
		self.cm_mof_flags_hdr = 'cm_mof_flags_' + str(suf)
                # For region file #
                self.ra_hdr = 'ra_' + str(suf)
                self.dec_hdr = 'dec_' + str(suf)
                self.a_hdr, self.b_hdr = None, None
                # For Eric Huff (EH) quality cuts #
                self.cm_s2n_r_hdr = None
                self.psfrec_t_hdr = None
                self.eh_cuts = False



	





class SOFCat():
        """Declare headers and axis labels for sof catalog in /data/des71.a/data/kuropat/des2247-4414_sof/y3v02/balrog_images/${re}/${t}/sof/${t}_sof.fits and /data/des71.a/data/kuropat/sof_stars/y3v02/balrog_images/{realization}/{tile}/mof/{tile}_mof.fits"""

	# Once matched, headers will have form 'hdr_1' or 'hdr_2' with a suffix (suf) #
        def __init__(self, inj, suf):
		"""Declare constants.

                Args:
                        inj (bool) -- Balrog injected catalog?
			suf (int) -- Refers to order in which catalog was matched in ms_matcher. Allowed values: 1 2
                """

                if inj:
                        self.title_piece = 'Inj SOF Cat'
                if inj is False:
                        self.title_piece = 'SOF Cat'
		self.axlabel = 'meas'
                # Headers are the same as MOFCat class with the exception of cm_mof_flags_hdr. Reproduced below for clarity #
		# Magnitude, is string of form '(mag_g, mag_r, mag_i, mag_z)' #
		self.mag_hdr = 'cm_mag_' + str(suf)
		self.mag_axlabel = 'cm_mag_meas'
                # For error calculation #
                self.mag_err_hdr = None
                self.cm_flux_hdr = 'cm_flux_' + str(suf)
                self.cm_flux_cov_hdr = 'cm_flux_cov_' + str(suf)
                # Size #
                self.cm_t_hdr = 'cm_T_'  + str(suf)
                self.cm_t_err_hdr = 'cm_T_err_'  + str(suf)
                self.cm_t_s2n_axlabel = 'cm_T_s2n_meas'
		# Flags #
                self.flags_hdr = 'flags_' + str(suf)
                self.obj_flags_hdr = 'obj_flags_' + str(suf)
                self.psf_flags_hdr = 'psf_flags_' + str(suf)
                self.cm_flags_hdr = 'cm_flags_' + str(suf)
                self.cm_max_flags_hdr = 'cm_max_flags_' + str(suf)
                self.cm_flags_r_hdr = 'cm_flags_r_' + str(suf)
                self.cm_mof_flags_hdr = None
                # For region file #
                self.ra_hdr = 'ra_' + str(suf)
                self.dec_hdr = 'dec_' + str(suf)
                self.a_hdr, self.b_hdr = None, None
                # For Eric Huff (EH) quality cuts #
                self.cm_s2n_r_hdr = None
                self.psfrec_t_hdr = None
                self.eh_cuts = False









class MOFCat():
	"""Declare headers and axes labels for MOF catalog. Currently, the galaxy truth catalogs are created using MOF and have the same headers. Works (mostly) with /data/des71.a/data/kuropat/sof_stars/y3v02/balrog_images/{realization}/{tile}/{tile}_{realization}_balrog_truth_cat_gals.fits, /data/des71.a/data/kuropat/sof_stars/y3v02/balrog_images/{realization}/{tile}/mof/{tile}_mof.fits, ..."""

	# Once matched, headers will have form 'hdr_1' or 'hdr_2' with a suffix (suf) #
        def __init__(self, inj, suf):
		"""Declare constants.

		Args:
			inj (bool) -- Balrog injected catalog?
			suf (int) -- Refers to order in which catalog was matched in ms_matcher. Allowed values: 1 2
		"""

		# For plot title #
		if inj:
			self.title_piece = 'Inj MOF Cat'
		if inj is False:
			self.title_piece = 'MOF Cat'
		self.axlabel = 'meas'
		# Magnitude, is string of form (mag_g, mag_r, mag_i, mag_z)  #
		self.mag_hdr = 'cm_mag_' + str(suf)
		self.mag_axlabel = 'cm_mag_meas'
		# For error calculation #
		self.mag_err_hdr = None 
		self.cm_flux_hdr = 'cm_flux_' + str(suf)
		self.cm_flux_cov_hdr = 'cm_flux_cov_' + str(suf)
		# Size #
		self.cm_t_hdr = 'cm_T_'  + str(suf)
		self.cm_t_err_hdr = 'cm_T_err_'  + str(suf)
		self.cm_t_s2n_axlabel = 'cm_T_s2n_meas'
		# Flags #
		self.flags_hdr = 'flags_' + str(suf)
                self.obj_flags_hdr = 'obj_flags_' + str(suf)
                self.psf_flags_hdr = 'psf_flags_' + str(suf)
                self.cm_flags_hdr = 'cm_flags_' + str(suf)
                self.cm_max_flags_hdr = 'cm_max_flags_' + str(suf)
                self.cm_flags_r_hdr = 'cm_flags_r_' + str(suf)
		self.cm_mof_flags_hdr = 'cm_mof_flags_' + str(suf)
		# For region file #
		self.ra_hdr = 'ra_' + str(suf)
		self.dec_hdr = 'dec_' + str(suf)
		self.a_hdr, self.b_hdr = None, None
		# For Eric Huff (EH) quality cuts #
                self.cm_s2n_r_hdr = None
                self.psfrec_t_hdr = None
                self.eh_cuts = False









class StarTruthCat(): #are there sep headers for MOFStarTruthCat and SOFStarTruthCat?
        """Declare headers and axes labels for star truth catalogs in /data/des71.a/data/kuropat/sof_stars/."""

	# Once matched, headers will have form 'hdr_1' or 'hdr_2' with a suffix (suf) #
	def __init__(self, inj, suf):
		"""Declare constants.
	
		Args:
			inj (bool) -- Balrog injected catalog?
			suf (int) -- Refers to order in which catalog was matched in ms_matcher. Allowed values: 1 2
		"""
		
		if inj:
			self.title_piece = 'Inj Star Truth Cat'
		if inj is False:
			self.title_piece = 'Star Truth Cat'
		self.axlabel = 'true'
		# Magnitude #
		self.mag_hdr = 'g_Corr_' + str(suf) 
		self.mag_axlabel = 'mag_true' #FIXME 
		# For error calculation #
                # Is of form 'PSF_MAG_ERR_{filter}_' + str(suf) #
		self.mag_err_hdr = 'PSF_MAG_ERR_' + str(suf)
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
		self.ra_hdr = 'RA_new_' + str(suf)
		self.dec_hdr = 'DEC_new_' + str(suf)
		self.a_hdr, self.b_hdr = None, None # Use cm_T
		# For Eric Huff (EH) quality cuts #
		self.cm_s2n_r_hdr = None
		self.psfrec_t_hdr = None
		self.eh_cuts = False
		








def get_class(cat_type, inj, suf):
        """Get the appropriate class for the catalog type.

        Args:
                cat_type -- Catalog type. Allowed values: 'gal_truth', 'mof', 'star_truth', 'sof', 'coadd'.
		inj (bool) -- Balrog injected catalog?
		suf (int) -- Refers to order in which catalog was matched in ms_matcher. Allowed values: 1 2
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









def get_catalog(basepath, cat_type, inj, realization_number, tile_name):
        """Get catalog to analyze.
	
        Args:
                basepath
                cat_type -- Catalog type. Allowed values: 'gal_truth', 'mof', 'star_truth', 'sof', 'coadd'.
                inj (bool)
                realization_number (str) -- Allowed values: '0' '1' '2' ...
                tile_name -- Different allowed values depending on catalog.
        Returns:
                fn -- Filename
        """

        if cat_type == 'gal_truth' and inj:
                fn = os.path.join(basepath, 'balrog_images', realization_number, tile_name, tile_name+'_'+realization_number+'_balrog_truth_cat_gals.fits')
        if cat_type == 'gal_truth' and inj is False:
                print 'No non-injected truth catalog exists.'

        if cat_type == 'star_truth' and inj:
                fn = os.path.join(basepath, 'balrog_images', realization_number, tile_name, tile_name+'_'+realization_number+'_balrog_truth_cat_stars.fits')
        if cat_type == 'star_truth' and inj is False:
                print 'No non-injected truth catalog exists.'

        if cat_type == 'sof' and inj:
                fn = os.path.join(basepath, 'balrog_images', realization_number, tile_name, 'sof', tile_name+'_sof.fits')
        if cat_type == 'sof' and inj is False:
                fn = os.path.join(basepath, tile_name, 'sof', tile_name+'_sof.fits')

        if cat_type == 'mof' and inj:
                fn = os.path.join(basepath, 'balrog_images', realization_number, tile_name, 'mof', tile_name+'_mof.fits')
        if cat_type == 'mof' and inj is False:
                fn = os.path.join(basepath, tile_name, 'mof', tile_name+'_mof.fits')

        #FIXME fn for coadds is dependent on filter
        if cat_type == 'coadd' and inj:
                fn = os.path.join(basepath, realization_number, tile_name, 'coadd', tile_name+'i_cat.fits')
        if cat_type == 'coadd' and inj is False:
                fn = os.path.join(basepath, tile_name, 'coadd', tile_name+'i_cat.fits')

        return fn









################################################################### Declare necessary variables ###################################################################

### For data analysis ###
# class1 refers to in1 in ms_matcher. in1 appends _1 to all the headers, hence suf=1. #
class1 = get_class(cat_type=MATCH_CAT1, inj=INJ1, suf=1)
class2 = get_class(cat_type=MATCH_CAT2, inj=INJ2, suf=2)

# Get arguments to pass to ms_matcher. Need to transform header of form 'ra_1' to 'ra', hence [:-2] #
ra1, ra2 = class1.ra_hdr[:-2], class2.ra_hdr[:-2]
dec1, dec2 = class1.dec_hdr[:-2], class2.dec_hdr[:-2]

# For plot labels #
axlabel1, axlabel2 = class1.axlabel, class2.axlabel

# Magnitudes #
m1_hdr, m2_hdr = class1.mag_hdr, class2.mag_hdr
mag_axlabel1, mag_axlabel2 = class1.mag_axlabel, class2.mag_axlabel

# For error calculation #
mag_err_hdr1, mag_err_hdr2 = class1.mag_err_hdr, class2.mag_err_hdr
cm_flux_hdr1, cm_flux_hdr2 = class1.cm_flux_hdr, class2.cm_flux_hdr
cm_flux_cov_hdr1, cm_flux_cov_hdr2 = class1.cm_flux_cov_hdr, class2.cm_flux_cov_hdr

# For signal to noise calculation #
cm_t_hdr1, cm_t_hdr2 = class1.cm_t_hdr, class2.cm_t_hdr
cm_t_err_hdr1, cm_t_err_hdr2 = class1.cm_t_err_hdr, class2.cm_t_err_hdr
cm_t_s2n_axlabel1, cm_t_s2n_axlabel2 = class1.cm_t_s2n_axlabel, class2.cm_t_s2n_axlabel

# Flags #
flags_1_hdr, flags_2_hdr = class1.flags_hdr, class2.flags_hdr
obj_flags_1_hdr, obj_flags_2_hdr = class1.obj_flags_hdr, class2.obj_flags_hdr
# psf_flags is a string of form '(0,0,0,0)'; must pass through get_floats_from_string() if used. #
psf_flags_1_hdr, psf_flags_2_hdr = class1.psf_flags_hdr, class2.psf_flags_hdr
cm_flags_1_hdr, cm_flags_2_hdr = class1.cm_flags_hdr, class2.cm_flags_hdr
cm_max_flags_1_hdr, cm_max_flags_2_hdr = class1.cm_max_flags_hdr, class2.cm_max_flags_hdr
cm_flags_r_1_hdr, cm_flags_r_2_hdr = class1.cm_flags_r_hdr, class2.cm_flags_r_hdr
cm_mof_flags_1_hdr, cm_mof_flags_2_hdr = class1.cm_mof_flags_hdr, class2.cm_mof_flags_hdr

# For quality cuts introduced by Eric Huff #
cm_s2n_r_hdr1, cm_s2n_r_hdr2 = class1.cm_s2n_r_hdr, class2.cm_s2n_r_hdr
psfrec_t_hdr1, psfrec_t_hdr2 = class1.psfrec_t_hdr, class2.psfrec_t_hdr

flaglist = [ flags_1_hdr, flags_2_hdr, cm_flags_1_hdr, cm_flags_2_hdr, cm_mof_flags_1_hdr, cm_mof_flags_2_hdr, obj_flags_1_hdr, obj_flags_2_hdr, psf_flags_1_hdr, psf_flags_2_hdr, cm_max_flags_1_hdr, cm_max_flags_2_hdr, cm_flags_r_1_hdr, cm_flags_r_2_hdr ]

# Used if LOG_FLAGS is True #
flag_idx = []


### For plot names, plot titles, log file names ###
title_piece1 = class1.title_piece
title_piece2 = class2.title_piece
MATCH_TYPE = match_type = get_match_type(title_piece1, title_piece2)


### Names for file directors (fd) ###
fd_flag_name, fd_mag_bins_name, fd_nop_name, fd_1sig_name = get_fd_names(outdir=outdir, title_piece1=title_piece1, title_piece2=title_piece2)

# Create log file for number of objects plotted (nop) #
fd_nop = open(fd_nop_name, 'w')
# Create log file for number of objects within 1-sigma #
fd_1sig = open(fd_1sig_name, 'w')
# Create log file for magnitude bins #
fd_mag_bins = open(fd_mag_bins_name, 'w')
fd_mag_bins.write('NUM_OBJS_IN_BIN, BIN_LHS, BIN_RHS, MEDIAN_HAXIS_MAG, MEDIAN_ERROR \n') #TODO will have logs for all griz bands - make this clear in file?
# Create log file for flags #
fd_flag = open(fd_flag_name, 'w')
fd_flag.write('TILE, FILTER, TYPE, REALIZATION, FLAG1_HEADER, FLAG2_HEADER, FLAG1_VALUE, FLAG2_VALUE, MAG1, MAG2 \n')
if LOG_FLAGS is False:
        fd_flag.write('Flags not logged because LOG_FLAGS is False.')









def matcher(basepath, outdir, realization_number, tile_name):
        """Match two catalogs on RA and DEC with a tolerance of 1 arcsecond via STILTS.

        Args:
                basepath (str) -- Path to catalogs to match. Ex: /data/des71.a/data/kuropat/sof_stars/y3v02/.
                outdir (str) -- Path to where matched catalogs are saved.
                realization_number (str or int) -- Currently allowed values: 0 1 2 3 4 5 6 7 8 9 depending on the basepath.
                tile_name (str) -- Currently allowed values: DES0347-5540  DES2329-5622  DES2357-6456 DES2247-4414 depending on the basepath.
        Returns:
                gal_outname (str) -- Name of catalog matched between galaxy truth catalog and mof catalog. Headers will have _1 appended for truth catalog and _2 appended for mof catalog.
                star_outname (str) -- Name of catalog matched between star truth catalog and mof catalog. Headers will have _1 appended for truth catalog and _2 appended for mo
    f catalog.
        """

        ### Get arguments to pass to ms_matcher ###

        # Input catalogs for STILTS #
        in1 = get_catalog(basepath=basepath, cat_type=MATCH_CAT1, inj=INJ1, realization_number=realization_number, tile_name=tile_name)
        in2 = get_catalog(basepath=basepath, cat_type=MATCH_CAT2, inj=INJ2, realization_number=realization_number, tile_name=tile_name)

        # Output catalog name for STILTS #
        outname = os.path.join(outdir, tile_name+'_'+realization_number+'_'+str(MATCH_TYPE)+'_match1and2.csv')

        # Overwrite matched catalog if it already exists? This must be a str. Allowed values: 'False' 'True'. #
        OVERWRITE = 'False'

	print '\nMatching ', in1, in2, '...\n'

	### Matching done in ms_matcher. Args: in1, in2, out, ra1, dec1, ra2, dec2, OVERWRITE ###
	# !!!!! Ensure that path to ms_matcher is correct #
        subprocess.call(['/data/des71.a/data/mspletts/balrog_validation_tests/scripts/ms_matcher', in1, in2, outname, ra1, dec1, ra2, dec2, OVERWRITE])

        return outname









def make_plots(m1_hdr, m2_hdr, ylow, yhigh):
	"""Makes plots.

	Args:
		m1_hdr, m2_hdr (str) -- Headers for magnitudes.
		ylow, yhigh (int or float) -- Limits for vertical axis. None is an allowed value.
	Returns:
		0
	"""

	for t in ALL_TILES:

		### For plotting all realizations at once, stacked ###
		if STACK_REALIZATIONS:

			# Filename #
			fn_stack = os.path.join(outdir, t+'_stacked_'+str(MATCH_TYPE)+'_match1and2.csv')
			
			# Check if stacked realization file already exists #
			if os.path.isfile(fn_stack):
				print 'Stacked realization catalog exists. Not overwriting ... \n'
				df = pd.read_csv(fn_stack)
			
			# Combine all realizations for one tile into a single catalog. Catalogs combined AFTER matching. #
			if os.path.isfile(fn_stack) is False:
				all_fn = []
				for r in ALL_REALIZATIONS:
					fn = matcher(basepath=basepath, outdir=outdir, realization_number=r, tile_name=t)
					all_fn.append(fn)

				print 'Stacking realizations. ', len(all_fn), 'files ...'
				df = pd.concat((pd.read_csv(fn) for fn in all_fn))
				print 'Stacking complete ... \n'

				# Save stacked catalog as DataFrame #
				df.to_csv(fn_stack, sep=',')
				print '-----> Saving as ', fn_stack

			# Name for plt.savefig() #
			fn = get_plot_save_name(outdir=outdir, realization_number='stacked_realizations', tile_name=t, title_piece1=title_piece1, title_piece2=title_piece2, ylow=ylow, yhigh=yhigh)
			# Title for plot #
			title = get_plot_suptitle(realization_number='stacked '+str(len(ALL_REALIZATIONS)), tile_name=t, title_piece1=title_piece1, title_piece2=title_piece2)


			### Handle star truth catalogs ###
			if MATCH_CAT1 == 'star_truth' or MATCH_CAT2 == 'star_truth':
				print 'Adding new column to matched csv ...'
				star_mag = get_star_mag(df=df)
				# 'mag_a' short for mag_all. New header must be of the form {base}_x where x is a single character because of the way m_axlabel is created from m_hdr #
				df.insert(len(df.columns), 'mag_a', star_mag)

			if MATCH_CAT1 == 'star_truth':
				m1_hdr = 'mag_a'
			if MATCH_CAT2 == 'star_truth':
				m2_hdr = 'mag_a'

			subplotter(cm_flux_hdr1=cm_flux_hdr1, cm_flux_hdr2=cm_flux_hdr2, cm_flux_cov_hdr1=cm_flux_cov_hdr1, cm_flux_cov_hdr2=cm_flux_cov_hdr2, cm_t_hdr1=cm_t_hdr1, cm_t_hdr2=cm_t_hdr2, cm_t_err_hdr1=cm_t_err_hdr1, cm_t_err_hdr2=cm_t_err_hdr2, cm_s2n_r_hdr1=cm_s2n_r_hdr1, cm_s2n_r_hdr2=cm_s2n_r_hdr2, psfrec_t_hdr1=psfrec_t_hdr1, psfrec_t_hdr2=psfrec_t_hdr2, df=df, fd_1sig=fd_1sig, fd_flag=fd_flag, fd_nop=fd_nop, flag_idx=flag_idx, flag_list=flaglist, mag_hdr1=m1_hdr, mag_hdr2=m2_hdr, mag_err_hdr1=mag_err_hdr1, mag_err_hdr2=mag_err_hdr2, plot_name=fn, plot_title=title, realization_number='stacked', run_type=None, tile_name=t, swap_hax=SWAP_HAX, stack_realizations=STACK_REALIZATIONS, cm_t_s2n_axlabel1=cm_t_s2n_axlabel1, cm_t_s2n_axlabel2=cm_t_s2n_axlabel2, axlabel1=axlabel1, axlabel2=axlabel2, fd_mag_bins=fd_mag_bins, ylow=ylow, yhigh=yhigh)




	### For plotting one realization at a time ###
        if STACK_REALIZATIONS is False:

		print 'Not stacking realizations...'

                for r in ALL_REALIZATIONS:
			
			# Filename #
                        fn = matcher(basepath=basepath, outdir=outdir, realization_number=r, tile_name=t)
			# DataFrame #
                        df = pd.read_csv(fn)

                        # Name for plt.savefig() #
                        fn = get_plot_save_name(outdir=outdir, realization_number=r, tile_name=t, title_piece1=title_piece1, title_piece2=title_piece2, ylow=ylow, yhigh=yhigh)
                        # Title for plot #
                        title = get_plot_suptitle(realization_number=r, tile_name=t, title_piece1=title_piece1, title_piece2=title_piece2)


                        if MATCH_CAT1 == 'star_truth' or MATCH_CAT2 == 'star_truth':
                                print 'Adding new column to matched csv ... \n'
                                star_mag = get_star_mag(df=df)
                                df.insert(len(df.columns), 'mag_a', star_mag)

                        if MATCH_CAT1 == 'star_truth':
                                m1_hdr = 'mag_a'
                        if MATCH_CAT2 == 'star_truth':
                                m2_hdr = 'mag_a'


			subplotter(cm_flux_hdr1=cm_flux_hdr1, cm_flux_hdr2=cm_flux_hdr2, cm_flux_cov_hdr1=cm_flux_cov_hdr1, cm_flux_cov_hdr2=cm_flux_cov_hdr2, cm_t_hdr1=cm_t_hdr1, cm_t_hdr2=cm_t_hdr2, cm_t_err_hdr1=cm_t_err_hdr1, cm_t_err_hdr2=cm_t_err_hdr2, cm_s2n_r_hdr1=cm_s2n_r_hdr1, cm_s2n_r_hdr2=cm_s2n_r_hdr2, psfrec_t_hdr1=psfrec_t_hdr1, psfrec_t_hdr2=psfrec_t_hdr2, df=df, fd_1sig=fd_1sig, fd_flag=fd_flag, fd_nop=fd_nop, flag_idx=flag_idx, flag_list=flaglist, mag_hdr1=m1_hdr, mag_hdr2=m2_hdr, mag_err_hdr1=mag_err_hdr1, mag_err_hdr2=mag_err_hdr2, plot_name=fn, plot_title=title, realization_number=r, run_type=None, tile_name=t, swap_hax=SWAP_HAX, stack_realizations=STACK_REALIZATIONS, cm_t_s2n_axlabel1=cm_t_s2n_axlabel1, cm_t_s2n_axlabel2=cm_t_s2n_axlabel2, axlabel1=axlabel1, axlabel2=axlabel2, fd_mag_bins=fd_mag_bins, ylow=ylow, yhigh=yhigh)


	return 0









################################################################### Run script. 0 returned when complete. ###################################################################

### Run once ###
print make_plots(m1_hdr=m1_hdr, m2_hdr=m2_hdr, ylow=YLOW, yhigh=YHIGH)



### Loop over vertical axis limits. Suggestions: for normalized plot with star truth catalog use y=[3, 10], for normalized plot with galaxy truth catalog use y=[3, 20]. For non-normalized plot with star truth catalog or galaxy truth catalog use y=[0.5, None]. ###
YLOOP = False

if YLOOP:
	for y in [0.5, None]:
		if y is None:
			YLOW, YHIGH = None, None
		if y is not None:
			YLOW, YHIGH = -1.0*y, y

		# Must pass constants as parameters here because used for plot name #
		print make_plots(m1_hdr=m1_hdr, m2_hdr=m2_hdr, ylow=YLOW, yhigh=YHIGH)



### Loop over possible colorbars ###
CBAR_LOOP = False

if CBAR_LOOP:
	# Possible combinations for colorbars. Only one True allowed at a time and NORMALIZE and HEXBIN must both be False. #
	cbar_bools_list =[[True, False, False, False, False], [False, True, False, False, False], [False, False, True, False, False]]
	for cbar_bools in cbar_bools_list:
		# Reset constants #
		CM_T_COLORBAR, CM_T_S2N_COLORBAR, CM_T_ERR_COLORBAR, NORMALIZE, HEXBIN = cbar_bools
		make_plots(m1_hdr=m1_hdr, m2_hdr=m2_hdr, ylow=YLOW, yhigh=YHIGH)



### Close log files ###
fd_1sig.close()
fd_flag.close()
fd_mag_bins.close()
fd_nop.close()

