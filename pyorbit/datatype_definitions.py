accepted_extensions = ['.yaml', '.yml', '.conf', '.config', '.input', ]
extended_dataset_keylist = ['expanded_dataset', 'extended_dataset', 'dataset_2d', 'dataset_expanded', 'dataset_extended', '2d_dataset']
external_dataset_keylist = ['external_dataset', 'dataset_external', 'external']
# Trying to guess all the possible mistakes....
datatype_definition = {
    'radial_velocity': ['RV', 'RVs', 'rv', 'rvs'],
    'transit_time': ['transit_time', 'transit_times', 'Tcent', 'TCent', 'Tc', 'TC', 'T0', 'TT'],
    'transit_duration': ['transit_duration', 'transit_durations', 'Tdur', 'TDur', 'Td', 'TD', 'T_dur', 'T_Dur'],
    'astrometry': ['astrometry', 'Astrometry', 'AstroMetry', 'Astro', 'astro', 'AM', 'Gaia', 'gaia'],
    'H-alpha': ['H', 'HA', 'h', 'ha', 'Halpha', 'H-alpha', 'halpha', 'h-alpha'],
    'photometry': ['P', 'Ph', 'p', 'ph', 'PHOT', 'Phot', 'phot', 'Photometry', 'photometry'],
    'FWHM': ['FWHM', 'fwhm'],
    'BIS': ['BIS', 'bis'],
    'EWs': ['EWs', 'EW', 'ew', 'ews'],
    'activity': ['Activity', 'activity', 'act_index', 'Act_index', 'activity_index', 'Activity_index'],
    'Ca_HK': ['Ca', 'Cahk', 'CaHK', 'Ca_hk', 'Ca_HK',
              'logR', 'logRhk', 'logRHK', 'logR_hk', 'logR_HK',
              'log(R)', 'log(Rhk)', 'log(RHK)', 'log(R_hk)', 'log(R_HK)'],
    'S_index': ['S', 'S_index', 'Shk', 'SHK', 'S_HK', 'S_hk'],
    'CCF': ['CCF', 'CCFs', 'ccf', 'ccfs'],
    'orbitize': ['orbitize', 'Orbitize', 'ORBITIZE'],
}

skip_plot = ['orbtize']
skip_model_plot = ['astrometry']

activity_datatype = ['H-alpha', 'FWHM', 'BIS', 'EWs', 'activity', 'Ca_HK', 'S_index']
activity_noderivative = ['H-alpha', 'FWHM', 'EWs', 'Ca_HK', 'S_index']