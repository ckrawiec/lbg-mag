import ConfigParser
import json

def parseconfigs(config_file):
    config = ConfigParser.SafeConfigParser()
    config.read(config_file)

    params = {}

    #Input & Output Files
    params['zp_files'] = config.get('I/O','zp_files')
    params['source_file'] = config.get('I/O','source_file')
    params['lens_file'] = config.get('I/O','lens_file')
    params['random_lens_file'] = config.get('I/O','random_lens_file')
    
    params['balrog_zp_files'] = config.get('I/O','balrog_zp_files')
    params['balrog_truth_files'] = config.get('I/O','balrog_truth_files')
    params['balrog_sim_file_format'] = config.get('I/O','balrog_sim_file_format')
    
    params['test_zp_file'] = config.get('I/O','test_zp_file')
    params['test_data_file'] = config.get('I/O','test_data_file')
    params['deep_file'] = config.get('I/O','deep_file')    

    params['deep_output'] = config.get('I/O','deep_output')
    params['balrog_output'] = config.get('I/O','balrog_output')
    params['output'] = config.get('I/O','output')
    
    if config.has_option('I/O', 'overwrite_type_files'):
        params['overwrite_type_files'] = config.getboolean('I/O','overwrite_type_files')
    else:
        params['overwrite_type_files'] = True
    if config.has_option('I/O', 'overwrite_corr_files'):
        params['overwrite_corr_files'] = config.getboolean('I/O','overwrite_corr_files')
    else:
        params['overwrite_corr_files'] = True
    
    #Table Column Names
    params['zp_id_column'] = config.get('Columns','zp_id_column')
    params['data_id_column'] = config.get('Columns','data_id_column')    
    params['test_data_id_column'] = config.get('Columns','test_data_id_column')
    params['balrog_id_column'] = config.get('Columns','balrog_id_column')
    
    params['test_z_column'] = config.get('Columns','test_z_column')
    params['deep_z_column'] = config.get('Columns','deep_z_column')
    
    params['deep_type_column'] = config.get('Columns','deep_type_column')
    params['balrog_type_column'] = config.get('Columns','balrog_type_column')
    
    params['deep_flux_column'] = config.get('Columns','deep_flux_column')
    params['balrog_flux_column'] = config.get('Columns','balrog_flux_column')
    
    params['deep_size_column'] = config.get('Columns','deep_size_column')
    params['balrog_size_column'] = config.get('Columns','balrog_size_column')

    params['test_mag_column'] = config.get('Columns','test_mag_column')
    
    #Assignment Criteria
    params['redshift_index'] = json.loads(config.get('Assignments','redshift_index'))
    params['below'] = json.loads(config.get('Assignments','below'))
    params['above'] = json.loads(config.get('Assignments','above'))
    params['true_ranges'] = json.loads(config.get('Assignments','true_ranges'))

    #dNdMu
    params['mu'] = config.getfloat('dNdMu','mu')
    params['flux_min'] = config.getfloat('dNdMu','flux_min')
    params['flux_max'] = config.getfloat('dNdMu','flux_max')

    return params
