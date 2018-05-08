import ConfigParser
import json

def parseconfigs(config_file):
    config = ConfigParser.SafeConfigParser()
    config.read(config_file)

    params = {}

    #Input & Output Files
    params['zp_files'] = config.get('I/O','zp_files')
    params['data_file'] = config.get('I/O','data_file')
    params['lens_file'] = config.get('I/O','lens_file')
    params['random_lens_file'] = config.get('I/O','random_lens_file')
    params['balrog_zp_files'] = config.get('I/O','balrog_zp_files')
    params['sim_file_format'] = config.get('I/O','sim_file_format')
    params['test_zp_file'] = config.get('I/O','test_zp_file')
    params['test_data_file'] = config.get('I/O','test_data_file')
    params['output'] = config.get('I/O','output')
    params['deep_file'] = config.get('I/O','deep_file')
    params['balrog_files'] = config.get('I/O','balrog_files')
    params['deep_output'] = config.get('I/O','deep_output')
    params['balrog_output'] = config.get('I/O','balrog_output')
    params['overwrite_fits'] = config.getboolean('I/O','overwrite_fits')
    
    #Table Column Names
    params['zp_id_col'] = config.get('Columns','zp_id_column')
    params['data_id_col'] = config.get('Columns','data_id_column')
    params['test_data_id_col'] = config.get('Columns','test_data_id_column')
    params['test_z_col'] = config.get('Columns','test_z_column')
    params['deep_z_column'] = config.get('Columns','deep_z_column')
    params['balrog_id_column'] = config.get('Columns','balrog_id_column')
    params['deep_type_column'] = config.get('Columns','deep_type_column')
    params['balrog_type_column'] = config.get('Columns','balrog_type_column')
    params['deep_flux_column'] = config.get('Columns','deep_flux_column')
    params['deep_size_column'] = config.get('Columns','deep_size_column')
    params['balrog_flux_column'] = config.get('Columns','balrog_flux_column')
    params['balrog_size_column'] = config.get('Columns','balrog_size_column')
    
    #Assignment Criteria
    params['redshift_index'] = json.loads(config.get('Assignments','redshift_index'))
    params['below'] = json.loads(config.get('Assignments','below'))
    params['above'] = json.loads(config.get('Assignments','above'))
    params['true_ranges'] = json.loads(config.get('Assignments','true_ranges'))

    #dNdMu
    params['mu'] = config.get('dNdMu','mu')
    params['flux_min'] = config.get('dNdMu','flux_min')
    params['flux_max'] = config.get('dNdMu','flux_max')

    return params
