
from .bin import bin
from .file_io import (load_dataset_from_files, create_file_header, validate_file_exists, validate_path_exists,
                      save_interpolated_df_as_file, save_binned_df_as_file, find_e0)
from .interpolate import interpolate

from .xas_logger import get_logger



def process_interpolate_bin(doc, db, draw_func_interp = None, draw_func_bin = None):
    logger = get_logger()
    if 'experiment' in db[doc['run_start']].start.keys():
        uid = doc['run_start']
        if db[uid].start['experiment'] == 'fly_energy_scan':
            path_to_file = db[uid].start['interp_filename']
            e0 = find_e0(db,uid)
            comments = create_file_header(db,uid)
            validate_path_exists(db,uid)

            path_to_file = validate_file_exists(path_to_file, file_type = 'interp')
            #print(f'>>>Path to file {path_to_file}')
            try:
                raw_df = load_dataset_from_files(db, uid)
                logger.info(f'Loading file successful for UID {uid}/{path_to_file}')
            except:
                logger.info(f'Loading file failed for UID {uid}/{path_to_file}')
            try:
                interpolated_df = interpolate(raw_df)
                logger.info(f'Interpolation successful for {path_to_file}')
                save_interpolated_df_as_file(path_to_file, interpolated_df, comments)
            except:
                logger.info(f'Interpolation failed for {path_to_file}')

            try:
                if e0 >0:
                    binned_df = bin(interpolated_df,e0)
                    logger.info(f'Binning successful for {path_to_file}')
                    save_binned_df_as_file(path_to_file, binned_df, comments)
                    if draw_func_interp is not None:
                        draw_func_interp(interpolated_df)
                    if draw_func_bin is not None:
                        draw_func_bin(binned_df, path_to_file)
                else:
                    print('Energy E0 is not defined')
            except:
                logger.info(f'Binning failed for {path_to_file}')

        elif db[uid].start['experiment'] == 'fly_energy_scan_em':
            logger.info('HAHAHAHHAHAH')





def process_interpolate_only(doc, db):
    if 'experiment' in db[doc['run_start']].start.keys():
        if db[doc['run_start']].start['experiment'] == 'fly_energy_scan':
            raw_df = load_dataset_from_files(db, doc['run_start'])
            interpolated_df = interpolate(raw_df)
            return interpolated_df


def process_interpolate_unsorted(uid, db):
     raw_df = load_dataset_from_files(db, uid)
     interpolated_df = interpolate(raw_df, sort=False)
     return interpolated_df