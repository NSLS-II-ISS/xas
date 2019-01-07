
from .bin import bin
from .file_io import (load_dataset_from_files, create_file_header, validate_file_exists, validate_path_exists,
                      save_interpolated_df_as_file, save_binned_df_as_file)
from .interpolate import interpolate


def process(doc, db, draw_func = None):
    if 'experiment' in db[doc['run_start']].start.keys():
        uid = doc['run_start']
        if db[uid].start['experiment'] == 'fly_energy_scan':
            raw_df = load_dataset_from_files(db, uid)
            interpolated_df = interpolate(raw_df)
            e0 = find_e0(db,uid)
            comments = create_file_header(db,uid)
            validate_path_exists(db,uid)
            path_to_file = validate_file_exists(db, uid )
            print(f'>>>Path to file {path_to_file}')
            save_interpolated_df_as_file(path_to_file, interpolated_df, comments)
            if e0 >0:
                binned_df = bin(interpolated_df,e0)
                save_binned_df_as_file(path_to_file, binned_df, comments)
            else:
                print('Energy E0 is not defined')

            if draw_func is not None:
                draw_func(interpolated_df)


def interpolate(doc, db):
    if 'experiment' in db[doc['run_start']].start.keys():
        if db[doc['run_start']].start['experiment'] == 'fly_energy_scan':
            raw_df = load_dataset_from_files(db, doc['run_start'])
            interpolated_df = interpolate(raw_df)
            return interpolated_df