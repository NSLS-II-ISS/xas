
from .bin import bin
from .file_io import (load_dataset_from_files, create_file_header, validate_file_exists, validate_path_exists,
                      save_interpolated_df_as_file, save_binned_df_as_file, find_e0, save_stepscan_as_file,
                      stepscan_remove_offsets, stepscan_normalize_xs, combine_xspress3_channels, filter_df_by_valid_keys,
                      save_primary_df_as_file, save_extended_data_as_file, dump_tiff_images)
from .db_io import load_apb_dataset_from_db, translate_apb_dataset, load_apb_trig_dataset_from_db, load_xs3_dataset_from_db, load_pil100k_dataset_from_db
from .interpolate import interpolate

from .xas_logger import get_logger
# import matplotlib.pyplot as plt
# import numpy as np
import time as ttime
import os
# from isscloudtools.slack import slack_upload_image
# from isscloudtools.cloud_dispatcher import generate_output_figures

from .vonhamos import process_von_hamos_scan #, save_vh_scan_to_file
import gc

def process_interpolate_bin(doc, db, draw_func_interp = None, draw_func_bin = None, cloud_dispatcher = None, print_func=None, dump_to_tiff=True):
    # logger = get_logger()
    if 'experiment' in db[doc['run_start']].start.keys():
        uid = doc['run_start']
        process_interpolate_bin_from_uid(uid, db, draw_func_interp=draw_func_interp, draw_func_bin=draw_func_bin, cloud_dispatcher=cloud_dispatcher, print_func=print_func, dump_to_tiff=dump_to_tiff)


def process_interpolate_bin_from_uid(uid, db, draw_func_interp = None, draw_func_bin = None, cloud_dispatcher = None, print_func=None, dump_to_tiff=False):

    logger = get_logger()
    hdr, primary_df, extended_data, comments, path_to_file, file_list, data_kind = get_processed_df_from_uid(uid, db,
                                                                          logger=logger,
                                                                          draw_func_interp=draw_func_interp,
                                                                          draw_func_bin=draw_func_bin,
                                                                          print_func=print_func)

    save_primary_df_as_file(path_to_file, primary_df, comments)

    try:
        save_extended_data_as_file(path_to_file, extended_data, data_kind=data_kind)
    except Exception  as e:
        print(e)
        pass

    if dump_to_tiff:
        if extended_data is not None:
            tiff_files = dump_tiff_images(path_to_file, primary_df, extended_data)
            print(f' >>>>>>>>>>Tiff file {tiff_files}')
            file_list += tiff_files


    # dfgd
    try:
        if cloud_dispatcher is not None:
            year, cycle, proposal = hdr.start['year'], hdr.start['cycle'], hdr.start['PROPOSAL']
            for f in file_list:
                cloud_dispatcher.load_to_dropbox(f, year=year, cycle=cycle, proposal=proposal)
                logger.info(f'({ttime.ctime()}) Sending data to the cloud successful for {path_to_file}')
    except Exception as e:
        logger.info(f'({ttime.ctime()}) Sending data to the cloud failed for {path_to_file}')
        raise e



    clear_db_cache(db)


def get_processed_df_from_uid(uid, db, logger=None, draw_func_interp=None, draw_func_bin = None, print_func=None, save_interpolated_file=True):
    if print_func is None:
        print_func = print
    if logger is None:
        logger = get_logger()

        # logger = get_logger(print_func=print_func)
    hdr = db[uid]
    experiment = hdr.start['experiment']
    comments = create_file_header(hdr)
    path_to_file = hdr.start['interp_filename']
    validate_path_exists(path_to_file)
    path_to_file = validate_file_exists(path_to_file, file_type='interp')
    e0 = find_e0(hdr)
    data_kind = 'default'
    file_list = []

    if experiment == 'fly_scan':

        # path_to_file = validate_file_exists(path_to_file, file_type='interp')
        stream_names = hdr.stream_names
        try:
            # default detectors
            apb_df, energy_df, energy_offset = load_apb_dataset_from_db(db, uid)
            raw_dict = translate_apb_dataset(apb_df, energy_df, energy_offset)

            for stream_name in stream_names:
                if stream_name == 'pil100k_stream':
                    apb_trigger_pil100k_timestamps = load_apb_trig_dataset_from_db(db, uid, use_fall=True,
                                                                                   stream_name='apb_trigger_pil100k')
                    pil100k_dict = load_pil100k_dataset_from_db(db, uid, apb_trigger_pil100k_timestamps)
                    raw_dict = {**raw_dict, **pil100k_dict}

                elif stream_name == 'xs_stream':
                    apb_trigger_xs_timestamps = load_apb_trig_dataset_from_db(db, uid, stream_name='apb_trigger_xs')
                    xs3_dict = load_xs3_dataset_from_db(db, uid, apb_trigger_xs_timestamps)
                    raw_dict = {**raw_dict, **xs3_dict}

            logger.info(f'({ttime.ctime()}) Loading file successful for UID {uid}/{path_to_file}')
        except Exception as e:
            logger.info(f'({ttime.ctime()}) Loading file failed for UID {uid}/{path_to_file}')
            raise e
        try:
            interpolated_df = interpolate(raw_dict)
            logger.info(f'({ttime.ctime()}) Interpolation successful for {path_to_file}')
            if save_interpolated_file:
                save_interpolated_df_as_file(path_to_file, interpolated_df, comments)
        except Exception as e:
            logger.info(f'({ttime.ctime()}) Interpolation failed for {path_to_file}')
            raise e

        try:
            if e0 > 0:
                processed_df = bin(interpolated_df, e0)
                (path, extension) = os.path.splitext(path_to_file)
                path_to_file = path + '.dat'
                logger.info(f'({ttime.ctime()}) Binning successful for {path_to_file}')

                if draw_func_interp is not None:
                    draw_func_interp(interpolated_df)
                if draw_func_bin is not None:
                    draw_func_bin(processed_df, path_to_file)
            else:
                print(f'({ttime.ctime()}) Energy E0 is not defined')
        except Exception as e:
            logger.info(f'({ttime.ctime()}) Binning failed for {path_to_file}')
            raise e

        # save_binned_df_as_file(path_to_file, processed_df, comments)


    elif (experiment == 'step_scan') or (experiment == 'collect_n_exposures'):
        # path_to_file = validate_file_exists(path_to_file, file_type='interp')
        df = stepscan_remove_offsets(hdr)
        df = stepscan_normalize_xs(df)
        processed_df = filter_df_by_valid_keys(df)
        # df_processed = combine_xspress3_channels(df)

    else:
        return

    processed_df = combine_xspress3_channels(processed_df)

    primary_df, extended_data = split_df_data_into_primary_and_extended(processed_df)

    ### WIP
    if 'spectrometer' in hdr.start.keys():
        if hdr.start['spectrometer'] == 'von_hamos':
            extended_data, comments, file_paths = process_von_hamos_scan(primary_df, extended_data, comments, hdr, path_to_file, db=db)
            data_kind = 'von_hamos'
            file_list = file_paths
            # save_vh_scan_to_file(path_to_file, vh_scan, comments)
    file_list.append(path_to_file)
    return hdr, primary_df, extended_data, comments, path_to_file, file_list, data_kind

import numpy as np
def split_df_data_into_primary_and_extended(df_orig):
    sec_cols = []
    for c in df_orig.columns:
        if len(df_orig[c][0].shape) > 0:
            sec_cols.append(c)

    if len(sec_cols) > 0:
        extended_data = {}
        df = df_orig.copy()
        for col in sec_cols:
            v = df.pop(col)
            extended_data[col] = np.array([i for i in v])

    else:
        df = df_orig
        extended_data = None
    return df, extended_data


def clear_db_cache(db):
    db._catalog._entries.cache_clear()
    gc.collect()



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



################


# ###########################################################
# # LOCKARD PROCESSING
# ###########################################################
# # storing table as hdf
# dump_to_hdf = True
# if dump_to_hdf:
#     hdf_fn = os.path.splitext(path_to_file)[0]+'.h5'
#     # hdr = db[uid]
#     t = db[uid].table(fill=True)
#     t.to_hdf(hdf_fn, 'data')
# dump_to_tiff_pushkar = True
# if dump_to_tiff_pushkar:
#     # deal with paths
#     tiff_storage_path = os.path.dirname(path_to_file) + '/tiff_storage/'
#     scan_name = os.path.splitext(os.path.basename(path_to_file))[0]
#     dat_file_fpath = tiff_storage_path + scan_name + '.dat'
#     tiff_images_path = tiff_storage_path + scan_name + '/'
# #
#     try:
#         os.mkdir(tiff_images_path)
#     except FileExistsError:
#         print('Warning Folder exists')
#         return
# #
# #     # deal with images
#     filename_list = []
#     for i, im in enumerate(t['pil100k_image']):
#         image_data = Image.fromarray(im)
# #
#         tiff_filename = '{}{:04d}{}'.format('image', i + 1, '.tiff')
#         tiff_path = tiff_images_path + tiff_filename
#         print(f'TIFF STORAGE: tiff will be saved in {tiff_path}')
#         image_data.save(tiff_path)
#         filename_list.append(tiff_filename)
# #
# #     # deal with table file
#     table_red = df[['hhm_energy', 'apb_ave_ch1_mean', 'apb_ave_ch2_mean', 'apb_ave_ch3_mean', 'apb_ave_ch4_mean']]
#
#     table_red = table_red.rename(
#         columns={'hhm_energy': '# energy', 'apb_ave_ch1': 'i0', 'apb_ave_ch2': 'it', 'apb_ave_ch3': 'ir',
#                  'apb_ave_ch4': 'iff'})
#     # table_red = df
#     table_red['filenames'] = filename_list
#     print(f'TIFF STORAGE: dat will be saved in {dat_file_fpath}')
#     table_red.to_csv(dat_file_fpath, sep='\t', index=False)

###########################################################
# PUSHKAR PROCESSING
###########################################################
# storing table as hdf
# dump_to_hdf = True
# if dump_to_hdf:
#     hdf_fn = os.path.splitext(path_to_file)[0]+'.h5'
#     hdr = db[uid]
#     t = db[uid].table(fill=True)
#     t.to_hdf(hdf_fn, 'data')
# dump_to_tiff_pushkar = True
# if dump_to_tiff_pushkar:
#     # deal with paths
#     tiff_storage_path = os.path.dirname(path_to_file) + '/tiff_storage/'
#     scan_name = os.path.splitext(os.path.basename(path_to_file))[0]
#     dat_file_fpath = tiff_storage_path + scan_name + '.dat'
#     tiff_images_path = tiff_storage_path + scan_name + '/'
#
#     try:
#         os.mkdir(tiff_images_path)
#     except FileExistsError:
#         print('Warning Folder exists')
#         return
#
#     # deal with images
#     # filename_list = []
#     for i, im in enumerate(t['pil100k_image']):
#         image_data = Image.fromarray(im)
#
#         tiff_filename = '{}{:04d}{}'.format('image', i + 1, '.tif')
#         tiff_path = tiff_images_path + tiff_filename
#         print(f'TIFF STORAGE: tiff will be saved in {tiff_path}')
#         image_data.save(tiff_path)
#         # filename_list.append(tiff_filename)
#
#     # deal with table file
#     table_red = t[['hhm_energy', 'apb_ave_ch1', 'apb_ave_ch2', 'apb_ave_ch3', 'apb_ave_ch4']]
#
#
#     table_red = table_red.rename(
#         columns={'hhm_energy': '# energy', 'apb_ave_ch1': 'i0', 'apb_ave_ch2': 'it', 'apb_ave_ch3': 'ir',
#                  'apb_ave_ch4': 'iff'})
#     # table_red['filenames'] = filename_list
#     print(f'TIFF STORAGE: dat will be saved in {dat_file_fpath}')
#     table_red.to_csv(dat_file_fpath, sep='\t', index=False)


# scratch

# import math
# from xas.math import gauss
# from scipy.optimize import curve_fit
#
#
# def get_image_data():
#     image = bpm_es.image.array_data.read()['bpm_es_image_array_data']['value'].reshape((960,1280))
#     image = image.astype(np.int16)
#     sum_lines = sum(image[:, [i for i in range(line - math.floor(n_lines/2), line + math.ceil(n_lines/2))]].transpose())
#     sum_lines = sum_lines - np.mean(sum_lines[:200])
#     index_max = sum_lines.argmax()
#     max_value = sum_lines.max()
#     min_value = sum_lines.min()
#     sum_lines_t = sum_lines.copy()
#     mask = sum_lines_t>max_value/2
#     idx_to_fit = np.where(mask)
#     sum_lines_t[~mask] = 0
#     x = np.arange(960)
#     coeff, var_matrix = curve_fit(gauss, x[idx_to_fit], sum_lines_t[idx_to_fit], p0=[1, index_max, 5])
#     return sum_lines, sum_lines_t, gauss(x, *coeff)
#
# data1, data1_t, fit1 = get_image_data()
#
# plt.figure()
# plt.plot(data1, label = 'Data 1')
# plt.plot(data1_t, label = 'Data 1 t')
# plt.plot(fit1, label = 'Fit 1')
# plt.legend()
# plt.grid()
#
#
# data2, data1_2, fit2 = get_image_data()
#
# plt.plot(data1, 'k-', label = 'Data 1')
# plt.plot(data1_t, 'k--', label = 'Data 1 t')
# plt.plot(fit1, 'r-', label = 'Fit 1')