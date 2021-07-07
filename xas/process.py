
from .bin import bin
from .file_io import (load_dataset_from_files, create_file_header, validate_file_exists, validate_path_exists,
                      save_interpolated_df_as_file, save_binned_df_as_file, find_e0, save_stepscan_as_file,
                      stepscan_remove_offsets, stepscan_normalize_xs, combine_xspress3_channels)
from .db_io import load_apb_dataset_from_db, translate_apb_dataset, load_apb_trig_dataset_from_db, load_xs3_dataset_from_db, load_pil100k_dataset_from_db
from .interpolate import interpolate

from .xas_logger import get_logger
import matplotlib.pyplot as plt
import numpy as np
import os
from isscloudtools.slack import slack_upload_image
from isstools.elements.cloud_plotting import generate_output_figures
from PIL import Image

def process_interpolate_bin(doc, db, draw_func_interp = None, draw_func_bin = None, cloud_dispatcher = None):
    logger = get_logger()
    if 'experiment' in db[doc['run_start']].start.keys():
        uid = doc['run_start']
        experiment =  db[uid].start['experiment']
        if  experiment.startswith('fly'):
            path_to_file = db[uid].start['interp_filename']
            e0 = find_e0(db,uid)
            comments = create_file_header(db,uid)
            validate_path_exists(db,uid)
            path_to_file = validate_file_exists(path_to_file, file_type = 'interp')

            # try:
            if experiment == 'fly_energy_scan':
                raw_dict = load_dataset_from_files(db, uid)
                key_base = 'i0'
            elif experiment == 'fly_energy_scan_apb':
                apb_df, energy_df, energy_offset = load_apb_dataset_from_db(db, uid)
                raw_dict = translate_apb_dataset(apb_df, energy_df, energy_offset)
                key_base = 'i0'
            elif experiment == 'fly_energy_scan_xs3':
                apb_df, energy_df, energy_offset = load_apb_dataset_from_db(db, uid)
                raw_dict = translate_apb_dataset(apb_df, energy_df, energy_offset)

                apb_trig_timestamps = load_apb_trig_dataset_from_db(db, uid)
                xs3_dict = load_xs3_dataset_from_db(db, uid, apb_trig_timestamps)

                raw_dict = {**raw_dict, **xs3_dict}
                key_base = 'CHAN1ROI1'
            elif experiment == 'fly_energy_scan_pil100k':
                # PLACEHOLDER !!!!
                apb_df, energy_df, energy_offset = load_apb_dataset_from_db(db, uid)
                raw_dict = translate_apb_dataset(apb_df, energy_df, energy_offset)
                apb_trig_timestamps = load_apb_trig_dataset_from_db(db, uid, use_fall=True, stream_name='apb_trigger_pil100k')
                pil100k_dict = load_pil100k_dataset_from_db(db, uid, apb_trig_timestamps)
                raw_dict = {**raw_dict, **pil100k_dict}
                key_base = 'pil100k_ROI1'

            logger.info(f'Loading file successful for UID {uid}/{path_to_file}')
        # except:
            logger.info(f'Loading file failed for UID {uid}/{path_to_file}')
            try:
                interpolated_df = interpolate(raw_dict, key_base = key_base)
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
            (path, extension) = os.path.splitext(path_to_file)
            path_to_binned = path + '.dat'

            try:
                if cloud_dispatcher is not None:
                    cloud_dispatcher.load_to_dropbox(path_to_binned)
                    logger.info(f'Sending data to the cloud successful for {path_to_binned}')
            #         #WIP workaround
            #         #channel = db[uid].start['slack_channel']
            #         #slack_service = cloud_dispatcher.slack_service
            #         #image_path = os.path.splitext(path_to_binned)[0] + '.png'
            #         #generate_output_figures(path_to_binned, image_path)
            #         #label = os.path.basename(path).split('.')[0]
            #         #slack_upload_image(slack_service,channel,image_path,label)
            #         #cloud_dispatcher.post_to_slack(path_to_binned ,db[uid].start['slack_channel'])
            #         logger.info(f'Sending data to the cloud successful for {path_to_binned}')
            except:
                logger.info(f'Sending data to the cloud failed for {path_to_binned}')

            # # if cloud_dispatcher is not None:
            # #     cloud_dispatcher.load_to_dropbox(path_to_binned)
            # #     cloud_dispatcher.post_to_slack(path_to_binned ,db[uid].start['slack_channel'])
            # #     logger.info(f'Sending data to the cloud successful for {path_to_binned}')
        elif experiment.startswith('step'):
            path_to_file = db[uid].start['interp_filename']
            validate_path_exists(db,uid)
            path_to_file = validate_file_exists(path_to_file, file_type = 'interp')
            comments = create_file_header(db, uid)
            df = stepscan_remove_offsets(db[uid])
            df = stepscan_normalize_xs(df)
            df = combine_xspress3_channels(df)
            # ghnfg
            save_stepscan_as_file(path_to_file, df, comments)


            ###########################################################
            # PUSHKAR PROCESSING
            ###########################################################
            #storing table as hdf
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


















def process_interpolate_bin_offline_pilatus_fly(uid, db, draw_func_interp = None, draw_func_bin = None, cloud_dispatcher = None):
    logger = get_logger()
    if 'experiment' in db[uid].start.keys():
        path_to_file = db[uid].start['interp_filename']
        path_to_file = path_to_file[:-4] + '_pil100k_processed' + path_to_file[-4:]

        e0 = find_e0(db,uid)
        comments = create_file_header(db,uid)
        validate_path_exists(db,uid)
        path_to_file = validate_file_exists(path_to_file, file_type = 'interp')
        try:
            apb_df, energy_df, energy_offset = load_apb_dataset_from_db(db, uid)
            raw_dict = translate_apb_dataset(apb_df, energy_df, energy_offset)

            apb_trig_timestamps = load_apb_trig_dataset_from_db(db, uid, use_fall=False)
            pil100k_dict = load_pil100k_dataset_from_db(db, uid, apb_trig_timestamps)

            raw_dict = {**raw_dict, **pil100k_dict}
            key_base = 'pil100k_ROI1'

            logger.info(f'Loading file successful for UID {uid}/{path_to_file}')
        except:
            logger.info(f'Loading file failed for UID {uid}/{path_to_file}')
        try:
        # return raw_dict
            interpolated_df = interpolate(raw_dict, key_base = key_base)

            logger.info(f'Interpolation successful for {path_to_file}')
            save_interpolated_df_as_file(path_to_file, interpolated_df, comments)
        except:
            logger.info(f'Interpolation failed for {path_to_file}')

        try:
            if e0 >0:
                binned_df = bin(interpolated_df, e0, skip_binning=True)
                logger.info(f'Binning successful for {path_to_file}')
                # return interpolated_df, binned_df
                save_binned_df_as_file(path_to_file, binned_df, comments)
                if draw_func_interp is not None:
                    draw_func_interp(interpolated_df)
                if draw_func_bin is not None:
                    draw_func_bin(binned_df, path_to_file)
            else:
                print('Energy E0 is not defined')
        except:
            logger.info(f'Binning failed for {path_to_file}')







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