path = os.path.dirname(path_to_binned)
tiff_storage_path = os.path.join(path, 'tiff_storage')
filename = os.path.basename(path_to_binned)
rixs_processing_file=  os.path.join( tiff_storage_path, filename)
#cloud_dispatcher.load_to_dropbox(path_to_binned)

name = filename.split('.')[0]

tiff_directory = os.path.join( tiff_storage_path, name)
zip_file =
os.system(f"zip eli.zip '{tiff_directory}/*.tif'")


from xas.file_io import load_extended_data_from_file
from PIL import Image

year = '2022'
cycle = '2'
proposal = '309855'
cloud_dispatcher = xlive_gui.cloud_dispatcher

for uid in range(265677, 265858+1):
     hdr = db[uid]
     if 'spectrometer' in hdr.start.keys():
         path_to_file = hdr.start['interp_filename']
         path, f_ext = os.path.split(hdr.start['interp_filename'])
         f, _ = os.path.splitext(f_ext)

         tiff_storage_path = os.path.join(os.path.dirname(path_to_file), 'tiff_storage') + '/'
         scan_name, _ = os.path.splitext(os.path.basename(path_to_file))
         dat_file_fpath = tiff_storage_path + scan_name + '.dat'
         tiff_images_path = tiff_storage_path + scan_name + '/'
         extended_data = load_extended_data_from_file(os.path.join(os.path.dirname(path_to_file), 'extended_data', f + '.h5'))
         if extended_data is not None:
             filepath_list = []
             filename_list = []

             for i, im in enumerate(extended_data['pil100k_image']):
                 image_data = Image.fromarray(im)
                 #
                 tiff_filename = '{}{:04d}{}'.format('image', i + 1, '.tif')
                 tiff_path = tiff_images_path + tiff_filename
                 print(f'TIFF STORAGE: tiff will be saved in {tiff_path}')
                 image_data.save(tiff_path)
                 filepath_list.append(tiff_path)
                 filename_list.append(tiff_filename)

             if zip:
                 # os.system(f"cd '{tiff_images_path}'")
                 zip_file = os.path.splitext(dat_file_fpath)[0] + '-1.zip'
                 # os.system(f"zip '{zip_file}' ./*.tif")
                 if os.path.exists(zip_file):
                     os.remove(zip_file)
                 os.system(f"cd '{tiff_images_path}'; zip '{zip_file}' ./*.tif")
                 filepath_list = [zip_file]

                 for f in filepath_list:
                     cloud_dispatcher.load_to_dropbox(f, year=year, cycle=cycle, proposal=proposal)
                     # logger.info(f'({ttime.ctime()}) Sending data to the cloud successful for {path_to_file}')


             # print(os.path.join(path, 'tiff_storage', f))
             # ext_data = load_extended_data_from_file





