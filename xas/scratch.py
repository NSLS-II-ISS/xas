import numpy as np

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

###########################
# Yet another attempt at combining the rixs planes

x = xview_gui.widget_data

rixs_dict = {}


for item in x.list_data.selectedItems():
    fname = os.path.join(x.working_folder, item.text())
    df, header = load_binned_df_from_file(fname)
    print(fname, df.shape)
    uid_idx1 = header.find('Scan.uid:') + 10
    uid_idx2 = header.find('\n', header.find('Scan.uid:'))
    uid = header[uid_idx1: uid_idx2]
    md = db[uid]['start']
    emission_energy = round(md['spectrometer_energy'], 2)
    if emission_energy not in rixs_dict.keys():
        rixs_dict[emission_energy] = [{'energy': df.energy, 'mu' : -df.pil100k_roi1 / df.i0}]
    else:
        rixs_dict[emission_energy].append({'energy': df.energy, 'mu' : -df.pil100k_roi1 / df.i0})



some_emission_energy = list(rixs_dict.keys())[0]
energy_in = rixs_dict[some_emission_energy][0]['energy'].values
energy_out = np.sort(list(rixs_dict.keys()))
rixs_data = np.zeros((energy_in.size, energy_out.size))

for i, each_energy_out in enumerate(energy_out):
    each_mu_set = []
    for ds in rixs_dict[each_energy_out]:

        each_mu_set.append(np.interp(energy_in, ds['energy'], ds['mu']))
    rixs_data[:, i] = np.mean(np.array(each_mu_set), axis=0)


ds_norm = xview_gui.project[0]
ds_norm_ka2 = xview_gui.project[11]

xes_norm = np.hstack((ds_norm_ka2.mu[:-2] / np.mean(ds_norm_ka2.mu[-2:]) * np.mean(ds_norm.mu[:2]), ds_norm.mu))

# plt.figure()
# plt.plot(ds_norm_ka2.energy[-2:], ds_norm_ka2.mu[-2:])
# plt.plot(ds_norm.energy[:2], ds_norm.mu[:2])
#
# plt.plot(ds_norm_ka2.energy, ds_norm_ka2.mu / np.mean(ds_norm_ka2.mu[-2:]) * np.mean(ds_norm.mu[:2]))
# plt.plot(ds_norm.energy, ds_norm.mu)

# xes_norm = ds_norm.mu
xes_norm /= xes_norm.max()
rixs_norm = np.mean(rixs_data[energy_in > 8350, :], axis=0)

# xes_8347 = xview_gui.project[1].mu

rixs_data_norm = rixs_data / rixs_norm[None, :] * xes_norm[None, :]
rixs_data_norm /= rixs_data_norm.max()

from xas.spectrometer import convert_rixs_to_energy_transfer

energy_et, rixs_et = convert_rixs_to_energy_transfer(energy_in, energy_out, rixs_data_norm.T)

##
plt.figure(1, clear=True)
plt.plot(energy_in, rixs_data_norm[:, 22]/0.823*1.2, label='RIXS')
plt.plot(herfd_data_energy, herfd_data_ni4, label='HERFD')
plt.legend()
plt.xlim(8325, 8355)

##

plt.figure(1, clear=True)
# plt.contourf(energy_in, energy_out, rixs_data_norm.T, 150, vmin=0, vmax=0.3)
plt.contourf(energy_in, energy_out, np.log10(rixs_data_norm.T), 250, vmin=np.log10(1e-3), vmax=np.log10(1))


plt.figure(2, clear=True)
# plt.contourf(energy_in, energy_out, rixs_data_norm.T, 150, vmin=0, vmax=0.3)
plt.contourf(energy_in, energy_et, np.log10(rixs_et.T), 250, vmin=np.log10(1e-3), vmax=np.log10(1))
# plt.contourf(energy_in, energy_et, rixs_et.T, 250, vmin=0, vmax=0.5)

plt.axis('square')

plt.figure(3, clear=True)
# plt.imshow(np.log10(rixs_et.T), vmin=np.log10(3e-3), vmax=np.log10(1))
plt.plot(energy_et, np.mean(rixs_et[48:52, :], axis=0))
plt.plot(energy_et, np.mean(rixs_et[62:66, :], axis=0))
plt.plot(energy_et, np.mean(rixs_et[78:82, :], axis=0))

plt.figure(4, clear=True)
# plt.plot(energy_out, xes_8347 / xes_8347.max())
plt.plot(energy_out, rixs_data_norm[115, :] / rixs_data_norm[115, :].max())

def save_data_to_files(energy_in, energy_out, rixs_norm, energy_herfd, mu_herfd, filename):
    rixs_matrix = np.block([[0, energy_out[None, :]], [energy_in[:, None], rixs_norm]])
    herfd_matrix = np.vstack((energy_herfd, mu_herfd)).T
    rixs_filename = os.path.join(x.working_folder, f'{filename}_rixs.dat')
    herfd_filename = os.path.join(x.working_folder, f'{filename}_herfd.dat')
    np.savetxt(rixs_filename, rixs_matrix, delimiter='\t')
    np.savetxt(herfd_filename, herfd_matrix, delimiter='\t', header='# energy (eV)\tmu (a.u.)')
    xview_gui.cloud_dispatcher.load_to_dropbox(rixs_filename, year='2022', cycle='2', proposal='310156')
    xview_gui.cloud_dispatcher.load_to_dropbox(herfd_filename, year='2022', cycle='2', proposal='310156')

save_data_to_files(energy_in, energy_out, rixs_data_norm, xview_gui.project[3].energy, xview_gui.project[3].mu, 'DW_NiIII')



herfd_data_energy = xview_gui.project[2].energy
herfd_data_ni4 = xview_gui.project[2].flat
herfd_data_ni3 = xview_gui.project[3].flat
herfd_preedge_ni4 = xview_gui.project[2].pre_edge
herfd_preedge_ni3 = xview_gui.project[3].pre_edge
herfd_postedge_ni4 = xview_gui.project[2].post_edge
herfd_postedge_ni3 = xview_gui.project[3].post_edge

n_out = 22
n_in = 48

# basis_ni4 = np.vstack((rixs_data_norm_ni4[:, n_out], np.ones(energy_in.size))).T
# basis_ni3 = np.vstack((rixs_data_norm_ni3[:, n_out], np.ones(energy_in.size))).T

rixs_data_norm_ni4_bkg_rem = rixs_data_norm_ni4 - np.mean(rixs_data_norm_ni4[energy_in<8325, :], axis=0)[None, :]
rixs_data_norm_ni3_bkg_rem = rixs_data_norm_ni3 - np.mean(rixs_data_norm_ni3[energy_in<8325, :], axis=0)[None, :]

rixs_data_norm_ni4_bkg_rem /= np.interp(energy_in, herfd_data_energy, herfd_postedge_ni4)[:, None]
rixs_data_norm_ni3_bkg_rem /= np.interp(energy_in, herfd_data_energy, herfd_postedge_ni3)[:, None]

scale_ni4, _, _, _ = np.linalg.lstsq(rixs_data_norm_ni4_bkg_rem[:, n_out][:, None],
                                     np.interp(energy_in, herfd_data_energy, herfd_data_ni4), rcond=-1)
scale_ni3, _, _, _ = np.linalg.lstsq(rixs_data_norm_ni3_bkg_rem[:, n_out][:, None],
                                     np.interp(energy_in, herfd_data_energy, herfd_data_ni3), rcond=-1)

rixs_data_norm_scaled_ni4 = rixs_data_norm_ni4_bkg_rem * scale_ni4
rixs_data_norm_scaled_ni3 = rixs_data_norm_ni3_bkg_rem * scale_ni3


plt.figure(5, clear=True)
plt.subplot(221)
plt.contourf(energy_in, energy_out, np.log10(rixs_data_norm_scaled_ni4.T + 1e-3), 250, vmin=np.log10(1e-3), vmax=np.log10(1))
# plt.contourf(energy_in, energy_out, rixs_data_norm_scaled_ni4.T, 250, vmin=0, vmax=1.05, cmap='cubehelix')
plt.title('Ni 4+')
plt.hlines([energy_out[n_out]], energy_in.min(), energy_in.max(), color='w')
plt.vlines([energy_in[n_in]], energy_out.min(), energy_out.max(), color='r')

plt.subplot(222)
plt.contourf(energy_in, energy_out, np.log10(rixs_data_norm_scaled_ni3.T + 1e-3), 250, vmin=np.log10(1e-3), vmax=np.log10(1))
plt.title('Ni 3+')
plt.hlines([energy_out[n_out]], energy_in.min(), energy_in.max(), color='w')
plt.vlines([energy_in[n_in]], energy_out.min(), energy_out.max(), color='r')

plt.subplot(223)
plt.plot(energy_in, rixs_data_norm_scaled_ni4[:, n_out], label='Ni 4+')
plt.plot(herfd_data_energy, herfd_data_ni4, label='Ni 4+ HERFD')

plt.plot(energy_in, rixs_data_norm_scaled_ni3[:, n_out], label='Ni 3+')
plt.plot(herfd_data_energy, herfd_data_ni3, label='Ni 3+ HERFD')
plt.legend()

plt.xlim(8330, 8360)

plt.subplot(224)
plt.plot(energy_out, rixs_data_norm_ni4[n_in, :] * scale_ni4, label='Ni 4+')
plt.plot(energy_out, rixs_data_norm_ni3[n_in, :] * scale_ni3, label='Ni 3+')
plt.legend()


__basis = np.vstack((rixs_data_norm_scaled_ni3[:, n_out],
                     np.interp(energy_in, herfd_data_energy, herfd_data_ni4))).T
__c, _, _, _ = np.linalg.lstsq(__basis, rixs_data_norm_scaled_ni4[:, n_out], rcond=-1)

plt.figure(6, clear=True)

plt.plot(energy_in, rixs_data_norm_scaled_ni4[:, n_out], label='Ni 4+')
plt.plot(energy_in, __basis @ __c, 'r:', label='Ni 4+ fit')

plt.plot(energy_in, rixs_data_norm_scaled_ni3[:, n_out], label='Ni 3+')
plt.plot(herfd_data_energy, herfd_data_ni4, label='Ni 4+ HERFD')


# __alpha = 0.35
# plt.plot(energy_in, (rixs_data_norm_scaled_ni4[:, n_out] - __alpha * rixs_data_norm_scaled_ni3[:, n_out]) / (1 - __alpha), 'r:', label='Ni 4+ recovered')
# plt.plot(herfd_data_energy, herfd_data_ni3, label='Ni 3+ HERFD')
plt.legend()

plt.xlim(8330, 8360)


###################

plt.figure(1)
plt.clf()

plt.subplot(1,2,1)
plt.plot(bla.dataset.x, bla.dataset.data - np.arange(bla.dataset.t.size) * 0.1)
plt.xlabel('Energy, eV')
plt.ylabel('mu norm')

plt.xlim(24300, 24550)

plt.subplot(2,2,2)
plt.plot(bla.dataset.x, bla.data_ref_fit)
plt.xlabel('Energy, eV')
plt.ylabel('mu norm')

plt.xlim(24300, 24550)

plt.subplot(2,2,4)
plt.plot(bla.dataset.t * 45 / 60, bla.c_fit.T, '.-')
plt.xlabel('time, min')
plt.ylabel('fraction')

plt.xlim(0, 19*45/60)