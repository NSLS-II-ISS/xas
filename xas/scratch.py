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

x = np.arange(0, 1292)
y = np.arange(0, 1292)

_cc = []
for _x, _y in zip(x, y):
    __cc = camera_sp1.calibration.C_dmotdpix(np.array((_x, _y)))
    _cc.append(__cc)

cc = np.array(_cc)

plt.figure(1, clear=True)
plt.plot(cc[:, 0, 0])
plt.plot(cc[:, 1, 0])
plt.plot(cc[:, 0, 1])
plt.plot(cc[:, 1, 1])

########
import copy
from xas.analysis import check_scan
from xas.process import get_processed_df_from_uid
from xas.outliers import outlier_rejection

def get_processed_df_from_uid_advanced(uid, db):
    hdr, df, _, _, _, _, _ = get_processed_df_from_uid(uid, db, save_interpolated_file=False)
    metadata = copy.deepcopy({**hdr.start})


    metadata['mu_good'] = check_scan(df, metadata)
    df['mut'] = -np.log(df["it"] / df["i0"])
    df['mur'] = -np.log(df["ir"] / df["it"])
    df['muf'] = df["iff"] / df["i0"]


    return df, metadata

# pick a sample_uid from scan with uid='5a4ba893-7336-4f9f-b1c6-cf35f3ae9854':
sample_uid = '369ea231-6fbf-4847-a'
uids = list(db.v2.search({'sample_uid': sample_uid}))

dfs, metadatas = [], []

for uid in uids:
    _df, _metadata = get_processed_df_from_uid_advanced(uid, db)
    dfs.append(_df)
    metadatas.append(_metadata)

df_mean, metadata_mean = outlier_rejection(dfs, uids, plot_diagnostics=True)

# ch = 'mut'
# method = 'trimmed_lof'
# for df, uid in zip(dfs, uids):
#     if uid in metadata_mean[ch][method]['inliers']:
#         plt.plot(df['energy'], df[ch], 'k-')
#     else:
#         plt.plot(df['energy'], df[ch], 'b-')
# plt.plot(df_mean['energy'], df_mean[ch][method], 'r-', lw=3)


uid = 'ac830c3a-f227-44c9-ab30-5bf43978e054'
e_nom, e_act, energy_ref_roi, mu_ref_roi, mu_fit = get_energy_offset(uid, db, db_proc, full_return=True)

df, metadata = get_processed_df_from_uid_advanced(uid, db)
plt.figure(2, clear=True)
plt.plot(energy_ref_roi, mu_ref_roi)

plt.plot(energy_ref_roi, mu_fit)

###########
from xas.file_io import load_binned_df_from_file
from xas.energy_calibration import compute_energy_shift_and_broadening_between_spectra
path = '/home/charles/Desktop/test_data/bender_scans'
files = \
[
'Bender scan at Pd-K edge - -264 N - 155.61 um.dat',
'Bender scan at Pd-K edge - -289 N - 150.52 um.dat',
'Bender scan at Pd-K edge - -314 N - 145.58 um.dat',
'Bender scan at Pd-K edge - -338 N - 140.62 um.dat',
'Bender scan at Pd-K edge - -363 N - 135.6 um.dat',
'Bender scan at Pd-K edge - -387 N - 130.6 um.dat',
'Bender scan at Pd-K edge - -413 N - 125.67 um.dat'
]

# energy_ref, mu_ref = db_proc.foil_spectrum('Pd', 'K')
df_ref, _ = load_binned_df_from_file(path + '/' + 'Pd K-edge foil energy calibration.dat')
energy_ref = df_ref["energy"]
mu_ref = -np.log(df_ref["ir"] / df_ref["it"])

dfs = []
for f in files:
    _df, _ = load_binned_df_from_file(path + '/' + f)
    _df['mur'] = -np.log(_df['ir'] / _df['it'])
    dfs.append(_df)

plt.figure(1, clear=True)
for i, df in enumerate(dfs):
    shift, sigma, energy_fit, mu_fit = compute_energy_shift_and_broadening_between_spectra(df['energy'].values, df['mur'].values,
                                                                                           energy_ref, mu_ref, e0=24350, de=200)
    print(i, shift, sigma)
    plt.plot(df['energy'].values, df['mur'].values - i, 'k-')
    plt.plot(energy_fit, mu_fit - i, 'r-')

plt.plot(energy_ref, mu_ref * 1.5 + 1, 'k-', lw=1)


##############

# RUN THIS IN NON-TILED ENVIRONMENT
from xas.process import get_processed_df_from_uid
import json
outpath = '/nsls2/data/iss/legacy/Sandbox/database/databroker_output'

def save_xas_data_as_json(uid, db, db_tag='new'):
    hdr = db[uid]
    if 'experiment' in hdr.start.keys():

        hdr, primary_df, _, _, _, _, _ = get_processed_df_from_uid(uid, db, save_interpolated_file=False)
        true_uid = hdr.start.uid
        _start = {**hdr.start}
        _stop = {'time_stop' : hdr.stop.time,
                 'exit_status' : hdr.stop.exit_status}
        metadata ={**_start, **_stop, 'db_tag': db_tag}
        output = {'metadata': metadata, 'data': primary_df.to_dict()}
        print(true_uid)
        fpath = os.path.join(outpath, true_uid)

        with open(fpath, 'w') as f:
            json.dump(output, f)

n1 = 317987
n2 = 321288
for i in range(n1, n1 + 50):
    save_xas_data_as_json(i, db, db_tag='new')


##############
# RUN THIS IN TILED ENVIRONMENT
import pandas as pd
import json
import os
from tiled.client import from_profile
outpath = '/nsls2/data/iss/legacy/Sandbox/database/databroker_output'

c_nsls2 = from_profile('nsls2', username='MYUSERNAME')
c = c_nsls2['iss']['sandbox']

def read_xas_data_from_json(fpath):
    with open(fpath, 'r') as f:
        data = json.loads(f.read())
    return data['metadata'], pd.DataFrame(data['data'])

files = [f for f in os.listdir(outpath)]
for f in files:
    md, df = read_xas_data_from_json(os.path.join(outpath, f))
    md['tag'] = 'test'
    md['tag_version'] = '0.0'
    md['tag_comment'] = 'initial data upload testing'
    x = c.write_dataframe(df, md)


#####################
from xas.file_io import load_binned_df_from_file

# foil
files = [
'/nsls2/data/iss/legacy/processed/2023/1/310956/foil Mn-K XES EXAFS 0001 vh_roi1.dat',
'/nsls2/data/iss/legacy/processed/2023/1/310956/foil Mn-K XES EXAFS 0003 vh_roi1.dat',
'/nsls2/data/iss/legacy/processed/2023/1/310956/foil Mn-K XES EXAFS 0004 vh_roi1.dat',
'/nsls2/data/iss/legacy/processed/2023/1/310956/foil Mn-K XES EXAFS 0006 vh_roi1.dat',
'/nsls2/data/iss/legacy/processed/2023/1/310956/foil Mn-K XES EXAFS 0007 vh_roi1.dat',
'/nsls2/data/iss/legacy/processed/2023/1/310956/foil Mn-K XES EXAFS 0008 vh_roi1.dat',
'/nsls2/data/iss/legacy/processed/2023/1/310956/foil Mn-K XES EXAFS 0009 vh_roi1.dat',
'/nsls2/data/iss/legacy/processed/2023/1/310956/foil Mn-K XES EXAFS 0010 vh_roi1.dat',
'/nsls2/data/iss/legacy/processed/2023/1/310956/foil Mn-K XES EXAFS 0011 vh_roi1.dat',
'/nsls2/data/iss/legacy/processed/2023/1/310956/foil Mn-K XES EXAFS 0012 vh_roi1.dat',
'/nsls2/data/iss/legacy/processed/2023/1/310956/foil Mn-K XES EXAFS 0001-r0002 vh_roi1.dat',
'/nsls2/data/iss/legacy/processed/2023/1/310956/foil Mn-K XES EXAFS 0002 vh_roi1.dat',
'/nsls2/data/iss/legacy/processed/2023/1/310956/foil Mn-K XES EXAFS 0001-r0003 vh_roi1.dat',
'/nsls2/data/iss/legacy/processed/2023/1/310956/foil Mn-K XES EXAFS 0002-r0002 vh_roi1.dat',
'/nsls2/data/iss/legacy/processed/2023/1/310956/foil Mn-K XES EXAFS 0003-r0002 vh_roi1.dat',
'/nsls2/data/iss/legacy/processed/2023/1/310956/foil 2 Mn-K XES EXAFS 0001 vh_roi1.dat',
'/nsls2/data/iss/legacy/processed/2023/1/310956/foil 2 Mn-K XES EXAFS 0002 vh_roi1.dat',
'/nsls2/data/iss/legacy/processed/2023/1/310956/foil 2 Mn-K XES EXAFS 0003 vh_roi1.dat',
'/nsls2/data/iss/legacy/processed/2023/1/310956/foil 2 Mn-K XES EXAFS 0004 vh_roi1.dat',
'/nsls2/data/iss/legacy/processed/2023/1/310956/foil 2 Mn-K XES EXAFS 0005 vh_roi1.dat',
'/nsls2/data/iss/legacy/processed/2023/1/310956/foil 2 Mn-K XES EXAFS 0006 vh_roi1.dat',
'/nsls2/data/iss/legacy/processed/2023/1/310956/foil 2 Mn-K XES EXAFS 0007 vh_roi1.dat',
'/nsls2/data/iss/legacy/processed/2023/1/310956/foil 2 Mn-K XES EXAFS 0008 vh_roi1.dat',
    ]

# LVMO
files = [
'/nsls2/data/iss/legacy/processed/2023/1/310956/LVMO (pos 047) Mn-K XES EXAFS 0001 vh_roi1.dat',
'/nsls2/data/iss/legacy/processed/2023/1/310956/LVMO (pos 046) Mn-K XES EXAFS 0001 vh_roi1.dat',
'/nsls2/data/iss/legacy/processed/2023/1/310956/LVMO (pos 045) Mn-K XES EXAFS 0001 vh_roi1.dat',
'/nsls2/data/iss/legacy/processed/2023/1/310956/LVMO (pos 044) Mn-K XES EXAFS 0001 vh_roi1.dat',
'/nsls2/data/iss/legacy/processed/2023/1/310956/LVMO (pos 043) Mn-K XES EXAFS 0001 vh_roi1.dat',
'/nsls2/data/iss/legacy/processed/2023/1/310956/LVMO (pos 042) Mn-K XES EXAFS 0001 vh_roi1.dat',
'/nsls2/data/iss/legacy/processed/2023/1/310956/LVMO (pos 041) Mn-K XES EXAFS 0001 vh_roi1.dat',
'/nsls2/data/iss/legacy/processed/2023/1/310956/LVMO (pos 040) Mn-K XES EXAFS 0001 vh_roi1.dat',
'/nsls2/data/iss/legacy/processed/2023/1/310956/LVMO (pos 039) Mn-K XES EXAFS 0001 vh_roi1.dat',
'/nsls2/data/iss/legacy/processed/2023/1/310956/LVMO (pos 038) Mn-K XES EXAFS 0001 vh_roi1.dat',
'/nsls2/data/iss/legacy/processed/2023/1/310956/LVMO (pos 037) Mn-K XES EXAFS 0001 vh_roi1.dat',
'/nsls2/data/iss/legacy/processed/2023/1/310956/LVMO (pos 036) Mn-K XES EXAFS 0001 vh_roi1.dat',
'/nsls2/data/iss/legacy/processed/2023/1/310956/LVMO (pos 035) Mn-K XES EXAFS 0001 vh_roi1.dat',
'/nsls2/data/iss/legacy/processed/2023/1/310956/LVMO (pos 034) Mn-K XES EXAFS 0001 vh_roi1.dat',
'/nsls2/data/iss/legacy/processed/2023/1/310956/LVMO (pos 033) Mn-K XES EXAFS 0001 vh_roi1.dat',
'/nsls2/data/iss/legacy/processed/2023/1/310956/LVMO (pos 032) Mn-K XES EXAFS 0001 vh_roi1.dat',
'/nsls2/data/iss/legacy/processed/2023/1/310956/LVMO (pos 031) Mn-K XES EXAFS 0001 vh_roi1.dat',
'/nsls2/data/iss/legacy/processed/2023/1/310956/LVMO (pos 030) Mn-K XES EXAFS 0001 vh_roi1.dat',
'/nsls2/data/iss/legacy/processed/2023/1/310956/LVMO (pos 029) Mn-K XES EXAFS 0001 vh_roi1.dat',
'/nsls2/data/iss/legacy/processed/2023/1/310956/LVMO (pos 028) Mn-K XES EXAFS 0001 vh_roi1.dat',
'/nsls2/data/iss/legacy/processed/2023/1/310956/LVMO (pos 027) Mn-K XES EXAFS 0001 vh_roi1.dat',
'/nsls2/data/iss/legacy/processed/2023/1/310956/LVMO (pos 026) Mn-K XES EXAFS 0001 vh_roi1.dat',
'/nsls2/data/iss/legacy/processed/2023/1/310956/LVMO (pos 025) Mn-K XES EXAFS 0001 vh_roi1.dat',
'/nsls2/data/iss/legacy/processed/2023/1/310956/LVMO (pos 024) Mn-K XES EXAFS 0001 vh_roi1.dat',
'/nsls2/data/iss/legacy/processed/2023/1/310956/LVMO (pos 023) Mn-K XES EXAFS 0001 vh_roi1.dat',
'/nsls2/data/iss/legacy/processed/2023/1/310956/LVMO (pos 022) Mn-K XES EXAFS 0001 vh_roi1.dat',
'/nsls2/data/iss/legacy/processed/2023/1/310956/LVMO (pos 021) Mn-K XES EXAFS 0001 vh_roi1.dat',
'/nsls2/data/iss/legacy/processed/2023/1/310956/LVMO (pos 020) Mn-K XES EXAFS 0001 vh_roi1.dat',
'/nsls2/data/iss/legacy/processed/2023/1/310956/LVMO (pos 019) Mn-K XES EXAFS 0001 vh_roi1.dat',
'/nsls2/data/iss/legacy/processed/2023/1/310956/LVMO (pos 018) Mn-K XES EXAFS 0001 vh_roi1.dat',
'/nsls2/data/iss/legacy/processed/2023/1/310956/LVMO (pos 017) Mn-K XES EXAFS 0001 vh_roi1.dat',
'/nsls2/data/iss/legacy/processed/2023/1/310956/LVMO (pos 016) Mn-K XES EXAFS 0001 vh_roi1.dat',
'/nsls2/data/iss/legacy/processed/2023/1/310956/LVMO (pos 015) Mn-K XES EXAFS 0001 vh_roi1.dat',
'/nsls2/data/iss/legacy/processed/2023/1/310956/LVMO (pos 014) Mn-K XES EXAFS 0001 vh_roi1.dat',
'/nsls2/data/iss/legacy/processed/2023/1/310956/LVMO (pos 013) Mn-K XES EXAFS 0001 vh_roi1.dat',
'/nsls2/data/iss/legacy/processed/2023/1/310956/LVMO (pos 012) Mn-K XES EXAFS 0001 vh_roi1.dat',
'/nsls2/data/iss/legacy/processed/2023/1/310956/LVMO (pos 011) Mn-K XES EXAFS 0001 vh_roi1.dat',
'/nsls2/data/iss/legacy/processed/2023/1/310956/LVMO (pos 010) Mn-K XES EXAFS 0001 vh_roi1.dat',
'/nsls2/data/iss/legacy/processed/2023/1/310956/LVMO (pos 009) Mn-K XES EXAFS 0001 vh_roi1.dat',
'/nsls2/data/iss/legacy/processed/2023/1/310956/LVMO (pos 008) Mn-K XES EXAFS 0001 vh_roi1.dat',
'/nsls2/data/iss/legacy/processed/2023/1/310956/LVMO (pos 007) Mn-K XES EXAFS 0001 vh_roi1.dat',
'/nsls2/data/iss/legacy/processed/2023/1/310956/LVMO (pos 006) Mn-K XES EXAFS 0001 vh_roi1.dat',
'/nsls2/data/iss/legacy/processed/2023/1/310956/LVMO (pos 005) Mn-K XES EXAFS 0001 vh_roi1.dat',
'/nsls2/data/iss/legacy/processed/2023/1/310956/LVMO (pos 004) Mn-K XES EXAFS 0001 vh_roi1.dat',
'/nsls2/data/iss/legacy/processed/2023/1/310956/LVMO (pos 003) Mn-K XES EXAFS 0001 vh_roi1.dat',
'/nsls2/data/iss/legacy/processed/2023/1/310956/LVMO (pos 002) Mn-K XES EXAFS 0001 vh_roi1.dat',

]

# LMO

files = [
 '/nsls2/data/iss/legacy/processed/2023/1/310956/LMO (pos 412) Mn-K XES EXAFS 0001 vh_roi1.dat',
 '/nsls2/data/iss/legacy/processed/2023/1/310956/LMO (pos 411) Mn-K XES EXAFS 0001 vh_roi1.dat',
 '/nsls2/data/iss/legacy/processed/2023/1/310956/LMO (pos 410) Mn-K XES EXAFS 0001 vh_roi1.dat',
 '/nsls2/data/iss/legacy/processed/2023/1/310956/LMO (pos 409) Mn-K XES EXAFS 0001 vh_roi1.dat',
 '/nsls2/data/iss/legacy/processed/2023/1/310956/LMO (pos 408) Mn-K XES EXAFS 0001 vh_roi1.dat',
 '/nsls2/data/iss/legacy/processed/2023/1/310956/LMO (pos 407) Mn-K XES EXAFS 0001 vh_roi1.dat',
 '/nsls2/data/iss/legacy/processed/2023/1/310956/LMO (pos 406) Mn-K XES EXAFS 0001 vh_roi1.dat',
 '/nsls2/data/iss/legacy/processed/2023/1/310956/LMO (pos 405) Mn-K XES EXAFS 0001 vh_roi1.dat',
 '/nsls2/data/iss/legacy/processed/2023/1/310956/LMO (pos 404) Mn-K XES EXAFS 0001 vh_roi1.dat',
 '/nsls2/data/iss/legacy/processed/2023/1/310956/LMO (pos 403) Mn-K XES EXAFS 0001 vh_roi1.dat',
 '/nsls2/data/iss/legacy/processed/2023/1/310956/LMO (pos 402) Mn-K XES EXAFS 0001 vh_roi1.dat',
 '/nsls2/data/iss/legacy/processed/2023/1/310956/LMO (pos 400) Mn-K XES EXAFS 0001 vh_roi1.dat',
 '/nsls2/data/iss/legacy/processed/2023/1/310956/LMO (pos 399) Mn-K XES EXAFS 0001 vh_roi1.dat',
 '/nsls2/data/iss/legacy/processed/2023/1/310956/LMO (pos 398) Mn-K XES EXAFS 0001 vh_roi1.dat',
 '/nsls2/data/iss/legacy/processed/2023/1/310956/LMO (pos 397) Mn-K XES EXAFS 0001 vh_roi1.dat',
 '/nsls2/data/iss/legacy/processed/2023/1/310956/LMO (pos 396) Mn-K XES EXAFS 0001 vh_roi1.dat',
 '/nsls2/data/iss/legacy/processed/2023/1/310956/LMO (pos 395) Mn-K XES EXAFS 0001 vh_roi1.dat',
 '/nsls2/data/iss/legacy/processed/2023/1/310956/LMO (pos 394) Mn-K XES EXAFS 0001 vh_roi1.dat',
 '/nsls2/data/iss/legacy/processed/2023/1/310956/LMO (pos 393) Mn-K XES EXAFS 0001 vh_roi1.dat',
 '/nsls2/data/iss/legacy/processed/2023/1/310956/LMO (pos 392) Mn-K XES EXAFS 0001 vh_roi1.dat',
 '/nsls2/data/iss/legacy/processed/2023/1/310956/LMO (pos 391) Mn-K XES EXAFS 0001 vh_roi1.dat',
 '/nsls2/data/iss/legacy/processed/2023/1/310956/LMO (pos 390) Mn-K XES EXAFS 0001 vh_roi1.dat',
 '/nsls2/data/iss/legacy/processed/2023/1/310956/LMO (pos 389) Mn-K XES EXAFS 0001 vh_roi1.dat',
 '/nsls2/data/iss/legacy/processed/2023/1/310956/LMO (pos 388) Mn-K XES EXAFS 0001 vh_roi1.dat',
 '/nsls2/data/iss/legacy/processed/2023/1/310956/LMO (pos 387) Mn-K XES EXAFS 0001 vh_roi1.dat',
 '/nsls2/data/iss/legacy/processed/2023/1/310956/LMO (pos 386) Mn-K XES EXAFS 0001 vh_roi1.dat',
 '/nsls2/data/iss/legacy/processed/2023/1/310956/LMO (pos 385) Mn-K XES EXAFS 0001 vh_roi1.dat',
 '/nsls2/data/iss/legacy/processed/2023/1/310956/LMO (pos 384) Mn-K XES EXAFS 0001 vh_roi1.dat',
 '/nsls2/data/iss/legacy/processed/2023/1/310956/LMO (pos 383) Mn-K XES EXAFS 0001 vh_roi1.dat',
 '/nsls2/data/iss/legacy/processed/2023/1/310956/LMO (pos 382) Mn-K XES EXAFS 0001 vh_roi1.dat',
 '/nsls2/data/iss/legacy/processed/2023/1/310956/LMO (pos 381) Mn-K XES EXAFS 0001 vh_roi1.dat',
 '/nsls2/data/iss/legacy/processed/2023/1/310956/LMO (pos 380) Mn-K XES EXAFS 0001 vh_roi1.dat',
 '/nsls2/data/iss/legacy/processed/2023/1/310956/LMO (pos 379) Mn-K XES EXAFS 0001 vh_roi1.dat',
 '/nsls2/data/iss/legacy/processed/2023/1/310956/LMO (pos 378) Mn-K XES EXAFS 0001 vh_roi1.dat',
 '/nsls2/data/iss/legacy/processed/2023/1/310956/LMO (pos 377) Mn-K XES EXAFS 0001 vh_roi1.dat',
 '/nsls2/data/iss/legacy/processed/2023/1/310956/LMO (pos 376) Mn-K XES EXAFS 0001 vh_roi1.dat',
 '/nsls2/data/iss/legacy/processed/2023/1/310956/LMO (pos 375) Mn-K XES EXAFS 0001 vh_roi1.dat',
 '/nsls2/data/iss/legacy/processed/2023/1/310956/LMO (pos 374) Mn-K XES EXAFS 0001 vh_roi1.dat',
 '/nsls2/data/iss/legacy/processed/2023/1/310956/LMO (pos 373) Mn-K XES EXAFS 0001 vh_roi1.dat',
 '/nsls2/data/iss/legacy/processed/2023/1/310956/LMO (pos 372) Mn-K XES EXAFS 0001 vh_roi1.dat',
 '/nsls2/data/iss/legacy/processed/2023/1/310956/LMO (pos 371) Mn-K XES EXAFS 0001 vh_roi1.dat',
 '/nsls2/data/iss/legacy/processed/2023/1/310956/LMO (pos 370) Mn-K XES EXAFS 0001 vh_roi1.dat',
 '/nsls2/data/iss/legacy/processed/2023/1/310956/LMO (pos 369) Mn-K XES EXAFS 0001 vh_roi1.dat',
 '/nsls2/data/iss/legacy/processed/2023/1/310956/LMO (pos 368) Mn-K XES EXAFS 0001 vh_roi1.dat',
 '/nsls2/data/iss/legacy/processed/2023/1/310956/LMO (pos 367) Mn-K XES EXAFS 0001 vh_roi1.dat',
 '/nsls2/data/iss/legacy/processed/2023/1/310956/LMO (pos 366) Mn-K XES EXAFS 0001 vh_roi1.dat',
 '/nsls2/data/iss/legacy/processed/2023/1/310956/LMO (pos 365) Mn-K XES EXAFS 0001 vh_roi1.dat',
 '/nsls2/data/iss/legacy/processed/2023/1/310956/LMO (pos 364) Mn-K XES EXAFS 0001 vh_roi1.dat',
 '/nsls2/data/iss/legacy/processed/2023/1/310956/LMO (pos 363) Mn-K XES EXAFS 0001 vh_roi1.dat',
 '/nsls2/data/iss/legacy/processed/2023/1/310956/LMO (pos 362) Mn-K XES EXAFS 0001 vh_roi1.dat',
 '/nsls2/data/iss/legacy/processed/2023/1/310956/LMO (pos 361) Mn-K XES EXAFS 0001 vh_roi1.dat',
 '/nsls2/data/iss/legacy/processed/2023/1/310956/LMO (pos 360) Mn-K XES EXAFS 0001 vh_roi1.dat',
 '/nsls2/data/iss/legacy/processed/2023/1/310956/LMO (pos 359) Mn-K XES EXAFS 0001 vh_roi1.dat',
 '/nsls2/data/iss/legacy/processed/2023/1/310956/LMO (pos 358) Mn-K XES EXAFS 0001 vh_roi1.dat',
 '/nsls2/data/iss/legacy/processed/2023/1/310956/LMO (pos 357) Mn-K XES EXAFS 0001 vh_roi1.dat',
 '/nsls2/data/iss/legacy/processed/2023/1/310956/LMO (pos 356) Mn-K XES EXAFS 0001 vh_roi1.dat',
 '/nsls2/data/iss/legacy/processed/2023/1/310956/LMO (pos 355) Mn-K XES EXAFS 0001 vh_roi1.dat',
 '/nsls2/data/iss/legacy/processed/2023/1/310956/LMO (pos 354) Mn-K XES EXAFS 0001 vh_roi1.dat',
 '/nsls2/data/iss/legacy/processed/2023/1/310956/LMO (pos 353) Mn-K XES EXAFS 0001 vh_roi1.dat',
 '/nsls2/data/iss/legacy/processed/2023/1/310956/LMO (pos 352) Mn-K XES EXAFS 0001 vh_roi1.dat',
 '/nsls2/data/iss/legacy/processed/2023/1/310956/LMO (pos 351) Mn-K XES EXAFS 0001 vh_roi1.dat',
 '/nsls2/data/iss/legacy/processed/2023/1/310956/LMO (pos 350) Mn-K XES EXAFS 0001 vh_roi1.dat',
 '/nsls2/data/iss/legacy/processed/2023/1/310956/LMO (pos 349) Mn-K XES EXAFS 0001 vh_roi1.dat',
 '/nsls2/data/iss/legacy/processed/2023/1/310956/LMO (pos 348) Mn-K XES EXAFS 0001 vh_roi1.dat',
 '/nsls2/data/iss/legacy/processed/2023/1/310956/LMO (pos 347) Mn-K XES EXAFS 0001 vh_roi1.dat',
 '/nsls2/data/iss/legacy/processed/2023/1/310956/LMO (pos 346) Mn-K XES EXAFS 0001 vh_roi1.dat',
 '/nsls2/data/iss/legacy/processed/2023/1/310956/LMO (pos 345) Mn-K XES EXAFS 0001 vh_roi1.dat',
 '/nsls2/data/iss/legacy/processed/2023/1/310956/LMO (pos 344) Mn-K XES EXAFS 0001 vh_roi1.dat',
 '/nsls2/data/iss/legacy/processed/2023/1/310956/LMO (pos 343) Mn-K XES EXAFS 0001 vh_roi1.dat',
 '/nsls2/data/iss/legacy/processed/2023/1/310956/LMO (pos 342) Mn-K XES EXAFS 0001 vh_roi1.dat',
 '/nsls2/data/iss/legacy/processed/2023/1/310956/LMO (pos 341) Mn-K XES EXAFS 0001 vh_roi1.dat',
 '/nsls2/data/iss/legacy/processed/2023/1/310956/LMO (pos 340) Mn-K XES EXAFS 0001 vh_roi1.dat',
 '/nsls2/data/iss/legacy/processed/2023/1/310956/LMO (pos 339) Mn-K XES EXAFS 0001 vh_roi1.dat',
 '/nsls2/data/iss/legacy/processed/2023/1/310956/LMO (pos 338) Mn-K XES EXAFS 0001 vh_roi1.dat',
 '/nsls2/data/iss/legacy/processed/2023/1/310956/LMO (pos 337) Mn-K XES EXAFS 0001 vh_roi1.dat',
 '/nsls2/data/iss/legacy/processed/2023/1/310956/LMO (pos 336) Mn-K XES EXAFS 0001 vh_roi1.dat',
 '/nsls2/data/iss/legacy/processed/2023/1/310956/LMO (pos 335) Mn-K XES EXAFS 0001 vh_roi1.dat',
 '/nsls2/data/iss/legacy/processed/2023/1/310956/LMO (pos 334) Mn-K XES EXAFS 0001 vh_roi1.dat',
 '/nsls2/data/iss/legacy/processed/2023/1/310956/LMO (pos 333) Mn-K XES EXAFS 0001 vh_roi1.dat',
 '/nsls2/data/iss/legacy/processed/2023/1/310956/LMO (pos 332) Mn-K XES EXAFS 0001 vh_roi1.dat',
 '/nsls2/data/iss/legacy/processed/2023/1/310956/LMO (pos 331) Mn-K XES EXAFS 0001 vh_roi1.dat',
 '/nsls2/data/iss/legacy/processed/2023/1/310956/LMO (pos 330) Mn-K XES EXAFS 0001 vh_roi1.dat',
 '/nsls2/data/iss/legacy/processed/2023/1/310956/LMO (pos 329) Mn-K XES EXAFS 0001 vh_roi1.dat',
 '/nsls2/data/iss/legacy/processed/2023/1/310956/LMO (pos 328) Mn-K XES EXAFS 0001 vh_roi1.dat',
 '/nsls2/data/iss/legacy/processed/2023/1/310956/LMO (pos 327) Mn-K XES EXAFS 0001 vh_roi1.dat',
 '/nsls2/data/iss/legacy/processed/2023/1/310956/LMO (pos 326) Mn-K XES EXAFS 0001 vh_roi1.dat',
 '/nsls2/data/iss/legacy/processed/2023/1/310956/LMO (pos 324) Mn-K XES EXAFS 0001 vh_roi1.dat',
 '/nsls2/data/iss/legacy/processed/2023/1/310956/LMO (pos 323) Mn-K XES EXAFS 0001 vh_roi1.dat',
 '/nsls2/data/iss/legacy/processed/2023/1/310956/LMO (pos 322) Mn-K XES EXAFS 0001 vh_roi1.dat',
 '/nsls2/data/iss/legacy/processed/2023/1/310956/LMO (pos 321) Mn-K XES EXAFS 0001 vh_roi1.dat',
 '/nsls2/data/iss/legacy/processed/2023/1/310956/LMO (pos 320) Mn-K XES EXAFS 0001 vh_roi1.dat',
 '/nsls2/data/iss/legacy/processed/2023/1/310956/LMO (pos 319) Mn-K XES EXAFS 0001 vh_roi1.dat',
 '/nsls2/data/iss/legacy/processed/2023/1/310956/LMO (pos 318) Mn-K XES EXAFS 0001 vh_roi1.dat',
 '/nsls2/data/iss/legacy/processed/2023/1/310956/LMO (pos 317) Mn-K XES EXAFS 0001 vh_roi1.dat',
 '/nsls2/data/iss/legacy/processed/2023/1/310956/LMO (pos 316) Mn-K XES EXAFS 0001 vh_roi1.dat',
 '/nsls2/data/iss/legacy/processed/2023/1/310956/LMO (pos 315) Mn-K XES EXAFS 0001 vh_roi1.dat',
 '/nsls2/data/iss/legacy/processed/2023/1/310956/LMO (pos 314) Mn-K XES EXAFS 0001 vh_roi1.dat',
 '/nsls2/data/iss/legacy/processed/2023/1/310956/LMO (pos 313) Mn-K XES EXAFS 0001 vh_roi1.dat',
 '/nsls2/data/iss/legacy/processed/2023/1/310956/LMO (pos 312) Mn-K XES EXAFS 0001 vh_roi1.dat',
 '/nsls2/data/iss/legacy/processed/2023/1/310956/LMO (pos 311) Mn-K XES EXAFS 0001 vh_roi1.dat',
 '/nsls2/data/iss/legacy/processed/2023/1/310956/LMO (pos 310) Mn-K XES EXAFS 0001 vh_roi1.dat',
 '/nsls2/data/iss/legacy/processed/2023/1/310956/LMO (pos 309) Mn-K XES EXAFS 0001 vh_roi1.dat',
 '/nsls2/data/iss/legacy/processed/2023/1/310956/LMO (pos 308) Mn-K XES EXAFS 0001 vh_roi1.dat',
 '/nsls2/data/iss/legacy/processed/2023/1/310956/LMO (pos 307) Mn-K XES EXAFS 0001 vh_roi1.dat',
 '/nsls2/data/iss/legacy/processed/2023/1/310956/LMO (pos 306) Mn-K XES EXAFS 0001 vh_roi1.dat',
 '/nsls2/data/iss/legacy/processed/2023/1/310956/LMO (pos 305) Mn-K XES EXAFS 0001 vh_roi1.dat',
 '/nsls2/data/iss/legacy/processed/2023/1/310956/LMO (pos 304) Mn-K XES EXAFS 0001 vh_roi1.dat',
 '/nsls2/data/iss/legacy/processed/2023/1/310956/LMO (pos 303) Mn-K XES EXAFS 0001 vh_roi1.dat',
 '/nsls2/data/iss/legacy/processed/2023/1/310956/LMO (pos 302) Mn-K XES EXAFS 0001 vh_roi1.dat',
 '/nsls2/data/iss/legacy/processed/2023/1/310956/LMO (pos 300) Mn-K XES EXAFS 0001 vh_roi1.dat',
 '/nsls2/data/iss/legacy/processed/2023/1/310956/LMO (pos 299) Mn-K XES EXAFS 0001 vh_roi1.dat',
 '/nsls2/data/iss/legacy/processed/2023/1/310956/LMO (pos 298) Mn-K XES EXAFS 0001 vh_roi1.dat',
 '/nsls2/data/iss/legacy/processed/2023/1/310956/LMO (pos 297) Mn-K XES EXAFS 0001 vh_roi1.dat',
 '/nsls2/data/iss/legacy/processed/2023/1/310956/LMO (pos 296) Mn-K XES EXAFS 0001 vh_roi1.dat',
 '/nsls2/data/iss/legacy/processed/2023/1/310956/LMO (pos 295) Mn-K XES EXAFS 0001 vh_roi1.dat',
 '/nsls2/data/iss/legacy/processed/2023/1/310956/LMO (pos 294) Mn-K XES EXAFS 0001 vh_roi1.dat',
 '/nsls2/data/iss/legacy/processed/2023/1/310956/LMO (pos 293) Mn-K XES EXAFS 0001 vh_roi1.dat',
 '/nsls2/data/iss/legacy/processed/2023/1/310956/LMO (pos 292) Mn-K XES EXAFS 0001 vh_roi1.dat',
 '/nsls2/data/iss/legacy/processed/2023/1/310956/LMO (pos 291) Mn-K XES EXAFS 0001 vh_roi1.dat',
 '/nsls2/data/iss/legacy/processed/2023/1/310956/LMO (pos 290) Mn-K XES EXAFS 0001 vh_roi1.dat',
 '/nsls2/data/iss/legacy/processed/2023/1/310956/LMO (pos 289) Mn-K XES EXAFS 0001 vh_roi1.dat',
 '/nsls2/data/iss/legacy/processed/2023/1/310956/LMO (pos 288) Mn-K XES EXAFS 0001 vh_roi1.dat',
 '/nsls2/data/iss/legacy/processed/2023/1/310956/LMO (pos 287) Mn-K XES EXAFS 0001 vh_roi1.dat',
 '/nsls2/data/iss/legacy/processed/2023/1/310956/LMO (pos 286) Mn-K XES EXAFS 0001 vh_roi1.dat',
 '/nsls2/data/iss/legacy/processed/2023/1/310956/LMO (pos 285) Mn-K XES EXAFS 0001 vh_roi1.dat',
 '/nsls2/data/iss/legacy/processed/2023/1/310956/LMO (pos 284) Mn-K XES EXAFS 0001 vh_roi1.dat',
 '/nsls2/data/iss/legacy/processed/2023/1/310956/LMO (pos 283) Mn-K XES EXAFS 0001 vh_roi1.dat',
 '/nsls2/data/iss/legacy/processed/2023/1/310956/LMO (pos 282) Mn-K XES EXAFS 0001 vh_roi1.dat',
 '/nsls2/data/iss/legacy/processed/2023/1/310956/LMO (pos 281) Mn-K XES EXAFS 0001 vh_roi1.dat',
 '/nsls2/data/iss/legacy/processed/2023/1/310956/LMO (pos 280) Mn-K XES EXAFS 0001 vh_roi1.dat',
 '/nsls2/data/iss/legacy/processed/2023/1/310956/LMO (pos 279) Mn-K XES EXAFS 0001 vh_roi1.dat',
 '/nsls2/data/iss/legacy/processed/2023/1/310956/LMO (pos 278) Mn-K XES EXAFS 0001 vh_roi1.dat',
 '/nsls2/data/iss/legacy/processed/2023/1/310956/LMO (pos 277) Mn-K XES EXAFS 0001 vh_roi1.dat',
 '/nsls2/data/iss/legacy/processed/2023/1/310956/LMO (pos 276) Mn-K XES EXAFS 0001 vh_roi1.dat',
 '/nsls2/data/iss/legacy/processed/2023/1/310956/LMO (pos 275) Mn-K XES EXAFS 0001 vh_roi1.dat',
 '/nsls2/data/iss/legacy/processed/2023/1/310956/LMO (pos 274) Mn-K XES EXAFS 0001 vh_roi1.dat',
 '/nsls2/data/iss/legacy/processed/2023/1/310956/LMO (pos 273) Mn-K XES EXAFS 0001 vh_roi1.dat',
 '/nsls2/data/iss/legacy/processed/2023/1/310956/LMO (pos 272) Mn-K XES EXAFS 0001 vh_roi1.dat',
 '/nsls2/data/iss/legacy/processed/2023/1/310956/LMO (pos 271) Mn-K XES EXAFS 0001 vh_roi1.dat',
 '/nsls2/data/iss/legacy/processed/2023/1/310956/LMO (pos 270) Mn-K XES EXAFS 0001 vh_roi1.dat',
 '/nsls2/data/iss/legacy/processed/2023/1/310956/LMO (pos 269) Mn-K XES EXAFS 0001 vh_roi1.dat',
 '/nsls2/data/iss/legacy/processed/2023/1/310956/LMO (pos 268) Mn-K XES EXAFS 0001 vh_roi1.dat',
 '/nsls2/data/iss/legacy/processed/2023/1/310956/LMO (pos 267) Mn-K XES EXAFS 0001 vh_roi1.dat',
 '/nsls2/data/iss/legacy/processed/2023/1/310956/LMO (pos 266) Mn-K XES EXAFS 0001 vh_roi1.dat',
 '/nsls2/data/iss/legacy/processed/2023/1/310956/LMO (pos 265) Mn-K XES EXAFS 0001 vh_roi1.dat',
 '/nsls2/data/iss/legacy/processed/2023/1/310956/LMO (pos 264) Mn-K XES EXAFS 0001 vh_roi1.dat',
 '/nsls2/data/iss/legacy/processed/2023/1/310956/LMO (pos 263) Mn-K XES EXAFS 0001 vh_roi1.dat',
 '/nsls2/data/iss/legacy/processed/2023/1/310956/LMO (pos 262) Mn-K XES EXAFS 0001 vh_roi1.dat',
 '/nsls2/data/iss/legacy/processed/2023/1/310956/LMO (pos 261) Mn-K XES EXAFS 0001 vh_roi1.dat',
 '/nsls2/data/iss/legacy/processed/2023/1/310956/LMO (pos 260) Mn-K XES EXAFS 0001 vh_roi1.dat',
 '/nsls2/data/iss/legacy/processed/2023/1/310956/LMO (pos 259) Mn-K XES EXAFS 0001 vh_roi1.dat',
 '/nsls2/data/iss/legacy/processed/2023/1/310956/LMO (pos 258) Mn-K XES EXAFS 0001 vh_roi1.dat',
 '/nsls2/data/iss/legacy/processed/2023/1/310956/LMO (pos 257) Mn-K XES EXAFS 0001 vh_roi1.dat',
 '/nsls2/data/iss/legacy/processed/2023/1/310956/LMO (pos 256) Mn-K XES EXAFS 0001 vh_roi1.dat',
 '/nsls2/data/iss/legacy/processed/2023/1/310956/LMO (pos 254) Mn-K XES EXAFS 0001 vh_roi1.dat',
 '/nsls2/data/iss/legacy/processed/2023/1/310956/LMO (pos 253) Mn-K XES EXAFS 0001 vh_roi1.dat',
 '/nsls2/data/iss/legacy/processed/2023/1/310956/LMO (pos 252) Mn-K XES EXAFS 0001 vh_roi1.dat',
 '/nsls2/data/iss/legacy/processed/2023/1/310956/LMO (pos 251) Mn-K XES EXAFS 0001 vh_roi1.dat',
 '/nsls2/data/iss/legacy/processed/2023/1/310956/LMO (pos 250) Mn-K XES EXAFS 0001 vh_roi1.dat',
 '/nsls2/data/iss/legacy/processed/2023/1/310956/LMO (pos 249) Mn-K XES EXAFS 0001 vh_roi1.dat',
 '/nsls2/data/iss/legacy/processed/2023/1/310956/LMO (pos 248) Mn-K XES EXAFS 0001 vh_roi1.dat',
 '/nsls2/data/iss/legacy/processed/2023/1/310956/LMO (pos 247) Mn-K XES EXAFS 0001 vh_roi1.dat',
 '/nsls2/data/iss/legacy/processed/2023/1/310956/LMO (pos 246) Mn-K XES EXAFS 0001 vh_roi1.dat',
 '/nsls2/data/iss/legacy/processed/2023/1/310956/LMO (pos 245) Mn-K XES EXAFS 0001 vh_roi1.dat',
 '/nsls2/data/iss/legacy/processed/2023/1/310956/LMO (pos 244) Mn-K XES EXAFS 0001 vh_roi1.dat',
 '/nsls2/data/iss/legacy/processed/2023/1/310956/LMO (pos 243) Mn-K XES EXAFS 0001 vh_roi1.dat',
 '/nsls2/data/iss/legacy/processed/2023/1/310956/LMO (pos 242) Mn-K XES EXAFS 0001 vh_roi1.dat',
 '/nsls2/data/iss/legacy/processed/2023/1/310956/LMO (pos 241) Mn-K XES EXAFS 0001 vh_roi1.dat',
 '/nsls2/data/iss/legacy/processed/2023/1/310956/LMO (pos 240) Mn-K XES EXAFS 0001 vh_roi1.dat',
 '/nsls2/data/iss/legacy/processed/2023/1/310956/LMO (pos 239) Mn-K XES EXAFS 0001 vh_roi1.dat',
 '/nsls2/data/iss/legacy/processed/2023/1/310956/LMO (pos 238) Mn-K XES EXAFS 0001 vh_roi1.dat',
 '/nsls2/data/iss/legacy/processed/2023/1/310956/LMO (pos 237) Mn-K XES EXAFS 0001 vh_roi1.dat',
 '/nsls2/data/iss/legacy/processed/2023/1/310956/LMO (pos 236) Mn-K XES EXAFS 0001 vh_roi1.dat',
 '/nsls2/data/iss/legacy/processed/2023/1/310956/LMO (pos 235) Mn-K XES EXAFS 0001 vh_roi1.dat',
 '/nsls2/data/iss/legacy/processed/2023/1/310956/LMO (pos 234) Mn-K XES EXAFS 0001 vh_roi1.dat',
 '/nsls2/data/iss/legacy/processed/2023/1/310956/LMO (pos 233) Mn-K XES EXAFS 0001 vh_roi1.dat',
 '/nsls2/data/iss/legacy/processed/2023/1/310956/LMO (pos 232) Mn-K XES EXAFS 0001 vh_roi1.dat',
 '/nsls2/data/iss/legacy/processed/2023/1/310956/LMO (pos 231) Mn-K XES EXAFS 0001 vh_roi1.dat',
 '/nsls2/data/iss/legacy/processed/2023/1/310956/LMO (pos 230) Mn-K XES EXAFS 0001 vh_roi1.dat',
 '/nsls2/data/iss/legacy/processed/2023/1/310956/LMO (pos 229) Mn-K XES EXAFS 0001 vh_roi1.dat',
 '/nsls2/data/iss/legacy/processed/2023/1/310956/LMO (pos 228) Mn-K XES EXAFS 0001 vh_roi1.dat',
 '/nsls2/data/iss/legacy/processed/2023/1/310956/LMO (pos 227) Mn-K XES EXAFS 0001 vh_roi1.dat',
 '/nsls2/data/iss/legacy/processed/2023/1/310956/LMO (pos 226) Mn-K XES EXAFS 0001 vh_roi1.dat',
 '/nsls2/data/iss/legacy/processed/2023/1/310956/LMO (pos 225) Mn-K XES EXAFS 0001 vh_roi1.dat',
 '/nsls2/data/iss/legacy/processed/2023/1/310956/LMO (pos 224) Mn-K XES EXAFS 0001 vh_roi1.dat',
 '/nsls2/data/iss/legacy/processed/2023/1/310956/LMO (pos 223) Mn-K XES EXAFS 0001 vh_roi1.dat',
 '/nsls2/data/iss/legacy/processed/2023/1/310956/LMO (pos 222) Mn-K XES EXAFS 0001 vh_roi1.dat',
 '/nsls2/data/iss/legacy/processed/2023/1/310956/LMO (pos 221) Mn-K XES EXAFS 0001 vh_roi1.dat',
 '/nsls2/data/iss/legacy/processed/2023/1/310956/LMO (pos 220) Mn-K XES EXAFS 0001 vh_roi1.dat',
 '/nsls2/data/iss/legacy/processed/2023/1/310956/LMO (pos 219) Mn-K XES EXAFS 0001 vh_roi1.dat',
 '/nsls2/data/iss/legacy/processed/2023/1/310956/LMO (pos 218) Mn-K XES EXAFS 0001 vh_roi1.dat',
 '/nsls2/data/iss/legacy/processed/2023/1/310956/LMO (pos 217) Mn-K XES EXAFS 0001 vh_roi1.dat',
 '/nsls2/data/iss/legacy/processed/2023/1/310956/LMO (pos 216) Mn-K XES EXAFS 0001 vh_roi1.dat',
 '/nsls2/data/iss/legacy/processed/2023/1/310956/LMO (pos 214) Mn-K XES EXAFS 0001-r0002 vh_roi1.dat',
 '/nsls2/data/iss/legacy/processed/2023/1/310956/LMO (pos 213) Mn-K XES EXAFS 0001 vh_roi1.dat',
 '/nsls2/data/iss/legacy/processed/2023/1/310956/LMO (pos 212) Mn-K XES EXAFS 0001 vh_roi1.dat',
 '/nsls2/data/iss/legacy/processed/2023/1/310956/LMO (pos 211) Mn-K XES EXAFS 0001 vh_roi1.dat',
 '/nsls2/data/iss/legacy/processed/2023/1/310956/LMO (pos 210) Mn-K XES EXAFS 0001 vh_roi1.dat',
 '/nsls2/data/iss/legacy/processed/2023/1/310956/LMO (pos 209) Mn-K XES EXAFS 0001 vh_roi1.dat',
 '/nsls2/data/iss/legacy/processed/2023/1/310956/LMO (pos 208) Mn-K XES EXAFS 0001 vh_roi1.dat',
 '/nsls2/data/iss/legacy/processed/2023/1/310956/LMO (pos 207) Mn-K XES EXAFS 0001 vh_roi1.dat',
 '/nsls2/data/iss/legacy/processed/2023/1/310956/LMO (pos 206) Mn-K XES EXAFS 0001 vh_roi1.dat',
 '/nsls2/data/iss/legacy/processed/2023/1/310956/LMO (pos 205) Mn-K XES EXAFS 0001 vh_roi1.dat',
 '/nsls2/data/iss/legacy/processed/2023/1/310956/LMO (pos 204) Mn-K XES EXAFS 0001 vh_roi1.dat',
 '/nsls2/data/iss/legacy/processed/2023/1/310956/LMO (pos 203) Mn-K XES EXAFS 0001 vh_roi1.dat',
 '/nsls2/data/iss/legacy/processed/2023/1/310956/LMO (pos 202) Mn-K XES EXAFS 0001 vh_roi1.dat',
 '/nsls2/data/iss/legacy/processed/2023/1/310956/LMO (pos 200) Mn-K XES EXAFS 0001 vh_roi1.dat',
 '/nsls2/data/iss/legacy/processed/2023/1/310956/LMO (pos 199) Mn-K XES EXAFS 0001 vh_roi1.dat',
 '/nsls2/data/iss/legacy/processed/2023/1/310956/LMO (pos 198) Mn-K XES EXAFS 0001 vh_roi1.dat',
 '/nsls2/data/iss/legacy/processed/2023/1/310956/LMO (pos 197) Mn-K XES EXAFS 0001 vh_roi1.dat',
 '/nsls2/data/iss/legacy/processed/2023/1/310956/LMO (pos 196) Mn-K XES EXAFS 0001 vh_roi1.dat',
 '/nsls2/data/iss/legacy/processed/2023/1/310956/LMO (pos 195) Mn-K XES EXAFS 0001 vh_roi1.dat',
 '/nsls2/data/iss/legacy/processed/2023/1/310956/LMO (pos 194) Mn-K XES EXAFS 0001 vh_roi1.dat',
 '/nsls2/data/iss/legacy/processed/2023/1/310956/LMO (pos 193) Mn-K XES EXAFS 0001 vh_roi1.dat',
 '/nsls2/data/iss/legacy/processed/2023/1/310956/LMO (pos 192) Mn-K XES EXAFS 0001 vh_roi1.dat',
 '/nsls2/data/iss/legacy/processed/2023/1/310956/LMO (pos 191) Mn-K XES EXAFS 0001 vh_roi1.dat',
 '/nsls2/data/iss/legacy/processed/2023/1/310956/LMO (pos 189) Mn-K XES EXAFS 0001 vh_roi1.dat',
 '/nsls2/data/iss/legacy/processed/2023/1/310956/LMO (pos 188) Mn-K XES EXAFS 0001 vh_roi1.dat',
 '/nsls2/data/iss/legacy/processed/2023/1/310956/LMO (pos 187) Mn-K XES EXAFS 0001 vh_roi1.dat',
 '/nsls2/data/iss/legacy/processed/2023/1/310956/LMO (pos 186) Mn-K XES EXAFS 0001 vh_roi1.dat',
 '/nsls2/data/iss/legacy/processed/2023/1/310956/LMO (pos 185) Mn-K XES EXAFS 0001 vh_roi1.dat',
 '/nsls2/data/iss/legacy/processed/2023/1/310956/LMO (pos 184) Mn-K XES EXAFS 0001 vh_roi1.dat',
 '/nsls2/data/iss/legacy/processed/2023/1/310956/LMO (pos 183) Mn-K XES EXAFS 0001 vh_roi1.dat',
 '/nsls2/data/iss/legacy/processed/2023/1/310956/LMO (pos 182) Mn-K XES EXAFS 0001 vh_roi1.dat',
 '/nsls2/data/iss/legacy/processed/2023/1/310956/LMO (pos 181) Mn-K XES EXAFS 0001 vh_roi1.dat',
 '/nsls2/data/iss/legacy/processed/2023/1/310956/LMO (pos 180) Mn-K XES EXAFS 0001 vh_roi1.dat',
 '/nsls2/data/iss/legacy/processed/2023/1/310956/LMO (pos 179) Mn-K XES EXAFS 0001 vh_roi1.dat',
 '/nsls2/data/iss/legacy/processed/2023/1/310956/LMO (pos 178) Mn-K XES EXAFS 0001 vh_roi1.dat',
 '/nsls2/data/iss/legacy/processed/2023/1/310956/LMO (pos 177) Mn-K XES EXAFS 0001 vh_roi1.dat',
 '/nsls2/data/iss/legacy/processed/2023/1/310956/LMO (pos 176) Mn-K XES EXAFS 0001 vh_roi1.dat',
 '/nsls2/data/iss/legacy/processed/2023/1/310956/LMO (pos 175) Mn-K XES EXAFS 0001 vh_roi1.dat',
 '/nsls2/data/iss/legacy/processed/2023/1/310956/LMO (pos 174) Mn-K XES EXAFS 0001 vh_roi1.dat',
 '/nsls2/data/iss/legacy/processed/2023/1/310956/LMO (pos 173) Mn-K XES EXAFS 0001 vh_roi1.dat',
 '/nsls2/data/iss/legacy/processed/2023/1/310956/LMO (pos 172) Mn-K XES EXAFS 0001 vh_roi1.dat',
 '/nsls2/data/iss/legacy/processed/2023/1/310956/LMO (pos 171) Mn-K XES EXAFS 0001 vh_roi1.dat',
 '/nsls2/data/iss/legacy/processed/2023/1/310956/LMO (pos 170) Mn-K XES EXAFS 0001 vh_roi1.dat',
 '/nsls2/data/iss/legacy/processed/2023/1/310956/LMO (pos 169) Mn-K XES EXAFS 0001 vh_roi1.dat',
 '/nsls2/data/iss/legacy/processed/2023/1/310956/LMO (pos 168) Mn-K XES EXAFS 0001 vh_roi1.dat',
 '/nsls2/data/iss/legacy/processed/2023/1/310956/LMO (pos 167) Mn-K XES EXAFS 0001 vh_roi1.dat',
 '/nsls2/data/iss/legacy/processed/2023/1/310956/LMO (pos 166) Mn-K XES EXAFS 0001 vh_roi1.dat',
 '/nsls2/data/iss/legacy/processed/2023/1/310956/LMO (pos 165) Mn-K XES EXAFS 0001 vh_roi1.dat',
 '/nsls2/data/iss/legacy/processed/2023/1/310956/LMO (pos 164) Mn-K XES EXAFS 0001 vh_roi1.dat',
 '/nsls2/data/iss/legacy/processed/2023/1/310956/LMO (pos 163) Mn-K XES EXAFS 0001 vh_roi1.dat',
 '/nsls2/data/iss/legacy/processed/2023/1/310956/LMO (pos 162) Mn-K XES EXAFS 0001 vh_roi1.dat',
 '/nsls2/data/iss/legacy/processed/2023/1/310956/LMO (pos 161) Mn-K XES EXAFS 0001 vh_roi1.dat',
 '/nsls2/data/iss/legacy/processed/2023/1/310956/LMO (pos 160) Mn-K XES EXAFS 0001 vh_roi1.dat',
 '/nsls2/data/iss/legacy/processed/2023/1/310956/LMO (pos 159) Mn-K XES EXAFS 0001 vh_roi1.dat',
 '/nsls2/data/iss/legacy/processed/2023/1/310956/LMO (pos 158) Mn-K XES EXAFS 0001 vh_roi1.dat',
 '/nsls2/data/iss/legacy/processed/2023/1/310956/LMO (pos 157) Mn-K XES EXAFS 0001 vh_roi1.dat',
 '/nsls2/data/iss/legacy/processed/2023/1/310956/LMO (pos 156) Mn-K XES EXAFS 0001 vh_roi1.dat',
 '/nsls2/data/iss/legacy/processed/2023/1/310956/LMO (pos 155) Mn-K XES EXAFS 0001 vh_roi1.dat',
 '/nsls2/data/iss/legacy/processed/2023/1/310956/LMO (pos 154) Mn-K XES EXAFS 0001 vh_roi1.dat',
 '/nsls2/data/iss/legacy/processed/2023/1/310956/LMO (pos 153) Mn-K XES EXAFS 0001 vh_roi1.dat',
 '/nsls2/data/iss/legacy/processed/2023/1/310956/LMO (pos 152) Mn-K XES EXAFS 0001 vh_roi1.dat',
 '/nsls2/data/iss/legacy/processed/2023/1/310956/LMO (pos 151) Mn-K XES EXAFS 0001 vh_roi1.dat',
 '/nsls2/data/iss/legacy/processed/2023/1/310956/LMO (pos 150) Mn-K XES EXAFS 0001 vh_roi1.dat',
 '/nsls2/data/iss/legacy/processed/2023/1/310956/LMO (pos 149) Mn-K XES EXAFS 0001 vh_roi1.dat',
 '/nsls2/data/iss/legacy/processed/2023/1/310956/LMO (pos 148) Mn-K XES EXAFS 0001 vh_roi1.dat',
 '/nsls2/data/iss/legacy/processed/2023/1/310956/LMO (pos 147) Mn-K XES EXAFS 0001 vh_roi1.dat',
 '/nsls2/data/iss/legacy/processed/2023/1/310956/LMO (pos 146) Mn-K XES EXAFS 0001 vh_roi1.dat',
 '/nsls2/data/iss/legacy/processed/2023/1/310956/LMO (pos 145) Mn-K XES EXAFS 0001 vh_roi1.dat',
 '/nsls2/data/iss/legacy/processed/2023/1/310956/LMO (pos 144) Mn-K XES EXAFS 0001 vh_roi1.dat',
 '/nsls2/data/iss/legacy/processed/2023/1/310956/LMO (pos 143) Mn-K XES EXAFS 0001 vh_roi1.dat',
 '/nsls2/data/iss/legacy/processed/2023/1/310956/LMO (pos 142) Mn-K XES EXAFS 0001 vh_roi1.dat',
 '/nsls2/data/iss/legacy/processed/2023/1/310956/LMO (pos 141) Mn-K XES EXAFS 0001 vh_roi1.dat',
 '/nsls2/data/iss/legacy/processed/2023/1/310956/LMO (pos 140) Mn-K XES EXAFS 0001 vh_roi1.dat',
 '/nsls2/data/iss/legacy/processed/2023/1/310956/LMO (pos 139) Mn-K XES EXAFS 0001 vh_roi1.dat',
 '/nsls2/data/iss/legacy/processed/2023/1/310956/LMO (pos 138) Mn-K XES EXAFS 0001 vh_roi1.dat',
 '/nsls2/data/iss/legacy/processed/2023/1/310956/LMO (pos 137) Mn-K XES EXAFS 0001 vh_roi1.dat',
 '/nsls2/data/iss/legacy/processed/2023/1/310956/LMO (pos 135) Mn-K XES EXAFS 0001 vh_roi1.dat',
 '/nsls2/data/iss/legacy/processed/2023/1/310956/LMO (pos 134) Mn-K XES EXAFS 0001 vh_roi1.dat',
 '/nsls2/data/iss/legacy/processed/2023/1/310956/LMO (pos 133) Mn-K XES EXAFS 0001 vh_roi1.dat',
 '/nsls2/data/iss/legacy/processed/2023/1/310956/LMO (pos 132) Mn-K XES EXAFS 0001 vh_roi1.dat',
 '/nsls2/data/iss/legacy/processed/2023/1/310956/LMO (pos 131) Mn-K XES EXAFS 0001 vh_roi1.dat',
 '/nsls2/data/iss/legacy/processed/2023/1/310956/LMO (pos 130) Mn-K XES EXAFS 0001 vh_roi1.dat',
 '/nsls2/data/iss/legacy/processed/2023/1/310956/LMO (pos 129) Mn-K XES EXAFS 0001 vh_roi1.dat',
 '/nsls2/data/iss/legacy/processed/2023/1/310956/LMO (pos 128) Mn-K XES EXAFS 0001 vh_roi1.dat',
 '/nsls2/data/iss/legacy/processed/2023/1/310956/LMO (pos 127) Mn-K XES EXAFS 0001 vh_roi1.dat',
 '/nsls2/data/iss/legacy/processed/2023/1/310956/LMO (pos 126) Mn-K XES EXAFS 0001 vh_roi1.dat',
 '/nsls2/data/iss/legacy/processed/2023/1/310956/LMO (pos 125) Mn-K XES EXAFS 0001 vh_roi1.dat',
 '/nsls2/data/iss/legacy/processed/2023/1/310956/LMO (pos 124) Mn-K XES EXAFS 0001 vh_roi1.dat',
 '/nsls2/data/iss/legacy/processed/2023/1/310956/LMO (pos 123) Mn-K XES EXAFS 0001 vh_roi1.dat',
 '/nsls2/data/iss/legacy/processed/2023/1/310956/LMO (pos 122) Mn-K XES EXAFS 0001 vh_roi1.dat',
 '/nsls2/data/iss/legacy/processed/2023/1/310956/LMO (pos 121) Mn-K XES EXAFS 0001 vh_roi1.dat',
 '/nsls2/data/iss/legacy/processed/2023/1/310956/LMO (pos 120) Mn-K XES EXAFS 0001 vh_roi1.dat',
 '/nsls2/data/iss/legacy/processed/2023/1/310956/LMO (pos 119) Mn-K XES EXAFS 0001 vh_roi1.dat',
 '/nsls2/data/iss/legacy/processed/2023/1/310956/LMO (pos 118) Mn-K XES EXAFS 0001 vh_roi1.dat',
 '/nsls2/data/iss/legacy/processed/2023/1/310956/LMO (pos 117) Mn-K XES EXAFS 0001 vh_roi1.dat',
 '/nsls2/data/iss/legacy/processed/2023/1/310956/LMO (pos 116) Mn-K XES EXAFS 0001 vh_roi1.dat',
 '/nsls2/data/iss/legacy/processed/2023/1/310956/LMO (pos 115) Mn-K XES EXAFS 0001 vh_roi1.dat',
 '/nsls2/data/iss/legacy/processed/2023/1/310956/LMO (pos 114) Mn-K XES EXAFS 0001 vh_roi1.dat',
 '/nsls2/data/iss/legacy/processed/2023/1/310956/LMO (pos 113) Mn-K XES EXAFS 0001 vh_roi1.dat',
 '/nsls2/data/iss/legacy/processed/2023/1/310956/LMO (pos 112) Mn-K XES EXAFS 0001 vh_roi1.dat',
 '/nsls2/data/iss/legacy/processed/2023/1/310956/LMO (pos 111) Mn-K XES EXAFS 0001 vh_roi1.dat',
 '/nsls2/data/iss/legacy/processed/2023/1/310956/LMO (pos 110) Mn-K XES EXAFS 0001 vh_roi1.dat',
 '/nsls2/data/iss/legacy/processed/2023/1/310956/LMO (pos 109) Mn-K XES EXAFS 0001 vh_roi1.dat',
 '/nsls2/data/iss/legacy/processed/2023/1/310956/LMO (pos 108) Mn-K XES EXAFS 0001 vh_roi1.dat',
 '/nsls2/data/iss/legacy/processed/2023/1/310956/LMO (pos 107) Mn-K XES EXAFS 0001 vh_roi1.dat',
 '/nsls2/data/iss/legacy/processed/2023/1/310956/LMO (pos 106) Mn-K XES EXAFS 0001 vh_roi1.dat',
 '/nsls2/data/iss/legacy/processed/2023/1/310956/LMO (pos 105) Mn-K XES EXAFS 0001 vh_roi1.dat',
 '/nsls2/data/iss/legacy/processed/2023/1/310956/LMO (pos 104) Mn-K XES EXAFS 0001 vh_roi1.dat',
 '/nsls2/data/iss/legacy/processed/2023/1/310956/LMO (pos 103) Mn-K XES EXAFS 0001 vh_roi1.dat',
 '/nsls2/data/iss/legacy/processed/2023/1/310956/LMO (pos 102) Mn-K XES EXAFS 0001 vh_roi1.dat',
 '/nsls2/data/iss/legacy/processed/2023/1/310956/LMO (pos 101) Mn-K XES EXAFS 0001 vh_roi1.dat',
 '/nsls2/data/iss/legacy/processed/2023/1/310956/LMO (pos 100) Mn-K XES EXAFS 0001 vh_roi1.dat',
 '/nsls2/data/iss/legacy/processed/2023/1/310956/LMO (pos 099) Mn-K XES EXAFS 0001 vh_roi1.dat',
 '/nsls2/data/iss/legacy/processed/2023/1/310956/LMO (pos 098) Mn-K XES EXAFS 0001 vh_roi1.dat',
 '/nsls2/data/iss/legacy/processed/2023/1/310956/LMO (pos 097) Mn-K XES EXAFS 0001 vh_roi1.dat',
 '/nsls2/data/iss/legacy/processed/2023/1/310956/LMO (pos 096) Mn-K XES EXAFS 0001 vh_roi1.dat',
 '/nsls2/data/iss/legacy/processed/2023/1/310956/LMO (pos 095) Mn-K XES EXAFS 0001 vh_roi1.dat',
 '/nsls2/data/iss/legacy/processed/2023/1/310956/LMO (pos 094) Mn-K XES EXAFS 0001 vh_roi1.dat',
 '/nsls2/data/iss/legacy/processed/2023/1/310956/LMO (pos 093) Mn-K XES EXAFS 0001 vh_roi1.dat',
 '/nsls2/data/iss/legacy/processed/2023/1/310956/LMO (pos 092) Mn-K XES EXAFS 0001 vh_roi1.dat',
 '/nsls2/data/iss/legacy/processed/2023/1/310956/LMO (pos 091) Mn-K XES EXAFS 0001 vh_roi1.dat',
 '/nsls2/data/iss/legacy/processed/2023/1/310956/LMO (pos 090) Mn-K XES EXAFS 0001 vh_roi1.dat',
 '/nsls2/data/iss/legacy/processed/2023/1/310956/LMO (pos 089) Mn-K XES EXAFS 0001 vh_roi1.dat',
 '/nsls2/data/iss/legacy/processed/2023/1/310956/LMO (pos 088) Mn-K XES EXAFS 0001 vh_roi1.dat',
 '/nsls2/data/iss/legacy/processed/2023/1/310956/LMO (pos 087) Mn-K XES EXAFS 0001 vh_roi1.dat',
 '/nsls2/data/iss/legacy/processed/2023/1/310956/LMO (pos 086) Mn-K XES EXAFS 0001 vh_roi1.dat',
 '/nsls2/data/iss/legacy/processed/2023/1/310956/LMO (pos 085) Mn-K XES EXAFS 0001 vh_roi1.dat',
 '/nsls2/data/iss/legacy/processed/2023/1/310956/LMO (pos 084) Mn-K XES EXAFS 0001 vh_roi1.dat',
 '/nsls2/data/iss/legacy/processed/2023/1/310956/LMO (pos 083) Mn-K XES EXAFS 0001 vh_roi1.dat',
 '/nsls2/data/iss/legacy/processed/2023/1/310956/LMO (pos 082) Mn-K XES EXAFS 0001 vh_roi1.dat',
 '/nsls2/data/iss/legacy/processed/2023/1/310956/LMO (pos 081) Mn-K XES EXAFS 0001 vh_roi1.dat',
 '/nsls2/data/iss/legacy/processed/2023/1/310956/LMO (pos 080) Mn-K XES EXAFS 0001 vh_roi1.dat',
 '/nsls2/data/iss/legacy/processed/2023/1/310956/LMO (pos 079) Mn-K XES EXAFS 0001 vh_roi1.dat',
 '/nsls2/data/iss/legacy/processed/2023/1/310956/LMO (pos 078) Mn-K XES EXAFS 0001 vh_roi1.dat',
 '/nsls2/data/iss/legacy/processed/2023/1/310956/LMO (pos 077) Mn-K XES EXAFS 0001 vh_roi1.dat',
 '/nsls2/data/iss/legacy/processed/2023/1/310956/LMO (pos 076) Mn-K XES EXAFS 0001 vh_roi1.dat',
 '/nsls2/data/iss/legacy/processed/2023/1/310956/LMO (pos 075) Mn-K XES EXAFS 0001 vh_roi1.dat',
 '/nsls2/data/iss/legacy/processed/2023/1/310956/LMO (pos 074) Mn-K XES EXAFS 0001 vh_roi1.dat',
 '/nsls2/data/iss/legacy/processed/2023/1/310956/LMO (pos 073) Mn-K XES EXAFS 0001 vh_roi1.dat',
 '/nsls2/data/iss/legacy/processed/2023/1/310956/LMO (pos 072) Mn-K XES EXAFS 0001 vh_roi1.dat',
 '/nsls2/data/iss/legacy/processed/2023/1/310956/LMO (pos 071) Mn-K XES EXAFS 0001 vh_roi1.dat',
 '/nsls2/data/iss/legacy/processed/2023/1/310956/LMO (pos 070) Mn-K XES EXAFS 0001 vh_roi1.dat',
 '/nsls2/data/iss/legacy/processed/2023/1/310956/LMO (pos 069) Mn-K XES EXAFS 0001 vh_roi1.dat',
 '/nsls2/data/iss/legacy/processed/2023/1/310956/LMO (pos 068) Mn-K XES EXAFS 0001 vh_roi1.dat',
 '/nsls2/data/iss/legacy/processed/2023/1/310956/LMO (pos 067) Mn-K XES EXAFS 0001 vh_roi1.dat',
 '/nsls2/data/iss/legacy/processed/2023/1/310956/LMO (pos 066) Mn-K XES EXAFS 0001 vh_roi1.dat',
 '/nsls2/data/iss/legacy/processed/2023/1/310956/LMO (pos 065) Mn-K XES EXAFS 0001 vh_roi1.dat',
 '/nsls2/data/iss/legacy/processed/2023/1/310956/LMO (pos 064) Mn-K XES EXAFS 0001 vh_roi1.dat',
 '/nsls2/data/iss/legacy/processed/2023/1/310956/LMO (pos 063) Mn-K XES EXAFS 0001 vh_roi1.dat',
 '/nsls2/data/iss/legacy/processed/2023/1/310956/LMO (pos 062) Mn-K XES EXAFS 0001 vh_roi1.dat',
 '/nsls2/data/iss/legacy/processed/2023/1/310956/LMO (pos 061) Mn-K XES EXAFS 0001 vh_roi1.dat',
 '/nsls2/data/iss/legacy/processed/2023/1/310956/LMO (pos 060) Mn-K XES EXAFS 0001 vh_roi1.dat',
 '/nsls2/data/iss/legacy/processed/2023/1/310956/LMO (pos 059) Mn-K XES EXAFS 0001 vh_roi1.dat',
 '/nsls2/data/iss/legacy/processed/2023/1/310956/LMO (pos 058) Mn-K XES EXAFS 0001 vh_roi1.dat',
 '/nsls2/data/iss/legacy/processed/2023/1/310956/LMO (pos 057) Mn-K XES EXAFS 0001 vh_roi1.dat',
 '/nsls2/data/iss/legacy/processed/2023/1/310956/LMO (pos 056) Mn-K XES EXAFS 0001 vh_roi1.dat',
 '/nsls2/data/iss/legacy/processed/2023/1/310956/LMO (pos 055) Mn-K XES EXAFS 0001 vh_roi1.dat',
 '/nsls2/data/iss/legacy/processed/2023/1/310956/LMO (pos 054) Mn-K XES EXAFS 0001 vh_roi1.dat',
 '/nsls2/data/iss/legacy/processed/2023/1/310956/LMO (pos 053) Mn-K XES EXAFS 0001 vh_roi1.dat',
 '/nsls2/data/iss/legacy/processed/2023/1/310956/LMO (pos 052) Mn-K XES EXAFS 0001 vh_roi1.dat',
 '/nsls2/data/iss/legacy/processed/2023/1/310956/LMO (pos 051) Mn-K XES EXAFS 0001 vh_roi1.dat',
 '/nsls2/data/iss/legacy/processed/2023/1/310956/LMO (pos 050) Mn-K XES EXAFS 0001 vh_roi1.dat',
 '/nsls2/data/iss/legacy/processed/2023/1/310956/LMO (pos 049) Mn-K XES EXAFS 0001 vh_roi1.dat',
 '/nsls2/data/iss/legacy/processed/2023/1/310956/LMO (pos 048) Mn-K XES EXAFS 0001 vh_roi1.dat',
 '/nsls2/data/iss/legacy/processed/2023/1/310956/LMO (pos 047) Mn-K XES EXAFS 0001 vh_roi1.dat',
 '/nsls2/data/iss/legacy/processed/2023/1/310956/LMO (pos 046) Mn-K XES EXAFS 0001 vh_roi1.dat',
 '/nsls2/data/iss/legacy/processed/2023/1/310956/LMO (pos 045) Mn-K XES EXAFS 0001 vh_roi1.dat',
 '/nsls2/data/iss/legacy/processed/2023/1/310956/LMO (pos 044) Mn-K XES EXAFS 0001 vh_roi1.dat',
 '/nsls2/data/iss/legacy/processed/2023/1/310956/LMO (pos 043) Mn-K XES EXAFS 0001 vh_roi1.dat',
 '/nsls2/data/iss/legacy/processed/2023/1/310956/LMO (pos 042) Mn-K XES EXAFS 0001 vh_roi1.dat',
 '/nsls2/data/iss/legacy/processed/2023/1/310956/LMO (pos 041) Mn-K XES EXAFS 0001 vh_roi1.dat',
 '/nsls2/data/iss/legacy/processed/2023/1/310956/LMO (pos 040) Mn-K XES EXAFS 0001 vh_roi1.dat',
 '/nsls2/data/iss/legacy/processed/2023/1/310956/LMO (pos 039) Mn-K XES EXAFS 0001 vh_roi1.dat',
 '/nsls2/data/iss/legacy/processed/2023/1/310956/LMO (pos 038) Mn-K XES EXAFS 0001 vh_roi1.dat',
 '/nsls2/data/iss/legacy/processed/2023/1/310956/LMO (pos 037) Mn-K XES EXAFS 0001 vh_roi1.dat',
 '/nsls2/data/iss/legacy/processed/2023/1/310956/LMO (pos 036) Mn-K XES EXAFS 0001 vh_roi1.dat',
 '/nsls2/data/iss/legacy/processed/2023/1/310956/LMO (pos 035) Mn-K XES EXAFS 0001 vh_roi1.dat',
 '/nsls2/data/iss/legacy/processed/2023/1/310956/LMO (pos 034) Mn-K XES EXAFS 0001 vh_roi1.dat',
 '/nsls2/data/iss/legacy/processed/2023/1/310956/LMO (pos 033) Mn-K XES EXAFS 0001 vh_roi1.dat',
 '/nsls2/data/iss/legacy/processed/2023/1/310956/LMO (pos 032) Mn-K XES EXAFS 0001 vh_roi1.dat',
 '/nsls2/data/iss/legacy/processed/2023/1/310956/LMO (pos 031) Mn-K XES EXAFS 0001 vh_roi1.dat',
 '/nsls2/data/iss/legacy/processed/2023/1/310956/LMO (pos 030) Mn-K XES EXAFS 0001 vh_roi1.dat',
 '/nsls2/data/iss/legacy/processed/2023/1/310956/LMO (pos 029) Mn-K XES EXAFS 0001 vh_roi1.dat',
 '/nsls2/data/iss/legacy/processed/2023/1/310956/LMO (pos 028) Mn-K XES EXAFS 0001 vh_roi1.dat',
 '/nsls2/data/iss/legacy/processed/2023/1/310956/LMO (pos 027) Mn-K XES EXAFS 0001 vh_roi1.dat',
 '/nsls2/data/iss/legacy/processed/2023/1/310956/LMO (pos 026) Mn-K XES EXAFS 0001 vh_roi1.dat',
 '/nsls2/data/iss/legacy/processed/2023/1/310956/LMO (pos 025) Mn-K XES EXAFS 0001 vh_roi1.dat',
 '/nsls2/data/iss/legacy/processed/2023/1/310956/LMO (pos 024) Mn-K XES EXAFS 0001 vh_roi1.dat',
 '/nsls2/data/iss/legacy/processed/2023/1/310956/LMO (pos 023) Mn-K XES EXAFS 0001 vh_roi1.dat',
 '/nsls2/data/iss/legacy/processed/2023/1/310956/LMO (pos 022) Mn-K XES EXAFS 0001 vh_roi1.dat',
 '/nsls2/data/iss/legacy/processed/2023/1/310956/LMO (pos 021) Mn-K XES EXAFS 0001 vh_roi1.dat',
 '/nsls2/data/iss/legacy/processed/2023/1/310956/LMO (pos 020) Mn-K XES EXAFS 0001 vh_roi1.dat',
 '/nsls2/data/iss/legacy/processed/2023/1/310956/LMO (pos 019) Mn-K XES EXAFS 0001 vh_roi1.dat',
 '/nsls2/data/iss/legacy/processed/2023/1/310956/LMO (pos 018) Mn-K XES EXAFS 0001 vh_roi1.dat',
 '/nsls2/data/iss/legacy/processed/2023/1/310956/LMO (pos 017) Mn-K XES EXAFS 0001 vh_roi1.dat',
 '/nsls2/data/iss/legacy/processed/2023/1/310956/LMO (pos 016) Mn-K XES EXAFS 0001 vh_roi1.dat',
 '/nsls2/data/iss/legacy/processed/2023/1/310956/LMO (pos 015) Mn-K XES EXAFS 0001 vh_roi1.dat',
 '/nsls2/data/iss/legacy/processed/2023/1/310956/LMO (pos 014) Mn-K XES EXAFS 0001 vh_roi1.dat',
 '/nsls2/data/iss/legacy/processed/2023/1/310956/LMO (pos 013) Mn-K XES EXAFS 0001 vh_roi1.dat',
 '/nsls2/data/iss/legacy/processed/2023/1/310956/LMO (pos 012) Mn-K XES EXAFS 0001 vh_roi1.dat',
 '/nsls2/data/iss/legacy/processed/2023/1/310956/LMO (pos 011) Mn-K XES EXAFS 0001 vh_roi1.dat',
 '/nsls2/data/iss/legacy/processed/2023/1/310956/LMO (pos 010) Mn-K XES EXAFS 0001 vh_roi1.dat',
 '/nsls2/data/iss/legacy/processed/2023/1/310956/LMO (pos 009) Mn-K XES EXAFS 0001 vh_roi1.dat',
 '/nsls2/data/iss/legacy/processed/2023/1/310956/LMO (pos 008) Mn-K XES EXAFS 0001 vh_roi1.dat',
 '/nsls2/data/iss/legacy/processed/2023/1/310956/LMO (pos 007) Mn-K XES EXAFS 0001 vh_roi1.dat',
 '/nsls2/data/iss/legacy/processed/2023/1/310956/LMO (pos 006) Mn-K XES EXAFS 0001 vh_roi1.dat',
 '/nsls2/data/iss/legacy/processed/2023/1/310956/LMO (pos 005) Mn-K XES EXAFS 0001 vh_roi1.dat',
 '/nsls2/data/iss/legacy/processed/2023/1/310956/LMO (pos 004) Mn-K XES EXAFS 0001 vh_roi1.dat',
 '/nsls2/data/iss/legacy/processed/2023/1/310956/LMO (pos 003) Mn-K XES EXAFS 0001 vh_roi1.dat',
 '/nsls2/data/iss/legacy/processed/2023/1/310956/LMO (pos 002) Mn-K XES EXAFS 0001 vh_roi1.dat',
 '/nsls2/data/iss/legacy/processed/2023/1/310956/LMO (pos 001) Mn-K XES EXAFS 0001 vh_roi1.dat'
]


rixs = []
for f in files:
    print(f)
    df, _ = load_binned_df_from_file(f)
    # print(energy_in.size)
    bla = df.columns.tolist()
    _rixs = -df[bla[1::2]].values / df[bla[2::2]].values
    if _rixs.shape[1] == 504:
        rixs.append(_rixs)

energy_out = df['energy'].values
energy_in = np.array([float('.'.join(i.split('_')[2:4])) for i in bla[1::2]])
rixs = np.array(rixs) / 1e6

rixs_av = np.mean(rixs, axis=0)

plt.figure(1, clear=True)
plt.contourf(energy_out, energy_in, rixs_av.T, 50)

n_out = [230, 245]
plt.vlines(energy_out[n_out], energy_in.min(), energy_in.max(), colors='r')


herfds = rixs_av[n_out, :].T
herfds -= np.mean(herfds[20:60, :], axis=0)[None, :]
herfds /= np.mean(herfds[400:, :], axis=0)[None, :]

plt.figure(2, clear=True)
plt.plot(energy_in, herfds)


e_in_min = 6570
e_out_min = 6472.5
e_out_max = 6494.5

e_out_glitch_mask = (energy_out <= 6477.25) | (energy_out >= 6477.75)

e_in_mask = energy_in >= e_in_min
e_out_mask = (energy_out >= e_out_min) & (energy_out <= e_out_max) & e_out_glitch_mask

rixs_av_roi = rixs_av[np.ix_(e_out_mask, e_in_mask)]
uu, ss, vvT = np.linalg.svd(rixs_av_roi)
vv = vvT.T


plt.figure(3, clear=True)
plt.subplot(221)
plt.contourf(energy_out[e_out_mask], energy_in[e_in_mask], rixs_av_roi.T, 50, vmin=rixs_av_roi.min(), vmax=rixs_av_roi.max()*0.5)

i_cmp = [0, 1, 2]

plt.subplot(222)
plt.semilogy(ss, 'ks-')

plt.subplot(223)
plt.plot(energy_out[e_out_mask], uu[:, i_cmp])

plt.subplot(224)
plt.plot(energy_in[e_in_mask], vv[:, i_cmp])


