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

uids = []
for item in x.list_data.selectedItems():
    fname = os.path.join(x.working_folder, item.text())
    df, header = load_binned_df_from_file(fname)
    print(fname, df.shape)
    uid_idx1 = header.find('Scan.uid:') + 10
    uid_idx2 = header.find('\n', header.find('Scan.uid:'))
    uid = header[uid_idx1: uid_idx2]
    uids.append(uid)



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
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
mpl.use("TkAgg")
from scipy.interpolate import CubicSpline
from xas.file_io import load_binned_df_from_file
from xas.energy_calibration import compute_energy_shift_and_broadening_between_spectra, fine_convolution_grid
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

# df = dfs[0]
# energy, mu = df["energy"], df["mur"]
# e0 = 24350
# de = 200

# cs = CubicSpline(energy_ref, mu_ref)

# roi_mask = (energy > (e0 - de / 2)) & (energy < (e0 + de / 2))
# energy_roi = energy[roi_mask]
# mu_roi = mu[roi_mask]

# fine_grid_energy_ref = fine_convolution_grid(energy_ref, sigma=0.01)
# fine_grid_mu_ref = cs(fine_grid_energy_ref)

plt.figure(1, clear=True)
for i, df in enumerate(dfs):
    shift, sigma, energy_fit, mu_fit = compute_energy_shift_and_broadening_between_spectra(df['energy'].values, df['mur'].values,
                                                                                           energy_ref - 14.85, mu_ref, e0=24350, de=200)
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
for i in range(n1+100, n1 + 300):
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
## _to_do: add muf/mut/mur to the df calculation; also check if mut mur muf are actually good names for the columns
for f in files:

    md, df = read_xas_data_from_json(os.path.join(outpath, f))
    md['mu_quality'] = check_scan(df, md)
    print(f, md['mu_quality'])

    md['tag'] = 'test'
    md['tag_version'] = '0.0'
    md['tag_comment'] = 'test quality check'
    x = c.write_dataframe(df, md)

###############
import xarray as xr
from tiled.client import from_profile
c_nsls2 = from_profile('nsls2', username='username')
c = c_nsls2['iss']['raw']
run = c[-1]
t = run['apb_stream']['data']['apb_stream'].read().squeeze()
data = np.array([list(i) for i in t])
channels = ['i0', 'itrans', 'iref', 'ifluo', 'aux1', 'aux2', 'aux3', 'aux4']
arr = xr.DataArray(data[:, 1:], coords=[data[:, 0], channels], dims=['time', 'channel'])

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
#########




print(f'with APB and SHUTTER: {np.mean([1.6482348442077637, 1.7039844989776611, 1.6537742614746094, 1.7817823886871338, 1.7271864414215088])}')
print(f'with APB: {np.mean([1.213865041732788, 1.268347978591919, 1.3224411010742188, 1.2643957138061523, 1.3515892028808594])}')
print(f'without APB: {np.mean([1.0586943626403809, 1.1041831970214844, 1.0835952758789062, 1.1153781414031982, 1.138744592666626])}')


######

from xas.db_io import load_apb_dataset_from_db, load_apb_trig_dataset_from_db




def plot_trace_for_uid(uid, offset=0):
    hdr = db[uid]
    t_apb = hdr.table(stream_name='apb_stream', fill=True).apb_stream[1]
    t_trig = hdr.table(stream_name='apb_trigger_pil100k', fill=True).apb_trigger_pil100k[1]

    t_min = t_apb[:, 0].min()

    apb_trace = t_apb[:, 6] / 1e6
    apb_trace2 = t_apb[:, 5] / 1e6
    # apb_trace -= np.percentile(apb_trace, 1)
    # apb_trace /= np.percentile(apb_trace, 99)


    # plt.plot(t_apb[:, 0], t_apb[:, 5] / 1e6)
    plt.plot(t_apb[:, 0] - t_min, apb_trace - offset, '-')
    plt.plot(t_apb[:, 0] - t_min, apb_trace2 - offset, '-')
    # plt.plot(t_trig[:, 0] - t_min, t_trig[:, 1] - offset, 'r-')



y_shift = 1.25
plt.figure(1, clear=True)
for i, uid in enumerate(uids):
    plot_trace_for_uid(uid, offset=i*y_shift)
    plt.text(0.2, -i*y_shift + y_shift*0.25, f'({i+1}) {uid}', ha='center', va='center')

# plt.legend()
plt.xlim(0, 1)

plot_trace_for_uid(-1, 0)

############

uids = ['452dc667-f81a-4daf-be60-e6e7d486211d',
        '2c5c0b8f-01e9-4e22-8748-d582dc75e8c0',
        'd178dd6b-531a-4780-8fa4-101a7b8569b8',
        '262acc5b-2941-4284-b01e-ba7bcb4a2122',
        '3572f520-1777-417e-a0b6-685ff8267761',
        '3a8af009-5213-40e7-9340-0de79b34551d',
        '87ec51b5-3848-4dfe-ba27-320418258ffe',
        '4e00831c-9163-4715-9f32-e413ec1b889a',
        'b4090950-b187-4416-b611-ca0d3ec8c8c7',
        'e23bcf19-d8c0-4fac-997f-d0f7dfe1296b'
        ]
# uid = '452dc667-f81a-4daf-be60-e6e7d486211d'

plt.figure(1, clear=True)

for i, uid in enumerate(uids):
    hdr = db[uid]

    stream_name='apb_trigger_pil100k'
    data_trig_raw = list(hdr.data(stream_name=stream_name, field=stream_name))[0]

    # from xas.db_io import load_apb_dataset_from_db, load_apb_trig_dataset_from_db
    df_apb, _, _ = load_apb_dataset_from_db(db, uid)

    derived_timestamps_trig = load_apb_trig_dataset_from_db(db, uid, stream_name=stream_name)

    t_apb = df_apb.timestamp.values
    s_apb = df_apb.aux2.values

    s_apb -= s_apb.min()
    s_apb /= s_apb.max()

    plt.plot(data_trig_raw[:, 0] - t_apb.min(), data_trig_raw[:, 1] - i * 1.15)
    plt.plot(t_apb - t_apb.min(), s_apb - i * 1.15)
    plt.vlines(derived_timestamps_trig.tolist() - t_apb.min(), 0 - i * 1.15, 1.1 - i * 1.15, colors='k')


##########


unique_samples = get_unique_values_for_metadata_key(node, 'sample_name')
unique_scans = get_unique_values_for_metadata_key(node, 'scan_name')
# unique_quality = get_unique_values_for_metadata_key(node, 'scan_quality')

denis_tree = {}

for _sample in unique_samples:
    denis_tree[_sample] = {}
    for _scan in unique_scans:
        denis_tree[_sample][_scan] = {}
        # for _quality_i, _quality in enumerate(unique_quality):
        print(_sample, _scan,)
        _node = node.search(Key('sample_name') == _sample).search(Key('scan_name') == _scan)
        # _uid_indexes = list(range(len(_uids)))
        denis_tree[_sample][_scan] = {'node': _node}#, 'index': _uid_indexes, 'checkbox_id': []}




#########
from matplotlib import pyplot as plt
plt.ion()



# plt.figure(1, clear=True)
# for i in range(10):
#     # df = load_apb_dataset_from_tiled(run)
#     # plt.plot(df.timestamp, df.i0)
#     df = read_hhm_encoder_dataset_from_tiled(run)
#     plt.plot(df.timestamp, df.energy)


apb_df = load_apb_dataset_from_tiled(run)



energy_df = load_hhm_encoder_dataset_from_tiled(run)

apb_dict = translate_dataset(apb_df)
energy_dict = translate_dataset(energy_df, columns=['energy'])

raw_dict = {**apb_dict, **energy_dict}

dataset = raw_dict

interpolated_dataset = interpolate(dataset)

rebinned_dataset = rebin(interpolated_dataset, run.metadata['start']['e0'])

# interpolated_dataset = {}
key= 'mufluor'
df = dataset[key]
# df['timestamp_bin'] = pd.cut(df['timestamp'], bins=timestamp_edges, labels=timestamp)
# df_mean = df.groupby('timestamp_bin').mean(numeric_only=False)
# interpolated_dataset[key] = df_mean[key]

plt.figure(1, clear=True)
# plt.plot(df['timestamp'], df[key], '.-')
# plt.plot(interpolated_dataset['timestamp'], interpolated_dataset[key], '.')
plt.plot(rebinned_dataset['timestamp'], rebinned_dataset[key] - rebinned_dataset['iff'] / rebinned_dataset['i0'], '.-')
# plt.plot(rebinned_dataset['timestamp'], , '.-')


####

plt.figure(1, clear=True)
plt.plot(apb_df['timestamp'], apb_df['mufluor'])
# plt.plot(apb_df['timestamp'], apb_df['iff'] / apb_df['i0'])

######

# plt.figure(1, clear=True)
# for i, (k, ds) in enumerate(raw_dict.items()):
#     plt.plot(ds.timestamp.values, np.ones(ds.timestamp.values.size) - i)

# imports

from ophyd import EpicsSignal
import time as ttime
import numpy as np
import matplotlib.pyplot as plt
plt.ion()

# ophyd signals
apb_fa_filename_bin = EpicsSignal('XF:08IDB-CT{PBA:1}:FA:Stream:Bin:File-SP', name='apb_filename_bin')
apb_fa_filename_txt = EpicsSignal('XF:08IDB-CT{PBA:1}:FA:Stream:Txt:File-SP', name='apb_filename_txt')
apb_fa_stream = EpicsSignal('XF:08IDB-CT{PBA:1}:FA:Stream-SP', name='apb_stream')

apb_trig_filename = EpicsSignal('XF:08IDB-CT{PBA:1}:Pulse:2:Filename-SP', name='apb_trig_filename')
apb_trig_stream = EpicsSignal('XF:08IDB-CT{PBA:1}:Pulse:2:Stream:Mode-SP', name='apb_trig_stream')
apb_trig_acquire = EpicsSignal('XF:08IDB-CT{PBA:1}:Pulse:2:Mode-SP', name='apb_trig_acquire')

def set_filenames(n):
    fname = f'/nsls2/data/iss/legacy/raw/apb/2023/04/05/trigger_test2_{n:02d}'
    fname_fa_bin = fname + '_fa.bin'
    fname_fa_txt = fname + '_fa.txt'
    apb_fa_filename_bin.put(fname_fa_bin)
    apb_fa_filename_txt.put(fname_fa_txt)

    fname_trig = fname + '_trig.bin'
    apb_trig_filename.put(fname_trig)
    return fname_fa_bin, fname_trig

def acquisition_sequence():
    ttime.sleep(0.2) # small delay before start to let the file PVs to write
    apb_fa_stream.put(1)
    ttime.sleep(1) # sleep a bit to ensure that streaming actually begins

    apb_trig_stream.put(1)
    ttime.sleep(0.5) # small delay to ensure that streaming of trigger is commenced before pulsing begins
    apb_trig_acquire.put(2)

    ttime.sleep(1) # let it run for a bit

    apb_fa_stream.put(0)
    apb_trig_stream.put(0)
    apb_trig_acquire.put(0)
    ttime.sleep(0.5) # small delay to ensure that files finish writing

# def acquisition_sequence():
#     ttime.sleep(0.2) # small delay before start to let the file PVs to write
#
#     apb_trig_acquire.put(1)  #Set Pulser acquire state to “FA”
#
#     ttime.sleep(0.2) # small delay before start set Pulser acquire state to “FA”
#
#     apb_trig_stream.put(1)   #Enable Pulser streaming to file.  Nothing will happen until FA streaming starts
#
#     ttime.sleep(0.2) # small delay before start set Pulser streaming to “Enable”
#     apb_fa_stream.put(1)  #Start the FA streaming.  This will immediately start the Pulser streaming as well.
#
#     ttime.sleep(3) # let it run for a bit
#
#     apb_fa_stream.put(0)
#     apb_trig_stream.put(0)
#     apb_trig_acquire.put(0)
#     ttime.sleep(0.5) # small delay to ensure that files finish writing

# def acquisition_sequence():
#     ttime.sleep(0.2) # small delay before start to let the file PVs to write
#
#     apb_trig_stream.put(1) #Enable Pulser Streaming.  Pulser is still OFF.
#     ttime.sleep(0.2) # small delay to ensure that streaming is Enabled.
#     apb_trig_acquire.put(2)  #Turn Pulser ON. Pulser transitions will be streamed to file.
#
#     ttime.sleep(0.2) # small delay to ensure that streaming is running.
#     apb_fa_stream.put(1) #Turn FA streaming ON
#
#     ttime.sleep(3) # let it run for a bit
#
#     apb_fa_stream.put(0)
#     apb_trig_stream.put(0)
#     apb_trig_acquire.put(0)
#     ttime.sleep(0.5) # small delay to ensure that files finish writing


# data readout
def read_fa_data(fpath):
    raw_data = np.fromfile(fpath, dtype=np.int32)

    columns = ['timestamp', 'i0', 'it', 'ir', 'iff', 'aux1', 'aux2', 'aux3', 'aux4']
    num_columns = len(columns) + 1
    raw_data = raw_data.reshape((raw_data.size // num_columns, num_columns))

    derived_data = np.zeros((raw_data.shape[0], raw_data.shape[1] - 1))
    derived_data[:, 0] = raw_data[:, -2] + raw_data[:, -1] * 8.0051232 * 1e-9  # Unix timestamp with nanoseconds
    for i in range(num_columns - 2):
        derived_data[:, i + 1] = raw_data[:, i]  # ((raw_data[:, i] ) - Offsets[i]) / Gains[i]

    return derived_data

def read_trig_data(fpath):
    raw_data = np.fromfile(fpath, dtype=np.int32)
    raw_data = raw_data.reshape((raw_data.size // 3, 3))

    derived_data = np.zeros((raw_data.shape[0], 2))
    derived_data[:, 0] = raw_data[:, 1] + raw_data[:, 2] * 8.0051232 * 1e-9  # Unix timestamp with nanoseconds
    derived_data[:, 1] = raw_data[:, 0]
    return derived_data

def get_data(fname_fa_bin, fname_trig):
    data_fa = read_fa_data(fname_fa_bin)
    data_trig = read_trig_data(fname_trig)

    time_fa = data_fa[:, 0]
    signal_fa = data_fa[:, 6]
    t_min = time_fa.min()

    time_trig = data_trig[:, 0]
    signal_trig = data_trig[:, 1]
    return (time_fa - t_min), signal_fa, (time_trig - t_min), signal_trig


# def acquire_data(n=3):
n = 1
plt.figure(1, clear=True)
for i in range(n):
    fname_fa_bin, fname_trig = set_filenames(i + 1)
    acquisition_sequence()
    offset = i * 1.15
    time_fa, signal_fa, time_trig, signal_trig = get_data(fname_fa_bin, fname_trig)
    signal_fa_ = (signal_fa - np.percentile(signal_fa, 1)) / (np.percentile(signal_fa, 99) - np.percentile(signal_fa, 1))
    plt.plot(time_fa, signal_fa_ - offset)
    plt.plot(time_trig, signal_trig - offset)


plt.figure(2)
plt.plot(time_fa, signal_fa*1e-6)

processed_scans = list(range(352382, 354132))
uids = list(db.v2.search({'proposal': '311347'}))

# scan_md = [[db[uid].start['time'], db[uid].start['scan_id']] for uid in uids]
# scan_ids = [db[uid].start['scan_id'] for uid in uids]

# for uid in uids:
for i in range(389939, 389963):
#     try:
    hdr = db[uid]
    if hdr.start['scan_id'] not in processed_scans:
        try:
            process_interpolate_bin(hdr.stop, db)
        except Exception as e:
            print(e)

############

offset = 0
def elasic_plot_func(x, y, x_fit, y_fit, Ecen, fwhm, roi_label='', roi_color='k'):
    global offset
    plt.plot(x, y - offset, '.', color=roi_color)
    plt.plot(x_fit, y_fit - offset, '-', color=roi_color)
    offset += 0.1

motor_pos_elastic = []
fwhm_elastic = []
max_elastic = []

plt.figure(1, clear=True)
offset = 0
for uid in range(389939, 389963 + 1):
    _fwhm, _max= analyze_linewidth_fly_scan(db, uid, x_key='energy', rois=[1], plot_func=elasic_plot_func)
    # _fwhm = analyze_elastic_fly_scan(db, uid, rois=[1], plot_func=elasic_plot_func)
    fwhm_elastic.append(_fwhm)
    max_elastic.append(_max)
    motor_pos_elastic.append(db[uid].start['tweak_motor_position'])

motor_pos_emission = []
fwhm_emission = []
max_emission = []
plt.figure(2, clear=True)
offset = 0
for uid in range(389461, 389485 + 1):
# for uid in range(389456, 389485 + 1):
    _fwhm, _max = analyze_linewidth_fly_scan(db, uid, rois=[1], plot_func=elasic_plot_func)
    fwhm_emission.append(_fwhm)
    max_emission.append(_max)
    motor_pos_emission.append(db[uid].start['tweak_motor_position'])



_, ax1 = plt.subplots(1, num=3, clear=True)
ax2 = ax1.twinx()
# ax3 = ax1.twinx()
# ax4 = ax1.twinx()

def my_plot_func(ax=ax1, marker='.', color='k', sign=1):
    def _bla(x, y, x_min, y_fit):
        _y = np.array(y) * sign
        _y_fit = np.array(y_fit) * sign
        ax.plot(x, _y, color + marker)
        ax.plot(x, _y_fit, color + '-')
        ax.vlines([x_min], np.min(_y), np.max(_y), colors=color)
    return _bla

get_optimal_crystal_alignment_position(motor_pos_elastic, fwhm_elastic, plot_func=my_plot_func(ax=ax1, marker='o', color='b', sign=1))
get_optimal_crystal_alignment_position(motor_pos_emission, fwhm_emission, plot_func=my_plot_func(ax=ax2, marker='o', color='r', sign=1))
# get_optimal_crystal_alignment_position(motor_pos_elastic, -np.array(max_elastic), plot_func=my_plot_func(ax=ax3, marker='s', color='b', sign=-1))
# get_optimal_crystal_alignment_position(motor_pos_emission, -np.array(max_emission), plot_func=my_plot_func(ax=ax4, marker='s', color='r', sign=-1))



total_img = []

elastic_scan = {}
for uid in range(389939, 389963):
    df = get_processed_df_from_uid(uid, db, save_interpolated_file=False, return_processed_df=True)

    hdr, primary_df, extended_data, comments, path_to_file, file_list, data_kind
    hdr = db[uid]
    x = hdr.start['tweak_motor_position']
    t = hdr.table(stream_name='pil100k_stream', fill=True)
    scan = []
    for img in t.pil100k_image[1]:
        scan.append(np.sum(img[85:120, 194:265]))
        if uid in [389939, 389950, 389962]:
            total_img.append(img)
    elastic_scan[x] = scan

bla = np.mean(np.array(total_img), axis=0)
plt.figure(2, clear=True)
plt.imshow(bla, vmin=0.04, vmax=0.07)


# data['Ecen'] = []
data['fwhm'] = []
# data['I_cor'] = []
# data['I_fit'] = []
# data['I_fit_raw'] = []
# data['E_scan'] = []
for uid in range(389939, 389963):
    Ecen0, fwhm0 = estimate_center_and_width_of_peak(E, I)
    Ecen, fwhm, I_cor, I_fit, I_fit_raw = fit_gaussian(E, I, Ecen0, fwhm0)
    # Ecen, fwhm, I_cor, I_fit, I_fit_raw, E_scan = analyze_elastic_fly_scan(db, uid=uid)
    # data['Ecen'].append(Ecen)
    data['fwhm'].append(fwhm)
    # data['I_cor'].append(I_cor)
    # data['I_fit'].append(I_fit)
    # data['I_fit_raw'].append(I_fit_raw)
    # data['E_scan'].append(E_scan)







#Elastic scan plot

data = {}
# data['x'] = []
for uid in range(389939, 389963):
    hdr = db[uid]
    x_val = hdr.start['tweak_motor_position']
    fwhm = analyze_elastic_fly_scan(db, uid=uid)
    data[x_val] = fwhm

fig, ax = plt.subplots(1,1, figsize=(4,6))

ax.plot(data.keys(), data.values(), '-ro')
ax.set_xlabel(r'R (mm)', size=12)
ax.set_ylabel(r'FWHM (eV)', size=12)
# ax.set_yticks(np.arange(1.52, 1.9, 0.05))

path = '/home/xf08id/Documents/Spectrometer_paper/'

plt.savefig(path + 'elastic.png', bbox_inches='tight', dpi=600)


#Emission cuk

#main = 389556, 389485
#aux2 = 389487, 389539
#aux3 = 389540, 389593

data_m = {}
for uid in range(389456, 389485):
    hdr = db[uid]
    x_val = hdr.start['tweak_motor_position']
    fwhm, _ = analyze_linewidth_fly_scan(db, uid=uid)
    data_m[x_val] = fwhm

data_2 = {}
for uid in range(389487, 389539):
    hdr = db[uid]
    if 'tweak_motor_position' in hdr.start.keys():
        _x = hdr.start['spectrometer_config']['bragg_registration']['pos_act']['motor_cr_assy_x'][0]
        x_val = hdr.start['tweak_motor_position']
        x_val = (x_val/1000) + _x
        x_key = 'johann_aux2_crystal_motor_cr_aux2_roll'
        fwhm, _ = analyze_linewidth_fly_scan(db, uid=uid, x_key=x_key)
        data_2[x_val] = fwhm
    else:
        pass

data_3 = {}
for uid in range(389540, 389593):
    hdr = db[uid]
    if 'tweak_motor_position' in hdr.start.keys():

        _x = hdr.start['spectrometer_config']['bragg_registration']['pos_act']['motor_cr_assy_x'][0]
        x_val = hdr.start['tweak_motor_position']
        x_val = (x_val/1000) + _x
        x_key = 'johann_aux3_crystal_motor_cr_aux3_roll'
        fwhm, _ = analyze_linewidth_fly_scan(db, uid=uid, x_key=x_key)
        data_3[x_val] = fwhm
    else:
        pass

fig, ax = plt.subplots(1,1, figsize=(6,4))

ax.plot(data_m.keys(), data_m.values(), '-ro', label='main')
ax.plot(data_2.keys(), data_2.values(), '-bo', label='aux2')
ax.plot(data_3.keys(), data_3.values(), '-go', label='aux3')

ax.set_xlabel(r'R (mm)', size=12)
ax.set_ylabel(r'FWHM (a.u.)', size=12)
# ax.set_yticks(np.arange(1.52, 1.9, 0.05))
ax.legend()

path = '/home/xf08id/Documents/Spectrometer_paper/'

plt.savefig(path + 'emission_cuk.png', bbox_inches='tight', dpi=600)

#herfd main 389700, 389724
#herfd aux2 389725, 389778
#herfd aux3 389779, 389832
#herfd aux4 389834, 389878
#herfd aux5 389880, 389938

X = []
for uid in range(389700, 389725):
    hdr = db[uid]
    x_val = hdr.start['tweak_motor_position']
    X.append(x_val)
    df = get_processed_df_from_uid(uid, db=db)
    np.savetxt(path + f'Main x {x_val:.3f}.dat', np.column_stack((df[1]['energy'], df[1]['i0'], df[1]['pil100k_roi1'])))

np.savetxt(path + f"Main x value.dat", np.column_stack((X)))

X = []
for uid in range(389725, 389779):
    hdr = db[uid]
    if 'tweak_motor_position' in hdr.start.keys():
        # _x = hdr.start['spectrometer_config']['bragg_registration']['pos_act']['motor_cr_assy_x'][0]
        _x = 986.700040041
        x_val = hdr.start['tweak_motor_position']
        x_val = (x_val / 1000) + _x
        X.append(x_val)
        df = get_processed_df_from_uid(uid, db=db)
        np.savetxt(path + f'Aux2 x {x_val:.3f}.dat', np.column_stack((df[1]['energy'], df[1]['i0'], df[1]['pil100k_roi1'])))
    else:
        pass
np.savetxt(path + f"Aux2 x value.dat", np.column_stack((X)))


X = []
for uid in range(389779, 389832):
    hdr = db[uid]
    if 'tweak_motor_position' in hdr.start.keys():
        # _x = hdr.start['spectrometer_config']['bragg_registration']['pos_act']['motor_cr_assy_x'][0]
        _x = 986.700040041
        x_val = hdr.start['tweak_motor_position']
        x_val = (x_val / 1000) + _x
        X.append(x_val)
        df = get_processed_df_from_uid(uid, db=db)
        np.savetxt(path + f'Aux3 x {x_val:.3f}.dat', np.column_stack((df[1]['energy'], df[1]['i0'], df[1]['pil100k_roi1'])))
    else:
        pass
np.savetxt(path + f"Aux3 x value.dat", np.column_stack((X)))

X = []
for uid in range(389834, 389878):
    hdr = db[uid]
    if 'tweak_motor_position' in hdr.start.keys():
        # _x = hdr.start['spectrometer_config']['bragg_registration']['pos_act']['motor_cr_assy_x'][0]
        _x = 986.700040041
        x_val = hdr.start['tweak_motor_position']
        x_val = (x_val / 1000) + _x
        X.append(x_val)
        df = get_processed_df_from_uid(uid, db=db)
        np.savetxt(path + f'Aux4 x {x_val:.3f}.dat', np.column_stack((df[1]['energy'], df[1]['i0'], df[1]['pil100k_roi1'])))
    else:
        pass
np.savetxt(path + f"Aux4 x value.dat", np.column_stack((X)))

X = []
for uid in range(389880, 389938):
    hdr = db[uid]
    if 'tweak_motor_position' in hdr.start.keys():
        # _x = hdr.start['spectrometer_config']['bragg_registration']['pos_act']['motor_cr_assy_x'][0]
        _x = 986.700040041
        x_val = hdr.start['tweak_motor_position']
        x_val = (x_val / 1000) + _x
        X.append(x_val)
        df = get_processed_df_from_uid(uid, db=db)
        np.savetxt(path + f'Aux5 x {x_val:.3f}.dat', np.column_stack((df[1]['energy'], df[1]['i0'], df[1]['pil100k_roi1'])))
    else:
        pass
np.savetxt(path + f"Aux5 x value.dat", np.column_stack((X)))

def get_timestamps_for_pbs(uid, db):
    hdr = db[uid]
    t_apb = hdr.table(stream_name='apb_stream', fill=True)
    t_enc = hdr.table(stream_name='pb9_enc1', fill=True)
    time_enc = t_enc['pb9_enc1'][1][:, 0] + t_enc['pb9_enc1'][1][:, 1]*1e-9
    time_apb = t_apb['apb_stream'][1][:, 0]
    return time_apb, time_enc


time_apb_sn, time_enc_sn =  get_timestamps_for_pbs('86011010-7b40-4987-abd0-106f8b31a86d', db)
time_apb_in, time_enc_in =  get_timestamps_for_pbs('a159e8ca-b4ed-4dde-b0ed-a2e6f4875067', db) # Sep 2023
# time_apb_in, time_enc_in =  get_timestamps_for_pbs('8d45c349', db) # June 2023

plt.figure(1, clear=True)
plt.plot(time_enc_sn - time_enc_sn[0], '.-')
plt.plot(time_enc_in - time_enc_in[0], '.-')


_, time_enc_in_now =  get_timestamps_for_pbs('a159e8ca-b4ed-4dde-b0ed-a2e6f4875067', db) # Sep 2023
_, time_enc_in_old = get_timestamps_for_pbs('8d45c349', db) # June 2023

plt.figure(1, clear=True)
plt.plot(time_enc_in_now - time_enc_in_now[0], '.-')
plt.plot(time_enc_in_old - time_enc_in_old[0], '.-')



# _, time_enc_cu_now =  get_timestamps_for_pbs('73761037-9d85-435b-a4a8-bb168baa73a0', db)
_, time_enc_cu_now =  get_timestamps_for_pbs('e6bbc310-aa06-4f87-8edc-bd1864088a38', db)
# _, time_enc_cu_now =  get_timestamps_for_pbs('e018c7bb-37ba-4a7d-8d9c-185ad579b7ab', db)
_, time_enc_cu_old = get_timestamps_for_pbs('3d265b3a', db)

plt.figure(1, clear=True)
plt.plot(time_enc_cu_now - time_enc_cu_now[0], '.-')
plt.plot(time_enc_cu_old - time_enc_cu_old[0], '.-')


R = 500

hkl1 = [1, 1, 1]
hkl2 = [1, 1, 0]
cos_alpha = np.dot(hkl1, hkl2)/np.linalg.norm(hkl1)/np.linalg.norm(hkl2)
alpha = np.arccos(cos_alpha)
alpha_deg = np.rad2deg(alpha)

bragg_deg = 90 - alpha_deg + 2.15
print(f'{bragg_deg=}, {alpha_deg=}')

bragg = np.deg2rad(bragg_deg)

# bragg = np.deg2rad(45)
# alpha = np.deg2rad(45) # 0*np.deg2rad(np.rad2deg(np.arccos(np.sqrt(2/3))))
crystal_x = R * np.sin(bragg + alpha) * np.sin(bragg - alpha)
crystal_y = R * np.sin(bragg + alpha) * np.cos(bragg - alpha)
det_x = 0
det_y = 2 * R * np.sin(bragg) * np.cos(bragg)
K_x = R * np.sin(bragg)**2 - R/2
K_y = R * np.sin(bragg)*np.cos(bragg)

phi = np.deg2rad(np.linspace(0, 360, 361))
KK_x = R/2 * np.cos(phi) + K_x
KK_y = R/2 * np.sin(phi) + K_y

plt.figure(1, clear=True)
plt.plot(0, 0, 'ko')

plt.plot(crystal_x, crystal_y, 'bo')
plt.plot(det_x, det_y, 'go')
plt.plot([0, crystal_x, det_x], [0, crystal_y, det_y], 'k-')

plt.plot(K_x, K_y, 'mo')
plt.plot(KK_x, KK_y, 'm-')

plt.xlim(-500, 1000)
plt.ylim(-500, 1000)
plt.axis('equal')


from xas.process import process_interpolate_bin_from_uid, get_processed_df_from_uid
from xas.file_io import save_extended_data_as_file

uids = [ 'd6b89a0d-1cf5-4b0b-8695-bc84cb5007e5',
         '318cec96-5213-4773-baf2-da50cae02f2d',
         '138e62df-c191-4831-9249-a5d252e8b6e8',
         '4989c59b-76a0-46dd-8723-51fd6f1ff6e5',
         'abaa4995-9bd1-4a5b-9b2d-4338ea30ab12',
         'b5d1704b-8c99-4f71-884f-14204166ab7b',
         '119523cd-0aae-4964-af3d-140458a51bcd',
         '4f94bb97-166e-4194-b145-bab39479bbf1',
         '23bfef24-7602-4a65-b761-cd9617152716',
         '8e71ee10-c9a3-4c49-9c59-11f169ca94d9' ]


# uids = list(range(42213, 42321 + 1))[4:]
uids = list(range(42304, 42321 + 1))
for i, uid in enumerate(uids):
    try:
        print('workding on', uid, db[uid].start.name, 'progress:', i / len(uids) * 100)
    except:
        pass
    if 'experiment' not in db[uid].start:
        continue
    hdr, primary_df, extended_data, comments, path_to_file, file_list, data_kind = get_processed_df_from_uid(uid, db,
                                                                              logger=None,
                                                                              draw_func_interp=None,
                                                                              draw_func_bin=None,
                                                                              print_func=None,
                                                                              save_interpolated_file=False,
                                                                              update_start=None,
                                                                              load_images=True)

    path_to_file = hdr.start['interp_filename'][:-3] + 'dat'
    # extended_data['pil100k_image'] = extended_data['pil100k_image'].astype(np.float32)
    save_extended_data_as_file(path_to_file, extended_data, data_kind=data_kind)



###
np.random.seed(1)

x = np.arange(100)
s = np.exp(-((x - 50)/ (np.sqrt(2) * 5))**2) * 7
noise1 = np.random.randn(100) * (np.sqrt(55)/np.sqrt(500))
noise2 = np.random.randn(100) * (np.sqrt(55/10)/np.sqrt(180))


plt.figure(1, clear=True)
plt.plot(x, s + noise1)
plt.plot(x, s + noise2)



####
from xas.file_io import load_binned_df_from_file


files = ['/nsls2/data/iss/legacy/processed/2024/1/309462/24-2254 reduced Pt-L3 RIXS 10s 0001 vh_roi1.dat',
         '/nsls2/data/iss/legacy/processed/2024/1/309462/24-2254 reduced Pt-L3 RIXS 10s 0002 vh_roi1.dat',
         '/nsls2/data/iss/legacy/processed/2024/1/309462/24-2254 reduced Pt-L3 RIXS 10s 0003 vh_roi1.dat',
         '/nsls2/data/iss/legacy/processed/2024/1/309462/24-2254 reduced Pt-L3 RIXS 10s 0004 vh_roi1.dat',
         '/nsls2/data/iss/legacy/processed/2024/1/309462/24-2254 reduced Pt-L3 RIXS 10s 0005 vh_roi1.dat',
         '/nsls2/data/iss/legacy/processed/2024/1/309462/24-2254 reduced Pt-L3 RIXS 10s 0006 vh_roi1.dat',
         '/nsls2/data/iss/legacy/processed/2024/1/309462/24-2254 reduced Pt-L3 RIXS 10s 0007 vh_roi1.dat',
         '/nsls2/data/iss/legacy/processed/2024/1/309462/24-2254 reduced Pt-L3 RIXS 10s 0008 vh_roi1.dat',
         '/nsls2/data/iss/legacy/processed/2024/1/309462/24-2254 reduced Pt-L3 RIXS 10s 0009 vh_roi1.dat',
         '/nsls2/data/iss/legacy/processed/2024/1/309462/24-2254 reduced Pt-L3 RIXS 10s 0010 vh_roi1.dat',
         '/nsls2/data/iss/legacy/processed/2024/1/309462/24-2254 reduced Pt-L3 RIXS 10s 0001-r0002 vh_roi1.dat',
         '/nsls2/data/iss/legacy/processed/2024/1/309462/24-2254 reduced Pt-L3 RIXS 10s 0002-r0002 vh_roi1.dat',
         '/nsls2/data/iss/legacy/processed/2024/1/309462/24-2254 reduced Pt-L3 RIXS 10s 0004-r0002 vh_roi1.dat',
         '/nsls2/data/iss/legacy/processed/2024/1/309462/24-2254 reduced Pt-L3 RIXS 10s 0003-r0002 vh_roi1.dat',
         '/nsls2/data/iss/legacy/processed/2024/1/309462/24-2254 reduced Pt-L3 RIXS 10s 0005-r0002 vh_roi1.dat',
         '/nsls2/data/iss/legacy/processed/2024/1/309462/24-2254 reduced Pt-L3 RIXS 10s 0006-r0002 vh_roi1.dat',
         '/nsls2/data/iss/legacy/processed/2024/1/309462/24-2254 reduced Pt-L3 RIXS 10s 0007-r0002 vh_roi1.dat',
         '/nsls2/data/iss/legacy/processed/2024/1/309462/24-2254 reduced Pt-L3 RIXS 10s 0008-r0002 vh_roi1.dat',
         '/nsls2/data/iss/legacy/processed/2024/1/309462/24-2254 reduced Pt-L3 RIXS 10s 0009-r0002 vh_roi1.dat',
         '/nsls2/data/iss/legacy/processed/2024/1/309462/24-2254 reduced Pt-L3 RIXS 10s 0010-r0002 vh_roi1.dat',
         '/nsls2/data/iss/legacy/processed/2024/1/309462/24-2254 reduced Pt-L3 RIXS 10s 0011 vh_roi1.dat',
         '/nsls2/data/iss/legacy/processed/2024/1/309462/24-2254 reduced Pt-L3 RIXS 10s 0012 vh_roi1.dat',
         '/nsls2/data/iss/legacy/processed/2024/1/309462/24-2254 reduced Pt-L3 RIXS 10s 0013 vh_roi1.dat',
         '/nsls2/data/iss/legacy/processed/2024/1/309462/24-2254 reduced Pt-L3 RIXS 10s 0014 vh_roi1.dat',
         '/nsls2/data/iss/legacy/processed/2024/1/309462/24-2254 reduced Pt-L3 RIXS 10s 0015 vh_roi1.dat',
         # '/nsls2/data/iss/legacy/processed/2024/1/309462/24-2254 reduced Pt-L3 RIXS 10s 0016 vh_roi1.dat',
         '/nsls2/data/iss/legacy/processed/2024/1/309462/24-2254 reduced Pt-L3 RIXS 10s 0017 vh_roi1.dat',
         '/nsls2/data/iss/legacy/processed/2024/1/309462/24-2254 reduced Pt-L3 RIXS 10s 0018 vh_roi1.dat',
         '/nsls2/data/iss/legacy/processed/2024/1/309462/24-2254 reduced Pt-L3 RIXS 10s 0019 vh_roi1.dat',
         '/nsls2/data/iss/legacy/processed/2024/1/309462/24-2254 reduced Pt-L3 RIXS 10s 0020 vh_roi1.dat',]

files = [
'/nsls2/data/iss/legacy/processed/2024/1/309462/23-24102 (pos 001) Pt-L3 RIXS 10s 0001 vh_roi1.dat',
'/nsls2/data/iss/legacy/processed/2024/1/309462/23-24102 (pos 001) Pt-L3 RIXS 10s 0001-r0002 vh_roi1.dat',
'/nsls2/data/iss/legacy/processed/2024/1/309462/23-24102 (pos 001) Pt-L3 RIXS 10s 0001-r0003 vh_roi1.dat',
'/nsls2/data/iss/legacy/processed/2024/1/309462/23-24102 (pos 001) Pt-L3 RIXS 10s 0001-r0004 vh_roi1.dat',
'/nsls2/data/iss/legacy/processed/2024/1/309462/23-24102 (pos 001) Pt-L3 RIXS 10s 0001-r0005 vh_roi1.dat',
'/nsls2/data/iss/legacy/processed/2024/1/309462/23-24102 (pos 001) Pt-L3 RIXS 10s 0001-r0006 vh_roi1.dat',
'/nsls2/data/iss/legacy/processed/2024/1/309462/23-24102 (pos 001) Pt-L3 RIXS 10s 0001-r0007 vh_roi1.dat',
'/nsls2/data/iss/legacy/processed/2024/1/309462/23-24102 (pos 001) Pt-L3 RIXS 10s 0001-r0008 vh_roi1.dat',
'/nsls2/data/iss/legacy/processed/2024/1/309462/23-24102 (pos 001) Pt-L3 RIXS 10s 0001-r0009 vh_roi1.dat',
'/nsls2/data/iss/legacy/processed/2024/1/309462/23-24102 (pos 001) Pt-L3 RIXS 10s 0001-r0010 vh_roi1.dat',
'/nsls2/data/iss/legacy/processed/2024/1/309462/23-24102 (pos 001) Pt-L3 RIXS 10s 0001-r0011 vh_roi1.dat',
'/nsls2/data/iss/legacy/processed/2024/1/309462/23-24102 (pos 001) Pt-L3 RIXS 10s 0001-r0012 vh_roi1.dat',
'/nsls2/data/iss/legacy/processed/2024/1/309462/23-24102 (pos 001) Pt-L3 RIXS 10s 0001-r0013 vh_roi1.dat',
'/nsls2/data/iss/legacy/processed/2024/1/309462/23-24102 (pos 001) Pt-L3 RIXS 10s 0001-r0017 vh_roi1.dat',
'/nsls2/data/iss/legacy/processed/2024/1/309462/23-24102 (pos 001) Pt-L3 RIXS 10s 0001-r0018 vh_roi1.dat',
'/nsls2/data/iss/legacy/processed/2024/1/309462/23-24102 (pos 001) Pt-L3 RIXS 10s 0001-r0019 vh_roi1.dat',
'/nsls2/data/iss/legacy/processed/2024/1/309462/23-24102 (pos 001) Pt-L3 RIXS 10s 0001-r0020 vh_roi1.dat',
'/nsls2/data/iss/legacy/processed/2024/1/309462/23-24102 (pos 001) Pt-L3 RIXS 10s 0001-r0021 vh_roi1.dat',
'/nsls2/data/iss/legacy/processed/2024/1/309462/23-24102 (pos 001) Pt-L3 RIXS 10s 0001-r0022 vh_roi1.dat',
'/nsls2/data/iss/legacy/processed/2024/1/309462/23-24102 (pos 001) Pt-L3 RIXS 10s 0001-r0023 vh_roi1.dat',
'/nsls2/data/iss/legacy/processed/2024/1/309462/23-24102 (pos 001) Pt-L3 RIXS 10s 0001-r0024 vh_roi1.dat',
'/nsls2/data/iss/legacy/processed/2024/1/309462/23-24102 (pos 001) Pt-L3 RIXS 10s 0001-r0025 vh_roi1.dat',
'/nsls2/data/iss/legacy/processed/2024/1/309462/23-24102 (pos 001) Pt-L3 RIXS 10s 0001-r0026 vh_roi1.dat',
'/nsls2/data/iss/legacy/processed/2024/1/309462/23-24102 (pos 001) Pt-L3 RIXS 10s 0001-r0027 vh_roi1.dat',
'/nsls2/data/iss/legacy/processed/2024/1/309462/23-24102 (pos 001) Pt-L3 RIXS 10s 0001-r0030 vh_roi1.dat',
'/nsls2/data/iss/legacy/processed/2024/1/309462/23-24102 (pos 001) Pt-L3 RIXS 10s 0001-r0031 vh_roi1.dat',
'/nsls2/data/iss/legacy/processed/2024/1/309462/23-24102 (pos 001) Pt-L3 RIXS 10s 0001-r0032 vh_roi1.dat',
'/nsls2/data/iss/legacy/processed/2024/1/309462/23-24102 (pos 001) Pt-L3 RIXS 10s 0001-r0033 vh_roi1.dat',
]


files = [
'/nsls2/data/iss/legacy/processed/2024/1/309462/23-63961 (pos 001) Pt-L3 RIXS 10s 0001 vh_roi1.dat',
'/nsls2/data/iss/legacy/processed/2024/1/309462/23-63961 (pos 001) Pt-L3 RIXS 10s 0001-r0002 vh_roi1.dat',
'/nsls2/data/iss/legacy/processed/2024/1/309462/23-63961 (pos 001) Pt-L3 RIXS 10s 0001-r0003 vh_roi1.dat',
'/nsls2/data/iss/legacy/processed/2024/1/309462/23-63961 (pos 001) Pt-L3 RIXS 10s 0001-r0004 vh_roi1.dat',
'/nsls2/data/iss/legacy/processed/2024/1/309462/23-63961 (pos 001) Pt-L3 RIXS 10s 0001-r0005 vh_roi1.dat',
'/nsls2/data/iss/legacy/processed/2024/1/309462/23-63961 (pos 001) Pt-L3 RIXS 10s 0001-r0006 vh_roi1.dat',
'/nsls2/data/iss/legacy/processed/2024/1/309462/23-63961 (pos 001) Pt-L3 RIXS 10s 0001-r0007 vh_roi1.dat',
'/nsls2/data/iss/legacy/processed/2024/1/309462/23-63961 (pos 001) Pt-L3 RIXS 10s 0001-r0008 vh_roi1.dat',
'/nsls2/data/iss/legacy/processed/2024/1/309462/23-63961 (pos 001) Pt-L3 RIXS 10s 0001-r0009 vh_roi1.dat',
'/nsls2/data/iss/legacy/processed/2024/1/309462/23-63961 (pos 001) Pt-L3 RIXS 10s 0001-r0010 vh_roi1.dat',
'/nsls2/data/iss/legacy/processed/2024/1/309462/23-63961 (pos 001) Pt-L3 RIXS 10s 0001-r0011 vh_roi1.dat',
'/nsls2/data/iss/legacy/processed/2024/1/309462/23-63961 (pos 001) Pt-L3 RIXS 10s 0001-r0012 vh_roi1.dat',
'/nsls2/data/iss/legacy/processed/2024/1/309462/23-63961 (pos 001) Pt-L3 RIXS 10s 0001-r0016 vh_roi1.dat',
'/nsls2/data/iss/legacy/processed/2024/1/309462/23-63961 (pos 001) Pt-L3 RIXS 10s 0001-r0017 vh_roi1.dat',
'/nsls2/data/iss/legacy/processed/2024/1/309462/23-63961 (pos 001) Pt-L3 RIXS 10s 0001-r0018 vh_roi1.dat',
'/nsls2/data/iss/legacy/processed/2024/1/309462/23-63961 (pos 001) Pt-L3 RIXS 10s 0001-r0019 vh_roi1.dat',
'/nsls2/data/iss/legacy/processed/2024/1/309462/23-63961 (pos 001) Pt-L3 RIXS 10s 0001-r0020 vh_roi1.dat',
'/nsls2/data/iss/legacy/processed/2024/1/309462/23-63961 (pos 001) Pt-L3 RIXS 10s 0001-r0021 vh_roi1.dat',
'/nsls2/data/iss/legacy/processed/2024/1/309462/23-63961 (pos 001) Pt-L3 RIXS 10s 0001-r0022 vh_roi1.dat',
'/nsls2/data/iss/legacy/processed/2024/1/309462/23-63961 (pos 001) Pt-L3 RIXS 10s 0001-r0023 vh_roi1.dat',
'/nsls2/data/iss/legacy/processed/2024/1/309462/23-63961 (pos 001) Pt-L3 RIXS 10s 0001-r0024 vh_roi1.dat',
'/nsls2/data/iss/legacy/processed/2024/1/309462/23-63961 (pos 001) Pt-L3 RIXS 10s 0001-r0025 vh_roi1.dat',
'/nsls2/data/iss/legacy/processed/2024/1/309462/23-63961 (pos 001) Pt-L3 RIXS 10s 0001-r0026 vh_roi1.dat',
# '/nsls2/data/iss/legacy/processed/2024/1/309462/23-63961 (pos 001) Pt-L3 RIXS 10s 0001-r0027 vh_roi1.dat',
'/nsls2/data/iss/legacy/processed/2024/1/309462/23-63961 (pos 001) Pt-L3 RIXS 10s 0001-r0030 vh_roi1.dat',
'/nsls2/data/iss/legacy/processed/2024/1/309462/23-63961 (pos 001) Pt-L3 RIXS 10s 0001-r0031 vh_roi1.dat',
'/nsls2/data/iss/legacy/processed/2024/1/309462/23-63961 (pos 001) Pt-L3 RIXS 10s 0001-r0032 vh_roi1.dat',
'/nsls2/data/iss/legacy/processed/2024/1/309462/23-63961 (pos 001) Pt-L3 RIXS 10s 0001-r0033 vh_roi1.dat',
]

rxes_list = []

for file in files:
    df, _ = load_binned_df_from_file(file)
    _energy_in = df.energy.values
    _energy_out = np.array([float('.'.join(i.split('_')[-2:])) for i in df.columns[1::2].tolist()])
    _counts = df.values[:, 1::2]
    _i0 = df.values[:, 2::2]
    rxes_list.append({'energy_in': _energy_in, 'energy_out': _energy_out, 'intensity': _counts / (_i0 / np.median(_i0, axis=1)[:, None])})

energy_in = _energy_in
energy_out = _energy_out
rxes = np.mean(np.array([_rxes['intensity'] for _rxes in rxes_list]), axis=0)

plt.figure(2, clear=True)
plt.contourf(energy_out, energy_in, rxes, 251, vmin=0, vmax=70)




# rxes_242254 = {'energy_out': energy_out, 'energy_in': energy_in, 'rxes': rxes}
# rxes_2324102 = {'energy_out': energy_out, 'energy_in': energy_in, 'rxes': rxes}
rxes_2363961 = {'energy_out': energy_out, 'energy_in': energy_in, 'rxes': rxes}


plt.figure(1, clear=True)
plt.contourf(rxes_242254['energy_out'], rxes_242254['energy_in'], rxes_242254['rxes'], 251, vmin=0, vmax=70)
plt.title('24-2254')

plt.figure(2, clear=True)
plt.contourf(rxes_2324102['energy_out'], rxes_2324102['energy_in'], rxes_2324102['rxes'], 251, vmin=0, vmax=70)
plt.title('23-24102')

plt.figure(3, clear=True)
plt.contourf(rxes_2363961['energy_out'], rxes_2363961['energy_in'], rxes_2363961['rxes'], 251, vmin=0, vmax=70)
plt.title('23-63961')


#########

mask = np.zeros(img.shape, dtype=bool)
mask[:, 59::61] = True
mask[:, 60::61] = True
mask[:, 61::61] = True
mask[96:99, :] = True

vsize, hsize = mask.shape
hpix, vpix = np.meshgrid(np.arange(hsize), np.arange(vsize))

hpix_arr = hpix.reshape(vsize * hsize)
vpix_arr = vpix.reshape(vsize * hsize)
mask_arr = mask.reshape(vsize * hsize)

dr = np.sqrt((hpix_arr[mask_arr][:, None] - hpix_arr[~mask_arr][None, :])**2 + (vpix_arr[mask_arr][:, None] - vpix_arr[~mask_arr][None, :])**2)
index_closest = np.argmin(dr, axis=1)

img_mask_arr = img.reshape(vsize * hsize)[mask_arr]
img_closeset_arr = img.reshape(vsize * hsize)[index_closest]

# mask_closest = np.zeros(img.shape, dtype=bool)
# for idx in index_closest:
#     vp = vpix_arr[~mask_arr][idx]
#     hp = hpix_arr[~mask_arr][idx]
#     # print(vp, hp)
#     mask_closest[vp, hp] = True


img_cor = img.copy()
img_cor[mask] *= 0.96
# img_cor[mask_closest] *= 10

plt.figure(1, clear=True)
plt.subplot(211)
plt.imshow(img, vmin=3, vmax=10)

plt.subplot(212)
plt.imshow(img_cor, vmin=3, vmax=10)

plt.figure(2, clear=True)
plt.plot(img_closeset_arr, img_closeset_arr / img_mask_arr, 'k.')
# plt.axis('square')