# Temperature-conversion program using PyQt
import numpy as np
import matplotlib.pyplot as plt
import pkg_resources
import scipy.integrate
from scipy import interpolate
import math
from PyQt5 import QtCore

import time as ttime
#import pexpect
from pexpect import pxssh
from ftplib import FTP
from subprocess import call

from . import xray
import pandas as pd
from isstools.dialogs.BasicDialogs import question_message_box
from xas.bin import xas_energy_grid

from collections import defaultdict


# Function to return a default
# values for keys that is not
# present
def _default_value():
    return "Not Present"


# Defining the dict




def stitch_two_points(t_orig, e1, v1, e2, v2, frac=0.5):
    t = t_orig.copy()
    t -= t.min()

    T1 = t.max()* frac
    T2 = t.max() * (1-frac)

    n = t.size
    f1 = (1 + np.cos( (t - T1/2) / (T1*2) * 4* np.pi))/2
    f2 = (1 + np.cos((t - T1 - T2/2) / (T2 * 2) * 4 * np.pi)) / 2

    f1[t > T1] = 0
    f2[t < T1] = 0

    F = np.vstack((f1,f2)).T

    q = t[1] - t[0]

    J = np.ones(n) * q
    J[0] = 0.5 * q
    J[-1] = 0.5 * q
    O = np.ones((n, n))*q
    A =  np.tril(O)  - np.diag(np.diag(O)) * 0.5
    A[:, 0] = 0.5 * q
    A[0, 0] = 0
    B = np.vstack((J @ A @ F, J @ F))
    b = np.hstack((e2 - e1 - J @ np.ones(n) * v1, v2 - v1))
    c,_,_,_ = np.linalg.lstsq(B, b, rcond=-1)
    return (A @ (A @ F @ c  + v1) + e1)



class TrajectoryCreator():
    def __init__(self, servocycle=16000, pulses_per_deg=360000):
        self.energy_grid = []
        self.encoder_grid = []
        self.servocycle = servocycle
        self.pulses_per_deg= pulses_per_deg

        self.e0 =  None
        self.direction = 'forward'
        self.edge_start = -30
        self.edge_end = 50

    def define(self, edge_energy = 11564, offsets = ([-200,-30,50,1000]), velocities = ([200, 20, 200]), stitching = ([75, 75, 10, 10, 100, 100]),
              padding_lo = 1, padding_hi=1, trajectory_type = 'Step', sine_duration = 20,
               dsine_preedge_duration = 10,
               dsine_edge_duration=8,
              dsine_postedge_duration = 20,
               dsine_preedge_frac=0.15,
               dsine_postedge_frac=0.85,
               pad_time=0):

        # self.servocycle=servocycle
        self.e0 = edge_energy
        self.direction = 'forward'
        self.edge_start = offsets[1]
        self.edge_end = offsets[2]
        self.gen_t_step = 1 / self.servocycle * 100

        if ((trajectory_type == 'Double Sine/Constant Edge') or
             (trajectory_type == 'standard')):
            t_pre = dsine_preedge_duration
            t_xanes = dsine_edge_duration
            t_exafs = dsine_postedge_duration
            print(t_pre, t_xanes, t_exafs)
            t_all = t_pre + t_xanes + t_exafs
            e_preedge_lo = edge_energy + offsets[0]
            e_preedge_hi = edge_energy + offsets[1]
            e_postedge_lo = edge_energy + offsets[2]
            e_postedge_hi = edge_energy + offsets[3]
            self.gen_t_step *= t_all/30

            t = np.arange(0, t_all, self.gen_t_step)
            e = np.zeros(t.shape)
            t_sel = (t >= t_pre) & (t <= (t_pre + t_xanes))
            e_velocity = (e_postedge_lo - e_preedge_hi) / (t[t_sel][-1] - t[t_sel][0])
            e_offset = e_preedge_hi - e_velocity * t_pre
            # e[t_sel] = e_velocity * t[t_sel] + e_offset
            e[t_sel] = np.linspace(e_preedge_hi, e_postedge_lo, int(np.sum(t_sel)))

            # t_sel = (t <= t_pre)
            # e[t_sel] = stitch_two_points(t[t_sel], e_preedge_lo, 0, e_preedge_hi,
            #                              e_velocity, frac=dsine_preedge_frac)
            # t_sel = (t >= (t_pre + t_xanes))
            # e[t_sel] = stitch_two_points(t[t_sel], e_postedge_lo, e_velocity,
            #                              e_postedge_hi, 0, frac=dsine_postedge_frac)

            e_preedge_hi_act = e[t_sel][0]
            e_postedge_lo_act = e[t_sel][-1]
            t_sel = (t < t_pre)
            e[t_sel] = stitch_two_points(t[t_sel], e_preedge_lo, 0, e_preedge_hi_act - e_velocity * self.gen_t_step,
                                         e_velocity, frac=dsine_preedge_frac)
            t_sel = (t > (t_pre + t_xanes))
            e[t_sel] = stitch_two_points(t[t_sel], e_postedge_lo_act + e_velocity * self.gen_t_step, e_velocity,
                                         e_postedge_hi, 0, frac=dsine_postedge_frac)

            self.time = t
            self.energy = e





        if trajectory_type == 'Step':
            preedge_lo = edge_energy+offsets[0]
            preedge_hi = edge_energy+offsets[1]
            edge_lo = preedge_hi
            edge_hi = edge_energy+offsets[2]
            postedge_lo = edge_hi
            postedge_hi = edge_energy+offsets[3]
            velocity_preedge = velocities[0]
            velocity_edge = velocities[1]
            velocity_postedge = velocities[2]

            preedge_stitch_lo = preedge_lo + stitching[0] #
            preedge_stitch_hi = preedge_hi - stitching[1] #
            edge_stitch_lo = edge_lo + stitching[2] #
            edge_stitch_hi = edge_hi - stitching[3] #
            postedge_stitch_lo = postedge_lo + stitching[4] #
            postedge_stitch_hi = postedge_hi - stitching[5] #

            # create padding points to ensure derivative is zero
            t_padding_lo = 0
            e_padding_lo = (1 - ((padding_lo/2.5) * 0.0125)) * preedge_lo#preedge_lo-20

            #create preedge
            t_current = padding_lo
            t_preedge_lo = t_current + (-preedge_lo + preedge_stitch_lo) / velocity_preedge
            e_preedge_lo = preedge_stitch_lo
            t_preedge_hi = t_current + (-preedge_lo + preedge_stitch_hi) / velocity_preedge
            e_preedge_hi = preedge_stitch_hi
            # e_step = 1*1e-2
            # e_preedge = np.arange(e_preedge_lo, e_preedge_hi + e_step, e_step)
            e_preedge = np.linspace(e_preedge_lo, e_preedge_hi, np.round(e_preedge_hi - e_preedge_lo)*10)
            # t_step = (t_preedge_hi - t_preedge_lo) / (len(e_preedge) - 1)
            # t_preedge = np.arange(t_preedge_lo, t_preedge_hi + (t_step / 2), t_step)
            t_preedge = np.linspace(t_preedge_lo, t_preedge_hi, e_preedge.size)
            t_step = t_preedge[1] - t_preedge[0]
            # dfgbdsfgbsd
            # # stich padding an preedge
            t_stich1 = np.linspace(t_padding_lo, t_preedge_lo, np.round((t_preedge_lo - t_padding_lo)/t_step))
            e_stich1 = stich_two_points(t_stich1,
                                        e_padding_lo, 0,
                                        e_preedge[0], (e_preedge[1] - e_preedge[0]) / (t_preedge[1] - t_preedge[0]))

            # create edge
            t_current = t_current + (-preedge_lo + preedge_hi) / velocity_preedge
            t_edge_lo = t_current + (-edge_lo + edge_stitch_lo) / velocity_edge
            e_edge_lo = edge_stitch_lo
            t_edge_hi = t_current + (-edge_lo + edge_stitch_hi) / velocity_edge
            e_edge_hi = edge_stitch_hi
            # e_step = 1*1e-2
            # e_edge = np.arange(e_edge_lo, e_edge_hi + e_step, e_step)
            e_edge = np.linspace(e_edge_lo, e_edge_hi, np.round(e_edge_hi - e_edge_lo)*10)
            # t_edge = np.arange(t_edge_lo, t_edge_hi + (t_step / 2), t_step)
            # t_step = (t_edge_hi - t_edge_lo) / (len(e_edge) - 1)
            t_edge = np.linspace(t_edge_lo, t_edge_hi, e_edge.size)
            t_step = t_edge[1] - t_edge[0]
            # t_edge = np.arange(t_edge_lo, t_edge_hi + (t_step / 2), t_step)

            # # stich preedge and xanes
            t_stich2 = np.linspace(t_preedge_hi, t_edge_lo, np.round((t_edge_lo - t_preedge_hi) / t_step*10))
            e_stich2 = stich_two_points(t_stich2, e_preedge_hi, velocity_preedge, e_edge_lo, velocity_edge)
            e_stich2 = stich_two_points(t_stich2,
                                        e_preedge[-1], (e_preedge[1] - e_preedge[0]) / (t_preedge[1] - t_preedge[0]),
                                        e_edge[0], (e_edge[1] - e_edge[0]) / (t_edge[1] - t_edge[0]))

            t_current = t_current + (-edge_lo + edge_hi) / velocity_edge
            t_postedge_lo = t_current + (-postedge_lo + postedge_stitch_lo) / velocity_postedge
            e_postedge_lo = postedge_stitch_lo
            t_postedge_hi = t_current + (-postedge_lo + postedge_stitch_hi) / velocity_postedge
            e_postedge_hi = postedge_stitch_hi
            # e_step = 1*1e-2
            # e_postedge = np.arange(e_postedge_lo, e_postedge_hi + e_step, e_step)
            e_postedge = np.linspace(e_postedge_lo, e_postedge_hi, np.round(e_postedge_hi - e_postedge_lo)*10)
            # t_step = (t_postedge_hi - t_postedge_lo) / (len(e_postedge) - 1)
            # t_postedge = np.arange(t_postedge_lo, t_postedge_hi + (t_step / 2), t_step)
            t_postedge = np.linspace(t_postedge_lo, t_postedge_hi, e_postedge.size)
            t_step = t_postedge[1] - t_postedge[0]

            # # stich preedge and xanes
            t_stich3 = np.linspace(t_edge_hi, t_postedge_lo, np.round((t_postedge_lo - t_edge_hi) / t_step) )
            e_stich3 = stich_two_points(t_stich3,
                                        e_edge[-1], (e_edge[1] - e_edge[0]) / (t_edge[1] - t_edge[0]),
                                        e_postedge[0], (e_postedge[1] - e_postedge[0]) / (t_postedge[1] - t_postedge[0]))

            t_current = t_current + (-postedge_lo + postedge_hi) / velocity_postedge
            t_padding_hi = t_current + padding_hi
            e_padding_hi = (1 + ((padding_hi/2.5) * 0.0125)) * postedge_hi #postedge_hi+20

            t_stich4 = np.linspace(t_postedge_hi, t_padding_hi, np.round((t_padding_hi - t_postedge_hi) / t_step))
            e_stich4 = stich_two_points(t_stich4,
                                        e_postedge[-1],  (e_postedge[1] - e_postedge[0]) / (t_postedge[1] - t_postedge[0]),
                                        e_padding_hi, 0)

            # concatenate the arrays
            #self.time = np.array([t_padding_lo, t_preedge_lo, t_preedge_hi, \
            #                      t_edge_lo, t_edge_hi, t_postedge_lo, t_postedge_hi, t_padding_hi])


            # _time = np.concatenate(([t_padding_lo], t_preedge, t_edge, t_postedge, [t_padding_hi]))
            # _energy = np.concatenate(([e_padding_lo], e_preedge, e_edge, e_postedge, [e_padding_hi]))


            self.time = np.concatenate(([t_padding_lo], t_stich1[1:-1], t_preedge, t_stich2[1:-1], t_edge, t_stich3[1:-1], t_postedge, t_stich4[1:-1], [t_padding_hi]))
            #self.energy = np.array([e_padding_lo, e_preedge_lo, e_preedge_hi,\
            #                        e_edge_lo, e_edge_hi, e_postedge_lo, e_postedge_hi,e_padding_hi])
            self.energy = np.concatenate(([e_padding_lo], e_stich1[1:-1], e_preedge, e_stich2[1:-1], e_edge, e_stich3[1:-1], e_postedge, e_stich4[1:-1], [e_padding_hi]))

            # _time_ = np.arange(_time[0], self.time[-1], 1 / self.servocycle)






        elif trajectory_type.lower() == 'sine':
            total_time = float(sine_duration)
            preedge_lo = edge_energy+offsets[0]
            postedge_hi = edge_energy+offsets[3]
            x = np.linspace(-np.pi / 2, np.pi / 2, 100)
            energy = (np.sin(x) * (postedge_hi - preedge_lo) / 2) + (postedge_hi + preedge_lo) / 2
            time = np.linspace(0, total_time, 100)
            self.energy = energy
            self.time = time
            self.gen_t_step *= total_time / 30

        # elif trajectory_type == 'Double Sine/Constant Edge':
        #     edge_duration = (offsets[2] - offsets[1]) / vel_edge
        #     total_time = float(dsine_preedge_duration) + float(edge_duration) + float(dsine_postedge_duration)
        #     preedge_dur = float(dsine_preedge_duration) / total_time
        #     edge_dur = float(edge_duration) / total_time
        #     postedge_dur = float(dsine_postedge_duration) / total_time
        #     preedge_lo = edge_energy + offsets[0]
        #     preedge_hi = edge_energy + offsets[1]
        #     postedge_lo = edge_energy + offsets[2]
        #     postedge_hi = edge_energy + offsets[3]
        #     edge = edge_energy
        #
        #     x1_num = int(np.round(preedge_dur * 1000))
        #     xedge_num = int(np.round(edge_dur * 1000))
        #     x2_num = int(np.round(postedge_dur * 1000))
        #     x1 = np.linspace(-np.pi / 2, (3 * np.pi / 2), x1_num)
        #     x2 = np.linspace(-np.pi / 2, (3 * np.pi / 2), x2_num)
        #
        #
        #     time1 = np.linspace(0, (preedge_dur * total_time), 2 * len(x1))
        #     time_int = np.linspace((preedge_dur * total_time) + (time1[1] - time1[0]), (preedge_dur * total_time) + (edge_dur * total_time), 2 * xedge_num)
        #     time2 = np.linspace(time_int[-1] + (time_int[1] - time_int[0]), total_time, 2 * len(x2))
        #     self.time = np.concatenate((time1, time_int, time2))
        #
        #     m_factor = 1
        #     m_factor_der = 0.005
        #     pos1 = [0, 1000 + abs(preedge_hi - preedge_lo)]
        #     last_error = abs((pos1[-1] - pos1[0]) - (preedge_hi - preedge_lo))
        #     while abs((pos1[-1] - pos1[0]) - (preedge_hi - preedge_lo)) > 0.5:
        #         accel1_1 = m_factor * (preedge_hi - preedge_lo) * (np.sin(x1) + 1) / 1000
        #         factor = (max(accel1_1) - (vel_edge * (time_int[1] - time_int[0]))) / max(accel1_1)
        #         accel1_2 = -factor * m_factor * (preedge_hi - preedge_lo) * (np.sin(x1) + 1) / 1000
        #
        #         vel1_1 = scipy.integrate.cumtrapz(accel1_1 * (1 / (x1_num/2)), initial = 0)
        #         vel1_2 = scipy.integrate.cumtrapz(accel1_2 * (1 / (x1_num/2)), initial = 0) + vel1_1[-1]
        #
        #         vel1 = np.concatenate((vel1_1, vel1_2))
        #         pos1 = scipy.integrate.cumtrapz(vel1, initial = 0) + preedge_lo
        #         #print("1:", abs((pos1[-1] - pos1[0]) - (preedge_hi - preedge_lo)), m_factor)
        #         error = abs((pos1[-1] - pos1[0]) - (preedge_hi - preedge_lo))
        #         if(error > last_error):
        #             m_factor_der = -(m_factor_der / 2)
        #         m_factor += m_factor_der
        #         last_error = error
        #
        #     vel_int = np.array([np.diff(time_int)[0]/np.diff(time1)[0] * vel1_2[-1]] * xedge_num * 2)
        #     pos_int = scipy.integrate.cumtrapz(vel_int, initial = 0) + pos1[-1] + (pos1[-1] - pos1[-2])
        #
        #     m_factor = 1
        #     m_factor_der = 0.005
        #     pos2 = [0, 1000 + abs(postedge_hi - postedge_lo)]
        #     last_error = abs((pos2[-1] - pos2[0]) - (postedge_hi - postedge_lo))
        #     while abs((pos2[-1] - pos2[0]) - (postedge_hi - postedge_lo)) > 0.5:
        #         accel2_1 = m_factor * (postedge_hi - postedge_lo) * (np.sin(x2) + 1) / 1000
        #         acc_factor = m_factor * (postedge_hi - postedge_lo) * (np.sin(x2) + 1) / 1000
        #         accel2_2 = - acc_factor * (((vel_int[0] * len(acc_factor) / 2) + sum(acc_factor)) / sum(acc_factor))
        #
        #         time_adj = np.diff(time2)[0]/np.diff(time_int)[0]
        #         vel2_1 = scipy.integrate.cumtrapz(accel2_1 * (1 / (x2_num/2)), initial = 0) + time_adj * vel_int[-1]
        #         vel2_2 = scipy.integrate.cumtrapz(accel2_2 * (1 / (x2_num/2)), initial = 0) + vel2_1[-1]
        #         vel2 = np.concatenate((vel2_1, vel2_2))
        #         pos2 = scipy.integrate.cumtrapz(vel2, initial = 0) + pos_int[-1] + (pos_int[-1] - pos_int[-2])
        #         error = abs((pos2[-1] - pos2[0]) - (postedge_hi - postedge_lo))
        #         if(error > last_error):
        #             m_factor_der = -(m_factor_der / 2)
        #         m_factor += m_factor_der
        #         last_error = error
        #
        #     self.energy = np.concatenate((pos1, pos_int, pos2))

        elif ((trajectory_type == 'Double Sine') or (trajectory_type == 'double_sine')):
            total_time = float(dsine_preedge_duration) + float(dsine_postedge_duration)
            half = float(dsine_preedge_duration) / total_time
            preedge_lo = edge_energy + offsets[0]
            postedge_hi = edge_energy + offsets[3]
            edge = edge_energy

            self.gen_t_step *= total_time / 30

            x1_num = int(np.round(half * 100))
            x2_num = int(np.round((1 - half) * 100))
            x1 = np.linspace(-np.pi / 2, (3 * np.pi / 2), x1_num) #2 * half / x_step1)
            x2 = np.linspace(-np.pi / 2, (3 * np.pi / 2), x2_num) #2 * (1 - half) / x_step1)

            accel1 = (edge_energy - preedge_lo) * (np.sin(x1) + 1) #39.584072231342851
            accel2 = (postedge_hi - edge_energy) * (np.sin(x2) + 1)
            accel1 = np.concatenate((accel1, -accel1))
            accel2 = np.concatenate((accel2, -accel2))

            vel1 = scipy.integrate.cumtrapz(accel1 * (1 / x1_num), initial = 0)
            vel2 = scipy.integrate.cumtrapz(accel2 * (1 / x2_num), initial = 0)

            pos1 = scipy.integrate.cumtrapz(vel1 * (1 / x1_num), initial = 0) + preedge_lo
            pos1 = (pos1 / pos1[len(pos1) - 1]) * edge
            pos2 = scipy.integrate.cumtrapz(vel2 * (1 / x2_num), initial = 0) + edge

            self.energy = np.concatenate((pos1, pos2))

            time = np.linspace(0, (half * total_time), 2 * len(x1))
            time2 = np.linspace((half * total_time) + (time[1] - time[0]), total_time, 2 * len(x2))
            self.time = np.concatenate((time, time2))

        if pad_time > 0:
            t_pad_pre = np.arange(-pad_time, 0, self.gen_t_step/100)
            t_pad_post = self.time[-1]-np.arange(-pad_time, 0, self.gen_t_step/100)[::-1]

            self.time = np.hstack((t_pad_pre, self.time, t_pad_post))
            self.time += pad_time
            self.energy = np.hstack((self.energy[0]*np.ones(t_pad_pre.shape),
                                     self.energy,
                                     self.energy[-1]*np.ones(t_pad_post.shape)))


    def define_from_dict(self, _scan_parameters):
        scan_parameters = defaultdict(_default_value)
        for k, v in _scan_parameters.items():
            scan_parameters[k] = v
        EXAFS_end =  xray.k2e(scan_parameters['EXAFS_end'], scan_parameters['e0']/1000)

        self.define(edge_energy=scan_parameters['e0'],
                    offsets=([scan_parameters['preedge_start'],
                              scan_parameters['XANES_start'],
                              scan_parameters['XANES_end'],
                              EXAFS_end]),
                    trajectory_type=scan_parameters['type'],
                    sine_duration=scan_parameters['duration'],
                    dsine_preedge_duration=scan_parameters['preedge_duration'],
                    dsine_edge_duration=scan_parameters['edge_duration'],
                    dsine_postedge_duration=scan_parameters['postedge_duration'],
                    dsine_preedge_frac=scan_parameters['preedge_flex'],
                    dsine_postedge_frac=scan_parameters['postedge_flex'],
                    pad_time=scan_parameters['pad'])


    def compute_time_per_bin(self):
        e0 = self.e0
        # edge_start = self.edge_start
        # edge_end = self.edge_end
        edge_start = np.max([-30, self.edge_start])
        edge_end = np.min([50, self.edge_end])
        preedge_spacing = 5
        if e0 < 14000:
            xanes_spacing = 0.2
        elif e0 >= 14000 and e0 < 21000:
            xanes_spacing = 0.3
        elif e0 >= 21000:
            xanes_spacing = 0.4
        else:
            xanes_spacing = 0.3
        exafs_k_spacing = 0.04

        e = None
        t = None

        if hasattr(self, 'energy_grid'):
            e = self.energy_grid
            t = self.time_grid


        if e is not None:
            try:
                e_bin = xas_energy_grid(e, e0, edge_start, edge_end, preedge_spacing, xanes_spacing, exafs_k_spacing)
                e_edges = np.hstack((e_bin[0], e_bin[:-1] + 0.5 * np.diff(e_bin), e_bin[-1]))
                npt, _ = np.histogram(e, e_edges)

                self.e_bin = e_bin
                self.time_per_bin = npt/self.servocycle
            except:
                self.e_bin = [0, 1]
                self.time_per_bin = [np.nan, np.nan]
        else:

            self.e_bin = [0, 1]
            self.time_per_bin = [np.nan, np.nan]




    def interpolate(self):
        cs = interpolate.CubicSpline(self.time, self.energy, bc_type='clamped')
        # cs = interpolate.CubicSpline(self.time, self.energy, bc_type='natural')
        self.time_grid = np.arange(self.time[0], self.time[-1], 1 / self.servocycle)
        self.energy_grid=cs(self.time_grid)
        # self.energy_grid_der=np.diff(self.energy_grid)/np.diff(self.time_grid)
        self.compute_time_per_bin()

    def revert_light(self):
        self.energy = self.energy[::-1]

    def revert(self):
        self.energy = self.energy[::-1]
        self.energy_grid = self.energy_grid[::-1]
        # self.energy_grid_der = self.energy_grid_der[::-1]
        self.direction = 'backward'

    def tile_light(self, reps = 1, single_direction = False):
        if not single_direction:
            self.time = np.append(self.time,
                                       (self.time + self.time[-1] + self.time[1] - self.time[0]))
            self.energy = np.append(self.energy, np.flipud(self.energy))
        # self.time_grid = np.tile(self.time_grid, reps)
        # self.time_grid = np.kron(np.arange(1, reps+1)/self.time_grid[-1], self.time_grid)
        self.time = np.hstack(
            [i * (self.time[-1] + self.time[1] - self.time[0]) + self.time for i in range(reps)])
        self.energy = np.tile(self.energy, reps)
        # self.compute_time_per_bin()
        if reps > 1:
            if self.direction == 'forward':
                self.direction += '-backward repeat'
            else:
                self.direction += '-forward repeat'

    def tile(self,reps = 1, single_direction = False):
        if not single_direction:
            self.time_grid = np.append(self.time_grid, (self.time_grid + self.time_grid[-1] + self.time_grid[1] - self.time_grid[0]))
            self.energy_grid = np.append(self.energy_grid, np.flipud(self.energy_grid))
        # self.time_grid = np.tile(self.time_grid, reps)
        # self.time_grid = np.kron(np.arange(1, reps+1)/self.time_grid[-1], self.time_grid)
        self.time_grid = np.hstack([i * (self.time_grid[-1] + self.time_grid[1] - self.time_grid[0]) + self.time_grid for i in range(reps)])
        self.energy_grid = np.tile(self.energy_grid, reps)
        self.compute_time_per_bin()
        if reps>1:
            if self.direction == 'forward':
                self.direction += '-backward repeat'
            else:
                self.direction += '-forward repeat'


    def define_complete(self, scan_parameters, lightweight=False):
        self.define_from_dict(scan_parameters)
        if lightweight:
            if scan_parameters['revert']: self.revert_light()
            self.tile_light(reps=scan_parameters['repeat'], single_direction=scan_parameters['single_direction'])
        else:
            self.interpolate()
            if scan_parameters['revert']: self.revert()
            self.tile(reps=scan_parameters['repeat'], single_direction=scan_parameters['single_direction'])
        self.element = scan_parameters['element']
        self.edge = scan_parameters['edge']
        self.E0 = scan_parameters['e0']



    def e2encoder(self, offset):
        self.encoder_grid = -xray.energy2encoder(self.energy_grid, self.pulses_per_deg, offset)

    def e2energy(self, offset):
        self.energy_grid = -xray.encoder2energy(self.encoder_grid, self.pulses_per_deg, offset)

    def plot(self):
        plt.plot(self.time, self.energy, 'r+')
        plt.plot(self.energy_grid,'b')
        plt.show()

    def load_trajectory_file(self, filename, offset, is_energy):
        array_out=[]
        with open(str(filename)) as f:
            first = f.readline()
            if first[0] != '#':
                array_out.append(float(first))
            else:
                pattern = 'E0: '
                idx =  first.find(pattern) + len(pattern)
                self.e0 = float(first[idx:])
            for line in f:
                array_out.append(float(line))
        array_out = np.array(array_out)
        if is_energy:
            # self.energy_grid_loaded = array_out
            self.energy_grid = array_out
        else:
            # self.energy_grid_loaded = -xray.encoder2energy(array_out, self.pulses_per_deg, offset)
            self.energy_grid = -xray.encoder2energy(array_out, self.pulses_per_deg, offset)
        self.time_grid = np.arange(self.energy_grid.size)/self.servocycle
        self.compute_time_per_bin()


    def save(self, filename):
        # if filename[-4:] == '.txt':
        #     filename = filename[:-4]
        # print(filename)
        # if len(filename):
        #     fileName, fileExtension = os.path.splitext(filename)
        #     if fileExtension is not '.txt':
        #         filename = fileName + '.txt'
        #     print(filename)
        #     if (os.path.isfile(filename)):
        #         ret = question_message_box('Save trajectory...', '{} already exists. Do you want to replace it?'.format(
        #             filename.rsplit('/', 1)[1]))
        #         if not ret:
        #             print('Aborted!')
        #             return
            np.savetxt(filename,
                       self.energy_grid, fmt='%.6f',
                       header = f'element: {self.element}, edge: {self.edge}, E0: {self.e0}, direction: {self.direction}')
            call(['chmod', '666', filename])
            # self.trajectory_path = filename[:filename.rfind('/')] + '/'
            # self.label_current_trajectory.setText(filename.rsplit('/', 1)[1])
            # self.push_plot_traj.setEnabled(True)
            print('Trajectory saved! [{}]'.format(filename))


# class trajectory_manager():
#     def __init__(self, hhm, **kwargs):
#         self.hhm = hhm
#         self.traj_info = {}
#
#     # Function used to count the number of lines in a file
#     def file_len(self, fname):
#         with open(fname) as f:
#             plusone = 0
#             if f.readline()[0] != '#':
#                 plusone = 1
#             for i, l in enumerate(f):
#                 pass
#         return i + 1 + plusone
#
#     def read_header(self, filename):
#         test = ''
#         line = '#'
#         with open(filename) as myfile:
#             while line[0] == '#':
#                 line = next(myfile)
#                 test += line
#         return test[:-len(line)-1]
#
#     ########## load ##########
#     # Transfer the trajectory file to the motor controller
#     # arg1 = orig_file_name             -> Filename of the trajectory: E.g.: 'traj_Cu_fast.txt'
#     # arg2 = new_file_path              -> LUT number where the new trajectory file will be stored
#     # arg3 (optional) = new_file_name     -> New name that will be used as filename in the controller. Currently, it MUST be 'hhm.txt'
#     # arg4 (optional) = orig_file_path     -> Path to look for the file that will be transfered. Default = '/GPFS/xf08id/trajectory/'
#     # arg5 (optional) = ip                 -> IP of the controller that will receive the file. Default = '10.8.2.86'
#     def load(self, orig_file_name, new_file_path, is_energy, offset, new_file_name = 'hhm.txt'):
#         ip = self.hhm.ip
#         orig_file_path = self.hhm.traj_filepath
#
#         print('[Load Trajectory] Starting...')
#         traj_fn = orig_file_name
#
#         # Check if new_file_path is between the possible values
#         if int(new_file_path) > 9 or int(new_file_path) < 1:
#             print("[Load Trajectory] Path '{}' not possible. Please use a value in the range 1 <= new_file_path <= 9.".format(new_file_path))
#             return False
#
#         # Get number of lines in file
#         file_size = self.file_len(orig_file_path + orig_file_name)
#         print('[Load Trajectory] Number of lines in file: {}'.format(file_size))
#
#         # Get min and max of trajectory in eV
#         if orig_file_path[-1] != '/':
#             fp += '/'
#
#         traj = pd.read_table('{}{}'.format(orig_file_path, orig_file_name), header=None, comment='#')
#         name = orig_file_name
#         header = self.read_header('{}{}'.format(orig_file_path, orig_file_name))
#         if is_energy:
#             min_energy = int(np.round(traj).min())
#             max_energy = int(np.round(traj).max())
#             enc = np.int64(np.round(xray.energy2encoder(-traj, self.hhm.pulses_per_deg, -offset)))
#             orig_file_name = '.energy_traj_aux.txt'
#             np.savetxt('{}{}'.format(orig_file_path, orig_file_name), enc, fmt='%d', header=header, comments='')
#         else:
#             min_energy = int(xray.encoder2energy((-traj, self.hhm.pulses_per_deg).min()))
#             max_energy = int(xray.encoder2energy((-traj, self.hhm.pulses_per_deg).max()))
#
#         print('[Load Trajectory] Min energy: {}'.format(min_energy))
#         print('[Load Trajectory] Max energy: {}'.format(max_energy))
#
#         # Create ftp connection with default credential
#         ftp = FTP(ip)
#         ftp.login()
#         s = pxssh.pxssh()
#         ssh_login = s.login (ip, 'root', 'deltatau')
#
#
#         if ssh_login:
#             # Check if the directory exists in /usrflash/lut/. If it does not, create it.
#             if str(new_file_path) != '':
#                 ftp.cwd('/usrflash/')
#                 dir_list = ftp.nlst()
#                 dir_exists = 0
#                 for dir_name in dir_list:
#                     if dir_name == 'lut':
#                         dir_exists = 1
#                 if not dir_exists:
#                     print('[Load Trajectory] mkdir: /usrflash/lut')
#                     ftp.mkd('/usrflash/lut')
#                     s.sendline ('chown ftp:root /var/ftp/usrflash/lut')
#                     s.sendline ('chmod a+wrx /var/ftp/usrflash/lut')
#
#                 ftp.cwd('/usrflash/lut/')
#                 dir_list = ftp.nlst()
#                 dir_exists = 0
#                 for dir_name in dir_list:
#                     if dir_name == str(new_file_path):
#                         dir_exists = 1
#                 if not dir_exists:
#                     print('[Load Trajectory] mkdir: /usrflash/lut/{}'.format(new_file_path))
#                     ftp.mkd('/usrflash/lut/{}'.format(new_file_path))
#                     s.sendline ('chown ftp:root /var/ftp/usrflash/lut/{}'.format(new_file_path))
#                     s.sendline ('chmod a+wrx /var/ftp/usrflash/lut/{}'.format(new_file_path))
#
#                 s.sendline ('chown ftp:root /var/ftp/usrflash/lut/{}/hhm.txt'.format(new_file_path))
#                 s.sendline ('chmod 777 /var/ftp/usrflash/lut/{}/hhm.txt'.format(new_file_path))
#
#             ftp_file_path = '/var/ftp/usrflash/lut/{}/{}'.format(new_file_path, new_file_name)
#
#         # Open file and transfer to the power pmac
#             f = open(orig_file_path + str(orig_file_name), 'rb')
#             if(f.readable()):
#                 line = f.readline().decode('utf-8')
#                 if line[0] == '#':
#                     element = line[line.find('element:') + 9: line.find(',')].lstrip()
#                     edge_value = line[line.find('edge:') + 6: line.find(',', line.find('edge:'))].lstrip()
#                     e0_value = line[line.find('E0:') + 4:].lstrip()
#                     curr_hhm_traj = getattr(self.hhm, 'traj{}'.format(new_file_path))
#                     curr_hhm_traj.filename.put(traj_fn)
#                     curr_hhm_traj.elem.put(element)
#                     curr_hhm_traj.edge.put(edge_value)
#                     curr_hhm_traj.e0.put(e0_value)
#                     curr_hhm_traj.bla = 'bla'
#                 else:
#                     curr_hhm_traj = getattr(self.hhm, 'traj{}'.format(new_file_path))
#                     curr_hhm_traj.filename.put(traj_fn)
#                     curr_hhm_traj.elem.put('')
#                     curr_hhm_traj.edge.put('')
#                     curr_hhm_traj.e0.put('')
#                     f.close()
#                     f = open(orig_file_path + str(orig_file_name), 'rb')
#                 result = ftp.storbinary('STOR ' + '/usrflash/lut/' + str(new_file_path) + '/' + new_file_name, f)
#                 if(result == '226 File receive OK.'):
#                     print('[Load Trajectory] File sent OK')
#                     s.sendline ('chown ftp:root /var/ftp/usrflash/lut/{}/{}'.format(new_file_path, new_file_name))
#                     s.sendline ('chmod a+wrx /var/ftp/usrflash/lut/{}/{}'.format(new_file_path, new_file_name))
#                     s.sendline ('echo "{}\n{}\n{}\n{}" > /var/ftp/usrflash/lut/{}/hhm-size.txt'.format(file_size, name, min_energy, max_energy, new_file_path))
#                     ttime.sleep(0.01)
#                     ftp.close()
#                     print('[Load Trajectory] Permissions OK')
#
#                 f.close()
#
#             s.logout()
#             s.pid = None
#             print('[Load Trajectory] Completed!')
#         else:
#             print('[Load Trajectory] Fail! Not able to ssh into the controller...')
#
#     ########## init ##########
#     # Transfer the trajectory from the flash to the ram memory in the controller
#     # It must be called everytime you decide to use a different trajectory
#     # arg1 = lut_number                -> lookup table number of the trajectory that will be used - must be a number between 1 and 9
#     # arg2 (optional) = ip            -> IP of the controller that will receive the file. Default = '10.8.2.86'
#     # arg3 (optional) = filename    -> Filename of the trajectory file in the controller. Currently, it MUST be 'hhm.txt'
#     def init(self, lut_number, filename = 'hhm.txt'):
#         ip = self.hhm.ip
#         print('[Init Trajectory] Starting...')
#
#         self.hhm.lut_number.put(lut_number)
#
#         ttime.sleep(0.5)
#         while (int(self.hhm.lut_number_rbv.get()) != int(lut_number)):
#             ttime.sleep(.001)
#             QtCore.QCoreApplication.processEvents()
#
#         self.hhm.lut_start_transfer.put("1")
#         while (self.hhm.lut_transfering.get() == 0):
#             ttime.sleep(.001)
#             QtCore.QCoreApplication.processEvents()
#         while (self.hhm.lut_transfering.get() == 1):
#             ttime.sleep(.001)
#             QtCore.QCoreApplication.processEvents()
#         ttime.sleep(.25)
#         #while (self.hhm.trajectory_loading.get() == 0):
#         #    ttime.sleep(.001)
#         #    QtCore.QCoreApplication.processEvents()
#         while (self.hhm.trajectory_loading.get() == 1):
#             ttime.sleep(.001)
#             QtCore.QCoreApplication.processEvents()
#
#         ftp = FTP(ip)
#         ftp.login()
#         ftp.cwd('/usrflash/lut/{}'.format(lut_number))
#
#         file_list = ftp.nlst()
#         file_exists = 0
#         for file_name in file_list:
#             if file_name == filename:
#                 file_exists = 1
#         if file_exists == 0:
#             print('[Init Trajectory] File not found. :(\nAre you sure \'{}\' is the correct lut number?'.format(lut_number))
#         else:
#             info = []
#             def handle_binary(more_data):
#                 info.append(more_data)
#
#             resp = ftp.retrlines('RETR hhm-size.txt', callback=handle_binary)
#             if(len(info) == 2):
#                 size = int(info[0])
#                 name = info[1]
#             elif(len(info) == 4):
#                 size = int(info[0])
#                 name = info[1]
#                 min_en = int(info[2])
#                 max_en = int(info[3])
#             else:
#                 print('[Init Trajectory] Could not find the size and name info in the controller. Please, try sending the trajectory file again using trajectory_load(...)')
#                 return False
#
#             if(size == 0):
#                 print('[Init Trajectory] Size seems to be equal to 0. Please, try sending the trajectory file again using trajectory_load(...)')
#                 return False
#             else:
#                 self.hhm.cycle_limit.put(size)
#                 while (self.hhm.cycle_limit_rbv.get() != size):
#                     ttime.sleep(.01)
#                 print('[Init Trajectory] New lut number: {}'.format(lut_number))
#                 print('[Init Trajectory] Trajectory name: {}'.format(name))
#                 print('[Init Trajectory] Number of points: {}'.format(size))
#                 print('[Init Trajectory] Completed!')
#                 self.hhm.trajectory_name.put(name)
#             ftp.close()
#
#
#
#     ########## read_info ##########
#     # Function that prints info about the trajectories currently stored in the controller
#     # arg1 (optional) = ip    -> IP of the controller. Default = '10.8.2.86's
#     def read_info(self, silent=False):
#         ip = self.hhm.ip
#         ftp = FTP(ip)
#         ftp.login()
#         ftp.cwd('/usrflash/lut/')
#         if not silent:
#             print('-'*62)
#             print('The trajectories found in the controller (ip: {}) are:'.format(ip))
#         self.traj_info.clear()
#
#         def handle_binary(more_data):
#             info.append(more_data)
#
#         ret_txt = ''
#         for i in range(1, 10):
#             ftp.cwd('/usrflash/lut/{}'.format(i))
#
#             file_list = ftp.nlst()
#             filename = 'hhm-size.txt'
#             file_exists = 0
#             for file_name in file_list:
#                 if file_name == filename:
#                     file_exists = 1
#             if file_exists == 1:
#
#                 info = []
#
#                 resp = ftp.retrlines('RETR hhm-size.txt', callback=handle_binary)
#                 if(len(info) == 2):
#                     size = int(info[0])
#                     name = info[1]
#                     self.traj_info[str(i)] = {'name':str(name), 'size':str(size)}
#                     if not silent:
#                         print('{}: {:<24} (Size: {})'.format(i, name, size))
#                 elif(len(info) == 4):
#                     size = int(info[0])
#                     name = info[1]
#                     min_en = int(info[2])
#                     max_en = int(info[3])
#                     self.traj_info[str(i)] = {'name':str(name), 'size':str(size), 'min':str(min_en), 'max':str(max_en)}
#                     if not silent:
#                         print('{}: {:<24} (Size: {}, min: {}, max: {})'.format(i, name, size, min_en, max_en))
#                 else:
#                     self.traj_info[str(i)] = {'name':'undefined', 'size':'undefined'}
#                     if not silent:
#                         print('{}: Could not find the size and name info'.format(i))
#             elif not silent:
#                 self.traj_info[str(i)] = {'name':'undefined', 'size':'undefined'}
#                 print('{}: Could not find the size and name info'.format(i))
#
#         if not silent:
#             print('-'*62)
#
#         return self.traj_info
#
#
#     def current_lut(self):
#         return self.hhm.lut_number_rbv.get()
#
# def read_trajectory_limits(hhm):
#     current_lut = int(hhm.lut_number_rbv.get())
#     traj_manager = trajectory_manager(hhm)
#     info = traj_manager.read_info(silent=True)
#     if 'max' not in info[str(current_lut)] or 'min' not in info[str(current_lut)]:
#         raise Exception(
#             'Could not find max or min information in the trajectory.'
#             ' Try sending it again to the controller.')
#
#     e_min = int(info[str(current_lut)]['min'])
#     e_max = int(info[str(current_lut)]['max'])
#     return e_min, e_max