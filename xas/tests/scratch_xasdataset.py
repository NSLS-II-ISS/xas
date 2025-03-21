from larch import Group as xafsgroup
from larch.xafs import pre_edge, autobk, mback, xftf
from larch import Interpreter
import numpy as np



class XASDataSetWithReference:
    _md = {}
    _filename = ''
    _larch = Interpreter()

    def __init__(self, name=None, md=None,
                 energy=None, mu=None,
                 energy_ref=None, mu_ref=None,
                 filename=None, datatype=None, process=True,
                 xasdataset=None, ext_data=None,
                 *args, **kwargs):
        self.larch = xafsgroup()
        self.larch_ref = xafsgroup()
        if md is not None:
            self._md = md
            if 'e0' in md:
                self.larch.e0 = float(md['e0'])
            elif 'edge' in md:
                edge = md['edge']
                self.larch.e0 = int(edge[edge.find('(') + 1: edge.find(')')])

        if mu is not None:
            self.larch.mu = np.array(mu)
            self._mu = np.array(mu)
        if energy is not None:
            self.larch.energy = np.array(energy)
            self.energy = self.larch.energy
        if mu_ref is not None:
            self.larch.mu = np.array(mu_ref)
            self._mu_ref = np.array(mu_ref)
        if energy_ref is not None:
            self.larch_ref.energy = np.array(energy)
            self.energy = self.larch.energy_ref

        if filename is not None:
            self._filename = filename
        if ext_data is not None:
            self.ext_data = ext_data

        if name is not None:
            self.name = name
        if datatype is not None:
            self.datatype = datatype
        if mu is not None and energy is not None:
            if process:
                self.nnorm = 2
                self.clamp_hi = 0
                self.clamp_lo = 0
                self.kweight = 2
                self.rbkg = 1
                self.kmax = 8
                self.kmin = 3
                self.kmin_ft = 3
                self.kmax_ft = 8
                try:
                    self.normalize()
                    self.deriv()
                    try:
                        self.extract_chi()
                        self.kmin_ft = 3
                        self.kmax_ft = self.kmax - 1

                    except:
                        pass

                except:
                    self.pre1 = None
                    self.pre2 = None
                    self.norm1 = None
                    self.norm2 = None
                    self.e0 = None
                    self.pre_edge = None
                    self.post_edge = None
                    self.edge_step = None

            if xasdataset is not None:
                self.clamp_hi = xasdataset.clamp_hi
                self.clamp_lo = xasdataset.clamp_lo
                self.kmin = xasdataset.kmin
                self.kmax = xasdataset.kmax
                self.kmin_ft = xasdataset.kmin_ft
                self.kmax_ft = xasdataset.kmax_ft
                self.kweight = xasdataset.kweight
                self.rbkg = xasdataset.rbkg
                self.pre1 = xasdataset.pre1
                self.pre2 = xasdataset.pre2
                self.norm1 = xasdataset.norm1
                self.norm2 = xasdataset.norm2
                self.nnorm = xasdataset.nnorm
                self.e0 = xasdataset.e0
                self.pre_edge = xasdataset.pre_edge
                self.post_edge = xasdataset.post_edge
                self.edge_step = xasdataset.edge_step
                # self.update_larch()

    def update_larch(self):
        if self.mu is not None:
            self.larch.mu = np.array(self.mu)
        if self.energy is not None:
            self.larch.energy = np.array(self.energy)

    def deriv(self):
        # mu_deriv=np.diff(np.transpose(self.mu.values))/np.diff(self.energy)
        mu_deriv = np.diff(self.mu) / np.diff(self.energy)
        self.mu_deriv = mu_deriv
        self.energy_deriv = (self.energy[1:] + self.energy[:-1]) / 2

    def flatten(self):
        step_index = int(np.argwhere(self.energy > self.e0)[0])
        zeros = np.zeros(step_index)
        ones = np.ones(self.energy.shape[0] - step_index)
        step = np.concatenate((zeros, ones), axis=0)
        diffline = (self.post_edge - self.pre_edge) / self.edge_step
        self.flat = self.norm + step * (1 - diffline)

    def normalize(self):
        pre_edge(self.larch, group=self.larch, _larch=self._larch)
        self.energy = self.larch.energy
        self.mu = self.larch.mu
        self.norm = self.larch.norm
        self.new_ds = False
        self.pre1 = self.larch.pre_edge_details.pre1
        self.pre2 = self.larch.pre_edge_details.pre2
        self.norm1 = self.larch.pre_edge_details.norm1
        self.norm2 = self.larch.pre_edge_details.norm2
        self.e0 = self.larch.e0
        self.pre_edge = self.larch.pre_edge
        self.post_edge = self.larch.post_edge
        self.edge_step = self.larch.edge_step
        self.flatten()

    def normalize_force(self):
        pre_edge(self.larch, group=self.larch, _larch=self._larch, e0=self.e0, pre1=self.pre1, pre2=self.pre2,
                 norm1=self.norm1, norm2=self.norm2, nnorm=self.nnorm)
        self.norm = self.larch.norm
        self.e0 = self.larch.e0
        self.pre_edge = self.larch.pre_edge
        self.post_edge = self.larch.post_edge
        self.edge_step = self.larch.edge_step
        self.flatten()

    def extract_chi(self):
        # print('chi ing')
        autobk(self.larch, group=self.larch, _larch=self._larch)

        self.chi = self.larch.chi
        self.bkg = self.larch.bkg
        self.kmin = self.larch.autobk_details.kmin
        self.kmax = self.larch.autobk_details.kmax
        # self.nclamp = 2
        # self.rbkg = 1

        # self.kmin_ft = self.kmin

    def extract_chi_force(self):
        # print('chi force reporting')
        autobk(self.larch, group=self.larch, _larch=self._larch,
               e0=self.e0, kmin=self.kmin, kmax=self.kmax, kweight=self.kweight,
               rbkg=self.rbkg, clamp_lo=self.clamp_lo, clamp_hi=self.clamp_hi)
        # autobk(self.larch, group=self.larch, _larch=self._larch, e0=self.e0, kmin=self.kmin, kmax=self.kmax,
        #        nclamp=2, clamp_hi=10)
        self.k = self.larch.k
        self.chi = self.larch.chi
        self.bkg = self.larch.bkg

    def extract_ft(self):
        # print('ft reporting')
        # print(self.kmin_ft)
        # xftf(self.larch, group=self.larch,  _larch=self._larch, kmin=self.kmin_ft, kmax=self.kmax)
        xftf(self.larch, group=self.larch, _larch=self._larch, kmin=self.kmin_ft, kmax=self.kmax_ft)

        self.r = self.larch.r
        self.chir = self.larch.chir
        self.chir_mag = self.larch.chir_mag
        self.chir_im = self.larch.chir_re
        self.chir_re = self.larch.chir_im
        # self.chir_pha = self.larch.chir_pha
        self.kmax_ft = self.kmax
        self.kwin = self.larch.kwin

    def extract_ft_force(self, window={}):
        # print('ft force reporting')
        if not window:
            xftf(self.larch, group=self.larch, _larch=self._larch,
                 kmin=self.kmin_ft, kmax=self.kmax_ft, kweight=self.kweight)
        else:
            window_type = window['window_type']
            tapering = window['tapering']
            r_weight = window['r_weight']
            # print('setting window')
            xftf(self.larch, group=self.larch, _larch=self._larch,
                 kmin=self.kmin_ft, kmax=self.kmax_ft, kweight=self.kweight,
                 window=window_type, dk=tapering, rweight=r_weight)
        self.r = self.larch.r
        self.chir = self.larch.chir
        self.chir_mag = self.larch.chir_mag
        self.chir_im = self.larch.chir_re
        self.chir_re = self.larch.chir_im
        # self.chir_pha = self.larch.chir_phas
        self.kwin = self.larch.kwin

    @property
    def md(self):
        return self._md

    @md.setter
    def md(self, md):
        self._md = md
        if 'e0' in md:
            # self.larch.e0 = int(md['e0'])
            self.larch.e0 = float(md['e0'])
            pass
        elif 'edge' in md:
            edge = md['edge']
            self.larch.e0 = int(edge[edge.find('(') + 1: edge.find(')')])

    @property
    def mu(self):
        return self._mu

    @mu.setter
    def mu(self, mu):
        if hasattr(mu, 'values'):
            values = mu.values
        else:
            values = mu
        # self._mu = pd.DataFrame(values, columns=['mu'])
        # self.larch.mu = self._mu
        self.larch.mu = values

    @property
    def filename(self):
        return self._filename

    @filename.setter
    def filename(self, filename):
        self._filename = filename
