import numpy as np
import pandas as pd
import larch
from larch.xafs import pre_edge, autobk

def calc_mus(df: pd.DataFrame):
    """Shorthand function to calculate mut, muf, and mur in dataframe using
    i0, it, ir, and iff"""
    df["mut"] = -np.log(df["it"] / df["i0"])
    df["mur"] = -np.log(df["ir"] / df["it"])
    df["muf"] = df["iff"] / df["i0"]


class LarchCalculator:
    """Container for wrappers around larch functionality"""
    def __init__(self) -> None:
        pass
    # def __init__(
    #     self, 
    #     input_group: larch.Group=None, 
    #     output_group: larch.Group=None,
    #     store_results=True,
    #     ) -> None:
    #     if store_results is True:
    #         if input_group is None:
    #             self.input_group = larch.Group()
    #         if output_group is None:
    #             self.output_group = larch.Group()

    @staticmethod
    def _intepret_params(params: dict, function_kwargs: list[str]):
        return {kwarg: (params[kwarg] if kwarg in params else None) for kwarg in function_kwargs}
        
    @staticmethod
    def custom_flatten(larch_group: larch.Group):
        step_index = int(np.argwhere(larch_group.energy > larch_group.e0)[0])
        zeros = np.zeros(step_index)
        ones = np.ones(larch_group.energy.shape[0] - step_index)
        step = np.concatenate((zeros, ones), axis=0)
        diffline = (larch_group.post_edge - larch_group.pre_edge) / larch_group.edge_step
        larch_group.flat = larch_group.norm + step * (1 - diffline)

    @staticmethod
    def normalize(
        energy,
        mu,
        flatten_output=True,
        return_norm_parameters=False,
        **params,
    ):
        """Wrapper around `larch.xafs.pre_edge`. By default returns mu after 
        normalization and flattening, in addition to the fitted pre-edge and
        post-edge curves."""

        larch_pre_edge_kwargs = LarchCalculator._intepret_params(params, [
            "e0", "step", "pre1", "pre2", "norm1", "norm2", "nnorm", "nvict"
        ])

        energy = np.array(energy)
        mu = np.array(mu)
        raw_larch_group = larch.Group(energy=energy, mu=mu)
        norm_larch_group = larch.Group(energy=energy)
        pre_edge(raw_larch_group, group=norm_larch_group, **larch_pre_edge_kwargs)
        LarchCalculator.custom_flatten(norm_larch_group)
        
        if flatten_output:
            mu_out = norm_larch_group.flat
        else:
            mu_out = norm_larch_group.norm
        
        if return_norm_parameters:
            norm_parameters = dict(
                e0=norm_larch_group.e0,
                step=norm_larch_group.edge_step,
                pre1=norm_larch_group.pre_edge_details.pre1,
                pre2=norm_larch_group.pre_edge_details.pre2,
                norm1=norm_larch_group.pre_edge_details.norm1,
                norm2=norm_larch_group.pre_edge_details.norm2,
                nnorm=norm_larch_group.pre_edge_details.nnorm,
                nvict=norm_larch_group.pre_edge_details.nvict,
            )
            return mu_out, norm_larch_group.pre_edge, norm_larch_group.post_edge, norm_parameters
        else:
            return mu_out, norm_larch_group.pre_edge, norm_larch_group.post_edge
        
    @staticmethod
    def auto_background(
        energy,
        mu,
        return_autobk_params,
        **params,
    ):
        """Wrapper around `larch.xafs.autobk` function"""

        # incomplete list of autobk kwargs, more can be added later
        larch_autobk_kwargs = LarchCalculator._intepret_params(params, [
            "kmin", "kmax", "clamp_lo", "clamp_hi", "rbkg", "kweight", "wing"
        ])

        energy = np.array(energy)
        mu = np.array(mu)
        larch_group_in = larch.Group(energy=energy, mu=mu)
        larch_group_out = larch.Group()
        autobk(larch_group_in, group=larch_group_out, **larch_autobk_kwargs)

        k_out = larch_group_out.k
        chi_out = larch_group_out.chi
        
        if return_autobk_params:
            autobk_params = dict(
                # there are other autobk params but these are all we care about for xdash gui
                kmin=larch_group_out.kmin,
                kmax=larch_group_out.kmax,
                clamp_lo=larch_group_out.clamp_lo,
                clamp_hi=larch_group_out.clamp_hi,
                rbkg=larch_group_out.rbkg,
                kweight=larch_group_out.kweight,
                win=larch_group_out.win,
            )
            return k_out, chi_out, autobk_params
        else:
            return k_out, chi_out