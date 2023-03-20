import numpy as np
import pandas as pd
import larch
from larch.xafs import pre_edge

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
        **larch_pre_edge_kwargs,
        ):
        """Wrapper around `larch.xafs.pre_edge`, by default returns mu after 
        normalization and flattening"""
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
                step=norm_larch_group.step,
                pre1=norm_larch_group.pre1,
                pre2=norm_larch_group.pre2,
                norm1=norm_larch_group.norm1,
                norm2=norm_larch_group.norm2,
                nnorm=norm_larch_group.nnorm,
                nvict=norm_larch_group.nvict,
            )
            return mu_out, norm_parameters
        else:
            return mu_out
