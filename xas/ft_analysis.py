import numpy as np
from matplotlib import pyplot as plt
from datetime import datetime


def data_ft(data, plotting=True):
    t = data[5:-5, 0]
    x = data[5:-5, 1:]
    x /= np.mean(x, axis=0)
    f = np.fft.fftfreq(t.size, t[1] - t[0])
    # f_norm = (f>0) & (f<250)

    x_fft = np.fft.fft(x, axis=0)
    x_fft_abs = np.abs(x_fft)
    # x_fft_abs /= np.max(x_fft_abs[f_norm, :], axis=0)
    f_sel = f >= 0
    if plotting:
        fpath = plot_ft(f[f_sel], x_fft_abs[f_sel, :])
        return fpath


def plot_ft(f, x_fft, channels=[0], f_max=250):
    if f_max is None:
        f_max = f.max()

    f_norm = (f > 5) & (f < 250)

    plt.ioff()
    fig = plt.figure()
    fig.set_tight_layout(True)
    plt.plot(f, x_fft[:, channels], lw=0.5)
    plt.xlim(0, f_max)
    plt.ylim(0, x_fft[f_norm, channels].max())
    plt.xlabel("frequency, Hz")
    plt.ylabel("amplitude")
    fpath = (
        "/nsls2/xf08id/log/diagnostics/"
        + "vibration_diagnostics_"
        + str(datetime.now()).replace(":", "-")[:-7]
        + ".png"
    )
    plt.savefig(fpath, dpi=300)
    plt.ion()
    plt.close(fig)
    return fpath
