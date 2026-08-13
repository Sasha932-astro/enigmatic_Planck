#!/usr/bin/env python3
import pathlib
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

HERE = pathlib.Path(__file__).resolve().parent
OUT = HERE / 'scenarios'
OUT.mkdir(exist_ok=True)

H = 6.62607015e-27      # erg s
C = 2.99792458e10       # cm/s
K = 1.380649e-16        # erg K^-1
LO, HI = 400.0, 1600.0


def read_clv(teff):
    path = HERE / f'mps_atlas/set1/MH0.0/teff{teff}/logg4.4/mpsa_intensity_spectra.dat'
    data = np.loadtxt(path, skiprows=2)
    wav_nm = data[:, 0]
    inu = data[:, 1:]
    with open(path) as f:
        f.readline()
        mu = np.array(f.readline().split()[2:], dtype=float)
    i_mu1 = int(np.argmin(np.abs(mu - 1.0)))
    assert abs(mu[i_mu1]-1.0) < 1e-8
    return wav_nm, inu[:, i_mu1]


def brightness_temperature_nu(wav_nm, inu):
    nu = C / (wav_nm * 1e-7)
    pref = 2.0 * H * nu**3 / C**2
    return (H * nu / K) / np.log1p(pref / inu)


def main():
    w57, i57 = read_clv(5700)
    w58, i58 = read_clv(5800)
    assert np.allclose(w57, w58)
    wav = w57
    tb57 = brightness_temperature_nu(wav, i57)
    tb58 = brightness_temperature_nu(wav, i58)
    dtb = tb58 - tb57
    m = (wav >= LO) & (wav <= HI) & np.isfinite(tb57) & np.isfinite(tb58)

    np.savez(OUT / '_brightness_temperature_5700_5800.npz',
             wavelength_nm=wav, Tb5700=tb57, Tb5800=tb58, deltaTb=dtb)

    fig, ax = plt.subplots(figsize=(10,6))
    ax.plot(wav[m], tb5700:=tb57[m], lw=1.6, label='MPS-ATLAS 5700 K')
    ax.plot(wav[m], tb5800:=tb58[m], lw=1.6, label='MPS-ATLAS 5800 K')
    ax.axhline(5700, ls='--', lw=1.2, alpha=.7)
    ax.axhline(5800, ls='--', lw=1.2, alpha=.7)
    ax.set(xlabel='Wavelength (nm)', ylabel='Brightness temperature at mu=1 (K)',
           title='MPS-ATLAS disk-centre brightness temperature')
    ax.grid(alpha=.25); ax.legend(); fig.tight_layout()
    fig.savefig(OUT / '_brightness_temperature_5700_5800.png', dpi=180)
    plt.close(fig)

    fig, ax = plt.subplots(figsize=(10,6))
    ax.plot(wav[m], dtb[m], lw=1.6)
    ax.axhline(100.0, ls='--', lw=1.4, label='Delta Teff = 100 K')
    ax.set(xlabel='Wavelength (nm)', ylabel='Tb(5800) - Tb(5700) (K)',
           title='Difference in disk-centre brightness temperature')
    ax.grid(alpha=.25); ax.legend(); fig.tight_layout()
    fig.savefig(OUT / '_brightness_temperature_delta_5700_5800.png', dpi=180)
    plt.close(fig)

    # modest diagnostics over 400-1600 nm
    print(f'Brightness-temperature experiment, mu=1, {LO:.0f}-{HI:.0f} nm')
    print(f'Tb5700 range: {tb57[m].min():.1f} .. {tb57[m].max():.1f} K')
    print(f'Tb5800 range: {tb58[m].min():.1f} .. {tb58[m].max():.1f} K')
    print(f'Delta Tb range: {dtb[m].min():.1f} .. {dtb[m].max():.1f} K')
    print(f'Delta Tb mean:  {dtb[m].mean():.1f} K')
    for lam in (400, 500, 700, 1000, 1300, 1600):
        j = int(np.argmin(abs(wav-lam)))
        print(f'{wav[j]:7.1f} nm: Tb57={tb57[j]:7.1f} K  Tb58={tb58[j]:7.1f} K  dTb={dtb[j]:6.1f} K')

if __name__ == '__main__':
    main()
