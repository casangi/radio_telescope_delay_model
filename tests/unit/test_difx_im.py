"""Full-precision check against a production DiFX correlator model.

``e21d15-1-b1_1036.{calc,im}`` are the delivered correlator input and model
files of one EHT 2021 job (3C273, MPIfR, DiFX 2.7.1; ALMA EC delivery
2019.1.01218.V). Running the embedded difxcalc11 on the .calc inputs with
the production conventions -- the delivered EOPs, an atmosphere-free phased
ALMA, and the packaged station catalogs -- must reproduce the .im
per-antenna u/v/w at the micrometre level for every station (measured
0.4-6.4 um rms). This guards the catalog contents (ocean loading and pole
tide coverage and code resolution for the EHT stations) as much as the
model: a missing or mis-resolved catalog entry shows up as millimetres.
"""

import os

import numpy as np
import pytest
from astropy.time import Time

from radio_telescope_delay_model.delay_model_cpp import calc_available
from radio_telescope_delay_model.difx_files import evaluate_im, read_calc, read_im

DATA = os.path.join(os.path.dirname(__file__), "data")


def test_read_calc_and_im():
    calc = read_calc(os.path.join(DATA, "e21d15-1-b1_1036.calc"))
    im = read_im(os.path.join(DATA, "e21d15-1-b1_1036.im"))
    assert calc["station_name"] == ["AA", "AX", "GL", "KT", "MG", "MM", "SW"]
    assert calc["source_name"] == "3C273"
    assert len(im) == 4 and im[0]["sec"] == 17760
    w = evaluate_im(im, "W", 0, [im[0]["mjd"]], [im[0]["sec"] + 60.0])[0]
    delay = evaluate_im(im, "DELAY", 0, [im[0]["mjd"]], [im[0]["sec"] + 60.0])[0]
    np.testing.assert_allclose(w, -299792458.0 * delay * 1e-6, rtol=1e-12)


@pytest.mark.skipif(not calc_available(), reason="CALC11 extension not built")
def test_model_reproduces_production_im_at_micrometres():
    from radio_telescope_delay_model.difxcalc11_core import (
        difxcalc11_samples,
        evaluate_samples,
    )

    calc = read_calc(os.path.join(DATA, "e21d15-1-b1_1036.calc"))
    im = read_im(os.path.join(DATA, "e21d15-1-b1_1036.im"))
    names = calc["station_name"]
    mjd = np.array([iv["mjd"] for iv in im for _ in range(3)])
    sec = np.array([iv["sec"] + t for iv in im for t in (12.0, 60.0, 108.0)])
    epochs = Time(
        mjd, sec / 86400.0, format="mjd", scale="utc"
    )  # two-part: exact seconds
    kw = dict(station_name=names, mount_type=calc["mount_type"], eop=calc["eop"])
    pos, pc = calc["antenna_position"], calc["phase_center_ra_dec"]
    on = difxcalc11_samples(pos, epochs, pc, **kw)
    off = difxcalc11_samples(pos, epochs, pc, add_atmosphere=False, **kw)
    for k, name in enumerate(names):
        samples = off if name == "AA" else on  # phased ALMA: no troposphere
        for q in "uvw":
            model = evaluate_samples(samples, q, epochs.unix)[:, k]
            truth = evaluate_im(im, q.upper(), k, mjd, sec)
            rms = np.sqrt(np.mean((model - truth) ** 2))
            assert rms < 10e-6, f"{name} {q}: {rms * 1e6:.1f} um"
