"""Compatibility constants for the legacy comparison scripts."""

from pathlib import Path

from py_vbc.config import load_config
from py_vbc.constants import *  # noqa: F403
from py_vbc import constants


PLANCK_CONFIG_FILE = Path(__file__).parents[1] / "planck2018_params.yaml"
TEST_CONFIG_FILE = Path(__file__).parents[1] / "test_params.yaml"
PLANCK_CONFIG = load_config(PLANCK_CONFIG_FILE)
TEST_CONFIG = load_config(TEST_CONFIG_FILE)

# Cosmology remains available under its YAML/runtime names for comparison code.
h = TEST_CONFIG.cosmology.h
omega_m = TEST_CONFIG.cosmology.omega_m
omega_b = TEST_CONFIG.cosmology.omega_b
sigma_8 = TEST_CONFIG.cosmology.sigma_8
ns = TEST_CONFIG.cosmology.ns
omega_r = TEST_CONFIG.cosmology.omega_r
costh = TEST_CONFIG.cosmology.costh


def test_definitions_and_constants_use_uppercase_names():
    assert constants.MPCTOKM == 3.08568e19
    assert constants.MPCTOCM == constants.MPCTOKM * 1e5
    assert (constants.DELTAC, constants.DELTAB, constants.DELTAT) == (0, 4, 8)
    assert (constants.REAL, constants.IMAG, constants.VELREAL, constants.VELIMAG) == (
        0,
        1,
        2,
        3,
    )
    assert (constants.DELC, constants.DELB, constants.VEL, constants.DELT) == (
        0,
        2,
        1,
        4,
    )
    assert constants.TCMB == 2.726
    assert constants.MPROTON == 1.6726e-24
    assert constants.MELECTRON == 9.11e-28
    assert constants.MEL_EV == 5.11e5
    assert constants.LIGHTSPEED == 3.0e10
    assert constants.CRITDENSITY == 1.8791e-29
    assert constants.SIGMAT == 0.665e-24
    assert constants.YHE == 0.25
    assert constants.BOLTZK == 1.3806e-16
    assert constants.MUB == 1.22
    assert constants.ASB == 7.56e-15
    assert constants.U_CMB == constants.ASB * constants.TCMB**4
    assert constants.T_GAMMA == (
        3
        * constants.LIGHTSPEED
        * constants.MELECTRON
        / (8 * constants.SIGMAT * constants.U_CMB)
    )


def test_legacy_lowercase_constant_names_are_not_exported():
    legacy_names = {
        "mpctokm",
        "mpctocm",
        "deltac",
        "deltab",
        "deltat",
        "real",
        "imag",
        "velreal",
        "velimag",
        "delc",
        "delb",
        "vel",
        "delt",
        "Tcmb",
        "mproton",
        "melectron",
        "mel_ev",
        "lightspeed",
        "critdensity",
        "sigmat",
        "yhe",
        "boltzk",
        "mub",
        "aSB",
        "u_cmb",
        "t_gamma",
    }
    assert legacy_names.isdisjoint(vars(constants))
