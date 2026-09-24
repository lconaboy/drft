import os
import numpy as np
import pytest

import py_vbc

VBC_VALUES = [0.0, 30.0, 60.0]
COMPONENTS = ["c", "b", "vc", "vb"]


@pytest.fixture(scope="module")
def reference_data():
    """Load precomputed reference data for standard cosmology."""
    data_dir = os.path.join(os.path.dirname(__file__), "data")
    npz_path = os.path.join(data_dir, "reference_vbc.npz")
    assert os.path.exists(npz_path), f"Reference data not found at {npz_path}"
    data = np.load(npz_path)
    return data


@pytest.fixture(scope="module")
def computed_results(reference_data):
    """
    Run py_vbc once for each v_bc value using the explicit reference wavenumber grid k.
    Returns a dictionary mapping vbc -> dict of component power spectra.
    """
    k_ref = reference_data["k"]
    zstart = 1000.0
    zend = 200.0
    dz = 3.0

    results = {}
    for vbc in VBC_VALUES:
        k_out, (p_c, p_b, p_vc, p_vb) = py_vbc.run_pyvbc(
            vbc=vbc,
            zstart=zstart,
            zend=zend,
            dz=dz,
            k=k_ref,
            delta=False,
            verbose=False,
        )
        results[vbc] = {
            "k": k_out,
            "c": p_c,
            "b": p_b,
            "vc": p_vc,
            "vb": p_vb,
        }
    return results


@pytest.mark.parametrize("vbc", VBC_VALUES)
@pytest.mark.parametrize("comp", COMPONENTS)
def test_component_power_spectrum(vbc, comp, computed_results, reference_data):
    """
    Verifies that the power spectrum for each component (c, b, vc, vb) and each
    v_bc value matches the reference output to high precision.
    """
    actual = computed_results[vbc][comp]
    v_key = int(vbc)
    expected = reference_data[f"p_{comp}_{v_key}"]

    np.testing.assert_allclose(
        actual,
        expected,
        rtol=1e-2,
        atol=1e-5,
        err_msg=(
            f"Power spectrum mismatch for component='{comp}' at v_bc={vbc} km/s. "
            f"Max relative difference: "
            f"{np.max(np.abs((actual - expected) / expected)):.2e}"
        ),
    )
