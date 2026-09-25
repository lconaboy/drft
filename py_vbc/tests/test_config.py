from pathlib import Path

from py_vbc.config import DEFAULT_RF_BASE, load_config


PACKAGE_ROOT = Path(__file__).parents[1]


def write_config(path: Path) -> Path:
    path.write_text(
        "cosmology:\n"
        "  h: 0.7\n"
        "  omega_m: 0.3\n"
        "  omega_b: 0.04\n"
        "  sigma_8: 0.8\n"
        "  ns: 0.97\n"
        "  costh: 1.0\n"
        "filenames:\n"
        "  transfer_functions:\n"
        "    0: tfs/runtime_transfer_z000.dat\n"
        "    1000: elsewhere/custom-name.dat\n"
        "  rf_base: recfast/runtime_recfast.dat\n",
        encoding="utf-8",
    )
    return path


def test_load_bundled_yaml():
    config = load_config(PACKAGE_ROOT / "planck2018_params.yaml")

    assert config.cosmology.h == 0.673
    assert config.cosmology.omega_m == 0.314
    assert config.cosmology.omega_b == 0.049
    assert config.cosmology.sigma_8 == 0.812
    assert config.cosmology.ns == 0.965
    assert config.cosmology.omega_r == 4.15e-5 / (config.cosmology.h**2.0)
    assert config.cosmology.costh == 1.0
    assert config.transfer_function_path(0) == (
        PACKAGE_ROOT / "tfs" / "planck2018_transfer_out_z000.dat"
    )
    assert config.rf_base == DEFAULT_RF_BASE


def test_optional_rf_base_defaults_and_notifies_once(monkeypatch, capsys):
    import py_vbc.config as config_module

    monkeypatch.setattr(config_module, "_DEFAULT_RF_NOTICE_SHOWN", False)
    config_path = PACKAGE_ROOT / "planck2018_params.yaml"

    first = load_config(config_path)
    second = load_config(config_path)
    message = capsys.readouterr().err

    assert first.rf_base == DEFAULT_RF_BASE
    assert second.rf_base == DEFAULT_RF_BASE
    assert message.count("rf_base was not specified") == 1


def test_yaml_relative_paths_and_explicit_tf_names(tmp_path):
    config_path = write_config(tmp_path / "params.yaml")
    config = load_config(config_path)

    assert config.transfer_function_path(0) == (
        tmp_path / "tfs" / "runtime_transfer_z000.dat"
    )
    assert config.transfer_function_path(1000) == (
        tmp_path / "elsewhere" / "custom-name.dat"
    )
    assert config.rf_base == tmp_path / "recfast" / "runtime_recfast.dat"
