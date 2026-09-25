"""Runtime configuration for :mod:`py_vbc`."""

import math
import sys
import threading
from dataclasses import dataclass
from pathlib import Path
from collections.abc import Iterable, Mapping
from types import MappingProxyType
from typing import Any

import yaml


class ConfigError(ValueError):
    """Raised when a py_vbc runtime configuration is invalid."""


DEFAULT_RF_BASE = Path(__file__).resolve().parent / "recfast" / "planck2018_recfast.dat"
_DEFAULT_RF_NOTICE_LOCK = threading.Lock()
_DEFAULT_RF_NOTICE_SHOWN = False


def _use_default_rf_base() -> Path:
    """Select the bundled RECFAST file and notify once per process."""
    global _DEFAULT_RF_NOTICE_SHOWN
    with _DEFAULT_RF_NOTICE_LOCK:
        if not _DEFAULT_RF_NOTICE_SHOWN:
            print(
                "rf_base was not specified; using the bundled py_vbc RECFAST file: "
                f"{DEFAULT_RF_BASE}",
                file=sys.stderr,
                flush=True,
            )
            _DEFAULT_RF_NOTICE_SHOWN = True
    return DEFAULT_RF_BASE


@dataclass(frozen=True)
class Cosmology:
    """Cosmological parameters used by py_vbc."""

    h: float
    omega_m: float
    omega_b: float
    sigma_8: float
    ns: float
    costh: float

    def __post_init__(self) -> None:
        for name in (
            "h",
            "omega_m",
            "omega_b",
            "sigma_8",
            "ns",
            "costh",
        ):
            value = getattr(self, name)
            if isinstance(value, bool) or not isinstance(value, (int, float)):
                raise ConfigError(f"cosmology.{name} must be a number")
            value = float(value)
            if not math.isfinite(value):
                raise ConfigError(f"cosmology.{name} must be finite")
            object.__setattr__(self, name, value)
        self.validate()

    def validate(self) -> None:
        if self.h <= 0:
            raise ConfigError("cosmology.h must be greater than zero")
        if self.omega_m <= 0:
            raise ConfigError("cosmology.omega_m must be greater than zero")
        if not 0 <= self.omega_b <= self.omega_m:
            raise ConfigError(
                "cosmology.omega_b must be between zero and cosmology.omega_m"
            )
        if self.sigma_8 <= 0:
            raise ConfigError("cosmology.sigma_8 must be greater than zero")
        if not -1 <= self.costh <= 1:
            raise ConfigError("cosmology.costh must be between -1 and 1")

    @property
    def omega_r(self) -> float:
        """Radiation density parameter derived from the little-h parameter."""
        return 4.15e-5 / (self.h**2.0)  # Dodelson (2002) Eq. 2.86


@dataclass(frozen=True)
class RuntimeConfig:
    """A validated, immutable py_vbc runtime configuration."""

    cosmology: Cosmology
    transfer_functions: Mapping[int, Path]
    source: Path
    rf_base: Path | None = None

    def __post_init__(self) -> None:
        if not isinstance(self.cosmology, Cosmology):
            raise ConfigError("cosmology must be a Cosmology instance")

        if not isinstance(self.transfer_functions, Mapping):
            raise ConfigError("transfer_functions must be a mapping")
        if not self.transfer_functions:
            raise ConfigError("transfer_functions must not be empty")

        normalized = {}
        for redshift, value in self.transfer_functions.items():
            if isinstance(redshift, bool) or not isinstance(redshift, int):
                raise ConfigError("transfer-function redshift keys must be integers")
            if redshift < 0:
                raise ConfigError("transfer-function redshifts must not be negative")
            normalized[redshift] = self._normalize_runtime_path(
                value, f"transfer_functions.{redshift}"
            )

        object.__setattr__(
            self,
            "transfer_functions",
            MappingProxyType(normalized),
        )
        if self.rf_base is None:
            rf_base = _use_default_rf_base()
        else:
            rf_base = self._normalize_runtime_path(self.rf_base, "rf_base")
        object.__setattr__(self, "rf_base", rf_base)
        object.__setattr__(
            self,
            "source",
            self._normalize_runtime_path(self.source, "source"),
        )

    @staticmethod
    def _normalize_runtime_path(value: Any, key: str) -> Path:
        if isinstance(value, Path):
            path = value
        elif isinstance(value, str) and value.strip():
            path = Path(value)
        else:
            raise ConfigError(f"{key} must be a path string or Path")
        return path.expanduser().resolve()

    def transfer_function_path(self, redshift: int) -> Path:
        """Return the configured transfer-function path for *redshift*."""
        redshift = int(redshift)
        try:
            return self.transfer_functions[redshift]
        except KeyError as exc:
            raise ConfigError(
                f"no transfer function is configured for redshift z={redshift}"
            ) from exc

    def validate_files(self, redshifts: Iterable[int]) -> None:
        """Check that the RECFAST and required transfer-function files exist."""
        missing = []
        unconfigured = []
        for redshift in sorted({int(z) for z in redshifts}):
            if redshift < 0:
                raise ConfigError("transfer-function redshifts must not be negative")
            if redshift not in self.transfer_functions:
                unconfigured.append(redshift)
                continue
            path = self.transfer_functions[redshift]
            if not path.is_file():
                missing.append(path)

        if not self.rf_base.is_file():
            missing.append(self.rf_base)

        messages = []
        if unconfigured:
            redshifts_text = ", ".join(str(z) for z in unconfigured)
            messages.append(
                f"transfer functions are not configured for z={redshifts_text}"
            )
        if missing:
            paths = "\n".join(f"  - {path}" for path in missing)
            messages.append(f"configured data files do not exist:\n{paths}")
        if messages:
            raise ConfigError("\n".join(messages))


def _required_section(data: dict[str, Any], name: str) -> dict[str, Any]:
    try:
        section = data[name]
    except KeyError as exc:
        raise ConfigError(f"missing required YAML section: {name}") from exc

    if not isinstance(section, dict):
        raise ConfigError(f"YAML section '{name}' must be a mapping")
    return section


def _number(section: dict[str, Any], section_name: str, key: str) -> float:
    try:
        value = section[key]
    except KeyError as exc:
        raise ConfigError(
            f"missing required YAML parameter: {section_name}.{key}"
        ) from exc

    if isinstance(value, bool) or not isinstance(value, (int, float)):
        raise ConfigError(f"YAML parameter {section_name}.{key} must be a number")
    return float(value)


def _path(value: Any, key: str, source: Path) -> Path:
    if not isinstance(value, str) or not value.strip():
        raise ConfigError(f"YAML parameter filenames.{key} must be a path string")

    path = Path(value).expanduser()
    if not path.is_absolute():
        path = source.parent / path
    return path.resolve()


def _transfer_function_paths(
    filenames: dict[str, Any], source: Path
) -> Mapping[int, Path]:
    try:
        values = filenames["transfer_functions"]
    except KeyError as exc:
        raise ConfigError(
            "missing required YAML parameter: filenames.transfer_functions"
        ) from exc

    if not isinstance(values, dict) or not values:
        raise ConfigError(
            "YAML parameter filenames.transfer_functions must be a non-empty mapping"
        )

    paths = {}
    for redshift, value in values.items():
        if isinstance(redshift, bool) or not isinstance(redshift, int):
            raise ConfigError("transfer-function redshift keys must be integers")
        if redshift < 0:
            raise ConfigError("transfer-function redshifts must not be negative")
        paths[redshift] = _path(value, f"transfer_functions.{redshift}", source)

    return MappingProxyType(paths)


def load_config(config_file: str | Path) -> RuntimeConfig:
    """Load and validate a py_vbc YAML configuration.

    Relative ``filenames`` paths are resolved relative to the YAML file so a
    configuration and its data files can be moved together. If ``rf_base`` is
    omitted, the bundled Planck RECFAST file is selected and reported once.
    """
    source = Path(config_file).expanduser().resolve()
    if not source.is_file():
        raise ConfigError(f"configuration file does not exist: {source}")

    try:
        with source.open("r", encoding="utf-8") as stream:
            data = yaml.safe_load(stream)
    except yaml.YAMLError as exc:
        raise ConfigError(f"invalid YAML in {source}: {exc}") from exc

    if not isinstance(data, dict):
        raise ConfigError("the YAML configuration root must be a mapping")

    cosmology_data = _required_section(data, "cosmology")
    filenames = _required_section(data, "filenames")

    rf_value = filenames.get("rf_base")

    cosmology = Cosmology(
        h=_number(cosmology_data, "cosmology", "h"),
        omega_m=_number(cosmology_data, "cosmology", "omega_m"),
        omega_b=_number(cosmology_data, "cosmology", "omega_b"),
        sigma_8=_number(cosmology_data, "cosmology", "sigma_8"),
        ns=_number(cosmology_data, "cosmology", "ns"),
        costh=_number(cosmology_data, "cosmology", "costh"),
    )

    return RuntimeConfig(
        cosmology=cosmology,
        transfer_functions=_transfer_function_paths(filenames, source),
        source=source,
        rf_base=(_path(rf_value, "rf_base", source) if rf_value is not None else None),
    )
