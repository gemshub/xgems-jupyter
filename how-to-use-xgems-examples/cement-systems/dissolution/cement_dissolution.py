"""
cement_dissolution.py
=====================

Extensible framework for kinetically controlled cement/mortar dissolution
simulations with xGEMS.

Architecture
------------
- ``ReactiveSurface``  : ABC mapping current amount -> reactive surface area (m^2).
                        Concrete: ``ConstantSurface``, ``SpecificSurface``,
                        ``ShrinkingSphereSurface``.
- ``KineticMaterial``  : ABC for any material whose amount evolves on a kinetic
                        (not equilibrium) law. Concrete: ``QuartzTST``.
                        Stubs:     ``MetalCorrosion``, ``OrganicDegradation``.
- ``CouplingHook``     : pre/post-equilibration hooks; useful for read-only
                        engine properties like pH. Stub: ``FixedPHHook``.
- ``CementDissolutionSimulation`` : orchestrates the time loop.

Conventions
-----------
- Time is in seconds throughout; helpers convert to/from years.
- Rates are in mol/s.
- Surface areas are in m^2; specific surface in m^2/mol.
- Failed equilibrations are recorded as NaN, never raised, so parameter sweeps
  complete instead of crashing midway.
- The baseline bulk composition (element vector ``b``) is captured ONCE, after
  the initial equilibration and BEFORE the time loop starts; coupling hooks
  receive a fresh copy on every step so their adjustments never compound.

Authors: refactor of original notebook by G. Dan Miron, Georg Kosakowski.
"""

from __future__ import annotations

from abc import ABC, abstractmethod
from dataclasses import dataclass, field
from typing import Callable, Mapping, Optional, Sequence

import math
import warnings

import numpy as np

try:
    import xgems as xg
except ImportError:  # pragma: no cover - allows inspecting the module without xgems
    xg = None


# ---------------------------------------------------------------------------
# Constants & helpers
# ---------------------------------------------------------------------------

YEAR_S: float = 365.25 * 24.0 * 3600.0  # Julian year, seconds
CELSIUS_TO_KELVIN: float = 273.15
R_GAS: float = 8.314462618  # J / (mol K)

# GEMS return codes we treat as "good"
_GEMS_OK_CODES = (2, 6)


def years_to_seconds(years: float) -> float:
    return years * YEAR_S


def seconds_to_years(seconds: float) -> float:
    return seconds / YEAR_S


def celsius_to_kelvin(t_celsius: float) -> float:
    return t_celsius + CELSIUS_TO_KELVIN


# ===========================================================================
# Reactive surface area models
# ===========================================================================

class ReactiveSurface(ABC):
    """Maps the current amount of a material to its reactive surface area (m^2).

    ``initial_amount_mol`` is provided so geometric models (shrinking-sphere,
    shrinking-core) have a reference point.
    """

    @abstractmethod
    def surface_area(self, current_amount_mol: float, initial_amount_mol: float) -> float:
        ...


@dataclass
class ConstantSurface(ReactiveSurface):
    """Fixed surface area, independent of remaining amount."""

    area_m2: float

    def surface_area(self, current_amount_mol, initial_amount_mol):
        return self.area_m2


@dataclass
class SpecificSurface(ReactiveSurface):
    """Surface area proportional to current amount: ``A = sa_per_mol * n``.

    Matches the original notebook's "rounded grains that don't fully dissolve"
    assumption (SA = qtz_amount * 9.065e-2 m^2/mol for ~1.5 mm spherical
    grains, Table 4.8 of NAB 18-05).
    """

    sa_per_mol_m2: float  # m^2/mol

    def surface_area(self, current_amount_mol, initial_amount_mol):
        return self.sa_per_mol_m2 * max(current_amount_mol, 0.0)


@dataclass
class ShrinkingSphereSurface(ReactiveSurface):
    """Surface area shrinks geometrically as ``(n/n0)^(2/3)``.

    Appropriate for monodisperse spherical grains that fully dissolve. Useful
    for metals and glasses where the original notebook's "don't fully dissolve"
    assumption doesn't apply.
    """

    initial_area_m2: float

    def surface_area(self, current_amount_mol, initial_amount_mol):
        if initial_amount_mol <= 0:
            return 0.0
        frac = max(current_amount_mol / initial_amount_mol, 0.0)
        return self.initial_area_m2 * (frac ** (2.0 / 3.0))


# ===========================================================================
# Engine state snapshot
# ===========================================================================

@dataclass
class EngineState:
    """Read-only snapshot of engine state, refreshed once per timestep.

    Holds the expensive-to-compute quantities so each material doesn't have to
    query the engine separately. Materials needing something not present here
    can still reach back to the raw engine via ``engine``.
    """

    time_s: float
    step: int
    temperature_K: float
    pressure_Pa: float
    pH: float
    aHplus: float
    aOHminus: float
    saturation_indices: Mapping[str, float] = field(default_factory=dict)  # phase name -> omega
    engine: object = None  # raw xgems.ChemicalEngine (escape hatch)


# ===========================================================================
# Kinetic materials
# ===========================================================================

class KineticMaterial(ABC):
    """A material whose amount evolves on a kinetic rate law.

    Subclasses must provide:
      - ``species_name``      : the xGEMS *species* (DC) name whose amount is
                                constrained each timestep.
      - ``compute_rate(...)`` : dissolution/reaction rate in mol/s.

    Optional:
      - ``saturation_phase``  : phase name whose omega the rate law needs.
                                Defaults to None.
    """

    species_name: str
    saturation_phase: Optional[str] = None

    @abstractmethod
    def compute_rate(self, state: EngineState, current_amount_mol: float) -> float:
        """Return the dissolution/reaction rate in mol/s.

        Positive value = material is being consumed (dissolution, corrosion,
        degradation). Return 0.0 to signal "no kinetic constraint this step"
        (e.g. supersaturated: let equilibrium handle it).
        """

    def validate(self, engine) -> None:
        """Raise if this material's species/phase is not in the loaded system."""
        try:
            engine.indexSpecies(self.species_name)
        except Exception as e:
            raise ValueError(
                f"{type(self).__name__}: species '{self.species_name}' not in GEMS system "
                f"(engine.indexSpecies raised: {e})"
            )
        if self.saturation_phase is not None:
            try:
                engine.indexPhase(self.saturation_phase)
            except Exception as e:
                raise ValueError(
                    f"{type(self).__name__}: phase '{self.saturation_phase}' not in GEMS "
                    f"system (engine.indexPhase raised: {e})"
                )


# ---------------------------------------------------------------------------
# Concrete: Quartz TST dissolution (faithful port of the original notebook)
# ---------------------------------------------------------------------------

@dataclass
class QuartzTST(KineticMaterial):
    """Transition-state-theory dissolution rate for quartz / silicates.

        r = ( k_n*f_n(T) + k_a*f_a(T)*a_H+^n_a + k_b*f_b(T)*a_OH-^n_b )
            * S * (1 - omega)

    where ``f_i(T) = exp(-Ea_i/R * (1/T - 1/T_ref))``.

    Default parameters are the original notebook's values (quartz, with the
    acid mechanism switched off). Override for other silicates via the
    dataclass fields.
    """

    species_name: str = "Qtz"
    saturation_phase: Optional[str] = "Quartz"
    surface: ReactiveSurface = field(default_factory=lambda: SpecificSurface(9.065e-2))

    # neutral mechanism
    k_neutral: float = 6.4e-14    # mol m^-2 s^-1 at T_ref
    Ea_neutral: float = 77_000.0  # J/mol

    # acid mechanism (k_acid = 0.0 => switched off)
    k_acid: float = 0.0
    Ea_acid: float = 0.0
    n_acid: float = 0.0

    # base mechanism
    k_base: float = 1.9e-10
    Ea_base: float = 80_000.0
    n_base: float = 0.34

    T_ref: float = 298.15

    # initial amount, captured on first compute_rate call; only used by
    # geometric surface models. Not part of the dataclass init.
    _initial_amount: Optional[float] = field(default=None, init=False, repr=False)

    def compute_rate(self, state, current_amount_mol):
        if self._initial_amount is None:
            self._initial_amount = current_amount_mol

        omega = state.saturation_indices.get(self.saturation_phase)
        if omega is None or not np.isfinite(omega):
            return 0.0

        # Dissolution-only kinetics: at omega >= 1 we let the equilibrium
        # solver decide (matches the original "only change lower limit if
        # dissolution happens" guard).
        if omega >= 1.0:
            return 0.0

        S = self.surface.surface_area(current_amount_mol, self._initial_amount)
        if S <= 0.0:
            return 0.0

        inv_T_diff = (1.0 / state.temperature_K) - (1.0 / self.T_ref)

        f_neutral = math.exp(-self.Ea_neutral / R_GAS * inv_T_diff)
        f_acid    = math.exp(-self.Ea_acid    / R_GAS * inv_T_diff)
        f_base    = math.exp(-self.Ea_base    / R_GAS * inv_T_diff)

        r_n = self.k_neutral * f_neutral
        r_a = self.k_acid    * f_acid * (state.aHplus  ** self.n_acid) if self.k_acid else 0.0
        r_b = self.k_base    * f_base * (state.aOHminus ** self.n_base) if self.k_base else 0.0

        return (r_n + r_a + r_b) * S * (1.0 - omega)


# ---------------------------------------------------------------------------
# STUB: metal corrosion
# ---------------------------------------------------------------------------

@dataclass
class MetalCorrosion(KineticMaterial):
    """STUB: simple constant-rate corrosion of a metallic species.

        r = k_corr * f_T(T) * S

    The Simulation will treat this exactly like QuartzTST: lower the species
    lower-limit by ``r * dt`` each step.

    Future extensions (write your own subclass or extend this one):
      - O2 dependence:  multiply by f(p_O2) or f([O2(aq)])
      - pH dependence:  piecewise k(pH) for passivation regimes
      - Butler-Volmer:  electrochemical control with Eh from the engine
      - Coupling to a growing corrosion-product layer (diffusion limitation)
        via a ``ShrinkingCoreSurface`` model
    """

    species_name: str = ""               # e.g. "Fe(metal)" or "Cu(metal)"
    saturation_phase: Optional[str] = None
    surface: ReactiveSurface = field(default_factory=lambda: ConstantSurface(1.0))

    k_corr: float = 0.0    # mol m^-2 s^-1 at T_ref
    Ea: float = 0.0        # J/mol
    T_ref: float = 298.15

    _initial_amount: Optional[float] = field(default=None, init=False, repr=False)

    def compute_rate(self, state, current_amount_mol):
        if self._initial_amount is None:
            self._initial_amount = current_amount_mol
        if current_amount_mol <= 0.0 or self.k_corr <= 0.0:
            return 0.0
        S = self.surface.surface_area(current_amount_mol, self._initial_amount)
        inv_T_diff = (1.0 / state.temperature_K) - (1.0 / self.T_ref)
        f_T = math.exp(-self.Ea / R_GAS * inv_T_diff)
        return self.k_corr * f_T * S


# ---------------------------------------------------------------------------
# STUB: organic degradation
# ---------------------------------------------------------------------------

@dataclass
class OrganicDegradation(KineticMaterial):
    """STUB: first-order degradation with optional Q10 temperature scaling.

        dn/dt = -k_decay * Q10^((T - T_ref)/10) * n

    The Simulation interprets the returned positive rate as "amount lost per
    second" and lowers the species' lower-limit accordingly.

    Future extensions:
      - Monod kinetics:    r = k_max * n / (K_half + n)
      - Electron-donor/acceptor coupling with redox state
      - Biomass coupling:  r depends on a separate biomass species
      - Multi-component organics (cellulose -> sugars -> CO2 chain)
    """

    species_name: str = ""
    saturation_phase: Optional[str] = None
    k_decay: float = 0.0   # 1/s
    Q10: float = 2.0
    T_ref: float = 298.15

    def compute_rate(self, state, current_amount_mol):
        if current_amount_mol <= 0.0 or self.k_decay <= 0.0:
            return 0.0
        T_factor = self.Q10 ** ((state.temperature_K - self.T_ref) / 10.0)
        return self.k_decay * T_factor * current_amount_mol


# ===========================================================================
# Coupling hooks (pre/post equilibration)
# ===========================================================================

class CouplingHook(ABC):
    """Hook fired around equilibration.

    Use when something can't be set directly on the engine (e.g. pH is
    read-only in the xGEMS API). Default implementations are no-ops, so
    subclasses only override what they need.

    Hooks receive a FRESH copy of the baseline bulk-composition vector every
    step. Modifications never compound across timesteps.
    """

    def before_equilibrate(self, engine, state: EngineState, bulk_b: np.ndarray) -> np.ndarray:
        """May return a modified bulk-composition vector. Default: no change."""
        return bulk_b

    def after_equilibrate(self, engine, state: EngineState) -> None:
        """Optional post-equilibration adjustment. Default: no-op."""
        return None


@dataclass
class FixedPHHook(CouplingHook):
    """STUB: target a fixed pH by iterative adjustment of a chosen handle element.

    pH is read-only in xGEMS, so the only way to control it is to perturb the
    bulk composition before equilibration. A robust implementation needs an
    inner Newton or bisection loop:

        equilibrate -> measure pH -> nudge b[handle] -> re-equilibrate
        until |pH - target| < tolerance.

    The current implementation is intentionally minimal: it shows the wiring
    so you can drop in a real algorithm. Use ``before_equilibrate`` to mutate
    the bulk vector that the Simulation will pass into ``engine.equilibrate``.

    Parameters
    ----------
    target_pH : float
        Desired pH at the end of the equilibration.
    tolerance : float
        Absolute pH tolerance (default 0.05).
    max_iter : int
        Max inner iterations per outer timestep (default 5).
    pH_handle : str
        Element symbol used to push pH up/down (e.g. "Cl" for acid, "Na" for base).
    """

    target_pH: float
    tolerance: float = 0.05
    max_iter: int = 5
    pH_handle: str = "Cl"

    def before_equilibrate(self, engine, state, bulk_b):
        # IMPLEMENT-ME: iterative adjustment of bulk_b[handle_idx] until pH
        # converges. See module docstring for the recommended pattern.
        return bulk_b


# ===========================================================================
# Configuration dataclasses
# ===========================================================================

@dataclass
class TimeConfig:
    """Time-stepping configuration. All values in seconds."""

    t_end_s: float = 5000.0 * YEAR_S
    dt_initial_s: float = 1.0 * YEAR_S
    dt_max_s: Optional[float] = None
    accel_after_s: float = 1000.0 * YEAR_S
    accel_factor: float = 1.01


@dataclass
class OutputConfig:
    """What to track during the simulation.

    Parameters
    ----------
    aqueous_elements : sequence of str
        Element symbols whose total aqueous concentration is recorded
        (mmol/L of aqueous phase).
    phase_groups : mapping[str, sequence[str]]
        Display-group -> list of GEMS phase names. The volume of each group
        is the sum over its phases (cm^3). Lets multiple GEMS phases collapse
        to a single display category (e.g. "ettringite" group covers
        SO4_CO3_AFt, CO3_SO4_AFt, ettringite).
    extras : mapping[str, callable]
        Extra scalars to record each step: name -> ``f(engine) -> float``.
    aqueous_phase : str
        GEMS phase name for the aqueous solution (default "aq_gen").
    """

    aqueous_elements: Sequence[str] = ("Ca", "Si", "Al", "Mg", "S", "C")
    phase_groups: Mapping[str, Sequence[str]] = field(default_factory=dict)
    extras: Mapping[str, Callable] = field(default_factory=dict)
    aqueous_phase: str = "aq_gen"


@dataclass
class CompositionOverride:
    """Programmatically override initial composition before the first equilibration.

    Two layered mechanisms:
      - ``species_lower_limits`` : set a per-species lower limit (mol). Use to
        force a minimum amount of an initial mineral (matches the original
        notebook's ``setSpeciesLowerLimit("Qtz", 26.679)`` pattern).
      - ``species_upper_limits`` : set a per-species upper limit (mol). Use to
        suppress an unwanted phase by capping it near zero.
      - ``bulk_elements``        : directly override the element-amount (b)
        vector before the first equilibration. Use to change overall mortar
        composition (more Ca, less Si, etc.) without editing the GEMS .lst file.

    Order of application: lower limits, upper limits, bulk element overrides.
    """

    species_lower_limits: Mapping[str, float] = field(default_factory=dict)
    species_upper_limits: Mapping[str, float] = field(default_factory=dict)
    bulk_elements: Mapping[str, float] = field(default_factory=dict)


# ===========================================================================
# Results container
# ===========================================================================

@dataclass
class SimulationResults:
    """Tidy container. All time series are aligned to ``time_years``.

    Failed-equilibration timesteps appear as NaN. Use ``to_dataframes()`` for
    pandas DataFrames (requires pandas)."""

    time_years: list = field(default_factory=list)
    pH: list = field(default_factory=list)
    gems_code: list = field(default_factory=list)
    aqueous_mmol_per_L: dict = field(default_factory=dict)   # element -> list[float]
    solids_cm3: dict = field(default_factory=dict)           # group -> list[float]
    rates_mol_per_s: dict = field(default_factory=dict)      # species -> list[float]
    extras: dict = field(default_factory=dict)               # name -> list[float]

    def to_dataframes(self):
        """Return (properties, aqueous, solids) pandas DataFrames."""
        import pandas as pd
        props = pd.DataFrame(
            {
                "time_years": self.time_years,
                "pH": self.pH,
                "gems_code": self.gems_code,
                **{f"rate_{sp}_mol_s": v for sp, v in self.rates_mol_per_s.items()},
                **self.extras,
            }
        )
        aq = pd.DataFrame({"time_years": self.time_years, **self.aqueous_mmol_per_L})
        solids = pd.DataFrame({"time_years": self.time_years, **self.solids_cm3})
        return props, aq, solids

    # ---- Persistence ------------------------------------------------------

    def save(self,
             output_dir: str = "results",
             run_name: Optional[str] = None,
             format: str = "csv",
             metadata: Optional[dict] = None,
             overwrite: bool = False) -> str:
        """Persist results to disk.

        Parameters
        ----------
        output_dir : str
            Directory to write into. Created if it does not exist.
        run_name : str, optional
            Subfolder/file prefix. Defaults to ``"run_YYYYMMDD_HHMMSS"``.
        format : {"csv", "xlsx"}
            ``"csv"`` writes three files (properties.csv, aqueous.csv,
            solids.csv) into ``output_dir/run_name/``. ``"xlsx"`` writes one
            workbook ``output_dir/run_name.xlsx`` with three sheets (requires
            pandas + openpyxl).
        metadata : dict, optional
            Saved as ``metadata.json`` next to the data (or as a separate sheet
            for xlsx). Use to record the simulation config for reproducibility.
        overwrite : bool
            If False and the target already exists, raises ``FileExistsError``.

        Returns
        -------
        str
            Absolute path of the directory (csv) or file (xlsx) written.
        """
        import os
        import json
        from datetime import datetime

        if format not in ("csv", "xlsx"):
            raise ValueError(f"format must be 'csv' or 'xlsx', got {format!r}")

        if run_name is None:
            run_name = "run_" + datetime.now().strftime("%Y%m%d_%H%M%S")

        os.makedirs(output_dir, exist_ok=True)

        # Build the three logical tables as dicts of column->list (used by both
        # the pandas path and the csv-stdlib fallback).
        props_table = {
            "time_years": self.time_years,
            "pH": self.pH,
            "gems_code": self.gems_code,
        }
        for sp, vals in self.rates_mol_per_s.items():
            props_table[f"rate_{sp}_mol_s"] = vals
        for name, vals in self.extras.items():
            props_table[name] = vals

        aq_table = {"time_years": self.time_years}
        aq_table.update(self.aqueous_mmol_per_L)

        solids_table = {"time_years": self.time_years}
        solids_table.update(self.solids_cm3)

        meta_payload = dict(metadata) if metadata else {}
        meta_payload.setdefault("run_name", run_name)
        meta_payload.setdefault("saved_at", datetime.now().isoformat(timespec="seconds"))
        meta_payload.setdefault("n_timesteps", len(self.time_years))
        meta_payload.setdefault(
            "n_failed_steps",
            sum(1 for x in self.pH if isinstance(x, float) and math.isnan(x)),
        )

        # --- CSV path: three files in a per-run subdirectory ---------------
        if format == "csv":
            run_dir = os.path.join(output_dir, run_name)
            if os.path.exists(run_dir) and not overwrite:
                raise FileExistsError(
                    f"{run_dir} already exists. Pass overwrite=True or choose a "
                    "different run_name."
                )
            os.makedirs(run_dir, exist_ok=True)

            try:
                import pandas as pd
                pd.DataFrame(props_table).to_csv(
                    os.path.join(run_dir, "properties.csv"), index=False)
                pd.DataFrame(aq_table).to_csv(
                    os.path.join(run_dir, "aqueous.csv"), index=False)
                pd.DataFrame(solids_table).to_csv(
                    os.path.join(run_dir, "solids.csv"), index=False)
            except ImportError:
                # Stdlib fallback so saving works even without pandas.
                import csv
                def _write_csv(path, table):
                    cols = list(table.keys())
                    n = len(table[cols[0]]) if cols else 0
                    with open(path, "w", newline="") as f:
                        w = csv.writer(f)
                        w.writerow(cols)
                        for i in range(n):
                            w.writerow([table[c][i] for c in cols])
                _write_csv(os.path.join(run_dir, "properties.csv"), props_table)
                _write_csv(os.path.join(run_dir, "aqueous.csv"),    aq_table)
                _write_csv(os.path.join(run_dir, "solids.csv"),     solids_table)

            with open(os.path.join(run_dir, "metadata.json"), "w") as f:
                json.dump(meta_payload, f, indent=2, default=str)

            return os.path.abspath(run_dir)

        # --- XLSX path: single workbook with three (+1) sheets -------------
        try:
            import pandas as pd
        except ImportError as e:
            raise ImportError(
                "format='xlsx' requires pandas (and openpyxl). Install them or "
                "use format='csv'."
            ) from e

        xlsx_path = os.path.join(output_dir, f"{run_name}.xlsx")
        if os.path.exists(xlsx_path) and not overwrite:
            raise FileExistsError(
                f"{xlsx_path} already exists. Pass overwrite=True or choose a "
                "different run_name."
            )

        # Flatten metadata into a two-column sheet (key, value)
        meta_rows = []
        def _flatten(prefix, obj):
            if isinstance(obj, dict):
                for k, v in obj.items():
                    _flatten(f"{prefix}.{k}" if prefix else k, v)
            elif isinstance(obj, (list, tuple)):
                meta_rows.append((prefix, ", ".join(str(x) for x in obj)))
            else:
                meta_rows.append((prefix, obj))
        _flatten("", meta_payload)

        with pd.ExcelWriter(xlsx_path, engine="openpyxl") as writer:
            pd.DataFrame(props_table).to_excel(writer, sheet_name="properties", index=False)
            pd.DataFrame(aq_table).to_excel(writer,    sheet_name="aqueous",    index=False)
            pd.DataFrame(solids_table).to_excel(writer, sheet_name="solids",    index=False)
            pd.DataFrame(meta_rows, columns=["key", "value"]).to_excel(
                writer, sheet_name="metadata", index=False)

        return os.path.abspath(xlsx_path)

    # ---- Plotting helpers -------------------------------------------------

    def plot_aqueous_and_pH(self, title: str = "Aqueous concentrations and pH",
                            log_x: bool = False, figsize=(7, 5)):
        import matplotlib.pyplot as plt
        fig, ax1 = plt.subplots(figsize=figsize)
        for elem, series in self.aqueous_mmol_per_L.items():
            ax1.plot(self.time_years, series, label=elem)
        ax1.set_xlabel("time (years)")
        ax1.set_ylabel("total conc. in solution (mmol/L)")
        ax1.set_title(title)
        ax1.grid(True)
        if log_x:
            ax1.set_xscale("log")
        ax2 = ax1.twinx()
        ax2.plot(self.time_years, self.pH, "k--", label="pH")
        ax2.set_ylabel("pH")
        l1, lab1 = ax1.get_legend_handles_labels()
        l2, lab2 = ax2.get_legend_handles_labels()
        ax1.legend(l1 + l2, lab1 + lab2, loc="upper right")
        fig.tight_layout()
        return fig, (ax1, ax2)

    def plot_solid_volumes(self, title: str = "Solid phase volume evolution",
                           secondary_series: Optional[str] = None,
                           secondary_label: Optional[str] = None,
                           ylim=None, log_x: bool = False, figsize=(12, 6)):
        """Stacked area plot of all solid phase groups.

        ``secondary_series`` : name in ``extras`` (e.g. "Ca/Si in CSH") plotted
        on a twin y-axis. Set to None to skip.
        """
        import matplotlib.pyplot as plt
        x = self.time_years
        labels = list(self.solids_cm3.keys())
        y_values = [self.solids_cm3[g] for g in labels]
        colors = plt.cm.tab20c(np.linspace(0, 1, len(labels)))
        hatches = ["/", "\\", "|", "-", "+", "x", "o", "O", ".", "*",
                   "//", "xx", "\\\\", "||", "++", "oo", "**", ".."]
        fig, ax = plt.subplots(figsize=figsize)
        stacks = ax.stackplot(x, *y_values, labels=labels, colors=colors, alpha=0.8)
        for stack, hatch in zip(stacks, hatches):
            stack.set_hatch(hatch)
        ax.set_xlabel("Time (years)")
        ax.set_ylabel("Phase volume (cm³)")
        ax.set_title(title)
        if ylim is not None:
            ax.set_ylim(*ylim)
        if log_x:
            ax.set_xscale("log")
        legends_lines, legends_labels = ax.get_legend_handles_labels()
        if secondary_series and secondary_series in self.extras:
            ax2 = ax.twinx()
            ax2.plot(x, self.extras[secondary_series], "k--",
                     label=secondary_label or secondary_series)
            ax2.set_ylabel(secondary_label or secondary_series)
            l2, lab2 = ax2.get_legend_handles_labels()
            legends_lines += l2
            legends_labels += lab2
        ax.legend(legends_lines, legends_labels, loc="upper right",
                  bbox_to_anchor=(1.35, 1.0), fontsize=10)
        fig.tight_layout()
        return fig, ax


# ===========================================================================
# Main simulation
# ===========================================================================

class CementDissolutionSimulation:
    """Time-stepped, kinetically controlled cement/mortar dissolution.

    Usage
    -----
        sim = CementDissolutionSimulation(
            system_file="gems_files/mortar-dat.lst",
            materials=[QuartzTST()],
            time_config=TimeConfig(t_end_s=years_to_seconds(5000)),
            output_config=OutputConfig(phase_groups={...}),
            composition_override=CompositionOverride(
                species_lower_limits={"Qtz": 26.679},
            ),
            hooks=[],
        )
        sim.run()
        sim.results.plot_aqueous_and_pH()

    Time loop (per step)
    --------------------
      1. Build EngineState (T, P, pH, a_H+, a_OH-, omegas for relevant phases).
      2. For each material:
            r = material.compute_rate(state, current_amount)
            if r > 0:  setSpeciesLowerLimit(species, max(0, n - r*dt))
            store r in results.rates_mol_per_s
      3. Run pre-equilibration hooks on a FRESH COPY of baseline b.
      4. equilibrate(T, P, b'); on failure, try a soft retry with relaxed limit;
         on persistent failure, record NaN and continue.
      5. Run post-equilibration hooks.
      6. Append results.
      7. Advance time; accelerate dt after ``accel_after_s`` if configured.
    """

    def __init__(
        self,
        system_file: str,
        materials: Sequence[KineticMaterial],
        temperature_K: float = 298.15,
        time_config: Optional[TimeConfig] = None,
        output_config: Optional[OutputConfig] = None,
        composition_override: Optional[CompositionOverride] = None,
        hooks: Sequence[CouplingHook] = (),
        verbose: bool = True,
    ):
        if xg is None:
            raise ImportError(
                "xgems is not importable. Install it (e.g. `pip install xgems`) "
                "before constructing a simulation."
            )
        self.system_file = system_file
        self.materials = list(materials)
        self.temperature_K = float(temperature_K)
        self.time_config = time_config or TimeConfig()
        self.output_config = output_config or OutputConfig()
        self.composition_override = composition_override
        self.hooks = list(hooks)
        self.verbose = verbose

        self.engine = None
        self.baseline_b: Optional[np.ndarray] = None

        self.results = SimulationResults()
        self._species_idx: dict = {}
        self._phase_idx: dict = {}
        self._element_idx: dict = {}
        self._setup_done = False

    # ----- setup ----------------------------------------------------------

    def setup(self) -> None:
        eng = xg.ChemicalEngine(self.system_file)
        self.engine = eng

        # Validate output config phase groups against the loaded system.
        for group, phases in self.output_config.phase_groups.items():
            for phase in phases:
                try:
                    self._phase_idx[phase] = eng.indexPhase(phase)
                except Exception as e:
                    raise ValueError(
                        f"OutputConfig.phase_groups: phase '{phase}' (in group "
                        f"'{group}') is not in the loaded GEMS system: {e}"
                    )
        # Aqueous phase index
        try:
            self._phase_idx[self.output_config.aqueous_phase] = eng.indexPhase(
                self.output_config.aqueous_phase
            )
        except Exception as e:
            raise ValueError(
                f"OutputConfig.aqueous_phase: '{self.output_config.aqueous_phase}' "
                f"not found: {e}"
            )

        # Element indices for tracked aqueous elements
        for el in self.output_config.aqueous_elements:
            try:
                self._element_idx[el] = eng.indexElement(el)
            except Exception as e:
                raise ValueError(f"Aqueous element '{el}' not in system: {e}")

        # Validate materials and cache their indices
        for m in self.materials:
            m.validate(eng)
            self._species_idx[m.species_name] = eng.indexSpecies(m.species_name)
            if m.saturation_phase and m.saturation_phase not in self._phase_idx:
                self._phase_idx[m.saturation_phase] = eng.indexPhase(m.saturation_phase)

        # Apply composition overrides BEFORE the initial equilibration. These
        # values then become part of the captured baseline_b.
        self._apply_composition_override()

        # Initial equilibration
        T = self.temperature_K
        P = eng.pressure()
        b = eng.elementAmounts().copy()
        eng.setColdStart()
        code = eng.equilibrate(T, P, b)
        if code not in _GEMS_OK_CODES:
            raise RuntimeError(
                f"Initial equilibration failed with code {code}. "
                "Check system file, temperature, and composition overrides."
            )

        # Capture baseline bulk composition AFTER initial equilibrate and
        # BEFORE the time loop. Per established best practice: hooks receive a
        # fresh copy each step, so their adjustments never compound.
        self.baseline_b = eng.elementAmounts().copy()

        # Cache species we always need
        self._species_idx["H+"] = eng.indexSpecies("H+")
        self._species_idx["OH-"] = eng.indexSpecies("OH-")

        # Initialize result buffers
        for m in self.materials:
            self.results.rates_mol_per_s.setdefault(m.species_name, [])
        for el in self.output_config.aqueous_elements:
            self.results.aqueous_mmol_per_L[el] = []
        for group in self.output_config.phase_groups:
            self.results.solids_cm3[group] = []
        for name in self.output_config.extras:
            self.results.extras[name] = []

        self._setup_done = True

    def _apply_composition_override(self) -> None:
        ov = self.composition_override
        if ov is None:
            return
        eng = self.engine
        for sp, val in ov.species_lower_limits.items():
            try:
                eng.setSpeciesLowerLimit(sp, val)
            except Exception as e:
                raise ValueError(f"setSpeciesLowerLimit('{sp}', {val}) failed: {e}")
        for sp, val in ov.species_upper_limits.items():
            try:
                eng.setSpeciesUpperLimit(sp, val)
            except Exception as e:
                raise ValueError(f"setSpeciesUpperLimit('{sp}', {val}) failed: {e}")
        if ov.bulk_elements:
            b = eng.elementAmounts().copy()
            for el, val in ov.bulk_elements.items():
                try:
                    idx = eng.indexElement(el)
                except Exception as e:
                    raise ValueError(f"bulk_elements: '{el}' not in system: {e}")
                b[idx] = val
            # Overwrite the engine's b by feeding it into the next equilibrate.
            # We store it for setup() to pick up; cleanest path is to call
            # equilibrate ourselves here so the override is committed.
            eng.setColdStart()
            code = eng.equilibrate(self.temperature_K, eng.pressure(), b)
            if code not in _GEMS_OK_CODES:
                raise RuntimeError(
                    f"Equilibration with overridden bulk composition failed "
                    f"(code {code})."
                )

    # ----- state snapshot -------------------------------------------------

    def _snapshot(self, time_s: float, step: int) -> EngineState:
        eng = self.engine
        ln_a = eng.lnActivities()
        aHplus = math.exp(ln_a[self._species_idx["H+"]])
        aOHminus = math.exp(ln_a[self._species_idx["OH-"]])
        # Collect omegas only for phases that some material actually needs
        sat: dict = {}
        for m in self.materials:
            if m.saturation_phase and m.saturation_phase not in sat:
                try:
                    sat[m.saturation_phase] = 10.0 ** eng.phaseSatIndex(
                        self._phase_idx[m.saturation_phase]
                    )
                except Exception:
                    sat[m.saturation_phase] = float("nan")
        return EngineState(
            time_s=time_s,
            step=step,
            temperature_K=self.temperature_K,
            pressure_Pa=eng.pressure(),
            pH=eng.pH(),
            aHplus=aHplus,
            aOHminus=aOHminus,
            saturation_indices=sat,
            engine=eng,
        )

    # ----- result recording ----------------------------------------------

    def _record_step(self, time_s: float, code: int, rates: dict, ok: bool) -> None:
        eng = self.engine
        self.results.time_years.append(seconds_to_years(time_s))
        self.results.gems_code.append(code)
        for sp, r in rates.items():
            self.results.rates_mol_per_s[sp].append(r)

        if not ok:
            # record NaNs for everything else
            self.results.pH.append(float("nan"))
            for el in self.output_config.aqueous_elements:
                self.results.aqueous_mmol_per_L[el].append(float("nan"))
            for group in self.output_config.phase_groups:
                self.results.solids_cm3[group].append(float("nan"))
            for name in self.output_config.extras:
                self.results.extras[name].append(float("nan"))
            return

        # pH
        self.results.pH.append(eng.pH())

        # Aqueous concentrations (mmol/L of aqueous phase)
        aq_idx = self._phase_idx[self.output_config.aqueous_phase]
        aq_volume_L = eng.phaseVolume(aq_idx) * 1000.0  # m^3 -> L
        aq_elements = eng.elementAmountsInPhase(aq_idx)
        for el in self.output_config.aqueous_elements:
            idx = self._element_idx[el]
            if aq_volume_L > 0:
                self.results.aqueous_mmol_per_L[el].append(
                    1000.0 * aq_elements[idx] / aq_volume_L
                )
            else:
                self.results.aqueous_mmol_per_L[el].append(float("nan"))

        # Solid phase volumes (group totals, cm^3)
        for group, members in self.output_config.phase_groups.items():
            vol_m3 = sum(eng.phaseVolume(self._phase_idx[p]) for p in members)
            self.results.solids_cm3[group].append(vol_m3 * 1e6)

        # Extras
        for name, fn in self.output_config.extras.items():
            try:
                self.results.extras[name].append(float(fn(eng)))
            except Exception:
                self.results.extras[name].append(float("nan"))

    # ----- main loop ------------------------------------------------------

    def run(self) -> SimulationResults:
        if not self._setup_done:
            self.setup()

        eng = self.engine
        T = self.temperature_K
        P = eng.pressure()
        tc = self.time_config

        time_s = 0.0
        dt = tc.dt_initial_s
        step = 0
        code = 2

        # First step starts at t = dt (matches original notebook)
        time_s += dt

        while time_s <= tc.t_end_s:
            step += 1
            if self.verbose:
                print(
                    f"\rt = {seconds_to_years(time_s):.2f} yr   "
                    f"dt = {seconds_to_years(dt):.4g} yr   "
                    f"code = {code}    ",
                    end="",
                )

            # 1) snapshot
            state = self._snapshot(time_s, step)

            # 2) per-material rates + lower-limit updates
            rates: dict = {}
            for m in self.materials:
                idx = self._species_idx[m.species_name]
                amount = eng.speciesAmount(idx)
                r = m.compute_rate(state, amount)
                # NaN-safe: a rate that came out NaN is treated as "no update"
                if not np.isfinite(r) or r <= 0.0:
                    rates[m.species_name] = 0.0
                    continue
                rates[m.species_name] = r
                new_lower = max(0.0, amount - r * dt)
                try:
                    eng.setSpeciesLowerLimit(m.species_name, new_lower)
                except Exception as e:
                    warnings.warn(
                        f"setSpeciesLowerLimit failed for {m.species_name}: {e}"
                    )

            # 3) pre-equilibration hooks on a fresh copy of baseline b
            b = self.baseline_b.copy()
            for h in self.hooks:
                b = h.before_equilibrate(eng, state, b)

            # 4) equilibrate, with a soft retry
            eng.setColdStart()
            code = eng.equilibrate(T, P, b)
            ok = code in _GEMS_OK_CODES
            if not ok:
                # Soft retry: relax the most-recently-set lower limit by 1 ULP-ish
                for m in self.materials:
                    if rates.get(m.species_name, 0.0) > 0:
                        idx = self._species_idx[m.species_name]
                        amount = eng.speciesAmount(idx)
                        relaxed = max(0.0, amount * (1.0 - 1e-14))
                        try:
                            eng.setSpeciesLowerLimit(m.species_name, relaxed)
                        except Exception:
                            pass
                try:
                    code = eng.reequilibrate(False)
                    ok = code in _GEMS_OK_CODES
                except Exception:
                    ok = False

            if not ok:
                # Record NaN and continue rather than crash.
                warnings.warn(
                    f"Equilibration failed at t = {seconds_to_years(time_s):.2f} yr "
                    f"(code {code}). Recording NaN and continuing."
                )
                self._record_step(time_s, code, rates, ok=False)
            else:
                # 5) post-equilibration hooks
                for h in self.hooks:
                    h.after_equilibrate(eng, state)
                # 6) record
                self._record_step(time_s, code, rates, ok=True)

            # 7) advance time
            time_s += dt
            if time_s > tc.accel_after_s:
                dt = dt * tc.accel_factor
                if tc.dt_max_s is not None:
                    dt = min(dt, tc.dt_max_s)

        if self.verbose:
            print("\nSimulation finished.")
        return self.results


# ===========================================================================
# Convenience: a Ca/Si-in-CSH extra
# ===========================================================================

def ca_si_in_phase_factory(phase_name: str = "CSHQ",
                           numerator: str = "Ca",
                           denominator: str = "Si") -> Callable:
    """Build an OutputConfig.extras callable that returns Ca/Si (or any element
    ratio) in a named phase.

    Usage::

        output_config = OutputConfig(
            phase_groups={...},
            extras={"Ca/Si in CSH": ca_si_in_phase_factory("CSHQ", "Ca", "Si")},
        )
    """

    def _extra(engine) -> float:
        idx_phase = engine.indexPhase(phase_name)
        idx_num = engine.indexElement(numerator)
        idx_den = engine.indexElement(denominator)
        amounts = engine.elementAmountsInPhase(idx_phase)
        denom = amounts[idx_den]
        if denom == 0:
            return float("nan")
        return amounts[idx_num] / denom

    return _extra
