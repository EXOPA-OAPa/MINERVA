"""
Run context: configuration, folders, run options and reference data shared by
all the simulations of a run. It replaces the `MINERVA` class of the version
used for the paper (same attributes, same loading order).
"""
import os
from pathlib import Path

import numpy as np
import yaml

from . import utils as minerva_utils
from .fast_gls import make_gls

# Folder with the reference data shipped with MINERVA (tables, tracks, default cadences)
DEFAULT_DATA_DIR = Path(__file__).resolve().parents[2] / "data"   # <repo>/data (src layout)

# Analysis parameters. Defaults are the values of the version used for the paper.
DEFAULT_ANALYSIS = {
    "n_start": 10,                     # observations at the first step of the sequential recovery
    "step": 5,                         # observations added at each step
    "sampling_strategy": "sequential", # order in which observations are added: "sequential" or "random"
    "fap_star": 0.10,                  # false alarm probability for stellar (intruder) peaks
    "fap_planet": 0.50,                # false alarm probability for the planet peak
    "delta_bic_threshold": 6.0,        # minimum Delta BIC (planet model vs null hypothesis) for a detection
    "random_seed": 42,                 # seed of the random generator (seed + i for the i-th catalog target)
    "gls_method": "recurrence",        # GLS implementation: "recurrence" (fastest), "exact", or "reference" (paper)
}

# Run options. Defaults are the values of the version used for the paper.
DEFAULT_OPTIONS = {
    "debug_plot": True,           # debug figures of pre-whitening / reverse pre-whitening
    "debug_fit_plot": True,       # debug figures of the fits
    "export_debug_pdf": True,     # write the debug figures to Output/Plots/DEBUG_<target>.pdf
    "show_on_screen": False,      # show the figures interactively (blocks the run)
    "debug_pdf_all_iter": True,   # save GLS / fit diagnostics for every iteration (False: first 3 only)
    "bootstrap": False,           # resample the stellar parameters within their uncertainties
    "rv_series_pdf": True,        # RV time series figure: True (every simulation), False, or N (first N)
    "n_workers": 1,               # catalog targets run in parallel: number of processes, or "auto"
}


def _as_dir(path):
    """Directory path as a string ending with a separator (the code concatenates file names)."""
    return os.path.join(str(path), "")


class RunContext:
    """Configuration, folders, options and reference data of a MINERVA run."""

    def __init__(self, config_path, overrides=None):
        """overrides: optional values that replace the configuration, from the command line:
        {"n_iterations": int, "n_workers": int or "auto"}."""
        self.overrides = dict(overrides or {})
        config_path = Path(config_path).resolve()
        self.config_path = str(config_path)
        with open(config_path, "r") as f:
            raw = yaml.safe_load(f)
        self.config = minerva_utils.load_config(str(config_path))

        # --- folders -----------------------------------------------------------
        paths = raw.get("paths") or {}
        base = config_path.parent
        self.data_dir = Path(paths.get("data_dir") or DEFAULT_DATA_DIR)
        self.user_dir = Path(paths.get("user_dir") or base)
        out_dir = Path(paths.get("output_dir") or base / "Output")
        if not self.data_dir.is_absolute():
            self.data_dir = base / self.data_dir
        if not self.user_dir.is_absolute():
            self.user_dir = base / self.user_dir
        if not out_dir.is_absolute():
            out_dir = base / out_dir

        self.InDir = _as_dir(self.data_dir)
        self.OutDir = _as_dir(out_dir)
        self.FigDir = _as_dir(out_dir / "Plots")
        os.makedirs(self.OutDir, exist_ok=True)
        os.makedirs(self.FigDir, exist_ok=True)

        # --- run options -------------------------------------------------------
        self.options = dict(DEFAULT_OPTIONS)
        self.options.update(raw.get("run_options") or {})
        unknown = set(self.options) - set(DEFAULT_OPTIONS)
        if unknown:
            raise ValueError(f"Unknown run_options in {config_path.name}: {sorted(unknown)}")
        if self.overrides.get("n_workers") is not None:
            self.options["n_workers"] = self.overrides["n_workers"]
        self.SHOW_ON_SCREEN = self.options["show_on_screen"]
        self.EXPORT_DEBUG_PDF = self.options["export_debug_pdf"]

        # --- analysis parameters ----------------------------------------------
        self.analysis = dict(DEFAULT_ANALYSIS)
        self.analysis.update(raw.get("analysis") or {})
        unknown = set(self.analysis) - set(DEFAULT_ANALYSIS)
        if unknown:
            raise ValueError(f"Unknown analysis parameters in {config_path.name}: {sorted(unknown)}")
        # GLS implementation used by the sequential analysis
        self.Gls = make_gls(self.analysis["gls_method"])

        # --- general variables -------------------------------------------------
        self.inputfile_sampling = self.config['name_sampling']
        self.inputfile_sampling_usp = self.config['name_sampling_usp']
        self.inputfile_nsts = self.config['inputfile_nsts']
        self.inputfile_mass_age = self.config['inputfile_mass_age']
        self.sig_floor = self.config['sig_floor']
        self.sig_pl = self.config['sig_pl']
        self.sig_best = self.config['sig_best']
        self.inputfile_activity = self.config['inputfile_activity']

        # --- reference data (same order as the paper version) --------------------
        self.time_sampling, self.time_sampling_usp = minerva_utils.import_file_sampling(
            InDir="", inputfile_sampling=self.resolve(self.inputfile_sampling),
            inputfile_sampling_usp=self.resolve(self.inputfile_sampling_usp))
        (self.tmp_mag, self.t_exp, self.errRV_F, self.errRV_G,
         self.errRV_K, self.errRV_M) = minerva_utils.import_nsts_noise(
            InDir="", inputfile_nsts=self.resolve(self.inputfile_nsts))
        self.mass_tmp, self.max_age = minerva_utils.import_max_age_from_mass(
            InDir="", inputfile_mass_age=self.resolve(self.inputfile_mass_age))
        self.df_activity = minerva_utils.import_file_activity(
            InDir="", inputfilename=self.resolve(self.inputfile_activity))

        # --- observing seasons of the default cadences ---------------------------
        xs = minerva_utils.find_seasons_automatic(self.time_sampling)
        for i, s in enumerate(minerva_utils.define_seasons(self.time_sampling, xs)):
            setattr(self, f's{i+1}', s)
        xs_usp = minerva_utils.find_seasons_automatic(self.time_sampling_usp)
        for i, s in enumerate(minerva_utils.define_seasons(self.time_sampling_usp, xs_usp)):
            setattr(self, f's{i+1}_usp', s)

        # --- Monte Carlo ---------------------------------------------------------
        self.numero_iterazioni = int(self.config.get('n_iterations', 100))
        if self.overrides.get("n_iterations") is not None:
            self.numero_iterazioni = int(self.overrides["n_iterations"])
        self.numerorv = np.zeros(self.numero_iterazioni)

    def resolve(self, filename):
        """Full path of an input file: absolute paths are kept; otherwise the file is
        looked up first in the user folder (catalogs, cadences), then in data/."""
        if filename is None:
            return None
        filename = str(filename)
        if os.path.isabs(filename):
            return filename
        user_file = self.user_dir / filename
        if user_file.exists():
            return str(user_file)
        return str(self.data_dir / filename)
