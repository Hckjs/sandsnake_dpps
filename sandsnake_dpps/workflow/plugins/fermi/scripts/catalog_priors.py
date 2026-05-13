from dataclasses import dataclass

import numpy as np
from astropy.table import QTable

EXTRAGALACTIC_CLASSES = {
    "agn",
    "bcu",
    "bll",
    "fsrq",
    "rdg",
    "sey",
    "nlsy1",
    "ssrq",
    "sbg",
    "css",
    "gal",
}

GALACTIC_CLASSES = {
    "psr",
    "msp",
    "pwn",
    "snr",
    "spp",
    "glc",
    "hmb",
    "lmb",
    "bin",
    "nov",
    "sfr",
    "mc",
    "hii",
}

UNCERTAIN_CLASSES = {
    "",
    "--",
    "nan",
    "none",
    "unk",
}


@dataclass(slots=True)
class RedshiftPriorConfig:
    class_column: str = "CLASS1"
    sed_column: str = "SED_class"
    redshift_column: str = "Redshift"
    lower_quantile: float = 0.16
    upper_quantile: float = 0.84
    min_sources_per_group: int = 50


def normalize_label(value) -> str:
    if np.ma.is_masked(value):
        return ""

    if isinstance(value, bytes):
        value = value.decode()

    text = str(value).strip().lower()
    if text in {"", "--", "nan", "none", "masked"}:
        return ""

    return text


def get_z_class(row, config: RedshiftPriorConfig) -> str:
    """
    Classify source origin for redshift treatment.

    Returns
    -------
    z_class : {"gal", "egal", "unk"}
        gal  : Galactic source class, use z = 0
        egal : Extragalactic source class
        unk  : Unknown / uncertain origin
    """
    source_class = normalize_label(row[config.class_column])

    if source_class in GALACTIC_CLASSES:
        return "gal"

    if source_class in EXTRAGALACTIC_CLASSES:
        return "egal"

    if source_class in UNCERTAIN_CLASSES:
        return "unk"

    return "unk"


def get_prior_group(row, config: RedshiftPriorConfig) -> str:
    """
    Assign only BLL/FSRQ-like redshift-prior groups.

    This intentionally does not create priors for all extragalactic classes.
    It follows the AGNpop-style two-population treatment:

    - bll  -> bll
    - fsrq -> fsrq
    - bcu + no SED_class -> bll
    - bcu + HSP/ISP     -> bll
    - bcu + LSP         -> fsrq

    All other classes get no prior group.
    """
    source_class = normalize_label(row[config.class_column])

    sed_class = ""
    if config.sed_column in row.colnames:
        sed_class = normalize_label(row[config.sed_column])

    if source_class == "bll":
        return "bll"

    if source_class == "fsrq":
        return "fsrq"

    if source_class == "bcu":
        if sed_class in {"", "hsp", "isp"}:
            return "bll"

        if sed_class == "lsp":
            return "fsrq"

        return ""

    return ""


def finite_positive_redshift(table: QTable, column: str) -> np.ndarray:
    z = np.asarray(table[column], dtype=float)
    return np.isfinite(z) & (z > 0)


def empty_redshift_prior_table() -> QTable:
    return QTable(
        names=[
            "z_prior_group",
            "n_z",
            "z_min",
            "z_q_low",
            "z_q_med",
            "z_q_high",
            "z_max",
        ],
        dtype=[
            "U16",
            int,
            float,
            float,
            float,
            float,
            float,
        ],
    )


def build_redshift_prior_table(
    catalog_table: QTable,
    config: RedshiftPriorConfig = RedshiftPriorConfig(),
) -> QTable:
    """
    Build empirical redshift-prior table only from measured BLL and FSRQ redshifts.

    BCUs are not used to build the distributions. They are only assigned to one of
    the two prior groups later.
    """
    source_classes = np.array(
        [normalize_label(row[config.class_column]) for row in catalog_table]
    )
    z = np.asarray(catalog_table[config.redshift_column], dtype=float)

    valid_z = np.isfinite(z) & (z > 0)

    prior_table = empty_redshift_prior_table()

    for group in ["bll", "fsrq"]:
        mask = valid_z & (source_classes == group)
        n_group = int(mask.sum())

        if n_group < config.min_sources_per_group:
            continue

        z_group = z[mask]

        prior_table.add_row(
            {
                "z_prior_group": group,
                "n_z": n_group,
                "z_min": float(np.min(z_group)),
                "z_q_low": float(np.quantile(z_group, config.lower_quantile)),
                "z_q_med": float(np.quantile(z_group, 0.50)),
                "z_q_high": float(np.quantile(z_group, config.upper_quantile)),
                "z_max": float(np.max(z_group)),
            }
        )

    return prior_table


def lookup_prior(prior_table: QTable, group: str):
    if not group:
        return None

    match = prior_table[prior_table["z_prior_group"] == group]
    if len(match) == 1:
        return match[0]

    return None


def append_redshift_prior_columns(
    catalog_table: QTable,
    prior_table: QTable,
    config: RedshiftPriorConfig = RedshiftPriorConfig(),
) -> QTable:
    """
    Append redshift columns used for later spectral simulations.

    Rules
    -----
    Galactic sources:
        z = 0, not measured.

    Unknown-origin sources:
        z = NaN.

    Extragalactic sources with measured z:
        use measured z.

    Extragalactic bll/fsrq/bcu without measured z:
        use BLL/FSRQ quantiles if a prior group can be assigned.

    Other extragalactic sources without measured z:
        z = NaN, no prior.
    """
    z_classes = []
    prior_groups = []
    has_measured_z = []
    z_measured = []
    z_q_low = []
    z_q_med = []
    z_q_high = []
    z_prior_n = []
    z_source = []

    has_redshift_column = config.redshift_column in catalog_table.colnames

    for row in catalog_table:
        z_class = get_z_class(row, config)
        prior_group = get_prior_group(row, config)

        z_catalog = np.nan
        if has_redshift_column:
            try:
                z_catalog = float(row[config.redshift_column])
            except (TypeError, ValueError):
                z_catalog = np.nan

        measured = np.isfinite(z_catalog) and z_catalog > 0

        z_classes.append(z_class)
        prior_groups.append(prior_group)

        if z_class == "gal":
            has_measured_z.append(False)
            z_measured.append(np.nan)
            z_q_low.append(np.nan)
            z_q_med.append(np.nan)
            z_q_high.append(np.nan)
            z_prior_n.append(0)
            z_source.append("galactic_zero")
            continue

        if z_class == "unk":
            has_measured_z.append(False)
            z_measured.append(np.nan)
            z_q_low.append(np.nan)
            z_q_med.append(np.nan)
            z_q_high.append(np.nan)
            z_prior_n.append(0)
            z_source.append("unknown_origin")
            continue

        if measured:
            has_measured_z.append(True)
            z_measured.append(z_catalog)
            z_q_low.append(z_catalog)
            z_q_med.append(z_catalog)
            z_q_high.append(z_catalog)
            z_prior_n.append(1)
            z_source.append("measured")
            continue

        prior = lookup_prior(prior_table, prior_group)

        if prior is not None:
            has_measured_z.append(False)
            z_measured.append(np.nan)
            z_q_low.append(float(prior["z_q_low"]))
            z_q_med.append(float(prior["z_q_med"]))
            z_q_high.append(float(prior["z_q_high"]))
            z_prior_n.append(int(prior["n_z"]))
            z_source.append(f"prior_{prior_group}")
            continue

        has_measured_z.append(False)
        z_measured.append(np.nan)
        z_q_low.append(np.nan)
        z_q_med.append(np.nan)
        z_q_high.append(np.nan)
        z_prior_n.append(0)
        z_source.append("no_prior")

    catalog_table["z_class"] = z_classes
    catalog_table["z_prior_group"] = prior_groups
    catalog_table["has_measured_redshift"] = has_measured_z
    catalog_table["z_measured"] = z_measured
    catalog_table["z_q_low"] = z_q_low
    catalog_table["z_q_med"] = z_q_med
    catalog_table["z_q_high"] = z_q_high
    catalog_table["z_prior_n"] = z_prior_n
    catalog_table["z_source"] = z_source

    return catalog_table
