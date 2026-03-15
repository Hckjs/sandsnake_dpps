from __future__ import annotations

import argparse
from collections.abc import Sequence

from astropy.table import QTable


def merge_catalog_group(
    sources: Sequence[str],
    output_file: str,
    group_name: str,
    *,
    overwrite: bool,
) -> None:
    rows = []

    for source in sources:
        table = QTable.read(source, format="ascii.ecsv")
        if len(table) == 0:
            raise ValueError(f"Input table is empty: {source}")
        rows.append(table[0])

    if not rows:
        return

    catalog_table = QTable(rows=rows)
    catalog_table.write(
        output_file,
        path=f"/{group_name}",
        format="hdf5",
        overwrite=overwrite,
        append=not overwrite,
        serialize_meta=True,
    )


def main(
    output_file: str,
    fgl_sources: Sequence[str],
    fhl_sources: Sequence[str],
) -> None:
    catalogs = {
        "FGL": list(fgl_sources),
        "FHL": list(fhl_sources),
    }

    first = True
    for name, sources in catalogs.items():
        if not sources:
            continue

        merge_catalog_group(
            sources=sources,
            output_file=output_file,
            group_name=name,
            overwrite=first,
        )
        first = False


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Merge per-source ECSV tables into grouped HDF5 tables."
    )
    parser.add_argument(
        "--fgl-sources",
        nargs="*",
        default=[],
        help="Input ECSV files for the FGL group.",
    )
    parser.add_argument(
        "--fhl-sources",
        nargs="*",
        default=[],
        help="Input ECSV files for the FHL group.",
    )
    parser.add_argument(
        "-o",
        "--output",
        required=True,
        help="Output HDF5 file.",
    )
    return parser.parse_args()


if __name__ == "__main__":
    smk = globals().get("snakemake")

    if smk is not None:
        main(
            output_file=smk.output[0],
            fgl_sources=smk.input.fgl_sources,
            fhl_sources=smk.input.fhl_sources,
        )
    else:
        args = parse_args()
        main(
            output_file=args.output,
            fgl_sources=args.fgl_sources,
            fhl_sources=args.fhl_sources,
        )
