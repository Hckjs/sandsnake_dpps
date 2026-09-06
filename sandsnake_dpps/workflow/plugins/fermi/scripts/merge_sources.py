import argparse
from collections.abc import Sequence

from astropy.table import QTable


def merge_catalog_group(
    sources: Sequence[str],
    output_file: str,
    group_name: str,
    *,
    overwrite: bool,
    expected_catalog: str,
) -> None:
    rows = []

    for source in sources:
        table = QTable.read(source, format="ascii.ecsv")
        if len(table) == 0:
            raise ValueError(f"Input table is empty: {source}")
        if "catalog" not in table.colnames:
            raise ValueError(
                f"Input table is missing required column 'catalog': {source}"
            )
        if str(table[0]["catalog"]) != expected_catalog:
            raise ValueError(
                f"Input table {source} has catalog {table[0]['catalog']!r}; "
                f"expected {expected_catalog!r}"
            )
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
    fhl3_sources: Sequence[str],
    fhl4_sources: Sequence[str],
) -> None:
    catalogs = {
        "4FGL": ("4FGL_DR4", list(fgl_sources)),
        "3FHL": ("3FHL", list(fhl3_sources)),
        "4FHL": ("4FHL", list(fhl4_sources)),
    }

    first = True
    for name, (expected_catalog, sources) in catalogs.items():
        if not sources:
            continue

        merge_catalog_group(
            sources=sources,
            output_file=output_file,
            group_name=name,
            overwrite=first,
            expected_catalog=expected_catalog,
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
        "--fhl3-sources",
        nargs="*",
        default=[],
        help="Input ECSV files for the 3FHL group.",
    )
    parser.add_argument(
        "--fhl4-sources",
        nargs="*",
        default=[],
        help="Input ECSV files for the 4FHL group.",
    )
    parser.add_argument(
        "-o",
        "--output",
        required=True,
        help="Output HDF5 file.",
    )
    return parser.parse_args()


def main_from_snakemake(snakemake) -> None:
    main(
        output_file=snakemake.output[0],
        fgl_sources=snakemake.input.fgl_sources,
        fhl3_sources=snakemake.input.fhl3_sources,
        fhl4_sources=snakemake.input.fhl4_sources,
    )


def main_from_args(args: argparse.Namespace) -> None:
    main(
        output_file=args.output,
        fgl_sources=args.fgl_sources,
        fhl3_sources=args.fhl3_sources,
        fhl4_sources=args.fhl4_sources,
    )


if "snakemake" in globals():
    main_from_snakemake(snakemake)  # noqa: F821
elif __name__ == "__main__":
    args = parse_args()
    main_from_args(args)
