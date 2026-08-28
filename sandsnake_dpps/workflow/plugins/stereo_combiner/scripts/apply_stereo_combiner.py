import argparse
import logging

from tqdm.auto import tqdm

from ctapipe.io import DataWriter, HDF5EventSource
from ctapipe.io.tableio import TelListToMaskTransform
from ctapipe.reco.reconstructor import ReconstructionProperty

from plugins.stereo_combiner.scripts.stereo_combiner import StereoCombiner


log = logging.getLogger(__name__)
log.setLevel(logging.DEBUG)


# Should be identical to used Quality Criterias in main analysis (and HillasReconstructor
# in process.yml if oyu want to compare these)
QUALITY_CRITERIA = [
    ("enough intensity", "parameters.hillas.intensity > 80"),
    ("Positive width", "parameters.hillas.width.value > 0"),
    ("Positive length", "parameters.hillas.length.value > 0"),
    ("enough pixels", "parameters.morphology.n_pixels > 3"),
    ("not clipped", "parameters.leakage.intensity_width_2 < 0.2"),
]


def apply_quality_selection(event, quality_query, prefix="disp", min_multiplicity=2):
    """
    Temporarily invalidate telescope-wise reconstructions that do not pass
    the configured StereoQualityQuery.

    If fewer than ``min_multiplicity`` telescopes remain, invalidate all
    telescope-wise reconstructions for this event so that no valid stereo
    reconstruction can be produced.

    Returns
    -------
    original_validity : dict
        Original ``is_valid`` states for restoring them after stereo
        reconstruction.
    multiplicity : int
        Number of telescopes passing both the mono validity requirement
        and the quality query.
    """
    original_validity = {}
    valid_telescopes = []

    for tel_id, dl2 in event.dl2.tel.items():
        if prefix not in dl2.geometry:
            continue

        mono = dl2.geometry[prefix]
        original_validity[tel_id] = mono.is_valid

        if not mono.is_valid:
            continue

        if tel_id not in event.dl1.tel:
            mono.is_valid = False
            continue

        parameters = event.dl1.tel[tel_id].parameters

        if parameters is None:
            mono.is_valid = False
            continue

        passes_quality = all(quality_query(parameters=parameters))

        if passes_quality:
            valid_telescopes.append(tel_id)
        else:
            mono.is_valid = False

    multiplicity = len(valid_telescopes)

    # Enforce stereo multiplicity >= 2.
    if multiplicity < min_multiplicity:
        for tel_id in valid_telescopes:
            event.dl2.tel[tel_id].geometry[prefix].is_valid = False

    return original_validity, multiplicity


def restore_mono_validity(event, original_validity, prefix="disp"):
    """Restore the original mono reconstruction validity flags."""
    for tel_id, is_valid in original_validity.items():
        event.dl2.tel[tel_id].geometry[prefix].is_valid = is_valid


def main(input, output, combiner):
    if combiner == "StereoDispCombinerAngCut":
        stereo_combiner = StereoCombiner.from_name(
            "StereoDispCombiner",
            prefix="disp",
            property=ReconstructionProperty.GEOMETRY,
            weights="aspect-weighted-intensity",
            min_ang_diff=20,
        )
    else:
        stereo_combiner = StereoCombiner.from_name(
            combiner,
            prefix="disp",
            property=ReconstructionProperty.GEOMETRY,
            weights="aspect-weighted-intensity",
        )

    # Explicitly use the same telescope-quality criteria for all combiners.
    stereo_combiner.quality_query.quality_criteria = QUALITY_CRITERIA

    with (
        HDF5EventSource(input_url=input) as source,
        DataWriter(
            event_source=source,
            write_dl2=True,
            output_path=output,
        ) as writer,
    ):
        for event in tqdm(
            source,
            desc=source.__class__.__name__,
            total=source.max_events,
            unit="ev",
            disable=False,
        ):
            # Apply identical telescope-level quality cuts before running
            # the stereo reconstruction.
            original_validity, _ = apply_quality_selection(
                event,
                quality_query=stereo_combiner.quality_query,
                prefix="disp",
                min_multiplicity=2,
            )

            stereo_combiner(event)

            # Restore the original telescope-wise mono validity flags.
            # Only the newly calculated stereo reconstruction should contain
            # the result of the quality selection.
            restore_mono_validity(
                event,
                original_validity,
                prefix="disp",
            )

            # Transformation of "telescopes" field in ReconstructedContainer
            # doesn't work with the HDF5EventSource. But it works with the
            # DataWriter. So we have to transform them by hand before writing
            # with DataWriter.
            trafo = TelListToMaskTransform(source.subarray)

            telescopes_list = event.dl2.stereo.geometry["disp"].telescopes
            event.dl2.stereo.geometry["disp"].telescopes = trafo(telescopes_list)

            for reco in event.dl2.stereo.keys():
                for algo in event.dl2.stereo[reco].keys():
                    telescopes_list = event.dl2.stereo[reco][algo].telescopes
                    event.dl2.stereo[reco][algo].telescopes = trafo._inverse(
                        telescopes_list
                    )

            writer(event)

        # Write shower distributions for correct effective area later on
        shower_dists = source.simulated_shower_distributions
        writer.write_simulated_shower_distributions(shower_dists)


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("--input", required=True)
    parser.add_argument("--output", required=True)
    parser.add_argument("--combiner", required=True)
    return parser.parse_args()


def main_from_snakemake(snakemake):
    main(
        snakemake.input.data,
        snakemake.output[0],
        snakemake.params.combiner,
    )


if "snakemake" in globals():
    main_from_snakemake(snakemake)  # noqa: F821
elif __name__ == "__main__":
    args = parse_args()
    main(
        args.input,
        args.output,
        args.combiner,
    )
