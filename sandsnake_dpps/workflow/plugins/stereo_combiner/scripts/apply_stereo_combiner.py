import argparse
import logging
from tqdm.auto import tqdm

from ctapipe.reco.reconstructor import ReconstructionProperty
from ctapipe.io import HDF5EventSource, DataWriter
from ctapipe.io.tableio import TelListToMaskTransform

from .stereo_combiner import StereoCombiner

log = logging.getLogger(__name__)
log.setLevel(logging.DEBUG)


def main(input, output, combiner):
    if combiner == "StereoDispCombinerSS":
        stereo_combiner = StereoCombiner.from_name(
            "StereoDispCombiner",
            prefix="disp",
            property=ReconstructionProperty.GEOMETRY,
            weights="konrad",
            sign_score_limit=0.8,
        )
    else:
        stereo_combiner = StereoCombiner.from_name(
            combiner,
            prefix="disp",
            property=ReconstructionProperty.GEOMETRY,
            weights="konrad",
        )

    with (
        HDF5EventSource(input_url=input) as source,
        DataWriter(event_source=source, write_dl2=True, output_path=output) as writer,
    ):
        for event in tqdm(
            source,
            desc=source.__class__.__name__,
            total=source.max_events,
            unit="ev",
            disable=False,
        ):
            stereo_combiner(event)

            # Transformation of "telescopes" field in ReconstructedContainer doesn't
            # work with the HDF5EventSource. But its working with the DataWriter.
            # So we have to transform them by hand before writing with DataWriter.
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


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--input", required=True)
    parser.add_argument("--output", required=True)
    parser.add_argument("--combiner", required=True)
    args = parser.parse_args()

    main(
        args.input,
        args.output,
        args.combiner,
    )
