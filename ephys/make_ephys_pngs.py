import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import argschema as ags
import os

from ipfx.dataset.create import create_ephys_data_set
import ipfx.script_utils as su
import ipfx.data_set_features as dsf
import ipfx.stimulus_protocol_analysis as spa
import ipfx.time_series_utils as tsu


class MakeEphysPngsParameters(ags.ArgSchema):
    nwb_dir = ags.fields.InputDir()
    metadata_file = ags.fields.InputFile()
    output_dir = ags.fields.OutputDir()
    n_per_png = ags.fields.Integer(default=3)


def ephys_plot(data_set, sweep_numbers, ax_ap, ax_lsq, subthresh_min_amp=-200):
    lsq_sweeps = data_set.sweep_set(sweep_numbers, enforce_equal_length=False)

    # Correct sweep amplitudes & lengths
    for swp in lsq_sweeps.sweeps:
        swp._i = swp._i * 1e-12
        swp._v = swp._v * 1e-3

        min_len = min(len(swp._i), len(swp._v))
        swp._i = swp._i[:min_len]
        swp._v = swp._v[:min_len]

    # get timing from first sweep
    swp = lsq_sweeps.sweeps[0]

    # Determine the sweep amplitude
    min_i = np.min(swp.i)
    max_i = np.max(swp.i)
    if min_i < 0:
        stim_amp = min_i
    else:
        stim_amp = max_i

    # Determine the start and end of the sweep
    at_amp_ind = np.flatnonzero(swp.i == stim_amp)
    start_ind = at_amp_ind[0]
    lsq_start = swp.t[at_amp_ind[0]]
    lsq_end = swp.t[at_amp_ind[-1]]

    lsq_spx, lsq_spfx = dsf.extractors_for_sweeps(
        lsq_sweeps,
        start=lsq_start,
        end=lsq_end,
        min_peak=-25,
        **dsf.detection_parameters(data_set.LONG_SQUARE)
    )

    lsq_an = spa.LongSquareAnalysis(lsq_spx, lsq_spfx,
        subthresh_min_amp=subthresh_min_amp)
    lsq_features = lsq_an.analyze(lsq_sweeps)

    subthresh_numbers = lsq_features["subthreshold_membrane_property_sweeps"].index.tolist()

    rheo_number = lsq_features["rheobase_sweep"].name
    rheo_amp = lsq_features["rheobase_i"]

    spiking_df = lsq_features["spiking_sweeps"]
    rheo_plus_80_number = spiking_df.loc[spiking_df["stim_amp"] == rheo_amp + 80, :].index.tolist()[0]

    for sn in subthresh_numbers:
        swp = lsq_sweeps.sweeps[sn]
        ax_lsq.plot(swp.t, swp.v, c="black", lw=0.5)

    swp = lsq_sweeps.sweeps[rheo_number]
    ax_lsq.plot(swp.t, swp.v, c="black", lw=0.5)

    spike_start = lsq_features["spikes_set"][rheo_number]["threshold_t"][0]
    spike_end = spike_start + 0.005
    spike_start_ind = tsu.find_time_index(swp.t, spike_start)
    spike_end_ind = tsu.find_time_index(swp.t, spike_end)
    ax_ap.plot(swp.t[spike_start_ind:spike_end_ind] - spike_start, swp.v[spike_start_ind:spike_end_ind])


    swp = lsq_sweeps.sweeps[rheo_plus_80_number]
    ax_lsq.plot(swp.t, swp.v + 120, c="black", lw=0.5)

    ax_ap.set_ylim(-100, 100)
    ax_ap.set(yticks=[], xticks=[])
    ax_lsq.set_ylim(-150, 200)
    ax_lsq.set_xlim(0, 1.0)
    ax_lsq.set(yticks=[], xticks=[])

    sns.despine(ax=ax_ap, left=True, bottom=True)
    sns.despine(ax=ax_lsq, left=True, bottom=True)



def identify_good_sweeps(data_set):
    good_sweep_numbers = []
    for sn in data_set._data.sweep_numbers:
        try:
            d = data_set.get_sweep_data(sn)
            swp = data_set.sweep(sn, enforce_equal_length=False)
            good_sweep_numbers.append(sn)
        except:
            continue
    return good_sweep_numbers


def main(nwb_dir, metadata_file, output_dir, n_per_png, **kwargs):
    metadata_df = pd.read_csv(metadata_file)

    types = metadata_df["cell_type_label"].unique()

    for t in types:
        print(t)

        fig, axes = plt.subplots(n_per_png, 2, gridspec_kw={"width_ratios": (1, 2),  "wspace":0.05}, figsize=(6, 4 * n_per_png))

        # choose random samples
        sub_meta = metadata_df.loc[metadata_df["cell_type_label"] == t, :]
        ax_counter = 0
        for n, r in sub_meta.sample(n=n_per_png, replace=False, random_state=42, axis=0).iterrows():
            print(r["cell_id"])
            ax_ap = axes[ax_counter, 0]
            ax_lsq = axes[ax_counter, 1]

            nwb2_path = os.path.join(nwb_dir, r["nwb_filename"])
            output_filename = os.path.join(output_dir, r["cell_id"] + ".png")
            data_set = create_ephys_data_set(nwb_file=nwb2_path, load_into_memory=False)
            good_sweep_numbers = identify_good_sweeps(data_set)
            if len(good_sweep_numbers) == 0:
                print(f"Could not find loadable sweeps for cell {r['cell_id']}")
                continue
            ephys_plot(data_set, good_sweep_numbers, ax_ap, ax_lsq)

            ax_counter += 1
        plt.savefig(os.path.join(output_dir, t + ".png"))

#     for n, r in metadata_df.iterrows():
#         print(r["cell_id"])
#
#         nwb2_path = os.path.join(nwb_dir, r["nwb_filename"])
#         output_filename = os.path.join(output_dir, r["cell_id"] + ".png")
#         data_set = create_ephys_data_set(nwb_file=nwb2_path, load_into_memory=False)
#         good_sweep_numbers = identify_good_sweeps(data_set)
#         if len(good_sweep_numbers) == 0:
#             print(f"Could not find loadable sweeps for cell {r['cell_id']}")
#             continue
#         save_ephys_plot(data_set, good_sweep_numbers, output_filename)


if __name__ == "__main__":
    module = ags.ArgSchemaParser(schema_type=MakeEphysPngsParameters)

    main(**module.args)
