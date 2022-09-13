import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import seaborn as sns
import argschema as ags
import os

from ipfx.dataset.create import create_ephys_data_set
import ipfx.script_utils as su
import ipfx.data_set_features as dsf
import ipfx.stimulus_protocol_analysis as spa
import ipfx.time_series_utils as tsu


class MakeEphysPngsParameters(ags.ArgSchema):
    metadata_file = ags.fields.InputFile()
    cluster_color_file = ags.fields.InputFile()
    output_dir = ags.fields.OutputDir()
    overwrite = ags.fields.Boolean(default=False)



def ephys_plot(cell_name, data_set, sweep_numbers, trace_color, ax_ap, ax_lsq_up, ax_lsq_down, subthresh_min_amp=-200):
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

    if cell_name != "20180720_sample_2":
        rheo_number = lsq_features["rheobase_sweep"].name
        rheo_amp = lsq_features["rheobase_i"]
    else:
        # hard code to handle problem sweep
        rheo_number = 17
        rheo_amp = 160.0

    spiking_df = lsq_features["spiking_sweeps"]
    rheo_plus_80_number = spiking_df.loc[spiking_df["stim_amp"] == rheo_amp + 80, :].index.tolist()[0]

    for sn in subthresh_numbers:
        swp = lsq_sweeps.sweeps[sn]
        ax_lsq_down.plot(swp.t, swp.v, c=trace_color, lw=0.5)

    swp = lsq_sweeps.sweeps[rheo_number]
    ax_lsq_down.plot(swp.t, swp.v, c=trace_color, lw=0.5)

    spike_start = lsq_features["spikes_set"][rheo_number]["threshold_t"][0]
    spike_end = spike_start + 0.005
    spike_start_ind = tsu.find_time_index(swp.t, spike_start)
    spike_end_ind = tsu.find_time_index(swp.t, spike_end)
    ax_ap.plot(swp.t[spike_start_ind:spike_end_ind] - spike_start, swp.v[spike_start_ind:spike_end_ind], c=trace_color)


    swp = lsq_sweeps.sweeps[rheo_plus_80_number]
    ax_lsq_up.plot(swp.t, swp.v, c=trace_color, lw=0.5)

    ax_ap.set_ylim(-100, 100)
#     ax_ap.set(yticks=[], xticks=[])
    ax_lsq_up.set_ylim(-150, 80)
    ax_lsq_up.set_xlim(0, 1.0)
    ax_lsq_down.set_ylim(-150, 80)
    ax_lsq_down.set_xlim(0, 1.0)
#     ax_lsq.set(yticks=[], xticks=[])

    for ax in (ax_ap, ax_lsq_up, ax_lsq_down):
        ax.set_ylabel("Vm (mV)")

    for ax in (ax_ap, ax_lsq_down):
        ax.set_xlabel("time (s)")

    sns.despine()



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


def main(metadata_file, cluster_color_file, output_dir, overwrite, **kwargs):
    metadata_df = pd.read_csv(metadata_file)
    colors_df = pd.read_csv(cluster_color_file, sep="\t", index_col=0)
    print(colors_df.head())

    types = metadata_df["RNA type"].unique()

    for t in types:
        t_color = colors_df.at[t, "cluster_color"]
        print(t, t_color)


        sub_meta = metadata_df.loc[metadata_df["RNA type"] == t, :]
        for n, r in sub_meta.iterrows():
            print(r["Cell"])
            nwb2_path = r["nwb_filepath"]
            output_filename = os.path.join(output_dir, f'{t.replace("/", "-")}_{r["Cell"]}.svg')
            if not overwrite:
                if os.path.exists(output_filename):
                    continue

            fig = plt.figure(figsize=(6, 4))
            g = gridspec.GridSpec(2, 2, width_ratios=(1, 2), wspace=0.3, hspace=0.3)

            ax_ap = plt.subplot(g[:, 0])
            ax_lsq_up = plt.subplot(g[0, 1])
            ax_lsq_down = plt.subplot(g[1, 1])

            data_set = create_ephys_data_set(nwb_file=nwb2_path, load_into_memory=False)
            good_sweep_numbers = identify_good_sweeps(data_set)
            if len(good_sweep_numbers) == 0:
                print(f"Could not find loadable sweeps for cell {r['Cell']}")
                continue
            ephys_plot(r["Cell"], data_set, good_sweep_numbers, t_color, ax_ap, ax_lsq_up, ax_lsq_down)

            plt.savefig(output_filename)
            plt.close()



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
