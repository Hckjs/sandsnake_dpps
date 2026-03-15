import matplotlib
import matplotlib.pyplot as plt
import yaml
import joblib

from astropy import units as u
from astropy.table import vstack
import numpy as np
import ctaplot


## Origin Reconstruction ##
def direction_results(data, axes):
    # subplots(2, 2)
    ctaplot.plot_theta2(
        (data["true_alt"].value * u.deg).to(u.rad),
        (data["disp_alt"].value * u.deg).to(u.rad),
        (data["true_az"].value * u.deg).to(u.rad),
        (data["disp_az"].value * u.deg).to(u.rad),
        ax=axes[0, 0],
        bins=50,
        range=(0, 1),
    )
    axes[0, 0].grid()
    axes[0, 0].set_title(r"$\Theta^{2}$ Distribution")

    ctaplot.plot_angular_resolution_per_energy(
        (data["true_alt"].value * u.deg).to(u.rad),
        (data["disp_alt"].value * u.deg).to(u.rad),
        (data["true_az"].value * u.deg).to(u.rad),
        (data["disp_az"].value * u.deg).to(u.rad),
        data["RandomForestRegressor_energy"].value * u.TeV,
        ax=axes[0, 1],
    )

    ctaplot.plot_angular_resolution_cta_requirement(
        "north", ax=axes[0, 1], color="black"
    )
    axes[0, 1].grid()
    axes[0, 1].legend()

    ctaplot.plot_migration_matrix(
        (data["true_alt"].value * u.deg).to(u.rad),
        (data["disp_alt"].value * u.deg).to(u.rad),
        ax=axes[1, 0],
        colorbar=True,
        xy_line=True,
        hist2d_args=dict(norm=matplotlib.colors.LogNorm()),
        line_args=dict(color="black"),
    )
    axes[1, 0].set_xlabel(r"$Altitude_{true}$")
    axes[1, 0].set_ylabel(r"$Altitude_{reco}$")

    true_az_mod = (data["true_az"].value + 180) % 360 - 180
    az_mod = (data["disp_az"].value + 180) % 360 - 180
    ctaplot.plot_migration_matrix(
        (true_az_mod * u.deg).to(u.rad),
        (az_mod * u.deg).to(u.rad),
        ax=axes[1, 1],
        colorbar=True,
        xy_line=True,
        hist2d_args=dict(norm=matplotlib.colors.LogNorm()),
        line_args=dict(color="black"),
    )
    axes[1, 1].set_xlabel(r"$Azimuth_{true}$")
    axes[1, 1].set_ylabel(r"$Azimuth_{reco}$")


def roc_curve_sign():
    return None


def plot_alt_az_reco(data, axes):
    # just copied from lstchains lst_sensitivity_pyirf notebook - needs to be adapted
    fig, axs = plt.subplots(1, 5, figsize=(30, 4))
    emin_bins = [0.0, 0.1, 0.5, 1, 5] * u.TeV
    emax_bins = [0.1, 0.5, 1, 5, 10] * u.TeV

    for i, ax in enumerate(axs):
        events = data[
            (data["reco_energy"] > emin_bins[i]) & (data["reco_energy"] < emax_bins[i])
        ]
        pcm = ax.hist2d(
            events["reco_az"].to_value(u.deg),
            events["reco_alt"].to_value(u.deg),
            bins=50,
        )
        ax.title.set_text(
            "%.1f-%.1f TeV" % (emin_bins[i].to_value(), emax_bins[i].to_value())
        )
        ax.set_xlabel("Az (º)")
        ax.set_ylabel("Alt (º)")
        fig.colorbar(pcm[3], ax=ax)
    plt.show()


## Energy Reconstruction ##
def plot_energy_results(data, axes):
    # subplots(2, 2)
    ctaplot.plot_energy_resolution(
        data["RandomForestRegressor_energy"].value * u.TeV,
        data["true_energy"].value * u.TeV,
        ax=axes[0, 0],
        bias_correction=False,
    )
    ctaplot.plot_energy_resolution_cta_requirement(
        "north", ax=axes[0, 0], color="black"
    )

    ctaplot.plot_energy_bias(
        data["RandomForestRegressor_energy"].value * u.TeV,
        data["true_energy"].value * u.TeV,
        ax=axes[1, 0],
    )
    ctaplot.plot_migration_matrix(
        np.log10(data["RandomForestRegressor_energy"].value) * u.TeV,
        np.log10(data["true_energy"].value) * u.TeV,
        ax=axes[0, 1],
        colorbar=True,
        xy_line=True,
        hist2d_args=dict(norm=matplotlib.colors.LogNorm()),
        line_args=dict(color="black"),
    )
    axes[0, 0].legend()
    axes[0, 0].set_title("")
    axes[0, 0].label_outer()
    axes[0, 1].set_xlabel(r"$log(E_{true}/TeV)$")
    axes[0, 1].set_ylabel(r"$log(E_{reco}/TeV)$")
    axes[1, 0].set_title("")
    axes[1, 0].set_ylabel("Energy bias")
    for ax in axes.ravel():
        ax.grid(True, which="both")
    axes[1, 1].remove()


## Gamma/Hadron Separation ##
def plot_roc_gamma(
    gamma_cl, gamma_true_e, proton_cl, proton_true_e, ax, energy_bins=None, **kwargs
):
    ax[1, 0].hist(
        gamma_cl,
        bins=100,
        histtype="step",
        density=True,
        label="Gammas",
    )
    ax[1, 1].hist(
        proton_cl,
        bins=100,
        histtype="step",
        density=True,
        label="Protons",
    )
    merged_cl_pred = vstack([gamma_cl, proton_cl])
    del proton_cl
    merged_true_e = vstack([gamma_true_e, proton_true_e])
    del gamma_true_e
    del proton_true_e
    true_labels = np.zeros(len(merged_cl_pred))
    true_labels[len(gamma_cl) :] = 1
    del gamma_cl

    if energy_bins is None:
        ctaplot.plot_roc_curve_gammaness(
            true_labels,
            merged_cl_pred["col0"].data,
            ax=ax[0, 0],
            **kwargs,
        )
    else:
        ctaplot.plot_roc_curve_gammaness_per_energy(
            true_labels,
            merged_cl_pred["col0"].data,
            merged_true_e["col0"].data,
            energy_bins=energy_bins,
            ax=ax[0, 0],
            **kwargs,
        )
    ax[0, 1].remove()
    del true_labels
    del merged_cl_pred
    del merged_true_e


## Miscellaneous ##
def plot_importances(model, features_names, ax, n_feat=15, **kwargs):
    importances = model.feature_importances_
    std = np.std([tree.feature_importances_ for tree in model.estimators_], axis=0)
    indices = np.argsort(importances)

    ordered_features = [features_names[index] for index in indices]

    n_feat = min(n_feat, len(ordered_features))

    selected_features = ordered_features[-n_feat:]
    selected_importances = importances[indices][-n_feat:]
    selected_std = std[indices][-n_feat:]

    ax.barh(
        selected_features,
        selected_importances,
        xerr=selected_std,
        align="center",
        **kwargs,
    )
    ax.set_yticks(range(n_feat))
    ax.set_yticklabels(selected_features)
    ax.grid()


def plot_models_feature_importances(
    e_reg_path, disp_path, p_clf_path, config_path, **kwargs
):
    with open(config_path, "r", encoding="utf-8") as file:
        config_file = yaml.safe_load(file)

    with open(e_reg_path, "rb") as file:
        e_reg_model = joblib.load(file)
    with open(disp_path, "rb") as file:
        disp_model = joblib.load(file)
    with open(p_clf_path, "rb") as file:
        p_clf_model = joblib.load(file)

    features = {
        "e_reg": config_file["TrainEnergyRegressor"]["EnergyRegressor"]["features"],
        "disp": config_file["TrainDispReconstructor"]["DispReconstructor"]["features"],
        "p_clf": config_file["TrainParticleClassifier"]["ParticleClassifier"][
            "features"
        ],
    }

    fig, axes = plt.subplots(4, 3, figsize=(18, 20))
    axes = axes.flatten()

    model_groups = [
        ("e_reg", e_reg_model._models.values(), "Energy Regression"),
        ("disp_norm", [pair[0] for pair in disp_model._models.values()], "DISP Norm"),
        ("disp_sign", [pair[1] for pair in disp_model._models.values()], "DISP Sign"),
        ("p_clf", p_clf_model._models.values(), "Particle Clf"),
    ]

    for row_idx, (key, models, title_prefix) in enumerate(model_groups):
        for i, model in enumerate(models):
            ax_idx = row_idx * 3 + i
            ax = axes[ax_idx]

            feat_key = "disp" if key.startswith("disp") else key
            plot_importances(model, features[feat_key], ax=ax, **kwargs)
            ax.set_title(f"{title_prefix} - TelType {i}")

    fig.suptitle("Feature Importances of All RF Models", fontsize=18)
    fig.tight_layout(rect=[0, 0.03, 1, 0.95])
    return fig


def plot_features(data_gammas, data_protons, axes):
    # Energy distribution
    axes[0].hist(
        np.log10(data_gammas["true_energy"]),
        histtype="step",
        bins=100,
        label="Gammas",
    )
    axes[0].hist(
        np.log10(data_protons["true_energy"]),
        histtype="step",
        bins=100,
        label="Protons",
    )
    axes[0].set_ylabel(r"# of events", fontsize=15)
    axes[0].set_xlabel(r"$log_{10}E$(GeV)")
    axes[0].legend()

    # disp_ distribution
    # TODO: per telescope?
    axes[1].hist(
        data_gammas["disp_tel_parameter"],
        histtype="step",
        bins=100,
        label="Gammas",
    )
    axes[1].hist(
        data_protons["disp_tel_parameter"],
        histtype="step",
        bins=100,
        label="Gammas",
    )
    axes[1].set_ylabel(r"# of events", fontsize=15)
    axes[1].set_xlabel(r"disp_ (m)")
