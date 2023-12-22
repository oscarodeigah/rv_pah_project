import pulse
import yaml
import json

from rv_pah_project.postprocess_model import (
    postprocess_simulation,
    compute_time_varying_elastance,
)


def main():
    animal = "CNT"
    meshfile = "CNT.h5"
    pv_datafile = "CNT.yml"
    lv_active_datafile = "LV_active_data.json"
    lv_matparam_a = 1.42
    result_file = "results_CNT.json"

    # Load the results from the model optimization
    with open(result_file, "r") as f:
        RV_results = json.load(f)

    features = [
        "cauchy_stress:fiber",
        "cauchy_stress:circumferential",
        "cauchy_stress:longitudinal",
        "cauchy_stress:radial",
        "green_strain:fiber",
        "green_strain:circumferential",
        "green_strain:longitudinal",
        "green_strain:radial",
        "gamma",
        "material_parameter_a",
        "measured_pressure",
        "measured_volume",
        "simulated_volume",
    ]

    matparam_dict = {animal: [lv_matparam_a, RV_results["passive_parameter"]]}

    # Load the mesh geometry (unloaded reference geometry)
    unloaded_geo = pulse.Geometry.from_file(h5name=meshfile, h5group="-1/unloaded")

    # Load the PV data and format it to be in the required shape
    with open(pv_datafile) as f:
        pv_data = yaml.safe_load(f)

    # Load active data (add initial point for the unloaded volume at pressure = 0)
    with open(lv_active_datafile, "r") as f:
        LV_active_ = json.load(f)

    LV_active_data = [0.0] + LV_active_["LV_active_data"]
    RV_active_data = [0.0] + RV_results["active_parameter"]

    # Postprocess the model optimization results
    postprocess_simulation(
        animal,
        unloaded_geo,
        features,
        matparam_dict,
        pv_data,
        LV_active_data,
        RV_active_data,
    )

    # Time-varying elastance is not computed for the unloaded (reference) geometry
    # hence, the initial points added to LV_active_data & RV_active_data are ignored.
    compute_time_varying_elastance(
        animal,
        unloaded_geo,
        matparam_dict,
        pv_data,
        LV_active_data[1:],
        RV_active_data[1:],
        outdir="elastance_results",
    )


if __name__ == "__main__":
    main()
