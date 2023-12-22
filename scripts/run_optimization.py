import pulse
import yaml
import json
import numpy as np
from model_optimization import passive_optimization, active_optimization


def main():

    meshfile = "CNT.h5"
    pv_datafile = "CNT.yml"
    lv_active_datafile = "LV_active_data.json"
    result_file = "results_CNT.json"
    lv_matparam_a = 1.42
    initial_rv_matparam_a = 0.8

    # Load the PV data and format it to be in the required shape
    with open(pv_datafile) as f:
        pv_data = yaml.safe_load(f)

    LVP_ = pv_data['LVP']; LVV_ = pv_data['LVV']
    RVP_ = pv_data['RVP']; RVV_ = pv_data['RVV']

    v_data_ = []
    p_data_ = []
    for i in range(len(LVV_)):
        v_data_.append((LVV_[i], RVV_[i]))
        p_data_.append((LVP_[i], RVP_[i]))

    p_data = np.array(p_data_)
    v_data = np.array(v_data_)

    # Load the mesh geometry (unloaded reference geometry)
    # It is also possible to ignore the "h5group" if that was not specified when creating the h5 file.
    unloaded_geo = pulse.Geometry.from_file(
        h5name=meshfile,
        h5group="-1/unloaded"
    )

    ################### Passive optimization ###################

    opt_rv_matparam_a = passive_optimization(
        unloaded_geo,
        [lv_matparam_a, initial_rv_matparam_a], 
        p_data,
        v_data
    )

    ################### Active optimization ###################

    # Load LV active data
    with open(lv_active_datafile, "r") as f:
        LV_active_ = json.load(f)

    active_control = pulse.RegionalParameter(unloaded_geo.cfun)

    RV_active_ = active_optimization(
        unloaded_geo,
        [lv_matparam_a, opt_rv_matparam_a],
        active_control,
        p_data,
        v_data,
        LV_active_["LV_active_data"]
    )

    # Save passive & active optimization results to file
    final_results = {"passive_parameter": opt_rv_matparam_a, "active_parameter": RV_active_}

    with open(result_file, 'w') as f:
        json.dump(final_results, f, indent=4)


if __name__ == "__main__":
    main()