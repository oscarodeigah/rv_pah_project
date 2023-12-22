import os

import h5py
import json

def save_dict_to_hdf5(file, group, dictionary):
    for key, value in dictionary.items():
        if isinstance(value, dict):
            subgroup = group.create_group(key)
            save_dict_to_hdf5(file, subgroup, value)
        else:
            group[key] = value


def load_hdf5_to_dict(file):
    data_dict = {}
    for key in file.keys():
        if isinstance(file[key], h5py.Group):
            data_dict[key] = load_hdf5_to_dict(file[key])
        elif isinstance(file[key], h5py.Dataset):
            data_dict[key] = file[key][()]
    return data_dict


def gather_pv_data(patients, filename):

    PV_dict = {}
    for patient in patients:
        with h5py.File(f'{patient}/all_results_{patient}.h5', 'r') as hf:
            res = load_hdf5_to_dict(hf)

        rv_pressure = res[patient]['measured_pressure']['2']
        measured_rv_volume = res[patient]['measured_volume']['2']
        simulated_rv_volume = res[patient]['simulated_volume']['2']

        PV_dict[patient.split("_")[0]] = {'rv_pressure': rv_pressure.tolist(), 'measured_volume_rv': measured_rv_volume.tolist(), 'simulated_volume_rv': simulated_rv_volume.tolist()}

    with open(f'{filename}.json', 'w') as f:
        json.dump(PV_dict, f, indent=4)


def gather_active_parameter_results(patients, filename):
    active_dict = {}
    for patient in patients:

        with h5py.File(f'{patient}/all_results_{patient}.h5', 'r') as hf:
            res = load_hdf5_to_dict(hf)

        active_dict[patient.split("_")[0]] = res[patient]['gamma']['2'][1:].tolist() # exclude the unloading point

    with open(f'{filename}.json', 'w') as f:
        json.dump(active_dict, f, indent=4)


def gather_cauchy_stress_results(patients, filename):

    stress_dict = {}

    for patient in patients:

        if not os.path.exists(f'{patient}/all_results_{patient}.h5'):
            continue

        with h5py.File(f'{patient}/all_results_{patient}.h5', 'r') as hf:
            res = load_hdf5_to_dict(hf)

        fiber_stress = (res[patient]['cauchy_stress:fiber']['2'][1:].tolist())
        circumferential_stress = (res[patient]['cauchy_stress:circumferential']['2'][1:].tolist())
        longitudinal_stress = (res[patient]['cauchy_stress:longitudinal']['2'][1:].tolist())

        stress_dict[patient.split("_")[0]] = {
            'fiber': fiber_stress,
            'circumferential': circumferential_stress,
            'longitudinal': longitudinal_stress
        }

    with open(f'{filename}.json', 'w') as f:
        json.dump(stress_dict, f, indent=4)


def gather_green_strain_results(patients, filename):

    strain_dict = {}

    for patient in patients:

        with h5py.File(f'{patient}/all_results_{patient}.h5', 'r') as hf:
            res = load_hdf5_to_dict(hf)

        fiber_strain = (res[patient]['green_strain:fiber']['2'][1:].tolist())
        circumferential_strain = (res[patient]['green_strain:circumferential']['2'][1:].tolist())
        longitudinal_strain = (res[patient]['green_strain:longitudinal']['2'][1:].tolist())

        strain_dict[patient.split("_")[0]] = {
            'fiber': fiber_strain,
            'circumferential': circumferential_strain,
            'longitudinal': longitudinal_strain
        }

    with open(f'{filename}.json', 'w') as f:
        json.dump(strain_dict, f, indent=4)


if __name__ == "__main__":

    patients = [
        'CNT',
        'PAHw4',
        'PAHw8',
        'PAHw12'
    ]

    gather_pv_data(patients, 'all_PV_results')
    gather_active_parameter_results(patients, 'all_active_results')
    gather_cauchy_stress_results(patients, 'all_stress_results')
    gather_green_strain_results(patients, 'all_strain_results')