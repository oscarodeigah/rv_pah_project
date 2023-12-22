import os
import h5py
import dolfin as df
import pulse
from pulse.utils import getLogger
from copy import deepcopy
import numpy as np
import json

from model_optimization import define_material, create_forward_problem
from postprocessing_utils import save_dict_to_hdf5, load_hdf5_to_dict

logger = getLogger(__name__)


def extract_features(forward_problem, geometry, material, lvp, rvp, active_control, target_pressure, target_activation):

    logger.info(f"\nComputing results for Pressure={target_pressure} and Activation={target_activation}")

    pulse.iterate.iterate(forward_problem, (lvp,rvp), target_pressure, initial_number_of_steps=20)
    pulse.iterate.iterate(forward_problem, active_control, target_activation, initial_number_of_steps=150)
    u, p = forward_problem.state.split(deepcopy=True)

    # Simulated volume
    simulated_lv_volume = pulse.dolfin_utils.get_cavity_volume(geometry, chamber='lv', u=u)
    simulated_rv_volume = pulse.dolfin_utils.get_cavity_volume(geometry, chamber='rv', u=u)

    F = pulse.kinematics.DeformationGradient(u)
    W = df.FunctionSpace(geometry.mesh, "DG", 1) 

    # Get the different stress directions
    circ_dir = F * geometry.c0
    long_dir = F * geometry.l0
    rad_dir = F * geometry.r0
    fiber_dir = F * geometry.f0

    # Get a measure over the LV and RV regions
    dx_lv = df.Measure("dx", domain=geometry.mesh, subdomain_data=geometry.cfun, subdomain_id=1)
    dx_rv = df.Measure("dx", domain=geometry.mesh, subdomain_data=geometry.cfun, subdomain_id=2)
    
    # Get the volume of each mesh region
    meshvol_lv = df.assemble(df.Constant(1.0) * dx_lv)
    meshvol_rv = df.assemble(df.Constant(1.0) * dx_rv)

    # Cauchy stress
    T = material.CauchyStress(F, p)
    Tfiber = df.project(df.inner(T * fiber_dir, fiber_dir), W)
    Tcirc = df.project(df.inner(T * circ_dir, circ_dir), W)
    Tlong = df.project(df.inner(T * long_dir, long_dir), W)
    Trad = df.project(df.inner(T * rad_dir, rad_dir), W)

    # Avg RV Cauchy stress
    avg_rv_Tfiber = df.assemble(df.inner(T * fiber_dir, fiber_dir) * dx_rv) / meshvol_rv
    avg_rv_Tcirc = df.assemble(df.inner(T * circ_dir, circ_dir) * dx_rv) / meshvol_rv
    avg_rv_Tlong = df.assemble(df.inner(T * long_dir, long_dir) * dx_rv) / meshvol_rv
    avg_rv_Trad = df.assemble(df.inner(T * rad_dir, rad_dir) * dx_rv) / meshvol_rv

    # Avg LV Cauchy stress
    avg_lv_Tfiber = df.assemble(df.inner(T * fiber_dir, fiber_dir) * dx_lv) / meshvol_lv
    avg_lv_Tcirc = df.assemble(df.inner(T * circ_dir, circ_dir) * dx_lv) / meshvol_lv
    avg_lv_Tlong = df.assemble(df.inner(T * long_dir, long_dir) * dx_lv) / meshvol_lv
    avg_lv_Trad = df.assemble(df.inner(T * rad_dir, rad_dir) * dx_lv) / meshvol_lv

    # Green strain
    E = pulse.kinematics.GreenLagrangeStrain(F)
    Efiber = df.project(df.inner(E * fiber_dir, fiber_dir), W)
    Ecirc = df.project(df.inner(E * circ_dir, circ_dir), W)
    Elong = df.project(df.inner(E * long_dir, long_dir), W)
    Erad = df.project(df.inner(E * rad_dir, rad_dir), W)

    # Avg RV Green strain
    avg_rv_Efiber = df.assemble(df.inner(E * fiber_dir, fiber_dir) * dx_rv) / meshvol_rv
    avg_rv_Ecirc = df.assemble(df.inner(E * circ_dir, circ_dir) * dx_rv) / meshvol_rv
    avg_rv_Elong = df.assemble(df.inner(E * long_dir, long_dir) * dx_rv) / meshvol_rv
    avg_rv_Erad = df.assemble(df.inner(E * rad_dir, rad_dir) * dx_rv) / meshvol_rv

    # Avg LV Green strain
    avg_lv_Efiber = df.assemble(df.inner(E * fiber_dir, fiber_dir) * dx_lv) / meshvol_lv
    avg_lv_Ecirc = df.assemble(df.inner(E * circ_dir, circ_dir) * dx_lv) / meshvol_lv
    avg_lv_Elong = df.assemble(df.inner(E * long_dir, long_dir) * dx_lv) / meshvol_lv
    avg_lv_Erad = df.assemble(df.inner(E * rad_dir, rad_dir) * dx_lv) / meshvol_lv
    
    cauchy_fiber = [avg_lv_Tfiber, avg_rv_Tfiber, df.Vector(Tfiber.vector())]
    cauchy_circumferential = [avg_lv_Tcirc, avg_rv_Tcirc, df.Vector(Tcirc.vector())]
    cauchy_longitudinal = [avg_lv_Tlong, avg_rv_Tlong, df.Vector(Tlong.vector())]
    cauchy_radial = [avg_lv_Trad, avg_rv_Trad, df.Vector(Trad.vector())]
    green_fiber = [avg_lv_Efiber, avg_rv_Efiber, df.Vector(Efiber.vector())]
    green_circumferential = [avg_lv_Ecirc, avg_rv_Ecirc, df.Vector(Ecirc.vector())]
    green_longitudinal = [avg_lv_Elong, avg_rv_Elong, df.Vector(Elong.vector())]
    green_radial = [avg_lv_Erad, avg_rv_Erad, df.Vector(Erad.vector())]

    return (
        u, simulated_lv_volume, simulated_rv_volume, 
        cauchy_fiber, cauchy_circumferential, cauchy_longitudinal, cauchy_radial,
        green_fiber, green_circumferential, green_longitudinal, green_radial
    )


def compute_time_varying_elastance(
        patient,
        geometry,
        matparam_dict,
        pv_data,
        LV_active_data,
        RV_active_data,
        outdir,
        increment = 0.1,
    ):

    if not os.path.exists(outdir):
        os.makedirs(outdir)

    lvp = df.Constant(0.0, name="lvp")
    rvp = df.Constant(0.0, name="rvp")
    active_control = pulse.RegionalParameter(geometry.cfun)

    material = define_material(
        geometry=geometry,
        matparam_a=matparam_dict[patient],
        active_control=active_control
    )

    forward_problem = create_forward_problem(geometry, material, p_endo_LV=lvp, p_endo_RV=rvp)

    num_points = len(pv_data['LVP'])

    E_t_dict = {"TVE": []} # Dictionary to store the time-varying elastance (TVE) values

    for i in range(num_points):
        logger.info(f"\nComputing elastance for {patient} at timepoint {i+1}/{num_points}")
        logger.info(f"\nInitial target pressure: LV={pv_data['LVP'][i]}; RV={pv_data['RVP'][i]}")
        vs = []
        ps = []
    
        # Iterate to the target pressures (and activation) for LV and RV
        pulse.iterate.iterate(forward_problem, (lvp, rvp), (pv_data['LVP'][i], pv_data['RVP'][i]), initial_number_of_steps=20)
        pulse.iterate.iterate(forward_problem, active_control, (LV_active_data[i], RV_active_data[i]), initial_number_of_steps=50)
        u, p = forward_problem.state.split(deepcopy=True)

        # compute the cavity volume for the RV
        RV_cav_vol = pulse.dolfin_utils.get_cavity_volume(geometry, chamber="rv", u=u)
        logger.info(f"\nCavity volume at initial pressure {pv_data['RVP'][i]} kPa is {RV_cav_vol} uL")

        # append to vs & ps lists
        vs.append(RV_cav_vol)
        ps.append(float(rvp))

        logger.info(f"\nPerturbed target pressure: RV={float(rvp)+increment}")

        # perturb RV pressure
        pulse.iterate.iterate(forward_problem, rvp, float(rvp)+increment)
        u, p = forward_problem.state.split(deepcopy=True)

        # compute the new RV cavity volume
        new_RV_cav_vol = pulse.dolfin_utils.get_cavity_volume(geometry, chamber="rv", u=u)
        logger.info(f"\nCavity volume at perturbed pressure {float(rvp)} kPa is {new_RV_cav_vol} uL")

        # append to vs & ps lists
        vs.append(new_RV_cav_vol)
        ps.append(float(rvp))

        elastance = np.mean(np.divide(np.diff(ps), np.diff(vs)))
        E_t_dict["TVE"].append(elastance)

        logger.info(f"\nRV Elastance at timepoint {i+1}/{num_points}:{patient} = {elastance} kPa/uL")

    # Save time-varying elastance to file
    fname = "time_varying_elastance_{}.json"
    json_path = "/".join([outdir, fname])

    with open(json_path.format(patient), 'w') as f:
        json.dump(E_t_dict, f, indent=4)
    


def make_spatial_plots(patient, geometry, results, features, outdir):

    if not os.path.exists(outdir):
        os.makedirs(outdir)

    mesh = geometry.mesh
    
    # Create Function Spaces
    gamma_space = df.FunctionSpace(mesh, "DG", 0)
    marker_space = df.FunctionSpace(mesh, "DG", 0)
    mat_space = df.FunctionSpace(mesh, "DG", 0)
    displacement_space = df.VectorFunctionSpace(mesh, "CG", 2)
    #cg1_space = df.FunctionSpace(mesh, "CG", 1)
    dg1_space = df.FunctionSpace(mesh, "DG", 1)

    times = range(len(results[patient]['gamma']['1']))

    # Add regions to a function
    zones = df.Function(marker_space, name="Freewall-zones")
    zones.vector()[:] = geometry.cfun.array()

    # Add material parameter to a function
    rmat = pulse.RegionalParameter(geometry.cfun)
    matfunc = df.Function(mat_space, name="material_parameter_a")
    matvec = []
    for r in set(geometry.cfun.array()):
        matvec.append(results[patient]['material_parameter_a'][str(r)][0])
    rmat.vector()[:] = matvec
    m = df.project(rmat._sum(), mat_space)
    matfunc.vector()[:] = m.vector()

    functions = {}
    
    for f in features:
        if f in ['material_parameter_a', 'measured_pressure', 'measured_volume', 'simulated_volume']:
            pass
        elif f == "gamma":
            functions[f] = df.Function(gamma_space, name="gamma")
        else:
            functions[f] = df.Function(dg1_space, name=f)

    u = df.Function(displacement_space)

    fname = "simulation_{}.xdmf"
    xdmf_path = "/".join([outdir, fname])

    xdmf = df.XDMFFile(df.MPI.comm_world, xdmf_path.format(patient))

    print("Time")
    for i, t in enumerate(times):
        print(f"{t}/{times[-1]}")

        # Write the material parameter and AHA-zones to the xdmf file
        xdmf.write_checkpoint(matfunc, "material_parameter_a", t, df.XDMFFile.Encoding.HDF5, True,)
        xdmf.write_checkpoint(zones, "Freewall-zones", t, df.XDMFFile.Encoding.HDF5, True,)
        
        u.vector()[:] = results[patient]['displacement'][str(t)][0]
        xdmf.write_checkpoint(u, "displacement", t, df.XDMFFile.Encoding.HDF5, True,)

        for f in functions.keys():

            if f == "gamma":
                pass

            else:
                functions[f].vector()[:] = results[patient][f]['global'][t]   
                xdmf.write_checkpoint(functions[f], f, t, df.XDMFFile.Encoding.HDF5, True,)

    xdmf.close()


def postprocess_simulation(
        patient,
        unloaded_geometry,
        features,
        matparam_dict,
        pv_data,
        LV_active_data,
        RV_active_data
    ):
        
    regions = set(unloaded_geometry.cfun.array())

    lvp = df.Constant(0.0, name="lvp")
    rvp = df.Constant(0.0, name="rvp")
    active_control = pulse.RegionalParameter(unloaded_geometry.cfun)

    # Add initial point for the unloaded volume at pressure = 0. 
    # This is a hack, fix this later. 
    LVP_ = [0.0] + pv_data['LVP']
    LVV_ = [pulse.dolfin_utils.get_cavity_volume(unloaded_geometry, chamber='lv')] + pv_data['LVV']
    RVP_ = [0.0] + pv_data['RVP']
    RVV_ = [pulse.dolfin_utils.get_cavity_volume(unloaded_geometry, chamber='rv')] + pv_data['RVV']

    material = define_material(
        geometry=unloaded_geometry,
        matparam_a=matparam_dict[patient],
        active_control=active_control
    )

    forward_problem = create_forward_problem(unloaded_geometry, material, p_endo_LV=lvp, p_endo_RV=rvp)

    num_points = len(LVP_)

    regional_results = {str(r): [] for r in list(regions)+["global"]}
    results = {f: deepcopy(regional_results) for f in features}
    results['displacement'] = {str(i): [] for i in range(num_points)}

    # Clean up results dictionary (delete 'global' key for features that do not have global values)
    # Maybe find a neater way of making the results dictionary in the future
    for k in results.keys():
        if k in ['material_parameter_a', 'gamma', 'measured_pressure', 'measured_volume', 'simulated_volume']:
            del results[k]['global']

    for i, r in enumerate(regions):
        results['material_parameter_a'][str(r)].append(matparam_dict[patient][i])

    for i, (lvp_, rvp_, lv_active_, rv_active_) in enumerate(zip(LVP_, RVP_, LV_active_data, RV_active_data)):

        logger.info(f"\nPostprocessing results for {patient} at timepoint {i+1}/{num_points}")
        u, sim_lvv, sim_rvv, T_fiber, T_circ, T_long, T_rad, E_fiber, E_circ, E_long, E_rad = extract_features(
            forward_problem, unloaded_geometry, material, lvp, rvp, active_control, (lvp_, rvp_), (lv_active_, rv_active_)
        )

        # Save displacement vector for use in making moving mesh.
        results['displacement'][str(i)].append(df.Vector(u.vector()))

        # Save the optimized and measured parameters for completeness.
        gamma = [lv_active_, rv_active_]
        measured_p = [lvp_, rvp_]
        measured_v = [LVV_[i], RVV_[i]]
        simulated_v = [sim_lvv, sim_rvv]

        for i, r in enumerate(list(regions) + ["global"]):
            results['cauchy_stress:fiber'][str(r)].append(T_fiber[i])
            results['cauchy_stress:circumferential'][str(r)].append(T_circ[i])
            results['cauchy_stress:longitudinal'][str(r)].append(T_long[i])
            results['cauchy_stress:radial'][str(r)].append(T_rad[i])
            results['green_strain:fiber'][str(r)].append(E_fiber[i])
            results['green_strain:circumferential'][str(r)].append(E_circ[i])
            results['green_strain:longitudinal'][str(r)].append(E_long[i])
            results['green_strain:radial'][str(r)].append(E_rad[i])

        # Results that do not have a global value
        for i, r in enumerate(regions):
            results['gamma'][str(r)].append(gamma[i])
            results['measured_pressure'][str(r)].append(measured_p[i])
            results['measured_volume'][str(r)].append(measured_v[i])
            results['simulated_volume'][str(r)].append(simulated_v[i])

    with h5py.File(f'all_results_{patient}.h5', 'w') as hf:
        group = hf.create_group(f'{patient}')
        save_dict_to_hdf5(hf, group, results)

    # Make spatial plots
    with h5py.File(f'all_results_{patient}.h5', 'r') as hf:
        res = load_hdf5_to_dict(hf)
    
    make_spatial_plots(patient, unloaded_geometry, res, features, outdir="spatial_plots")


