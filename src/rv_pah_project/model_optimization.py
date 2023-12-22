from datetime import datetime
import numpy as np
from scipy.optimize import minimize
import dolfin as df
import pulse
from pulse.utils import getLogger
from pulse.dolfin_utils import get_cavity_volume
import warnings

from .custom_mechanics_problem import CustomMechanicsProblem
from ffc.quadrature.deprecation import QuadratureRepresentationDeprecationWarning

warnings.simplefilter("ignore", QuadratureRepresentationDeprecationWarning)

logger = getLogger(__name__)


def cost_func(rv_model, rv_obs):
    return ((rv_model - rv_obs) / rv_obs) ** 2


def define_material(
    geometry, matparam_a, active_control=None, active_model="active_stress"
):
    a = pulse.RegionalParameter(geometry.cfun)
    a.vector()[0] = matparam_a[0]  # LVFW
    a.vector()[1] = matparam_a[1]  # RVFW

    material_parameters = dict(
        a=a,
        a_f=2.582,
        b=5.0,
        b_f=5.0,
        a_s=0.0,
        b_s=0.0,
        a_fs=0.0,
        b_fs=0.0,
    )

    material = pulse.HolzapfelOgden(
        active_model=active_model,
        parameters=material_parameters,
        activation=active_control,
        T_ref=1.0,
        eta=0.2,
        f0=geometry.f0,
        s0=geometry.s0,
        n0=geometry.n0,
    )
    return material


def create_forward_problem(geometry, material, p_endo_LV, p_endo_RV):
    lv_marker = geometry.markers["ENDO_LV"][0]
    rv_marker = geometry.markers["ENDO_RV"][0]
    lv_pressure = pulse.NeumannBC(traction=p_endo_LV, marker=lv_marker, name="lv")
    rv_pressure = pulse.NeumannBC(traction=p_endo_RV, marker=rv_marker, name="rv")
    neumann_bc = [lv_pressure, rv_pressure]

    # Robin BC
    pericardium_spring = 0.5
    robin_bc = [
        pulse.RobinBC(
            value=df.Constant(pericardium_spring), marker=geometry.markers["EPI"][0]
        )
    ]

    # Fix the basal plane in the longitudinal direction
    # 0 in V.sub(0) refers to x-direction, which is the longitudinal direction
    def fix_basal_plane(W):
        V = W if W.sub(0).num_sub_spaces() == 0 else W.sub(0)
        bc = df.DirichletBC(
            V.sub(0),
            df.Constant(0.0, name="fix_base"),
            geometry.ffun,
            geometry.markers["BASE"][0],
        )
        return bc

    dirichlet_bc = [fix_basal_plane]

    # Collect boundary conditions
    bcs = pulse.BoundaryConditions(
        dirichlet=dirichlet_bc, neumann=neumann_bc, robin=robin_bc
    )

    forward_problem = CustomMechanicsProblem(geometry, material, bcs)

    return forward_problem


def passive_optimization(
    unloaded_geometry,
    lv_matparam_a,
    initial_rv_matparam_a,
    pressure_data,
    volume_data,
    no_passive_filling_points=4,
):
    def objective(x):
        logger.info(
            f"Evaluate objective function for passive optimization with x = {x}"
        )
        material = define_material(unloaded_geometry, [lv_matparam_a, x[0]])

        lvp = df.Constant(0.0, name="lvp")
        rvp = df.Constant(0.0, name="rvp")
        forward_problem = create_forward_problem(
            unloaded_geometry, material, p_endo_LV=lvp, p_endo_RV=rvp
        )

        try:
            for target_ps, target_vs in zip(
                pressure_data[:no_passive_filling_points],
                volume_data[:no_passive_filling_points],
            ):
                pulse.iterate.iterate(
                    forward_problem, (lvp, rvp), target_ps, initial_number_of_steps=20
                )
        except pulse.mechanicsproblem.SolverDidNotConverge:
            logger.info(f"\nPassive inflation failed for RV matparam, a = {x[0]} kPa")
            return 1e6
        u, _ = forward_problem.state.split(deepcopy=True)
        LVV = pulse.dolfin_utils.get_cavity_volume(unloaded_geometry, chamber="lv", u=u)
        RVV = pulse.dolfin_utils.get_cavity_volume(unloaded_geometry, chamber="rv", u=u)
        LVV_obs, RVV_obs = volume_data[no_passive_filling_points - 1]
        f = cost_func(RVV, RVV_obs)

        logger.info(f"\nLVV = {LVV}, obs = {LVV_obs}, RVV = {RVV}, obs = {RVV_obs}")
        logger.info(f"\nCost function: {f}")
        return f

    def print_control(x):
        logger.info(f"\nCurrent x: {x}")

    x0 = np.array([initial_rv_matparam_a])
    bnds = [(0.05, 10)]

    logger.info("\nStarting passive optimization")
    start = datetime.now()
    res = minimize(
        objective, x0=x0, bounds=bnds, options={"disp": 5}, callback=print_control
    )
    print(res)
    logger.info(
        f"\nPassive optimization done. Optimized RV matparam, a = {res.x[0]} kPa"
    )
    end = datetime.now()

    total_run_time = end - start
    mins, secs = divmod(total_run_time.seconds, 60)
    hours, mins = divmod(mins, 60)
    logger.info(
        f"\nPassive optimization run time: {hours} hours, {mins} minutes and {secs} seconds"
    )

    return res.x[0]


def active_optimization(
    unloaded_geometry,
    matparam_a,
    active_control,
    pressure_data,
    volume_data,
    LV_active_data,
    no_passive_filling_points=4,
):
    assert (
        pressure_data.shape == volume_data.shape
    ), "Pressure and volume arrays are of different sizes"
    assert (
        len(LV_active_data) == pressure_data.shape[0]
    ), "LV active data and pressure data should have the same length"

    material = define_material(
        geometry=unloaded_geometry, matparam_a=matparam_a, active_control=active_control
    )

    lvp = df.Constant(0.0, name="lvp")
    rvp = df.Constant(0.0, name="rvp")
    forward_problem = create_forward_problem(
        unloaded_geometry, material, p_endo_LV=lvp, p_endo_RV=rvp
    )

    # Passive inflation to ED pressures
    logger.info(
        f"\nStarting passive inflation to LV EDP = {pressure_data[no_passive_filling_points-1][0]} kPa; RV EDP = {pressure_data[no_passive_filling_points-1][1]} kPa"
    )

    for target_ps in pressure_data[:no_passive_filling_points]:
        pulse.iterate.iterate(
            forward_problem, (lvp, rvp), target_ps, initial_number_of_steps=20
        )

    u, p = forward_problem.state.split(deepcopy=True)
    LV_EDV = pulse.dolfin_utils.get_cavity_volume(unloaded_geometry, chamber="lv", u=u)
    RV_EDV = pulse.dolfin_utils.get_cavity_volume(unloaded_geometry, chamber="rv", u=u)

    logger.info(f"\nPassive inflation done. LV EDV = {LV_EDV} uL; RV EDV = {RV_EDV} uL")

    logger.info("\nStarting active optimization")

    RV_active_values = [
        0.0
    ] * no_passive_filling_points  # Fill RV_active_values list with zeros for the passive filling points.
    x0 = np.array([0.0])
    overall_start = datetime.now()
    for i, (target_ps, target_vs, target_LV_active) in enumerate(
        zip(
            pressure_data[no_passive_filling_points:],
            volume_data[no_passive_filling_points:],
            LV_active_data[no_passive_filling_points:],
        )
    ):
        logger.info(
            f"\nOptimizing pressure step, pressures = {target_ps}, volumes = {target_vs}"
        )
        logger.info(f"\nLV active parameter {target_LV_active}")
        pulse.iterate.iterate(
            forward_problem, (lvp, rvp), target_ps, initial_number_of_steps=20
        )

        def objective(x):
            pulse.iterate.iterate(
                forward_problem,
                active_control,
                (target_LV_active + 1e-4, x[0] + 1e-4),
                initial_number_of_steps=150,
            )

            u, _ = forward_problem.state.split(deepcopy=True)
            LVV = pulse.dolfin_utils.get_cavity_volume(
                unloaded_geometry, chamber="lv", u=u
            )
            RVV = pulse.dolfin_utils.get_cavity_volume(
                unloaded_geometry, chamber="rv", u=u
            )
            LVV_obs, RVV_obs = target_vs
            f = cost_func(RVV, RVV_obs)
            logger.info(f"\nLVV = {LVV}, obs = {LVV_obs}, RVV = {RVV}, obs = {RVV_obs}")
            logger.info(f"\nCost function: {f}")
            return f

        def jac(x, eps=1e-4):
            u, _ = forward_problem.state.split(deepcopy=True)
            LVV0 = pulse.dolfin_utils.get_cavity_volume(
                unloaded_geometry, chamber="lv", u=u
            )
            RVV0 = pulse.dolfin_utils.get_cavity_volume(
                unloaded_geometry, chamber="rv", u=u
            )
            _, RVV_obs = target_vs
            active_control.vector()[:] += eps
            _, RVV1 = forward_problem.solve_volumes()
            active_control.vector()[:] -= eps
            j = [
                2 * ((RVV0 - RVV_obs) / RVV_obs) * (1 / RVV_obs) * ((RVV1 - RVV0) / eps)
            ]
            logger.info(f"\nComputed jacobian: {j}")
            return j

        bnds = [(0, 1000.0)]

        start = datetime.now()
        res = minimize(objective, x0=x0, jac=jac, bounds=bnds, options={"disp": 5})
        print(res)
        end = datetime.now()

        total_run_time = end - start
        mins, secs = divmod(total_run_time.seconds, 60)
        hours, mins = divmod(mins, 60)
        logger.info(f"\nRun_time: {hours} hours, {mins} minutes and {secs} seconds")

        RV_active_values.append(res.x[0])
        x0 = np.array([res.x])

    logger.info(f"\nActive optimization done.")
    overall_end = datetime.now()
    total_run_time = overall_end - overall_start
    mins, secs = divmod(total_run_time.seconds, 60)
    hours, mins = divmod(mins, 60)
    logger.info(
        f"\nActive optimization run time: {hours} hours, {mins} minutes and {secs} seconds"
    )

    return RV_active_values
