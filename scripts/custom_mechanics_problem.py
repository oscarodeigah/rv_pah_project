import dolfin as df
import pulse
from pulse.dolfin_utils import get_cavity_volume


class CustomMechanicsProblem(pulse.MechanicsProblem):

    def _init_solver(self):

        self._problem = df.NonlinearVariationalProblem(
            J=self._jacobian,
            F=self._virtual_work,
            u=self.state,
            bcs=self._dirichlet_bc,
        )
        self.solver = df.NonlinearVariationalSolver(self._problem)
        prm = self.solver.parameters
        prm['newton_solver']['absolute_tolerance'] = 1E-8
        prm['newton_solver']['relative_tolerance'] = 1E-7
        prm['newton_solver']['maximum_iterations']=15
        prm['newton_solver']['relaxation_parameter'] = 1.0
        prm['newton_solver']['linear_solver'] = 'mumps'

    def solve_volumes(self):
        self.solve()
        u, _ = self.state.split(deepcopy=True)
        LVV = get_cavity_volume(self.geometry, chamber="lv",u=u)
        RVV = get_cavity_volume(self.geometry, chamber="rv",u=u)
        
        return LVV, RVV





    
