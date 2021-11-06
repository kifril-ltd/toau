import numpy as np
from laptimize import Solver
from scipy.stats import norm
from scipy.optimize import linprog
from scipy.optimize import minimize

class StpOptimizer:

    def __init__(self, problem):
        self.resNumber = int(problem['resources_number'])
        self.prodNumber = int(problem['products_number'])
        self.a = np.reshape(np.array(problem['a']), (self.prodNumber, self.resNumber)).astype(float)
        self.b = np.array(problem['b']).astype(float)
        self.c = np.array(problem['c']).astype(float)
        self.p = np.array(problem['p']).astype(float)
        self.s = np.reshape(np.array(problem['s']), (self.prodNumber, self.resNumber)).astype(float)
        self.t = np.array(problem['t']).astype(float)
        self.dl = np.array(problem['dl']).astype(float)
        self.du = np.array(problem['du']).astype(float)

        self.stp_solution = None
        self.nonstp_solution = None

    def optimizeStp(self, alpha = None):
        objective = self.getObjective()
        if alpha:
            constraints = self.getConstraints(alpha)
        else:
            constraints = self.getConstraints(self.p)
        capacities = self.getCapacities()
        x0 = np.array([0] * self.prodNumber)
        sol = minimize(objective, x0, method='trust-constr', bounds=capacities, constraints=constraints)

        if not alpha:
            self.stp_solution = {}
            self.stp_solution['func_value'] = round(-sol['fun'], 2)
            self.stp_solution['x'] = np.round(sol['x'], 2)

        return {'func_value': round(-sol['fun'], 2), 'x': np.round(sol['x'], 2)}


    # def optimizeStp(self, alpha = None):
    #     objective = self.getObjective()
    #     if alpha:
    #         linear_constraints = self.getLinearConstraints(alpha)
    #     else:
    #         linear_constraints = self.getLinearConstraints(self.p)
    #     nonlinear_constraints = self.getNonLinearConstraints()
    #     capacities = self.getCapacities()
    #
    #     problem = {
    #         'objective': objective,
    #         'capacity' : capacities
    #     }
    #
    #     i = 0
    #     for constraint in linear_constraints + nonlinear_constraints:
    #         problem['constraints_' + str(i + 1)] = constraint
    #         i += 1
    #
    #     solution = Solver(problem, partition_len=0.1).solve()
    #     if not alpha:
    #         self.stp_solution = {}
    #         self.stp_solution['solution'] = []
    #         for key, value in solution[0].items():
    #             if 'x' in key:
    #                 self.stp_solution['solution'].append(value)
    #             if key == 'obj_value':
    #                 self.stp_solution['func_value'] = -value
    #
    #     return solution

    def optimizeNonStp(self):
        bounds = list(zip(self.dl, self.du))
        solution = linprog(-self.c, A_ub=self.a, b_ub=self.b, bounds=bounds)

        self.nonstp_solution = {}
        self.nonstp_solution['x'] = np.round(solution['x'], 2)
        self.nonstp_solution['func_value'] = round(-solution['fun'], 2)
        return self.nonstp_solution

    def calcRIF(self):
        if not self.stp_solution:
            self.optimizeStp()

        if not self.nonstp_solution:
            self.optimizeNonStp()

        e_nonstp = self.nonstp_solution['func_value']
        e_stp = self.stp_solution['func_value']

        return round(abs(e_nonstp - e_stp) * 100 / e_nonstp, 2)

    def calcRICR(self):
        if not self.stp_solution:
            self.optimizeStp()

        ricr = np.array([])
        for i in range(self.prodNumber):
            t_a = norm.ppf(self.p[i])
            ksi = t_a * (sum(self.s[i] ** 2 * self.stp_solution['x'][i]**2) + self.t[i]) ** (1 / 2)
            ricr = np.append(ricr, ((ksi / (sum(self.a[i]) + ksi)) * 100))

        return np.round(ricr, 2)


    def pValueStat(self, step = 0.1):
        solutions = []
        alpha = 0.5
        while alpha < 1:
            alpha_list = [alpha] * self.prodNumber
            solution = self.optimizeStp(alpha_list)
            print(f'Решено при alpha={alpha}')
            solutions.append(solution)
            alpha += step

        return solutions


    def getObjective(self):
        return lambda x: sum(-self.c * x)

    # def getObjective(self):
    #     objectiveDict = {}
    #     for i in range(self.prodNumber):
    #         objectiveDict['x' + str(i+1)] = lambda x: -self.c[i] * x
    #
    #     return objectiveDict

    def constraintsFunc(self, x, args):
        i = args[0]
        a = args[1]
        s = args[2]
        b = args[3]
        t = args[4]
        p = args[5]
        return - sum(a[i] * x) - norm.ppf(p) * (sum(s[i] ** 2 * x ** 2) + t[i] ** 2) ** (1 / 2) + b[i]

    def getConstraints(self, alpha):
        constraints = []

        for i in range(self.prodNumber):
            args = ((i, self.a, self.s, self.b, self.t, alpha,),)
            constraints.append({
                'type': 'ineq',
                'fun': self.constraintsFunc,
                'args': args
            })

        return constraints

    def getLinearConstraints(self, alpha):
        constraints = []
        for i in range(self.prodNumber):
            constraint = {}
            for j in range(self.resNumber):
                constraint['x' + str(j + 1)] = lambda x: self.a[i, j] * x
            t_a = norm.ppf(alpha[i])
            constraint['t' + str(i + 1)] = lambda t: t_a * t
            constraint['value'] = self.b[i]
            constraints.append(constraint)

        return constraints

    def getNonLinearConstraints(self):
        constraints = []
        for i in range(self.prodNumber):
            constraint = {}
            for j in range(self.resNumber):
                constraint['x' + str(j + 1)] = lambda x: -(self.s[i, j] ** 2) * x ** 2
            constraint['t' + str(i + 1)] = lambda t: t ** 2
            constraint['value'] = self.t[i]**2
            constraints.append(constraint)

        return constraints

    def getCapacities(self):
        return np.array(list(zip(self.dl, self.du)))

    # def getCapacities(self):
    #     capacities = {}
    #
    #     #ограничсения для x
    #     for i in range(self.prodNumber):
    #         capacities['x' + str(i + 1)] = [self.dl[i], self.du[i]]
    #
    #     #ограничения для t
    #     for i in range(self.prodNumber):
    #         tl = (sum(self.s[i] ** 2 * self.dl[i] ** 2) + self.t[i] ** 2) ** (1 / 2)
    #         tu = (sum(self.s[i] ** 2 * self.du[i] ** 2) + self.t[i] ** 2) ** (1 / 2)
    #         capacities['t' + str(i + 1)] = [tl, tu]
    #
    #     return capacities