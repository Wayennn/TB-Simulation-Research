import CurveFitting as cf
import Simulation as sm

# values important for the graphs
beginYear = 2011        # beginning year of curve-fitting
endYear = 2018          # end year of curve-fitting, beginning year of simulation
endYearSim = 2025       # end year of simulation, unimportant here
scale = 100000          # per *scale* people
ac = 8                  # maximum number of curves graphed

# values important for curve-fitting
iterations = 2          # number of iterating the trial and error
rrr = 17                 # should be odd, narrowness of search window for initial values, wideness of beta search scope
betatest = 20           # beta estimate
newData = True          # True if new dataset is used
bound = 1               # incidence fit: 0 if average, 1 if lower bound, 2 if upper bound

p0, beta = cf.curveFit(beginYear, endYear, endYearSim, scale, ac, iterations, rrr, betatest, newData, bound)
sm.simulate(beginYear, endYear, endYearSim, scale, newData, p0, beta)