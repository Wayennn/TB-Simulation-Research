import matplotlib.pyplot as plt
# from matplotlib.ticker import (MultipleLocator)
import numpy as np
import csv
from scipy.integrate import odeint

def curveFit(beginYear, endYear, endYearSim, scale, ac, iterations, rrr, betatest, newData, bound):
    
    def f(p0, t, consta, constb, beta):
        p = p0
        d_i = consta[4]-consta[8]
        d_t = consta[5]-consta[8]
        c1 = beta*consta[9]*(p[2] + consta[6]*p[3])
        c2 = p[6] - consta[8] - d_i*p[2] - d_t*p[3]
        dLAdt = c1*(p[4]+constb[0]*(p[5]+p[1])) - p[0]*(consta[0]+consta[1]+consta[8]+c2)
        dLBdt = consta[1]*p[0] + consta[3]*p[2] - p[1]*(constb[0]*c1+consta[2]+consta[8]+c2)
        dIdt = consta[0]*p[0] + consta[2]*p[1] + consta[10]*p[3] - p[2]*(consta[3]+constb[4]+consta[4]+c2)
        dTdt = constb[4]*p[2] - p[3]*(constb[3]*constb[1]+consta[10]+consta[5]+c2)
        dSAdt = (1-constb[2])*p[6] - p[4]*(c1+consta[8]+c2)
        dSBdt = constb[2]*p[6] + constb[3]*constb[1]*p[3] - p[5]*(constb[0]*c1+consta[8]+c2)
        dNdt = -p[6]*c2
        return np.array([dLAdt, dLBdt, dIdt, dTdt, dSAdt, dSBdt, dNdt])
    
    def ploti(_ax, yr, d=False):
        _ax.set_xlabel("year")
        _ax.set_ylabel("population proportion per " + str(scale) + " persons")
        _ax.set_xticks(np.arange(yr, yr+endYear-beginYear+1, 1))
        _ax.set_ylim(0, 1000)
        _ax.set_xlim(beginYear, endYear)
        # if yr==beginYear:
        #     _ax.set_yticks(np.arange(0, 0.00401*scale, 0.0005*scale))
        #     _ax.yaxis.set_minor_locator(MultipleLocator(0.00025*scale))
        # elif yr==endYear:
        #     _ax.set_yticks(np.arange(0.00140*scale, 0.00281*scale, 0.0002*scale))
        #     _ax.yaxis.set_minor_locator(MultipleLocator(0.0001*scale))
        # if d:
        #     _ax.set_yticks(np.arange(0, 0.00251*scale, 0.0005*scale))
        #     _ax.yaxis.set_minor_locator(MultipleLocator(0.0001*scale))
        _ax.legend()
    
    def plotp(_ax, yr):
        _ax.set_xlabel("year")
        _ax.set_ylabel("population proportion per " + str(scale) + " persons")
        _ax.set_xticks(np.arange(yr, yr+endYear-beginYear+1, 1))
        _ax.set_ylim(0, 1000)
        _ax.set_xlim(beginYear, endYear)
        # if yr==beginYear:
        #     _ax.set_yticks(np.arange(0, 0.00901*scale, 0.0010*scale))
        #     _ax.yaxis.set_minor_locator(MultipleLocator(0.0005*scale))
        # elif yr==endYear:
        #     _ax.set_yticks(np.arange(0.00250*scale, 0.00501*scale, 0.0005*scale))
        #     _ax.yaxis.set_minor_locator(MultipleLocator(0.0001*scale))
        _ax.legend()
    
    def fnd(_y1, par, delt):
        if par=="inci":
            return _y1[:, 3]/delt
        elif par=="preva":
            return _y1[:, 3]+_y1[:, 2]
        
    def fit(ct, rrr, ac, beta, p0, ax):
        bet = beta
        b = beta/np.floor(rrr/2)
        decimal = 0.1
        for count in range(ct):
            rsum = np.array([])
            scandeg = 9.95073125*decimal/(rrr-1)
            beta -= b*np.ceil(rrr/2)
            p0[1] -= scandeg*np.ceil(rrr/2)
            p0[5] += scandeg*np.ceil(rrr/2)
            p0[0] -= scandeg*np.ceil(rrr/2)
            p0[4] += scandeg*np.ceil(rrr/2)
            for l in range(rrr):
                p0[1] += scandeg
                p0[5] -= scandeg
                if p0[1]<=scandeg/10 or p0[5]<=scandeg/10:
                    for i in range(rrr**2):
                        rsum = np.append(rsum, 1)
                    continue
                for k in range(rrr):
                    p0[0] += scandeg
                    p0[4] -= scandeg
                    if p0[0]<=scandeg/10 or p0[4]<=scandeg/10:
                        for i in range(rrr):
                            rsum = np.append(rsum, 1)
                        continue
                    for j in range(rrr):
                        beta += b
                        y = odeint(f, p0, t, args=(consta, constb, beta))
                        rs = np.array([])
                        for i in range(len(incidence[bound])):
                            rs = np.append(rs, ((incidence[bound, i]-(y[i, 3]/constb[4])))**2)
                        for i in range(len(prevalence[1])):
                            rs = np.append(rs, ((prevalence[1, i]-(y[i, 3]+y[i, 2]))**2))
                        if count==ct-1 and k%(np.ceil(rrr/ac))==0 and j%(np.ceil(rrr/ac))==0 and l%(np.ceil(rrr/ac))==0:
                            ax[0].plot(t, (y[:, 3]/constb[4])*scale, linewidth=0.1, linestyle=":")
                            ax[1].plot(t, (y[:, 3]+y[:, 2])*scale, linewidth=0.1, linestyle=":")
                        rsum = np.append(rsum, np.sum(rs))
                    beta -= b*(rrr)   
                p0[0] -= scandeg*(rrr)
                p0[4] += scandeg*(rrr)
            p0[1] -= scandeg*(rrr)
            p0[5] += scandeg*(rrr)
            
            beta = beta + b
            p0[0] = p0[0] + scandeg
            p0[4] = p0[4] - scandeg
            p0[1] = p0[1] + scandeg
            p0[5] = p0[5] - scandeg
            
            beta = b*(np.argmin(rsum)%rrr) + beta
            p0[0] = scandeg*np.floor((np.argmin(rsum)%(rrr**2))/rrr) + p0[0]
            p0[4] = -scandeg*np.floor((np.argmin(rsum)%(rrr**2))/rrr) + p0[4]
            p0[1] = scandeg*np.floor((np.argmin(rsum))/(rrr**2)) + p0[1]
            p0[5] = -scandeg*np.floor((np.argmin(rsum))/(rrr**2)) + p0[5]
            if count==ct-1:
                print("iterations = " + str(ct) + "\nrrr = " + str(rrr) + "\ninitial beta = " + str(bet) + "\n-------------------------\n(p at 2011)")
                for i in range(len(p0)):
                    print("p0[" + str(i) + "] = " + str(p0[i]))
                print("beta = " + str(beta))
                y = odeint(f, p0, t, args=(consta, constb, beta))
                print("RSME = " + str(np.sqrt(np.amin(rsum))))
                print("-------------------------\n(p at 2018)")
                for i in range(len(p0)):
                    print("p[" + str(i) + "] = " + str(y[-1, i]))
            if beta==bet*2:
                raise Exception("WARNING! Parameter beta may be larger than estimate.")
            b = b/np.floor(rrr/2)
            decimal = decimal/(rrr-1)
    
        ax[0].set_title("Tuberculosis Incidence Rate in the Philippines")
        ax[1].set_title("Tuberculosis Prevalence Rate in the Philippines")
        ax[0].plot(t, (y[:, 3]/constb[4])*scale, linewidth=2, label="Incidence rate")
        # ax[0].plot(np.arange(beginYear, endYear+1, 1), incidence[0]*scale, linewidth=2, label="Estimated Incidence rate")
        # ax[0].plot(np.arange(beginYear, endYear+1, 1), incidence[1]*scale, linestyle=":", linewidth=2, label="Exp Inc rate lower bound")
        # ax[0].plot(np.arange(beginYear, endYear+1, 1), incidence[2]*scale, linestyle=":", linewidth=2, label="Exp Inc rate upper bound")
        ax[0].plot(np.arange(beginYear, endYear+1, 1), incidence[1]*scale, linewidth=2, label="Estimated Incidence rate")                   #commend this line and uncomment the previous three to reveal the truth
        ax[1].plot(t, (y[:, 3]+y[:, 2])*scale, linewidth=2, label="Prevalence rate")
        ax[1].plot(np.arange(beginYear, endYear+1, 1), prevalence[1]*scale, linewidth=2, label="Estimated Prevalence rate")
        ploti(ax[0], beginYear)
        plotp(ax[1], beginYear)
        
        return beta
    
    consta = []
    constb = []
    p0 = []
    incidence = []
    prevalence = []
    
    if newData:
        a = 3
    else:
        a = 1
    
    with open("Summary-of-model-parameters.csv", "r") as csvfile:
        csvreader = csv.reader(csvfile, delimiter = ",")
        i = 0
        for row in csvreader:
            if i<11:
                consta.append(float(eval(row[a])))          # change to 'row[1]' for old values
            elif i>11 and i<17:
                constb.append(float(eval(row[a])))          # change to 'row[1]' for old values
            elif i==18:
                beta = float(eval(row[a]))                  # change to 'row[1]' for old values
            elif i>19:
                try:
                    p0.append(float(eval(row[a])))          # change to 'row[1]' for old values
                except:
                    p0.append(float(eval(row[1])))
            i += 1
    
    with open("Philippines incidence rates.csv", "r") as csvfile:
        csvreader = csv.reader(csvfile, delimiter = ",")
        j = 0
        for row in csvreader:
            if j<3:
                incidence.append([])
                for i in row:
                    incidence[j].append(float(i)/100000)
            else:
                prevalence.append([])
                for i in row:
                    prevalence[j-3].append(float(i)/100000)
            j += 1
    
    incidence = np.array(incidence)
    prevalence = np.array(prevalence)
    consta = np.array(consta)
    constb = np.array(constb)
    p0 = np.array(p0)
    
    p0[3] = incidence[bound, 0]*constb[4]
    p0[2] = (prevalence[1, 0] - incidence[bound, 0]*constb[4])
    
    p0[0] = 0.2487682813
    p0[1] = 0.2487682813
    p0[4] = 0.2487682813
    p0[5] = 0.2487682813
    
    beta = betatest
    
    t = np.arange(beginYear, endYear+1, 1)
    
    fig, ax = plt.subplots(1, 2, figsize=(15, 6))
    fig.suptitle("Fitted Model for Years (" + str(beginYear) + "-" + str(endYear) +")", fontsize=20)
    
    beta = fit(iterations, rrr, ac, beta, p0, ax)
    
    # plt.show()
    
    return p0, beta