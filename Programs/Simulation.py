import matplotlib.pyplot as plt
from matplotlib.ticker import (MultipleLocator)
import numpy as np
import csv
from scipy.integrate import odeint

def simulate(beginYear, endYear, endYearSim, scale, newData, p01, beta1):

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
    
    def ploti(_ax, yr):
        _ax.set_xlabel("year")
        _ax.set_ylabel("population proportion per " + str(scale) + " persons")
        _ax.set_xticks(np.arange(yr, yr+endYear-beginYear+1, 1))
        if yr==beginYear:
            _ax.set_ylim(0, 0.00325*scale)
            _ax.set_yticks(np.arange(0, 0.00325*scale, 0.0005*scale))
            _ax.yaxis.set_minor_locator(MultipleLocator(0.0001*scale))
        elif yr==endYear:
            _ax.set_ylim(0, 0.00275*scale)
            _ax.set_yticks(np.arange(0, 0.00275*scale, 0.0005*scale))
            _ax.yaxis.set_minor_locator(MultipleLocator(0.0001*scale))
        _ax.legend()
    
    def plotp(_ax, yr):
        _ax.set_xlabel("year")
        _ax.set_ylabel("population proportion per " + str(scale) + " persons")
        _ax.set_xticks(np.arange(yr, yr+endYear-beginYear+1, 1))
        if yr==beginYear:
            _ax.set_ylim(0, 0.00650*scale)
            _ax.set_yticks(np.arange(0, 0.00650*scale, 0.001*scale))
            _ax.yaxis.set_minor_locator(MultipleLocator(0.00025*scale))
        elif yr==endYear:
            _ax.set_ylim(0, 0.00650*scale)
            _ax.set_yticks(np.arange(0, 0.00650*scale, 0.001*scale))
            _ax.yaxis.set_minor_locator(MultipleLocator(0.0002*scale))
        _ax.legend()
    
    def fnd(_y1, par, delt):
        if par=="inci":
            return _y1[:, 3]/delt
        elif par=="preva":
            return _y1[:, 3]+_y1[:, 2]
    
    consta = []
    constb = []
    p0 = []
    
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
    
    consta = np.array(consta)
    constb = np.array(constb)
    p0 = np.array(p0)
    t = np.arange(beginYear, endYear, 0.01)
    
    p0 = p01.copy()
    beta = beta1
    
    y = odeint(f, p0, t, args=(consta, constb, beta))
    
    fig, ax = plt.subplots(1, 2, figsize=(15, 6))
    fig.suptitle("Fitted Model for Years (" + str(beginYear) + "-" + str(endYear) +")", fontsize=20)
    ax[0].set_title("Tuberculosis Incidence Rate in the Philippines")
    ax[1].set_title("Tuberculosis Prevalence Rate in the Philippines")
    ax[0].plot(t, (y[:, 3]/constb[4])*scale, linewidth=2, label="Incidence rate")
    ax[1].plot(t, (y[:, 3]+y[:, 2])*scale, linewidth=2, label="Prevalence rate")
    ploti(ax[0], beginYear)
    plotp(ax[1], beginYear)
    
    t1 = np.arange(endYear, endYearSim, 0.01)
    
    fig1, axs1 = plt.subplots(2, 2, figsize=(15, 10))
    fig1.suptitle("Simulated Incidence Rates (Years " + str(endYear) + "-" + str(endYearSim) + ")", fontsize=20)
    axs1[0, 0].set_title("Varied Partial Immunity")
    axs1[0, 1].set_title("Varied Vaccine Coverage")
    axs1[1, 0].set_title("Varied Treatment Success")
    axs1[1, 1].set_title("Varied Treatment Duration")
    fig2, axs2 = plt.subplots(2, 2, figsize=(15, 10))
    fig2.suptitle("Simulated Prevalence Rates (Years " + str(endYear) + "-" + str(endYearSim) + ")", fontsize=20)
    axs2[0, 0].set_title("Varied Partial Immunity")
    axs2[0, 1].set_title("Varied Vaccine Coverage")
    axs2[1, 0].set_title("Varied Treatment Success")
    axs2[1, 1].set_title("Varied Treatment Duration")
    
    for i in range(2):
        ax[i].set_xlim([beginYear, endYear])
        for j in range(2):
            axs1[i, j].set_xlim([endYear, endYearSim])
            axs2[i, j].set_xlim([endYear, endYearSim])
    
    constbv = constb.copy()
    for i in range(6):
        y1 = odeint(f, y[-1], t1, args=(consta, constbv, beta))
        y1i = fnd(y1, "inci", constbv[4])
        y1p = fnd(y1, "preva", constbv[4])
        axs1[0, 0].plot(t1, y1i*scale, linewidth=2, label="Partial immunity = {:.2f}".format(constbv[0]))
        axs2[0, 0].plot(t1, y1p*scale, linewidth=2, label="Partial immunity = {:.2f}".format(constbv[0]))
        ploti(axs1[0, 0], endYear)
        plotp(axs2[0, 0], endYear)
        constbv[0] -= constb[0]/5
    
    constbv = constb.copy()
    for i in range(6):
        y1 = odeint(f, y[-1], t1, args=(consta, constbv, beta))
        y1i = fnd(y1, "inci", constbv[4])
        y1p = fnd(y1, "preva", constbv[4])
        axs1[0, 1].plot(t1, y1i*scale, linewidth=2, label="Vaccine Coverage = {:.2f}".format(constbv[2]))
        axs2[0, 1].plot(t1, y1p*scale, linewidth=2, label="Vaccine Coverage = {:.2f}".format(constbv[2]))
        ploti(axs1[0, 1], endYear)
        plotp(axs2[0, 1], endYear)
        constbv[2] += (1-constb[2])/5
    
    constbv = constb.copy()
    for i in range(7):
        y1 = odeint(f, y[-1], t1, args=(consta, constbv, beta))
        y1i = fnd(y1, "inci", constbv[4])
        y1p = fnd(y1, "preva", constbv[4])
        axs1[1, 0].plot(t1, y1i*scale, linewidth=2, label="Treatment Success = {:.2f}".format(constbv[3]))
        axs2[1, 0].plot(t1, y1p*scale, linewidth=2, label="Treatment Success = {:.2f}".format(constbv[3]))
        ploti(axs1[1, 0], endYear)
        plotp(axs2[1, 0], endYear)
        constbv[3] += (1-constb[3])/6
    
    constbv = constb.copy()
    for i in range(7):
        y1 = odeint(f, y[-1], t1, args=(consta, constbv, beta))
        y1i = fnd(y1, "inci", constbv[4])
        y1p = fnd(y1, "preva", constbv[4])
        axs1[1, 1].plot(t1, y1i*scale, linewidth=2, label="Treatment Duration = {:.2f} months".format(12/constbv[1]))
        axs2[1, 1].plot(t1, y1p*scale, linewidth=2, label="Treatment Duration = {:.2f} months".format(12/constbv[1]))
        ploti(axs1[1, 1], endYear)
        plotp(axs2[1, 1], endYear)
        if constbv[1]==12:
            constbv[1] = 24
        else:
            constbv[1] = (constbv[1]*12)/(12-constbv[1])
    
    plt.show()