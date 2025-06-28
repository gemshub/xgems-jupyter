# author R.A.PatelR.A.Patel 2021
# modified G.D. Miron

import numpy as np
from scipy.integrate import solve_ivp
from run.GEMSCalc import GEMS
import matplotlib.pylab as plt
#from bqplot import pyplot as plt

def bouges_composition(CaO,SiO2,Al2O3,Fe2O3,SO3):
    """
    computation of clinker phases from oxides using bouge's conversion
    """
    clink_phases={}
    clink_phases["C3S"]=C3S=4.07*CaO -7.6*SiO2-6.72*Al2O3-1.43*Fe2O3-2.85*SO3 
    clink_phases["C2S"]=C2S=2.87*SiO2-0.75*C3S 
    clink_phases["C3A"]=C3A=2.65*Al2O3-1.69*Fe2O3 
    clink_phases["C4AF"]=C4AF=3.04*Fe2O3
    return clink_phases
    
#parameters for parrot and killoh model source: lothenbach et al. 2008
K1 = {"C3S":1.5,"C2S":0.5,"C3A":1.,"C4AF":0.37}
N1 = {"C3S":0.7,"C2S":1.0,"C3A":0.85,"C4AF":0.7}
K2 = {"C3S":0.05,"C2S":0.02,"C3A":0.04,"C4AF":0.015}
K3 = {"C3S":1.1,"C2S":0.7,"C3A":1.0,"C4AF":0.4}
N3 = {"C3S":3.3,"C2S":5.0,"C3A":3.2,"C4AF":3.7}
H  = {"C3S":1.8,"C2S":1.35,"C3A":1.6,"C4AF":1.45}
Ea = {"C3S":41570,"C2S":20785,"C3A":54040,"C4AF":34087} #from barbara's excel sheets in paper only 2 leading digits are given
T0 = 20#c
ref_fineness = 385

output_times = [0.,0.1,1,3,7,28,30,60,90,365, 730,1000,2000] #range for output of time [days]. Always at 0 for initial conditions


"""
This cell provides implementation of parrot and killoh rate model implemented above
"""
def nuclandgrowth(alpha_t,K1,N1,fineness,ref_fineness):
    """
    Nucleation and growth rate equation of parrot and killoh model
    """
    ln = np.log
    rt = (K1/N1) * (1-alpha_t) * (-ln(1-alpha_t))**(1-N1) * (fineness/ref_fineness)
    return rt

def diffusion(alpha_t,K2):
    """
    Diffusion rate equation of parrot and killoh model
    """
    rt = K2 * (1-alpha_t)**(2/3)/(1-(1-alpha_t))**(1/3)
    return rt

def hydrationshell(alpha_t,K3,N3):
    """
    """
    rt = K3 * (1-alpha_t)**N3
    return rt

def f_wc_lothenbach(wc,H,alpha_t):
    """
    functional dependence on w/c in parrot and killoh model
    """
    alpha_cr  = H*wc
    if alpha_t > alpha_cr:
        return (1+3.333*(H*wc-alpha_t))**4
    else:
        return 1

def f_rh(RH):
    """
    functional dependence on relative humidity
    """
    return ((RH-0.55)/0.45)**4

def f_temp(T,T0,Ea):
    """
    functional dependence on temperature
    """
    exp = np.exp
    R = 8.314
    T0  = T + 273.15
    T  =  T + 273.15
    return  exp((-Ea/R) * ((1/T) - (1/T0)))

def overall_rate(t,alpha_t,K1,N1,K2,K3,N3,H,wc,RH,fineness,ref_fineness,T,T0,Ea):
    """
    overall rate of cement clinker phase based on parrot and killoh model
    """
    if alpha_t >= 1: alpha_t = 0.9999
    r1= nuclandgrowth(alpha_t,K1,N1,fineness,ref_fineness)
    r2= diffusion(alpha_t,K2)
    r3= hydrationshell(alpha_t,K3,N3)
    r= min(r1,r2,r3)
    r= r * f_wc_lothenbach(wc,H,alpha_t) * f_rh(RH) * f_temp(T,T0,Ea)
    return r
    
def parrot_killoh(wc, RH, T, fineness):
    DoH = {} # degree of hydration 
    for phase in ["C3S","C2S","C3A","C4AF"]:
        ivp_out = solve_ivp(overall_rate, [0,np.max(output_times)], [1e-15],t_eval=output_times
                           ,args=(K1[phase],N1[phase],K2[phase],K3[phase],N3[phase],
                           H[phase],wc,RH,fineness,ref_fineness,T,T0,Ea[phase]))
        DoH[phase] = np.squeeze(ivp_out.y[0])
    return DoH
    
def plot_bars(results):
    #fix color notation for each phase
    unique_phases = []
    for t in results.keys():
        for phase in results[t].keys():
            if results[t][phase] > 0:
                if phase not in unique_phases:
                    unique_phases.append(phase)
    colors = plt.get_cmap('tab20')(np.linspace(0., 1, len(unique_phases)+1))
    phase_colors = {}
    for phase,color in zip(unique_phases,colors):
        phase_colors[phase]=color
    i  = 0
    ytks= []
    ytksval = []
    for t in results.keys():
        prev =0
        for key,val in results[t].items():
            if key in unique_phases:
                plt.barh(i, val, left=prev, height=0.5,
                        label=key,color =phase_colors[key])
                prev += val
        ytks.append(i)
        ytksval.append(str(t))
        i+=1
        if t==0:
            plt.legend(ncol=3,bbox_to_anchor=(0, 1.5), loc='upper left')
    plt.yticks(ytks,ytksval)
    plt.xlabel('Volume fraction [m$^3$/m$^3$ of initial volume]')
    plt.ylabel("Time [days]")    
    return plt  
    
def to_phase_first_dict(results):
    """
    convert time first dict to phase first dict for phase diagrams
    """  
    for t in results.keys():
        if t== 0: 
            out = {}
            for phase in results[t].keys():
                out[phase] = []
        for phase in results[t].keys():
            out[phase].append(results[t][phase])
    for phase in out.keys():
        out[phase] = np.array(out[phase])
    return out

def phase_plot(results,datatype="vfrac"):
    """
    plot phase diagrams for visualization
    """
    unique_phases=[]
    for phase in results.keys():
        if np.count_nonzero(results[phase]) > 0: 
            unique_phases.append(phase)
    colors = plt.get_cmap('tab20')(np.linspace(0., 1, len(unique_phases)+1))
    phase_colors = {}
    for phase,color in zip(unique_phases,colors):
        phase_colors[phase]=color
    prev = 0
    for phase in unique_phases:
        plt.fill_between(output_times, results[phase]+prev, prev,color=phase_colors[phase],label=phase)
        prev += results[phase]
    plt.legend(ncol=3,bbox_to_anchor=(0, 1.5), loc='upper left')
    if datatype=="vfrac":
        plt.ylabel('Volume fraction [m$^3$/m$^3$ of initial volume]')
    elif datatype=="masses":
        plt.ylabel('mass [g/100 g of cement]')

    plt.xlabel("Time [days]") 
    plt.xscale("log")
    if datatype=="vfrac": plt.ylim(0,1)
    plt.xlim((output_times[2]+output_times[1])/2.,output_times[-1])
    return plt     
    
def run_hydration(clink_phases, wc, CSH2,T, DoH):
   input_file = 'gems_files/CemHyds-dat.lst'
   gemsk = GEMS(input_file) # initalize GEMS 
   all_species = clink_phases.copy()
   all_species["H2O@"] = wc * 100 # add water
   all_species["Gp"]=CSH2 # add gypsum
   print(all_species)
   for name in all_species.keys():
       all_species[name]*=1e-3
   gemsk.T = T + 273.15 # set T
   gemsk.add_multiple_species_amt(all_species,units = "kg") # add comosition in the form of clinkers 
   gemsk.add_species_amt("O2",1e-6) # to reduce stiffness related to redox
   gems_vol_frac = {}
   gems_phase_masses={}
   density = []
   for i in range(len(output_times)):
       if i > 0: gemsk.warm_start()
       for phase in clink_phases:
       #    set lower bound to limit dissolution for unreacted phases, limits the hydration of clinker phases to values calculated with the PK model 
           gemsk.species_lower_bound(phase,clink_phases[phase]*(1-DoH[phase][i])*1e-3,units="kg")
       print("Time-->",str(output_times[i]),":",gemsk.equilibrate())
       # collects the results
       gems_phase_masses[output_times[i]]=gemsk.phase_masses
       if output_times[i]==0: init_vol = gemsk.system_volume
       gems_vol_frac[output_times[i]]=gemsk.phase_volumes 
       for key,val in gems_vol_frac[output_times[i]].items():
           gems_vol_frac[output_times[i]][key]/=init_vol
       density.append(gemsk.system_mass/gemsk.system_volume)
   return gems_vol_frac, gems_phase_masses, density

