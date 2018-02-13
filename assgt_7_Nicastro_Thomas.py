import numpy as np
import matplotlib.pyplot as plot
from scipy.optimize import minimize as scip

#Loads info from text file
def loadTxt(fname):
    with open(fname) as f:
        x = []
        for row in f.readlines()[2:]:
            line = row.split()
            #Length of milky way is only 41 kpc I believe so will cut off after 41 kpc
            if(float(line[0]) < 41.0):
                x.append([float(r) for r in line])
            else:
                f.close()
                
    r_kpc = [i[0] for i in x]
    dr_kpc = [i[1] for i in x]
    v_kms = [i[2] for i in x]
    dv_kms = [i[3] for i in x]
    
    return r_kpc, dr_kpc, v_kms, dv_kms

#Takes the square root of all 3 velocities and returns it
def V(v_1, v_2, v_3):
    return np.sqrt(v_1+v_2+v_3)

#Calculates the velocity squared 
def v_circ(R, M, a):
    G = 4.49*(10**-6)
    v_c_sqrd = ((R*R)*G*M)/np.sqrt(((R*R)+(a*a))**3)
    return v_c_sqrd

#Funtion for optimizing use the chi-squared method
def v_opt(parm, v, dv, r):
    m_h = parm[0]
    m_d = parm[1]
    m_b = parm[2]
    a_h = parm[3]
    a_d = parm[4]
    a_b = parm[5]
    opt = 0

    #Sums up all values for optimization
    for i in range(56):
        V_opt = v[i] - V(v_circ(r[i], m_h, a_h), v_circ(r[i], m_b, a_b), v_circ(r[i], m_d, a_d))
        opt_i = ((V_opt)**2)/(dv[i]**2)
        opt += opt_i
        
    return opt
#Plots data
def plotData(x, y, y_model, y_model2=0):
    plot.scatter(x, y, label= "Data")
    plot.plot(x, y_model, label = "Model")
    #plot.plot(x, y_model2, label = "Guess")
    plot.title("Radius [kpc] vs rotation speed[kms]")
    plot.xlabel("Radius [kpc]")
    plot.ylabel("Rotation speed [kms]")
    plot.ylim([0,300])
    plot.legend()
    plot.show()

def main():
#****************Initial guesses***********************
    m_h_AU = 10**12
    a_h_kpc = 12

    m_b_AU = 5.0*(10**10)
    a_b_kpc = 0.5

    m_d_AU = 5*(10**11)
    a_d_kpc = 4
#****************Initial guesses***********************
    #Getting values from text file
    r_kpc, dr_kpc, v_kms, dv_kms = loadTxt("LogRC_data.dat")
    
    initial = [m_h_AU, m_d_AU, m_b_AU, a_h_kpc, a_d_kpc, a_b_kpc]
    
    #Optimizes function using nelder mead method
    res = scip(v_opt, initial, args=(v_kms, dv_kms, r_kpc,), method='Nelder-Mead')

    m_h_mod_AU = res.x[0]
    m_d_mod_AU = res.x[1]
    m_b_mod_AU = res.x[2] 
    a_h_mod_kpc = res.x[3]
    a_d_mod_kpc = res.x[4]
    a_b_mod_kpc = res.x[5]
    
    v_model_kms = []
    
    for R_kpc in r_kpc: 
        v_model_kms.append(V(v_circ(R_kpc, m_h_mod_AU, a_h_mod_kpc), v_circ(R_kpc, m_b_mod_AU, a_b_mod_kpc), v_circ(R_kpc, m_d_mod_AU, a_d_mod_kpc)))
    print(len(r_kpc))
    print(len(v_kms))
    print(res.x)
    plotData(r_kpc, v_kms, v_model_kms)
    
main()
