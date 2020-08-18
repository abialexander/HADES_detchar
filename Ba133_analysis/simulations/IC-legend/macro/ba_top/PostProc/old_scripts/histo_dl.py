import sys
import uproot
import math
import numpy as np
import sympy
import matplotlib.pyplot as plt    
import yaml
from pprint import pprint
import matplotlib.scale as scale
import random
import os

if(len(sys.argv) != 6):
    print('Usage: python3 histo_dl.py constants_IC160A.conf IC160A_m_q.conf [A value] [C value] output_name')
    sys.exit()

conf_path = sys.argv[1]
m_q_path = sys.argv[2]
A = float(sys.argv[3])    #where the transition layer ends
C = float(sys.argv[4])    #where the dead layer ends
output= str(sys.argv[5])


with open(conf_path) as fid:
    config = yaml.load(fid, Loader=yaml.Loader)
pprint(config)

r_c     = config['r_c']
R_b     = config['R_b']
R_u     = config['R_u']
h_c     = config['h_c']
H       = config['H']
H_u     = config['H_u']
offset  = config['offset']


with open(m_q_path) as m_q:
    m_q_config = yaml.load(m_q, Loader=yaml.Loader)
pprint(m_q_config)


r_E     = m_q_config['r_E']
z_F     = m_q_config['z_F']
m       = m_q_config['m']
q       = m_q_config['q']
q_A     = m_q_config['q_A']






def f_crystal(x):
    if x<C:
       return 0*x
    elif C<= x <=A:
       return 1/(A-C)*x-C/(A-C)
    else:
        return 1+0*x

def f_cavity(x):
    if x<C/2:
       return 0*x
    elif C/2<= x <=A/2:
       return 2/(A-C)*x-C/(A-C)
    else:
        return 1+0*x


def f_side_up(x):
    return m*x+q

def f_side_up_A(x):
    return m*x+q_A


def f_smear(x):
    a=0.35
    b=1.99e-3
    return math.sqrt(a+b*x)

def main():

    file=uproot.open("one_chain.root")
    tree=file["one_tree"] 
    Edep=tree.array("Edep")
    X=tree.array("x")
    Y=tree.array("y")
    Z=tree.array("z")


    total_energy=[]
    i=0
    print(len(X))

    script_dir = os.path.dirname(__file__) 
    rel_path = "energy_lists_new/"
    abs_file_path = os.path.join(script_dir, rel_path)

    with open(abs_file_path + 'energy_list_'+ output +'.txt', 'w+') as file_output:
        while i<len(X):
            sum_Edep=0
            j=1
            while j<len(X[i]):
                r=math.sqrt(math.pow(X[i][j],2)+math.pow(Y[i][j],2))
                z=Z[i][j]
                if R_b-A<=r<=R_b and z_F<=z<=H+offset:   #around crystal
                    side=R_b-r
                    evaluated=f_crystal(side)
                    Edep[i][j]=Edep[i][j]*evaluated
                elif r_E<=r<=R_b and offset<=z<=z_F and f_side_up(r)<= z<=f_side_up_A(r):                #up_side
                    side=R_b-r
                    up=z-offset
                    side_up=abs(z-f_side_up(r))/math.sqrt(1+math.pow(m,2))
                    minimum=min(side,up,side_up)
                    evaluated=f_crystal(minimum)
                    Edep[i][j]=Edep[i][j]*evaluated
                elif r_c+A/2<=r<r_E  and offset<=z<=A+offset:  #up
                    up=z-offset
                    evaluated=f_crystal(up)
                    Edep[i][j]=Edep[i][j]*evaluated
                elif r_c<=r<r_c+A/2 and offset<=z<=A+offset:  #corner cavity
                    side=r-r_c
                    up=z-offset
                    minimum=min(side,up)
                    if up==minimum:
                        evaluated=f_crystal(minimum)
                    else:
                        evaluated=f_cavity(minimum)
                    Edep[i][j]=Edep[i][j]*evaluated
                elif r_c<=r<=r_c+A/2 and offset+A/2<= z<=h_c+offset:  #around cavity
                    side_cavity=r-r_c
                    evaluated=f_cavity(side_cavity)
                    Edep[i][j]=Edep[i][j]*evaluated
                elif r<=A/2+r_c and h_c+offset<z<=h_c+A/2+offset:   #down cavity    
                    side_cavity=r-r_c
                    down_cavity=z-(offset+h_c)
                    minimum=min(side_cavity,down_cavity)
                    evaluated=f_cavity(minimum)
                    Edep[i][j]=Edep[i][j]*evaluated
                else:
                    evaluated=1.0
                    Edep[i][j]=Edep[i][j]*evaluated
                #print(evaluated)
                sum_Edep+=Edep[i][j]
                j=j+1   
                #print("sum_Edep",sum_Edep)
            if sum_Edep!=0:
                sum_Edep=sum_Edep*1000
                smear_sum_Edep=random.gauss(sum_Edep,(f_smear(sum_Edep))/2.355)
                total_energy.append(smear_sum_Edep)
                file_output.write("{}\n".format(smear_sum_Edep))
            i=i+1
    #print("total_energy",total_energy)

    print(len(total_energy))
    print(R_b)
    plt.hist(total_energy,200)
    plt.yscale('log')
    plt.show()


if __name__=="__main__":
    main()

