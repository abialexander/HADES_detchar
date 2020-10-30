import math as m
Qbb_keV = 2039.0
res_func_A = 0.039
res_func_C = 287.466
res_Qbb_keV = res_func_A*m.sqrt(Qbb_keV+res_func_C)
print("Energy resolution at Q_bb in keV:",res_Qbb_keV)
Qbb_MeV = 2.039
res_Qbb_MeV_Dave = (res_func_A/1000.)*m.sqrt(1000.*Qbb_MeV+res_func_C)
print("Energy resolution at Q_bb in MeV (Dave):",res_Qbb_MeV_Dave)
res_Qbb_MeV_Abi = 10*m.sqrt(10)*res_func_A*m.sqrt(Qbb_MeV+1000*res_func_C)
print("Energy resolution at Q_bb in MeV (Abi):",res_Qbb_MeV_Abi)