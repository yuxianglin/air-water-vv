#!/usr/bin/env python3
#Filename: makerun.py
#import numpy as np
#import os
#import shutil

#params = np.loadtxt("wave_parameters.csv", delimiter=',', skiprows=1)
#os.system("mkdir Solitary_runs")
for i in range(6):
    #filepath = "Solitary_runs/"+ "run" + `i`
    #os.system("mkdir " + filepath)
    #os.system("cp *.py " + filepath)
    #os.system("cp tank.stampede.slurm " + filepath)
    f = open("context.options_H="+str((i+1)*0.05),'w')
    #f.write("parallel=True ")
    #f.write("wave_type='single-peaked' ")
    #f.write("gauges=True ")
    #f.write("depth=" + `params[i][0]` + " ")
    f.write("wave_height=" + str((i+1)*0.05) + " ")
    #f.write("peak_period=" + `params[i][2]` + " ")
    #f.write("peak_wavelength=" + `params[i][3]` + " ")
    #f.write("tank_height=" + `params[i][4]`)
    f.close()
