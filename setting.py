import numpy as np
from numpy import random
np.random.seed()
#orbit settings
sigma_e=0.01
sigma_i=0.01
k_Hill=11.
P_inner=7.


#parallel settings
max_runs=1
num_proc = 1
nnodes=200

#default configuration
N_pl = np.random.random_integers(7,11)

#integration settings
integrator="ias15"
t_max=1.e8#*2.*np.pi
Noutputs=1000
doubling_time = 1.e6/4.5
acc_time=doubling_time/1000.


#path settings
basename='/home/shore/mnt/code/'
pythonpath='/home/shore/src/virtualenv-1.5.2/ve/'
rundir='/home/shore/mnt/data/Warm_Jupiter/round2/'
#rundir='/home/shore/mnt/data/Warm_Jupiter/'
subfile="qsubrebound"


#restart settings
restartinput="testeject.txt"
#restart from binary file
frombin=False #USE ONLY FOR RESTARTING RUNS
binfile=""
#restarting from a .txt file
fromtxt = False
txtfile =""

#output settings
checkpoint=True

#other flags
debug=False
addmass=True
gridsearch=True
check_distances = False
check_merger_time = False

#
statuscode={"eject":2,"star":3,"collision":1,"survive":0}
