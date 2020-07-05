import numpy as np
from math import e
import math
import  scipy as sy
import scipy.integrate as integrate
import scipy.stats as sp

f = open("data.csv", "r")

line_A = []
line_B = []
line_C=[]
start = f.readline()
for line in f:
    reading = line.rstrip().split(',')
    reading_as_float = [float(n) for n in reading]
    line_A.append(reading_as_float[0])
    line_B.append(reading_as_float[1])
    line_C.append(reading_as_float[2])


# checking column A
########################
#A as normal
#calculate miu
miu = sum(line_A)/len(line_A)
# calculate std
var = 0
for num in line_A:
    var+= (num-miu)**2
std = (var/len(line_A))**0.5
# calculate L (Normal distribution)
sum_L = 0
for num in line_A:
    sum_L -= (((num-miu)**2))

L_N = ((-len(line_A)/2)*np.log(2*np.pi))-((len(line_A)/2)*np.log(std**2)) +sum_L/(2*std**2)
#calculate AIC (Normal distribution)
K = 2 # we found miw and std
AIC_N = -2*(L_N)+2*K
#################################3
#A as Poisson
# A can not be a Poisson
##########################################
#A as Gamma
avg_A=sum(line_A)/len(line_A)
log=0
for num in line_A:
    log+=np.log(num)
avg_log=log/len(line_A)
a=(0.5)/(np.log(avg_A)-avg_log)
b=avg_A/a
n = len(line_A)
g = 16.3210242895836 # we found it from the excel

L_G = len(line_A)*(a-1)*(avg_log)-(len(line_A))*(g)-(len(line_A))*a*np.log(b)-((len(line_A))*avg_A)/b

#for num in line_A:
#    L_G*=((b**a)*(e**(-b*num))*(num**(a-1)))/ans_integral
K = 2  # we found a and b
AIC_G = -2*(L_G)+2*K
######################################################
#A as Triangular
line_A.sort()
a = min(line_A) #min number
b = max(line_A) #max number
c = 42.4
value_c = 0
max_likelihood = -1000000000000000
while c<b:
    f_x = 0
    for num in line_A[1:499]:
        if a<= num and num<c:
            f_x += (np.log(2)+np.log(num-a)-(np.log(c-a)+np.log(b-a)))
        if c<num and num<=b:
            f_x += (np.log(2)+np.log(b-num)-(np.log(b-c)+np.log(b-a)))
    if f_x>max_likelihood:
        max_likelihood = f_x
        value_c = c
    c+=0.1
L_T = max_likelihood
K = 3 # we found a, b and c
AIC_T = -2*(L_T)+2*K
#######################################
# checking who is the smallest
list_of_AIC_A = [AIC_N, AIC_G,AIC_T]
min_AIC = min(list_of_AIC_A)

print ('we got that the triangle has the smallest AIC',min(list_of_AIC_A))

print ('we got in the excel by KS test that it is triangle distribution with a =',a,'b =',b,'c =',value_c)
############################################################################
# checking column B and C
####################################
avg_B = sum(line_B)/len(line_B)
avg_C = sum(line_C)/len(line_C)
n_b = len(line_B)
n_c = len (line_C)

def std_equal():
    s2_B = 0
    s2_C = 0
    sum_B = 0
    sum_c = 0

    for num in line_B:
        sum_B += (num-avg_B)**2
    s2_B = (sum_B)/((len(line_B))-1)

    for num in line_C:
        sum_c += (num-avg_C)**2
    s2_C = (sum_c)/((len(line_C))-1)

    stat_F = s2_B/s2_C

    alpha = 0.05

    crit_F_1 = sp.f.ppf(1-alpha/2, len(line_B)-1, len(line_C)-1)
    crit_F_2 = 1/crit_F_1 # B and C have the same len

    if stat_F>crit_F_1 or stat_F<crit_F_2:
        return False
    else:
        return True

std_is_equal = 0
if std_equal()==True:
    std_is_equal = True
else:
    std_is_equal = False

def miu_equal():
    d0 = 0 # we want them to be equal
    sum_B = 0
    sum_c = 0
    alpha = 0.05

    for num in line_B:
        sum_B += (num-avg_B)**2
    s2_B = (sum_B)/((len(line_B))-1)

    for num in line_C:
        sum_c += (num-avg_C)**2
    s2_C = (sum_c)/((len(line_C))-1)

    sp2 = ((n_b-1)*s2_B+(n_c-1)*s2_C)/(n_b+n_c-2)

    stat_T = (avg_B-avg_C-d0)/((sp2*((1/n_b)+(1/n_c)))**0.5)

    critic_T_1 = sp.t.ppf(1-alpha/2, n_b+n_c-2)
    critic_T_2 = -critic_T_1

    if stat_T>critic_T_1 or stat_T<critic_T_2:
        return False
    else:
        return True

miu_is_equal = 0
if miu_equal()==True:
    miu_is_equal = True
else:
    miu_is_equal = False

# cheking if B and C are have the same miu amd std^2
if miu_is_equal == True and std_is_equal == True:
    print ('B and C are have the same miu and std^2')
else:
    print ('B and C are not have the same miu and std^2')

# join column B with C
join_list =  line_B+line_C

miu_join_list = sum(join_list)/len(join_list)
for num in join_list:
    var+= (num-miu)**2
std_join_list = (var/len(join_list))**0.5

print ('After we join both column we got that the miu =', miu_join_list, 'and the std = ', std_join_list)









f.close()