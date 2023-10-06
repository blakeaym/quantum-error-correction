from distutils.log import error
import random
# from matplotlib.pyplot import figure
import numpy as np
import itertools
from matplotlib import pyplot as plt
from sympy import pretty, pretty_print
import math

input_parity=np.array([[1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1, 0, 1, 0, 0, 0, 1, 1, 0, 0], 
[0, 1, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1, 0, 1, 1, 1, 0, 0, 1, 1, 0], 
[0, 0, 1, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1, 0, 0, 1, 0, 0, 0, 0, 1, 1], 
[0, 0, 0, 1, 0, 1, 0, 0, 0, 0, 0, 0, 1, 1, 1, 0, 1, 1, 0, 0, 1, 0], 
[0, 0, 0, 0, 1, 1, 0, 0, 0, 0, 0, 0, 0, 1, 1, 0, 0, 1, 1, 0, 1, 1], 
[0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 1, 1, 0, 0, 1, 1, 1, 0, 1], 
[0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 1, 0, 1, 0, 1, 1, 1, 1, 1, 0, 1, 1], 
[0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 1, 0, 0, 0, 1, 0, 1, 0, 1, 0, 1, 0], 
[0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 0, 1, 0, 1, 0, 0, 1, 0, 1, 0, 0], 
[0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0]])

def pauli_print(bin):
    str=""
    for i in bin:
        if i==1:
            str+="X"
        elif i==2:
            str+="Z"
        elif i==3:
            str+="Y"
        else:
            str+="I"
    print(str)


def kones(n,k):
    results=[]
    for bits in itertools.combinations(range(n),k):
        s=[0]*n
        for bit in bits:
            s[bit]=1
        results.append(s)
    return results

def make_parity():
    splitted=np.hsplit(input_parity,2)
    H_Z=np.concatenate((splitted[0], np.zeros(splitted[0].shape)),axis=1)
    H_X=np.concatenate((np.zeros(splitted[1].shape),splitted[1]),axis=1)
    return np.concatenate((H_X,H_Z),axis=0)

def pauli_ter_to_bin(ter):
    n=len(ter)
    res=[0]*2*n
    for i,t in enumerate(ter):
        if t==1:
            res[i]=1
        elif t==2:
            res[i+n]=1
        elif t==3:
            res[i]=1
            res[i+n]=1
    return res

def make_look_up_table(check_matrix,n,up_to_k):
    # syndromes=np.zeros(4)
    paulis = [1,2,3]
    errors=[]
    for k in range(1,up_to_k+1):
        for e in np.array(kones(n,k)):
            ones=np.where(e==1)[0]
            for ps in itertools.product(paulis, repeat=k):
                # print(ps)
                temp=[0]*n
                for ip in range(len(ps)):
                    temp[ones[ip]]=ps[ip]

                errors.append(temp)

    bin_errors=list(map(pauli_ter_to_bin, errors))
    bin_errors.insert(0,[0]*22)
    for i,e in enumerate(bin_errors):
        if i<=1000:
            print(e)
    # print(bin_errors)
    look_up_table=dict()
    for e in bin_errors:
        syndrome=np.dot(check_matrix,e).astype(int) %2
        if tuple(syndrome.tolist()) not in look_up_table:
            look_up_table[tuple(syndrome.tolist())]=tuple(e)
    return look_up_table
                

def make_errors(error, bit,rate,n):
    r=random.random()
    if r<rate:
        # print(r)
        r=random.random()
        if r<=1/3:
            # X error
            error[bit]=1
        elif r<=2/3:
            # Z error
            error[n+bit]=1
        else:
            # Y error
            error[bit]=1
            error[bit+n]=1
    return error

def simulate(noise_rate, trials,n,lookup,parity):
    num_decode_err=0
    for k in range(trials):
        error=[0]*2*n
        for bit in range(n):
            error=make_errors(error,bit,noise_rate,n)
        # print("error=",error)
        syndrome=np.dot(parity,error) %2
        # print(lookup.get(tuple(syndrome)))
        correction=lookup.get(tuple(syndrome))
        # print(syndrome,": ", correction)
        if correction==None:
            # print(error,": ", correction)
            num_decode_err+=1
            if(error.count(1)<=2):
                print("!!!")
                print(error,": ", correction)
        elif error!=list(correction):
            # print(error,": ", correction)
            num_decode_err+=1
            if(error.count(1)<=2):
                print("!!!")
                print(error,": ", correction)
    return num_decode_err

# print(np.hsplit(input_parity,2))
def trial():
    parity=make_parity()
    lookup=make_look_up_table(parity,11,4)
    # print(lookup)
    print(len(lookup))

    probs=[0.1,0.09,0.08,0.07,0.06,0.05,0.04,0.03,0.02]
    res=[]
    for i in range(2):
        for p in probs:
            act_p=p*(10**-i)
            trials=10**5
            num_decode_err=simulate(act_p,trials,11,lookup,parity)
            print(num_decode_err)
            print(num_decode_err/trials)
            res.append([act_p, num_decode_err/trials])
    print(res)
    for r in res:
        p3=r[0]**3
        if r[1]!=0:
            print(math.log10(r[1]/p3))

def test():
    parity=make_parity()
    lookup=make_look_up_table(parity,11,4)
    key=np.dot(parity,[0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0])
    print(lookup.get(tuple(key)))
    # print(lookup)
    # print(len(lookup))
    # num_decode_err=simulate(0.1,10000,11,lookup,parity)
    # print(num_decode_err)

def graph():
    # data=np.array([[0.1, 0.00603], [0.06, 0.00131], [0.04642, 0.00054], [0.02154, 6e-05]])
    # data=np.array([[0.1, 0.005914], [0.06, 0.001269], [0.04642, 0.000549], [0.02154, 5.5e-05], [0.015, 1.2e-05], [0.010000000000000002, 5e-06], [0.006, 2e-06], [0.004642, 0.0], [0.002154, 0.0], [0.0015, 0.0]])
    # data=np.array([[0.1, 0.005907], [0.06, 0.001259], [0.04642, 0.00055], [0.02154, 5.8e-05], [0.015, 1.5e-05], [0.010000000000000002, 1e-05], [0.006, 2e-06], [0.004642, 0.0], [0.002154, 1e-06], [0.0015, 0.0]])
    # data=np.array([[0.1, 0.09001], [0.09, 0.06999], [0.08, 0.05358], [0.07, 0.03768], [0.06, 0.02518], [0.05, 0.01542], [0.04, 0.00833], [0.03, 0.00362], [0.02, 0.00139], [0.010000000000000002, 0.00019], [0.009, 0.00014], [0.008, 7e-05], [0.007000000000000001, 3e-05], [0.006, 7e-05], [0.005000000000000001, 5e-05], [0.004, 1e-05], [0.003, 0.0], [0.002, 0.0], [0.001, 0.0], [0.0009, 0.0], [0.0008, 0.0], [0.0007000000000000001, 0.0], [0.0006, 0.0], [0.0005, 0.0], [0.0004, 0.0], [0.0003, 0.0], [0.0002, 0.0]])
    data=np.array([[0.1, 0.089552], [0.09, 0.069504], [0.08, 0.05176], [0.07, 0.036894], [0.06, 0.024688], [0.05, 0.01511], [0.04, 0.008212], [0.03, 0.00368], [0.02, 0.001218], [0.010000000000000002, 0.000165], [0.009, 0.000109], [0.008, 7.2e-05], [0.007000000000000001, 5.1e-05], [0.006, 3.1e-05], [0.005000000000000001, 1.9e-05], [0.004, 8e-06], [0.003, 4e-06], [0.002, 0.0]])
    for i in data:
        i[1]=i[1]/(i[0]**3)
    print(data)
    log_data=np.log10(data)
    # log_break_even=np.log10(break_even)
    print(log_data)

    # plt.figure(figsize=(10,10))
    # # plt.plot(log_data[:,0],log_data[:,1],'ko',log_data[:,0],log_data[:,1],'k-',log_break_even,log_break_even,'r--')
    plt.plot(log_data[:,0],log_data[:,1],'ko')
    plt.show()

# test()
# print(len((0, 1, 1, 1, 1, 0, 1, 1, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 0, 0)))

trial()
# data=np.array([[0.1, 0.225351], [0.04642, 0.063134], [0.02154, 0.015193], [0.01, 0.003453], [0.004642, 0.000753], [0.002154, 0.000171], [0.001, 3.9e-05], [0.0004642, 8e-06], [0.0002154, 1e-06]])
# # data=np.array([[0.1, 0.225351], [0.04642, 0.063134], [0.02154, 0.015193], [0.01, 0.003453], [0.004642, 0.000753], [0.002154, 0.000171], [0.001, 3.9e-05], [0.0004642, 8e-06], [0.0002154, 1e-06]])
# break_even=np.array([0.1,0.04642,0.02154,0.01,0.004642,0.002154,0.001,0.0004642,0.0002154])

graph()