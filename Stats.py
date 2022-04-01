from cmath import sqrt
import statistics as stats
import scipy.stats as st

############################################
def EpsilonC(p, pe, n):

    print(abs(p - pe)/sqrt((p*(1 - p)/n)))

def EpsilonH(pe1, pe2, n1, p2, p):
    
    print(abs(pe1 - pe2)/sqrt( (p*(1 - p)/(1/n1 + 1/n2) ) ) )

#############################################
def Min(n2, pe2):
    
    M1 = n2*pe2
    M2 = n2*(1 - pe2)
    
    if M1 < M2:
        print(M1)
    else:
        print(M2)

#############################################
def P(pe1, pe2, m1, m2):
    
    print((m1*pe1 + m2*pe2)/(m1+m2))
    
#############################################
def mysum(a, b):
    result = 0
    for i in a:
        result += i - b
        
    return result

def NormDis(ArrMean, SD):
    
    print(stats.NormalDist(ArrMean, SD))
    
def NormDisMan(N, Mean, X):
    
    print(sqrt( (mysum(X, Mean)**2)/N ))

#############################################
def Conformity(Ti, Oi):
    result = 0
    for t in Ti:
        for o in Oi:
            result += (t - o)**2 /t
    
    print(result)

##############################################
def HomInde(Oij):
        
    SumR = []
    SumC = []

    for r in Oij:
        for i in range(Rows):
            SumR.append(sum(r))
        for c in r:
            for j in range(len(Oij[0])):
                SumC.append(c)
            
    Tij = [[]*len(Oij[0])]*Rows
    
    
    for r in range(Rows):
        for c in range(len(Oij[0])):
            Tij[r].append(SumR[r]*SumC[c]/sum(SumC))
            
    Chi = 0
    for r in range(Rows):
        for c in range(len(Oij[0])):
            Chi += (Tij[r][c] - Oij[r][c])**2 /Tij[r][c]
            
    print(Chi)
    
##############################################
def SSquare(n, xarr, xarr2):

    xbar = xarr/n

    s2 = xarr2/n - xbar**2

    print("s² = " + s2)
    print(st.norm.ppf(s2))

    return s2

##############################################
    
##################################################################################################################################

print("1 = Epsilon Conformity, 2 = Epsilon Homogenité, 3 = Min, 4 = P, 5 = Normal Distribution, 6 = Test Conformity, 7 = Test Homogenité/Indépendance\nWhat formula?")
form = int(input())

#Epsilon Conformity
if form == 1 :
    print("p = ")
    p = float(input())
    print("n = ")
    n = float(input())
    print("pe = (if you don't know type in 0)")
    pe = float(input())
    if pe == 0:
        print("missing values? Lemme help..")
        pe = float(input())/n

    EpsilonC(p, pe, n)
    
#Epsilon Homogenity
if form == 2:
    
    print("pe1 = ")
    pe1 = float(input())
    print("pe2 = ")
    pe2 = float(input())
    print("n1 = ")
    n1 = float(input())
    print("n2 = ")
    n2 = float(input())
    print("p =  (if you don't know type in 0)")
    p = float(input())
    if p == 0:
        print("missing values? Lemme help..")
        p = float(input())/n

    EpsilonH(pe1, pe2, n1, n2, p)

#Min
if form == 3 :
    print("n2 =")
    n2 = float(input())
    print("pe2 = ")
    pe2 = float(input())
    
    Min(n2, pe2)
    
#P Po
if form == 4:
    print("pe1 =")
    pe1 = float(input())
    print("pe2 =")
    pe2 = float(input())
    print("m1 =")
    m1 = float(input())
    print("m2 =")
    m2= float(input())
    
    P(pe1, pe2, m1, m2)

#Stadard Deviation
if form == 5:
    
    print("If you don't have the Arithmetic mean and Stadard Deviation input 0 for either")
    
    print("Arithmetic Mean =")
    ArrMean = float(input())
    print("Standard deviation")
    SD = float(input())
    
    if ArrMean and SD != 0 :
        NormDis(ArrMean, SD)
        
    else:
        print("Population/Sample size =")
        N = int(input())
        print("Pop/Sample Mean =")
        Mean = float(input())
        print("Each value in the Pop/Sample (input -1 when done) =")
        X = []
        
        i = 0
        while True :
            X += float(input())
            
            if X[i] < 0 :
                X.pop(i)
                NormDisMan(N, Mean, X)
                break;
            
            else :
                i+=1
            
#Test Conformity. Ti = Effectif Theorique, Oi = Effectif Observé.            
if form == 6:
    print("Ti values (input -1 when done) =")
    Ti = []
    
    i = 0
    while True :
        Ti.append(float(input()))
        
        if Ti[i] < 0 :
            Ti.pop(i)
            break;
        
        else :
            i+=1
                
    print("Oi values (input -1 when done) =")
    Oi = []
    
    i = 0
    while True :
        Oi.append(float(input()))
        
        if Oi[i] < 0 :
            Oi.pop()
            break;
        
        else :
            i+=1
    
    Conformity(Ti, Oi)
    
#Test d'Homogénité ou d'indépendance
if form == 7:
    print("Number of Rows?")
    Rows = int(input())
    
    print("Oij Row values (input -1 when done) =")
    Oij = [[]]*Rows

    i = 0
    for j in range(Rows):
        print("Row", j+1)
        
        while True :

            Oij[j].append(float(input()))
            
            if Oij[j][i] < 0 :
                Oij[j].pop()
                break;
            
            else :
                i+=1
    HomInde(Oij)

#SSquare
if form == 8:
    print("n =")
    n = int(input())
    print("Sum xi =")
    xarr = float(input())
    print("Sum xi² =")
    xarr2 = float(input())

    SSquare(n, xarr, xarr2)