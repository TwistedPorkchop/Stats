from cmath import sqrt
import statistics as stats
import scipy.stats as st

############################################
def EpsilonCP(p, pe, n):

    print((abs(p - pe)/sqrt((p*(1 - p)/n))).real)
    
def EpsilonCSigma(xbar, mu, sigma, n):
    
    print((abs(xbar - mu)/(sigma/sqrt(n))).real)

def EpsilonHP(pe1, pe2, n1, n2, p):
    
    print((abs(pe1 - pe2)/sqrt( (p*(1 - p)/(1/n1 + 1/n2) ) )).real )
    
def EpsilonHSigma(xbar1, xbar2, n1, n2, sig1, sig2):
    
    print((abs(xbar1 - xbar2)/sqrt( (sig1/n1 + sig2/n2) ) ).real )
    
def Comparison(p, f, q, n):
    
    print( abs(p - f)/(sqrt( (p*q/n) )) )
    
def ComparisonMulti(p1, p2, p, n1, n2):
    
    print( abs(p1 - p2)/(sqrt( p*(1-p)*(1/n1 + 1/n2) ) ) )

#############################################
def CI(xbar, n, s, z):
    
    print("Confidence Interval [", xbar - z*s/sqrt(n), ",", xbar + z*s/sqrt(n), "]")
    
def Z_Value(xbar, CI, s, n):
    
    print("Z = ", (xbar - CI)*sqrt(n)/s)
    
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

#############################################
def NormDis(ArrMean, SD):
    
    print(stats.NormalDist(ArrMean, SD))
    
def NormDisMan(N, Mean, X):
    
    print(sqrt( (mysum(X, Mean)**2)/N ))

#############################################
def Conformity(Ti, Oi):
    result = 0
    for i in range(len(Ti)):
        result += (Ti[i] - Oi[i])**2 /Ti[i]
    
    print(result)

#############################################
def Ei(N, Pi):
    Ei = []

    for i in range(len(Pi)):
        Ei.append(N*Pi[i])
        
    print(Ei)

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
def Adjustment(Ei, Oi):

    result = 0

    for i in range(len(Oi)):
        result += ((Oi[i] - Ei[i])**2)/Ei
    
    print("Chi")
##############################################
def SSquare(n, xarr, xarr2, mu):

    xbar = xarr/n

    s2 = xarr2/n - xbar**2

    print("s² = ", s2)
    
    t = (abs(xbar - mu)/(sqrt(s2/(n-1)))).real
    
    print("t = ", t)

##############################################
def Fisher(n1, s1, n2, s2):
    ec1 = (n1*s1**2)/(n1 - 1)
    ec2 = (n2*s2**2)/(n2 - 1)
    
    print(ec2/ec1)
    
##############################################
def Testn(xbar1, xbar2, ss, n1, n2):
    
    print( abs(xbar1 - xbar2)/(sqrt(ss(1/n1 + 1/n2))) )
    
def ss(n1, s1, s2, n2):
    return (n1*s1**2 + n2*s2**2)/(n1+n2-2)
    
##################################################################################################################################

print("1 = Epsilon Conformity, 2 = Epsilon Homogenité, 3 = Min, 4 = P, 5 = Normal Distribution, 6 = Test Conformity, 7 = Test Homogenité/Indépendance\nWhat formula?")
form = int(input())

#Epsilon Conformity
if form == 1 :
    
    print("With p(1) or Stadard deviation(2)?")
    if int(input()) == 1:
        
        print("p = ")
        p = float(input())
        print("n = ")
        n = float(input())
        print("pe = (if you don't know type in 0)")
        pe = float(input())
        if pe == 0:
            print("missing values? Lemme help..")
            pe = float(input())/n

        EpsilonCP(p, pe, n)
        
    else:
        print("xbar =")
        xbar = float(input())
        print("mu = ")
        mu = float(input())
        print("sigma =")
        sigma = float(input())
        print("n = ")
        n = int(input())
        
        EpsilonCSigma(xbar, mu, sigma, n)
    
#Epsilon Homogenity
if form == 2:
    
    print("With p(1) or Stadard deviation(2)?")
    if int(input()) == 1:
    
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

        EpsilonHP(pe1, pe2, n1, n2, p)
    
    else:
        print("X bar 1 = ")
        xbar1 = float(input())
        print("X bar 2 = ")
        xbar2 = float(input())
        print("n1 = ")
        n1 = float(input())
        print("n2 = ")
        n2 = float(input())
        print("Sigma 1 = ")
        sig1 = float(input())
        print("Sigma 2 = ")
        sig2 = float(input())
        
        EpsilonHSigma(xbar1, xbar2, n1, n2, sig1, sig2)

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
    print("mu =")
    mu = float(input())

    SSquare(n, xarr, xarr2, mu)
    
#Test de Fisher
if form == 9:
    print("n1 =")
    n1 = int(input())
    print("s1 =")
    s1 = float(input())
    print("n2 =")
    n2 = int(input())
    print("s2 =")
    s2 = float(input())
    
    Fisher(n1, s1, n2, s2)
    
#Test n<30
if form == 10:
    
    print("xbar1 =")
    xbar1 = float(input())
    print("xbar2 =")
    xbar2 = float(input())
    print("n1 =")
    n1 = int(input())
    print("n2 =")
    n2 = int(input())
    
    print("s1 =")
    s1 = float(input())
    print("n2 =")
    s2 = float(input())
    ss = ss(n1, s1, n2, s2)
    
    Testn(xbar1, xbar2, ss, n1, n2)
    
#Comparison
if form == 11:
    
    print("p =")
    p = float(input())
    print("f =")
    f = float(input())
    print("q =")
    q = float(input())
    print("n =")
    n = float(input())
    
    
    Comparison(p, f, q, n)

#Comparison Multiple Population
if form == 12:
    
    print("p1 =")
    p1 = float(input())
    print("p2 =")
    p2 = float(input())
    print("p =")
    p = float(input())
    print("n1 =")
    n1 = float(input())
    print("n2 =")
    n2 = float(input())
    
    ComparisonMulti(p1, p2, p, n1, n2)
    
#Confidence Intervals
if form == 13:
    
    print("sample mean =")
    xbar = float(input())
    print("sample size =")
    n = int(input())
    print("sample standard deviation = ")
    s = float(input())
    print("confidence level? (0-9) 50% 60% 70% 80% 90% 95% 98% 99% 99.8% 99.9%")
    z = [0.674, 0.842, 1.036, 1.282, 1.645, 1.960, 2.326, 2.576, 3.090, 3.291,]
    i = int(input())
    
    CI(xbar, n, s, z[i])
    
#Z_Value
if form == 14:
    
    print("sample mean =")
    xbar = float(input())
    print("sample size =")
    n = int(input())
    print("sample standard deviation = ")
    s = float(input())
    print("Confidance interval (smallest value)")
    CI = float(input())
    
    Z_Value(xbar, CI, s, n)
    print("Look at the Z-Value table for the confidence interval")
    
#Adjustment Chi
if form == 15:

    print("Expected values (-1when done)")
    
    Ei = []
    
    i = 0
    while True :
        Ei.append(float(input()))
        
        if Ei[i] < 0 :
            Ei.pop()
            break;
        
        else :
            i+=1

    print("Observed values (-1when done)")
    
    Oi = []
    
    i = 0
    while True :
        Oi.append(float(input()))
        
        if Oi[i] < 0 :
            Oi.pop()
            break;
        
        else :
            i+=1

    Adjustment(Ei, Oi)

#Ei
if form == 16:
    
    print("N = ")
    n = int(input())
    print("Pi values as decimal (-1when done)")
    
    Pi = []
    
    i = 0
    while True :
        Pi.append(float(input()))
        
        if Pi[i] < 0 :
            Pi.pop()
            break;
        
        else :
            i+=1
    
    Ei(n, Pi)