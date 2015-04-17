#(C) Lee Carraher. All rights reserved.
#blah, blah, I'm partial to BSD distribution license

from pylab import *
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D



def parse(filename,fnc):
    #<-----------parsers---------------------->
    ptsmap = {}

    fl = file(filename,'r')
    l = fl.readline()
    l = fl.readline()
    time = 0
    while not l =="" :
        f = l.split(',')
        
        i=2
        team1=f[i]
        i=i+1
        while not str.isdigit(f[i]):
            team1=team1+' '+f[i]
            i=i+1    
        score1 = int(f[i])
        
        i=16+i
        team2 = f[i]
        i=i+1
        
        while not str.isdigit(f[i]):
            team2=team2+' '+f[i]
            i=i+1
        score2 = int(f[i])
        
        if ptsmap.has_key(team1):
            ptsmap[team1][team2] = fnc(score1,score2,time)
        else: ptsmap[team1] = {team2:fnc(score1,score2,time)}

        if ptsmap.has_key(team2):
            ptsmap[team2][team1] = fnc(score2,score1,time)
        else: ptsmap[team2] = {team1:fnc(score2,score1,time)}
        
        time = time+1

        l = fl.readline()
        
        
    teams = ptsmap.keys()

    m = len(teams)
    data = [[0]*m for i in xrange(m)]

    for i in xrange(m):
        team = teams[i]
        for j in xrange(m):
            if ptsmap[team].has_key(teams[j]):
                data[i][j] = ptsmap[team][teams[j]]

    data = array(data)
    return data,teams
    
#<-----------Maths-----------> 
# Directional distance metric
# just square the differences for euclidean distance
def distance(X,Y=[],d=0):
    if d == 0: d = len(X)
    if len(Y)==0:Y=[0.0]*d
    return sum ([ (X[j]-Y[j])**2. for j in range(d)] )    
    
def rank(data,teams,reduceDegree):

    teamAndWL = {}
    U, s, Vh = linalg.svd(data,full_matrices=True)
    #S = diag(s)
    #truncate the svd
    #U = U[:,0:reduceDegree] #orthogonal "Hanger" matrix
    #S = S[0:reduceDegree,0:reduceDegree] #the singular value matrix "stretcher"
    V =  array([  Vh[k,:]*s[k] for k in range(reduceDegree)]) #aligner matrix
    #dataHat = dot(dot(U,S),V) #estimate of the data

    S = []
    for i in xrange(len(teams)):
        w = 0
        l = 0
        for j in xrange(len(teams)):
            if data[i][j]>0: w=w+1
            if data[i][j]<0: l = l+1
        teamAndWL[teams[i]] = (w,l)
        
        # rank on w/l weighted stable distributions of the random-surfer model
        #useful for giving overall ranking a direction (otherwise just being 
        #abnormal is rewarded, uhmm 'Grambling St.')
        S.append((distance(V[:,i])*(float(w)/l),teams[i]+"("+str(w)+':'+str(l)+") "))
        
        # rank on unweighted stable distributions of the random-surfer model
        #S.append((distance(V[:,i]),teams[i]+"("+str(w)+':'+str(l)+") "))
        #rank just on win loss records
        #S.append((float(w)/float(l),teams[i]+"("+str(w)+':'+str(l)+") "))

    S = sorted(S)
    ret = {}
    
    
    for i in range(len(S)):
        l = S[i]
        ret[l[1]]=i
        #print l[1] + ": "+str(l[0])
    return ret,V,teamAndWL



#<---------------Main------------>
filename = '2013.csv'


#you can of course define your own distance algorithms
def fnc(x,y,time):
    
    #if x-y>0:return 1.0*(time+1)**.5
    #return -1.0*(time+1)**.5
    return float(x-y)*log(time+1)
        
data,teams = parse(filename,fnc)


ret,V,t = rank(data,teams,1)
#output ranking along the principle eigenvector
#for r in ret:
#    print ret


#since ranking is an ill-posed problem, lets look at some
#rankings among other eigenvectors, and maintain a list of
#team rankings along these vectors, then aggregate their
#ranking averages
fullData = {}
for team in ret.keys():
    fullData[team]=[ret[team]]

for i in range(2,30):
    print i
    ret,V,t = rank(data,teams,i)
    for team in ret.keys():
        fullData[team].append(ret[team])

S = []
for team in fullData.keys():
    avg = sum(fullData[team])/float(len(fullData[team]))
    S.append((avg,team))
    
S=sorted(S)

for i in range(len(S)):
    l = S[i]
    ret[l[1]]=i
    print l[1] + ": "+str(sum(fullData[l[1]])/len(fullData[l[1]]))
        


#<other distance and filtering methods---------->

def cosineDistance(X,Y,d=3):
    sumtop = 0.0
    sumA = 0.0
    sumB = 0.0
    for i in range(d):
        sumtop = (Y[i])*X[i] + sumtop
        sumA =   (Y[i])*(Y[i]) +sumA
        sumB = (X[i])*(X[i]) + sumB  
    return sumtop/(((sumA)**.5)*((sumB)**.5))
 

def filterData(X,d=3,k=.1):
    keep = []
    cutoffVector = [0.0]*d
    for i in range(d):
        s = float(sum(X[i]))/float(len(X[i]))
        var = sum( map(lambda x: (x-s)*(x-s), X[i]))/(len(X[i])-1)
        cutoffVector[i] = s + var*k
    for j in xrange(len(X)):
        if distance(X[j][0:d],cutoffVector,d)>0.0 :
            keep.append(j)
    return keep

def normalize(X,d=3):
    for i in range(d):
        s = float(sum(X[i]))/float(len(X[i]))
        var = sum( map(lambda x: (x-s)*(x-s), X[i]))/(len(X[i])-1)
        X[i] = map(lambda x: (x-s)/var, X[i])
    return X
    
    
#<-----------------3d Stuff--------------------->
#comment out if it crashes/you just want rankings
from mpl_toolkits.mplot3d import Axes3D
cm = plt.get_cmap("RdYlGn")
fig = plt.figure(1)
ax = Axes3D(fig)
col = []
    
#keep = filterData(data,8)
#data = [data[i] for i in keep]
#teams = [teams[i] for i in keep]
ret,V,t = rank(data,teams,10)
dims = [1,2,3]

plt.rcParams.update({'font.size': 8})

for i in xrange(len(teams)):
    col.append(log(1+float(t[teams[i]][0])/t[teams[i]][1])) 
    if float(t[teams[i]][0])/t[teams[i]][1]>3.:
        ax.text(V[dims[0]][i],V[dims[1]][i],V[dims[2]][i],teams[i])
        for d in xrange(len(data[i])):
            if not data[i][d]==0 and float(t[teams[d]][0])/t[teams[d]][1]>3.:
                plt.plot([V[dims[0]][i],V[dims[0]][d]],[V[dims[1]][i],V[dims[1]][d]],[V[dims[2]][i],V[dims[2]][d]],color="grey")             
        
ax.scatter(V[dims[0]],V[dims[1]],V[dims[2]],s=15, c=col, marker=',')

ax.set_xlabel('Component 1')
ax.set_ylabel('Component 2')
ax.set_zlabel('Component 3')
plt.show()


