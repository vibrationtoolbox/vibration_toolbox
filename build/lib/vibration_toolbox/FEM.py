import numpy as np
import matplotlib.pyplot as plt
import time
#import pdb
#pdb.set_trace()

__all__ = ['']

def FEM_pre():
    """Finite Element Method - Preprocessor.

        Returns t, x, v, zeta, omega, omega_d and A resulting from the
        free response of a second order linear ordinary differential
        equation defined by
        :math:`m\ddot{x} + c \dot{x} + k x = 0`
        given initial conditions :math:`x_0` and :math:`\dot{x}_0 = v_0` for
        :math:`0 < t < t_{max}`

        Parameters
        ----------
       node :  floats, optional
            element nodes
        x0, v0:  floats, optional
            initial printlacement, initial velocity
        max_time: float, optional
             time for :math:`x(t)`
    """

# Defining nodes and node numbers
x,y = input("Enter x and y location of node. (eg.: x y)").split()
x = float(x)
y = float(y)
loc = [x, y]
node = np.array(loc)
for i in range(1, 1000):
    x,y = input('Enter x and y location of node (ie. x y) or (-1 -1) to . ').split()
    x = float(x)
    y = float(y)
    loc = [x, y]
    if loc[0]<0.0 and loc[1]<0.0:
        break
    else:
        node = np.vstack((node,loc))
        #node = np.app([node], [np.array(loc)],axis=0)
        L = max([max(node[:,0])-min(node[:, 0]),max(node[:,1])-min(node[:, 1])])
        d = node
        xlo = np.min([node[:, 0], d[:, 0]])
        xho = np.max([node[:, 0], d[:, 0]])
        ylo = np.min([node[:, 1], d[:, 1]])
        yho = np.max([node[:, 1], d[:, 1]])
        xsp = xho - xlo
        ysp = yho - ylo
        if .618 * xsp > ysp:
            xh = xho + .1 * xsp
            xl = xlo - .1 * xsp
            yh = (yho + ylo) / 2 + (xh - xl) * .618 / 2
            yl = (yho + ylo) / 2 - (xh - xl) * .618 / 2
        else:
            yh = yho + .1 * ysp
            yl = ylo - .1 * ysp
            xh = (xho + xlo) / 2 + .5 * (yh - yl) / .618
            xl = (xho + xlo) / 2 - .5 * (yh - yl) / .618

        dx = .01 * (xh - xl)
        dy = .035 * (yh - yl)
        #fig = plt.figure()
        plt.scatter(node[:, 0], node[:, 1], marker='*', color='b')
        sz = np.shape(node)
        m = sz[0]
        n = sz[1]
    for j in range(0, m):
        print(j+1)
    an = plt.annotate(j+1,(node[j,0]+dx,node[j,1]+dy))
    #plt.annotate(1, (node[0, 0] + dx, node[0, 1] + dy))
    plt.axis('square')
    plt.axis([xl, xh, yl, yh])
    plt.axis('image')
# Connecting nodes
point = 'y'
print("Pick nodes to connect with elements.")
#if point=='y':
#    print("Use pointing device or arrow keys and return")
#else:
#    print("Enter node numbers one at a time")
time.sleep(1)

for i in range(0, 1000):
    if i != 0:
        print("Pick nodes to connect with elements.")

    if point == 'y':
        n1 = plt.ginput(1)
        pts = np.array(n1)
        x1 = pts[:, 0]
        y1 = pts[:, 1]
        dis = (node[:,0]-x1)**2+(node[:,1]-y1)**2
        nodenum1 = np.argmin(dis)
    else:
        nodenum1 = input("Enter node number 1: ")
        nodenum1 = int(nodenum1)-1

    plt.plot(node[nodenum1, 0], node[nodenum1, 1], marker='o', color='k')

    if point == 'y':
        n2 = plt.ginput(1)
        pts = np.array(n2)
        x2 = pts[:, 0]
        y2 = pts[:, 1]
        dis = (node[:,0]-x2)**2 + (node[:,1]-y2)** 2
        nodenum2 = np.argmin(dis)
    else:
        nodenum2 = input("Enter node number 2: ")
        nodenum2 = int(nodenum2)-1

    plt.plot(node[nodenum2, 0], node[nodenum2, 1], marker='o', color='k')
    plt.plot([node[nodenum1, 0],node[nodenum2, 0]], [node[nodenum1, 1],node[nodenum2, 1]], '-b')

    answer2 = 'n'
    if i > 0:
        answer2 = input("Same properties as previous element? (y/n)")

    ncon = [0,0,0,0,0,0,0]
    if answer2 == 'n':
        E = float(input(" Enter the modulus of elasticity of the member. "))
        G = float(input(" Enter the shear modulus of the member (zero for EB beam). "))
        I = float(input(" Enter the moment of area of the member. "))
        A = float(input(" Enter the cross sectional area of the member. "))
        Rho = float(input(" Enter the density per unit length of the member. "))

    ncon = np.vstack((ncon, [np.array([nodenum1 + 1, nodenum2 + 1, E, A, I, G, Rho])]))
    answer = input(" Enter another element? (y/n) ")
    if answer != 'y':
        break

ncon = np.delete(ncon, (0), axis=0) # Delete first rows of zeros

#Adding concentrated masses and inertias
conm = [0, 0, 0]
for i in range(0, 1000):
    if i == 0:
        answer = input(" Add concentrated masses and rotational inertias? (y/n) ")

    if i > 0:
        answer = input(" Add more concentrated masses and rotational inertias? (y/n) ")

    if answer != 'y':
        break

    print(" ")
    print(" Pick node to add mass/rotational inertia to.")
    time.sleep(.5)
    if point == 'y':
        n1 = plt.ginput(1)
        pts = np.array(n1)
        x1 = pts[:, 0]
        y1 = pts[:, 1]
        dis = (node[:, 0] - x1) ** 2 + (node[:, 1] - y1) ** 2
        nodenum = np.argmin(dis)
    else:
        nodenum = int(input(" Enter node number: "))-1

    plt.plot(node[nodenum, 0], node[nodenum, 1], '*w')
    plt.plot(node[nodenum, 0], node[nodenum, 1], 'xr')
    answer = input(" Add mass or rotational inertia?(m,i,n(one)) ")
    massval = int(input(" Enter magnitude of mass/inertia. "))
    coni = np.array([0, 0, 0])
    coni[0] = nodenum

    if answer == 'm':
        coni[1] = massval

    if answer == 'i':
        coni[2] = massval

    conm = np.vstack((conm, coni))

    plt.plot(node[nodenum, 0], node[nodenum, 1], 'xi')
    plt.plot(node[nodenum, 0], node[nodenum, 1], '*b')

conm =  np.delete(conm, (0), axis=0)

#Zeroing of displacements
zero = [0, 0]
ij = 0
for i in range(0, 1000):
    if i == 0:
        answer=input(" Add boundary conditions? (y/n) ")

    if i > 0:
        answer=input(" Zero another displacement? (y/n) ")

    if answer != 'y':
        break

    print(" ")
    print(" Pick node to zero ")
    time.sleep(1.)

    if point=='y':
        n1 = plt.ginput(1)
        pts = np.array(n1)
        x1 = pts[:, 0]
        y1 = pts[:, 1]
        dis = (node[:, 0] - x1) ** 2 + (node[:, 1] - y1) ** 2
        nodenum = np.argmin(dis)
    else:
        nodenum = int(input(" Enter node number: ")) - 1

    plt.plot(node[nodenum, 0], node[nodenum, 1],'*w')
    plt.plot(node[nodenum, 0], node[nodenum, 1],'xr')

    print(" ")
    print(" ")
    print(" Zero which displacement(s)?(x,y,t(heta),n(one)) ")
    print("(ie xt for x and theta)")
    answern=input(" ")
    sanswern=len(answern)
    for ii in range (0, sanswern):
        ij = ij+1
        answer=answern[ii]

    if answer == 'x':
        plt.plot(node[nodenum, 0], node[nodenum, 1],'*b')
        zero = np.vstack((zero, np.array([nodenum,0])))

    if answer=='y':
        plt.plot(node[nodenum, 0], node[nodenum, 1],'*b')
        zero = np.vstack((zero, np.array([nodenum, 1])))

    if answer=='t':
        plt.plot(node[nodenum, 0], node[nodenum, 1],'*b')
        zero = np.vstack((zero, np.array([nodenum, 2])))

    if answer=='n':
        plt.plot(node[nodenum, 0], node[nodenum, 1],'*b')

zero =  np.delete(zero, (0), axis=0)

force=[0, 0, 0]
#Adding loads
for i in range(1, 1000):
    answer=input(" Add loads? (y/n) ")
    if answer!='y':
        break

    print(" ")
    print(" Pick node to load")
    time.sleep(1.5)

    if point=='y':
        n1 = plt.ginput(1)
        pts = np.array(n1)
        x1 = pts[:, 0]
        y1 = pts[:, 1]
        dis = (node[:, 0] - x1) ** 2 + (node[:, 1] - y1) ** 2
        nodenum = np.argmin(dis)
    else:
        nodenum = int(input(" Enter node number: ")) - 1

    plt.plot(node[nodenum, 0], node[nodenum, 1],'*w')
    plt.plot(node[nodenum, 0], node[nodenum, 1],'xr')
    answer=input("Load which displacement?(x,y,t(heta),n(one)) ")

    if answer=='x':
        loadval=input(" Enter magnitude of load. ")
        force = np.vstack((force, np.array([nodenum, 1, loadval])))

    if answer=='y':
        loadval=input(" Enter magnitude of load. ")
        force = np.vstack((force, np.array([nodenum, 2, loadval])))

    if answer=='t':
        loadval=input(" Enter magnitude of load. ")
        force = np.vstack((force, np.array([nodenum, 3, loadval])))

    answer=input(" Load another node? (y/n) ")

    if answer != 'y':
        break

    plt.plot(node[nodenum, 0], node[nodenum, 1],'*b')

force =  np.delete(force, (0), axis=0)

answer=input("Save configuration file (Else all will have been in vain)? (y/n) ")

#if answer == 'y':
    # Create .con file in Pyton