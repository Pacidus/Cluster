import matplotlib.pyplot as plt
import tkinter as tk
from tkinter.filedialog import askopenfilename
import numpy as np
from Potential import CartoR, Ep

def ask():
	root = tk.Tk();
	root.withdraw();

	print("Donées des positions AIRSS :");
	Ref = askopenfilename();
	print(Ref)

	return(Ref);

Ref = ask();
Nam = Ref.split("/")[-1];

if(Nam[0] == 'A'):
	_,_,n,r0,_,E0,Np = Nam.split(" ");
else:
	_,n,r0,_,E0,Np = Nam.split(" ");

n = int(n.partition("n")[-1]);
r0 = float(r0.partition("r")[-1]);
E0 = float(E0.partition("E0")[-1]);
Np = Np.partition("Np")[-1];
Np = int(Np.partition(".")[0]);

Molec = np.loadtxt(Ref);
Molec = Molec.reshape((Np,n-1,3));



PIV = np.zeros((Np,int((n*(n-1))/2)));
E = np.zeros((Np));

for i in range(Np):
		R = CartoR(Molec[i]);
		np.sort(np.sqrt(R));
		PIV[i,:] = np.sort(np.sqrt(R));
		E[i] = Ep(R, E0, r0);

PIV = 1/(1+(PIV/(r0))**6);

Bool = (E == E.min());
PIV0 = PIV - PIV[Bool];

Enorm = E - E.min();

print(E.min())

D = np.sqrt((PIV0*PIV0).sum(1));

plt.grid(True);
plt.plot(D,Enorm,".");
#plt.plot(PIV0[:,0],PIV0[:,1]);


plt.show();
