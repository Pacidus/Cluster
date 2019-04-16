from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt
import tkinter as tk
from tkinter.filedialog import askopenfilename
import numpy as np

def ask():
	root = tk.Tk();
	root.withdraw();

	print("Donées des Énergies :");
	Ref = askopenfilename();
	print(Ref)

	return(Ref);

Ref = ask();
Nam = Ref.split("/")[-1];
Np = Nam.split(" ")[-1];
Np = Np.partition("Np")[-1];
Np = int(Np.partition(".")[0]);
n = Nam.split(" ")[1];
n = int(n.partition("n")[-1]);
Name = "./Cart" + Nam.partition("E")[-1];

E = np.loadtxt(Ref);
Cart = np.loadtxt(Name);

C = Cart.reshape((Np,n-1,3))[E==E.min()][0].T

fig = plt.figure();
ax = fig.add_subplot(111, projection='3d');

x = np.zeros((n));
y = np.zeros((n));
z = np.zeros((n));


x[1:] = C[0];
y[1:] = C[1];
z[1:] = C[2];

print(C);

ax.scatter(x,y,z,'r.', label = "E = %.9f" %(E.min()) );

ax.set_xlabel('X');
ax.set_ylabel('Y');
ax.set_zlabel('Z');

plt.legend();
plt.show();
