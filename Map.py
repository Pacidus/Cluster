import matplotlib.pyplot as plt
import numpy as np
from Potential import Map
from Mail import send

n = 10;
r0 = 2;
E0 = 25;
l = r0*(n**(1/3));
Np = 40;

PIV , E, Cart = Map(Np,n,E0,r0,l);

Bool = (E == E.min());
PIV0 = PIV - PIV[Bool];

Enorm = E - E.min();

D = np.sqrt((PIV0*PIV0).sum(1));

np.savetxt("Data/PIV n%d r0%f l%f E0%f Np%d.txt" %(n,r0,l,E0,Np),PIV);
np.savetxt("Data/E n%d r0%f l%f E0%f Np%d.txt" %(n,r0,l,E0,Np),E);

namefile = "Data/Cart n%d r0%f l%f E0%f Np%d.txt" %(n,r0,l,E0,Np)

Molec = open(namefile, 'wb');
Molec.write('# Array shape: {0}\n'.format(Cart.shape).encode('utf-8'));
for molecule in Cart:
	Molec.write('# New molecule\n'.encode('utf-8'));
	np.savetxt(Molec,molecule, fmt='%-12.8f');
Molec.close();



namefile = "Data/Cart n%d r0%f l%f E0%f Np%d.xyz" %(n,r0,l,E0,Np)

Molec = open(namefile, 'w');
for i in range(Np):
	Molec.write(str(n));
	Molec.write("\n");
	Molec.write("E = {0}\n".format(E[i]));
	Molec.write("C 0 0 0\n");
	for j in range(n-1):
		Molec.write("C {0} {1} {2}\n".format(*Cart[i,j]));
Molec.close();


plt.plot(D,Enorm,".");
plt.grid(True);

send(n,r0,l,E0,Np);

plt.show();
