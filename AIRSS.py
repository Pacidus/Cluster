import matplotlib.pyplot as plt
import numpy as np
from Potential import AIRSS, Ep, CartoR
from Mail import send

n = 13;
r0 = 2;
E0 = 25;
l = n/3;
Np = 1000;
delta = 5*10**-5

Cart = AIRSS(Np,n,E0,r0,l, delta);

namefile = "Data/ACart delta%f n%d r0%f l%f E0%f Np%d.txt" %(delta, n,r0,l,E0,Np)

Molec = open(namefile, 'wb');
Molec.write('# Array shape: {0}\n'.format(Cart.shape).encode('utf-8'));
for molecule in Cart:
	Molec.write('# New molecule\n'.encode('utf-8'));
	np.savetxt(Molec,molecule, fmt='%-12.8f');
Molec.close();



namefile = "Data/ACart n%d r0%f l%f E0%f Np%d.xyz" %(n,r0,l,E0,Np)

Molec = open(namefile, 'w');
for i in range(Np):
	R = CartoR(Cart[i]);
	Molec.write(str(n));
	Molec.write("\n");
	Molec.write("E = {0}\n".format(Ep(R, E0, r0)));
	Molec.write("C 0 0 0\n");
	for j in range(n-1):
		Molec.write("C {0} {1} {2}\n".format(*Cart[i,j]));
Molec.close();

send(n,r0,l,E0,Np);

