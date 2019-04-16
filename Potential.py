import numpy as np

def Ep(R, E0, r0):
	"""
	array*double*double -> double
	Ep: Calcul l'énergie potentielle du système
	hypothèse: R est un array de taille 1*((n*(n-1))/2) n étant le nombre d'atome
	dans le système & E0 > 0 & r0 > 0
	"""
	V = Vp(R, r0);
	Ei = 4*E0*V*(V-1);
	return(np.sum(Ei));

def dEp(R, E0, r0):
	"""
	array*double*double -> double
	Ep: Calcul de la dérivée de l'énergie potentielle du système
	hypothèse: R est un array de taille 1*((n*(n-1))/2) n étant le nombre d'atome
	dans le système & E0 > 0 & r0 > 0
	"""
	V = Vp(R, r0);
	Ei = 12*E0*V*(1-V)/(np.sqrt(R));
	return(Ei);

def Amov(n,E0,r0, CM):
	"""
	int*double*double -> array
	Amov: mouvement des atomes afin de minimiser l'énergie
	hypothèse: E0 > 0 & r0 > 0 & n > 0
	"""
	Mov = np.zeros((n-1,3));
	R = CartoR(CM);
	dE = dEp(R, E0, r0);
	C = np.append(np.array([[0,0,0]]),CM, axis = 0);
	count = range(n-1);
	d = 0;
	f = 0;
	for i in count:
		d = f;
		m = n - (i+1);
		f = d + m;
		sweep = range(d,f);
		for j in sweep:
			Vr = (CM[i+(j%m)]-C[i]);
			norm = np.sqrt((Vr*Vr).sum())
			Ur = Vr/norm;
			Mov[i] += Ur*dE[j];
	return(Mov);

def Vp(R, r0):
	"""
	array*double -> array
	Vp: Calcul intermédiaire de Ep
	hypothèse: R est un array de taille 1*((n*(n-1))/2) n étant le nombre d'atome
	dans le système & r0 > 0
	"""
	return((r0**6)/(2*(R**3)));



def CartoR(CM):
	"""
	array -> array
	CartoR: transforme une matrice de l'espace cartésien en matrice des distances carée
	hypothèse: np.shape(CM) == (n,3)
	"""
	R = np.sum(CM*CM,1);
	n,_ = CM.shape;
	N = n-1;

	for i in range(N):
		A = CM[1+i:,:] - CM[i,:];
		R = np.append(R,np.sum(A*A,1));

	return(R);

def Map(Np,n,E0,r0,l):
	"""
	int*int*double*double*double -> tuple(array,array)
	clustering: Trouve les etats stables dans l'espace PIV
	hypothèse: Np > 0 & n > 0 & E0 > 0 & r0 > 0 & l > 0
	"""
	PIV = np.zeros((Np,int((n*(n-1))/2)));
	E = np.zeros((Np));
	Cart = np.zeros((Np,n-1,3));
	for i in range(Np):

		CM = l*np.random.rand(n-1,3) - l/2;
		R = CartoR(CM);
		Ei = Ep(R, E0, r0);
		PIV[i,:] = R;
		Cart[i] = CM;
		E[i] = Ei/n;

	m = E.min();

	Bool = E<(m+10);

	while not Bool.all():
		for i in range(Np):
			while E[i] > (m+10):
				CM = l*np.random.rand(n-1,3) -l/2;
				Cart[i] = CM;
				R = CartoR(CM);
				Ei = Ep(R, E0, r0);
				PIV[i,:] = R;
				E[i] = Ei/n;
			m = E.min();
			Bool = E<(m+10);
			pc = (1*Bool.sum()/Np)*100;
			print(pc, m, E[i]);
			if pc == 100:
				break;
	print((1*Bool.sum()/Np)*100, m, "last");
	np.sort(R);
	PIV = (1 - (np.sqrt(PIV)/(r0))**6)/(1-(np.sqrt(PIV)/(r0))**12);
	return((PIV,E,Cart));

def AIRSS(Np,n,E0,r0,l):
	"""
	int*int*double*double*double -> tuple(array,array)
	clustering: Trouve les etats stables dans l'espace PIV
	hypothèse: Np > 0 & n > 0 & E0 > 0 & r0 > 0 & l > 0
	"""
	PIV = np.zeros((Np,int((n*(n-1))/2)));
	E = np.zeros((Np));
	Cart = np.zeros((Np,n-1,3));

	for i in range(Np):
		CM = l*np.random.rand(n-1,3) - l/2;


	print((1*Bool.sum()/Np)*100, m, "last");
	PIV = 1/(1+(np.sqrt(PIV)/(r0*1.75))**6);

	return((PIV,E,Cart));
