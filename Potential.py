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
	Ei = 12*E0*(V*2)*(1-(V*2))/(np.sqrt(R));
	return(Ei);

def Vp(R, r0):
	"""
	array*double -> array
	Vp: Calcul intermédiaire de Ep
	hypothèse: R est un array de taille 1*((n*(n-1))/2) n étant le nombre d'atome
	dans le système & r0 > 0
	"""
	return((r0**6)/((R**3)*2));

def Amov(n,E0,r0, CM):
	"""
	int*double*double -> array
	Amov: mouvement des atomes afin de minimiser l'énergie
	hypothèse: E0 > 0 & r0 > 0 & n > 0
	"""
	Mov = np.zeros((n-1,3));
	C = np.append(np.array([[0,0,0]]),CM, axis = 0);
	i = 1;
	while i < n:
		j = 0;
		while j < n:
			if i != j:
				Vr = C[i]-C[j];
				R = (Vr*Vr).sum();
				dE = dEp(R, E0, r0);
				Ur = Vr/R;
				Mov[i-1] += Ur*dE;
			j += 1;
		i += 1;
	return(Mov);

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
	PIV = 1/(1+(np.sqrt(PIV)/(r0*1.75))**6);
	return((PIV,E,Cart));

def AIRSS(Np,n,E0,r0,l,delta):
	"""
	int*int*double*double*double -> tuple(array,array)
	clustering: Trouve les etats stables dans l'espace PIV
	hypothèse: Np > 0 & n > 0 & E0 > 0 & r0 > 0 & l > 0
	"""
	PIV = np.zeros((Np,int((n*(n-1))/2)));
	E = np.zeros((Np));
	Cart = l*np.random.rand(Np,n-1,3) - l/2;
	Moln = range(Np);
	itt = 100000;
	atom = range(n-1);
	d2 = delta*delta;
	pas = 1;
	p = .985;
	for j in Moln:
		CM = Cart[j,::,::];
		i = 0;
		pas = 0.5;
		while i < itt:
			i += 1;
			mov = Amov(n,E0,r0,CM);
			norm = np.sqrt((mov*mov).sum(1));
			for k in atom:
				if(norm[k]*norm[k] > d2):
					U = mov[k]/norm[k];
					CM[k] = CM[k] - pas*U;
					pas = pas*p + 10**-6;
		Cart[j,::,::] = CM;
		print(100*(j+1)/Np, i);
	return(Cart);
