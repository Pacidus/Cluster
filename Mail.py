# Import smtplib for the actual sending function
import smtplib

# Import the email modules we'll need
from email.message import EmailMessage

def send(n,r0,l,E0,Np):
	"""
	int*int*double*double*double -> bool
	envoie un e-mail
	"""
	datafile = open("../ID.go", "r")
	username = datafile.readline();
	password = datafile.readline();


	fromaddr = username;
	toaddrs  = username;

	msg = """From: %s\nTo: %s\nSubject: %s\n\n%s\n
	 n = %d \n r0 = %f \n l = %f \n E0 = %f \n Np = %d
	""" % (fromaddr, toaddrs, "Simulations Stage", "Simulation terminee", n,r0,l,E0,Np)
	try:
		server = smtplib.SMTP('smtp.gmail.com:587')
		server.starttls()
		server.login(username,password)
		server.sendmail(fromaddr, toaddrs, msg)
		server.quit()
		return(True);
	except:
		return(False);
