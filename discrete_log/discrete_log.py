import primefac
import gmpy

def chinese_remainder(x_list,modules_list):
	"""
	just solves such system: x=x_i mod n
	"""
	M=1
	summ=0
	for i in range(len(modules_list)):
		M*=modules_list[i]
	for i in range(len(modules_list)):
		M_i=M/modules_list[i]
		M_i_1=gmpy.invert(M_i,modules_list[i])
		summ+=x_list[i]*M_i*M_i_1
	return summ%M
def baby_step_giant_step(g,b,p):
	"""
	This function implements Baby step - giant step algorithm for solving discrete logarithm problem ( works good only for prime module, that's the reason of implementing Pohlig_Hellman. g^x=b mod p
	"""
	m=gmpy.sqrt(p)+1
	table=dict()
	for j in range(m+1):
		table[pow(g,j,p)]=j
	mult=gmpy.invert(pow(g,m,p),p)
	gamma=b
	for i in range(m+1):
		if(gamma in table.keys()):
			return i*m+table[gamma]
		else:
			gamma*=mult
	print(g,b,p)
	return -1

def prime_power_groups(g,h,p,e):
	"""
	This function is a subroutine of Pohlig-Hellman algorithm for solving g^x= mod p^e where p is prime
	"""
	x=0
	n=pow(p,e)
	gamma=pow(g,pow(p,e-1),p)
	for k in range(0,e):
		gx=pow(g,x,n)
		h_k=pow(gmpy.invert(gx,p)*h,pow(p,e-1-k),p)
		d_k=baby_step_giant_step(gamma,h_k,p)
		x+=(d_k*pow(p,k))
	return x

def Pohlig_Hellman(divs,g,h,p):
	"""
	This function implements Pohlig-Hellman algorithm for solving discrete logarithm problem 
	divs must be a dict of divisors with their order [(p_i,e_i)]
	"""
	modules_list=[pow(div,divs[div]) for div in divs.keys()]
	x_list=[]
	for div in divs.keys():
		g_i=pow(g,p/pow(div,divs[div]),p)
		h_i=pow(h,p/pow(div,divs[div]),p)
		x_i=prime_power_groups(g,h,div,divs[div])
		x_list.append(x_i)
	x=chinese_remainder(x_list,modules_list)
	return x

def find_discrete_log(g,h,p):
	"""
	Tries to find discrete logarithm for g^x=h mod p
	chooses optimal algorithm
	"""
	d=primefac.factorint(p)
	if(len(d.keys())==1):
		for div in d.keys():
			if(d[div]==1):
				return baby_step_giant_step(g,h,p)
			else:
				return prime_power_groups(g,h,div,d[div])
			return -1
	else:
		return Pohlig_Hellman(d,g,h,p)
print find_discrete_log(5,7,18)
