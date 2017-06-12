"""
File name  : data_handler.py
Description:
Author     : Alioune Schurz
Lab        : Computational Molecular Design and Detection Lab (CMDD lab)
Institute  : National Taiwan University
"""

import sys,os,pybel,itertools
from openbabel import OBAtomAtomIter,OBBuilder
from numpy import *
from operator import mul

# list of reserved characters used in SMILE encoding of molecules (will be used to process the data)
smiles_char=['(',')','=','*','[',']','#','+','-',"@"]+[str(n) for n in range(0,9)]

def is_hydrogen_sidechain(smile):
	'''
	Will assert if the sidechain only contains hydrogenes.
	Args:
		string smile   : a SMILE string of a sidechain
	Returns:
		bool result    : the results
	'''
	if 'H' in smile:
		for char in smiles_char:
			#removing all reserved characters from the SMILE. Only atom names remain
			smile =smile.replace(char,'')
		# now we check that only H are present in the string, otherwise it is not a hydrogen side chain 
		ishydrogen = True
		for char in smile:
			if not char=='H':
				ishydrogen=False
		#--
		return ishydrogen


# * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * 
# * * * * * * * * * * * * * *   CSCCP   * * * * * * * * * * * * * * * * * * * * * * *

class TrivialCSCCP(Exception):
	def __init__(self,message):
		Exception.__init__(self,message)

class BadCSCCP(Exception):
	def __init__(self,message):
		Exception.__init__(self,message)

class CSCCP:
	''' This class represents the parameters to be given to a CSCCP solver. '''
	# The procedure to create an instance of this class is to first create a Scaffold object,
	# then load a scaffold structure into it with the method load_scaffold(), and finally use the method generate_csccp().This object is iterable.


	class Result:
		'''This class represents the parameters to be given to a CSCCP solver.'''
		def __init__(self,data):
			if not str(data.dtype)=="[('p', '<f8'), ('w', '<f8'), ('comp', 'O')]":
				print "CSCCP.Result.__init__(): Something is wrong with the dtype of the argument... terminating"
				exit()
			self.data=data
	
		def __iter__(self):
			for x in self.data:
				yield x
	
		def __len__(self):
			return len(self.data)
	
		def __str__(self):
			return str(self.data)
		
		def __getitem__(self,index):
			return self.data[index]
	
	def __init__(self):
		# structure; pybel.Molecule: pybel object representing a molecule (it is wrapper of open-babel.OBMol) 
		self.structure=0
		self.scaffold_name=""
		# n:
		self.n=0
		# R:
		self.R=10
		# ppm:
		self.ppm=5
		# mass_peak        : corresponds to w0
		mass_peak          =0
		# mass_peak_min    :
		self.mass_peak_min =0
		# mass_peak_max    :
		self.mass_peak_max =0
		
		# scaffold_weight_rel2config         : weight of the scaffold and hydrogen atoms at all hydrogen positions
		self.scaffold_weight_rel2config      =0 
		# scaffold_probability_rel2config    : cumulative probability of hydrogen atoms at all hydrogen positions
		self.scaffold_probability_rel2config =1.0
		# min_possible_weight, max_possible_weight: min and max weight the scaffold for all possible configurations and combinations  
		self.min_possible_weight       =0
		self.max_possible_weight       =0

		# wmin: minimum weight possible for any feasible solution (depends on the mass peak and resolution of the mass spectrometer)
		self.wmin           =0
		# wmax: maximum weight possible for any feasible solution (depends on the mass peak and resolution of the mass spectrometer)
		self.wmax           =0

		# W    : list of side chains weights 
		self.W =[] 
		# P    : list of non hydrogene side chains probabilities
		self.P =[]
		# K    : list of the number of sidechains at each non hydrogen position
		self.K =[]
		# sidechain_smiles    : list of SMILES of the non hydrogen side chains at each position
		self.sidechain_smiles =[]
		self.configIdx = 0
		# configuration    : position of the sidechains in the scaffold
		self.configuration=[]
		# number_of_compounds    : total number of compounds that can be generated if brute force is used
		self.number_of_compounds=0

	def __iter__(self):
		'''This object is made iterable so we can go trough all possible compounds for this scaffold.'''
		# The choice of having an iterator rather than generating a list of combinations is because of memory
		# limitations. There is no guaranty that the program won't stop because the memory isn't enough.
		# This will is used for all the brute force algorithms.

		#Example:
		#	for compound in csccp:
		#		print compound
		
		X=[range(0,k) for k in self.K]
		for compound in itertools.product(*X):
			yield list(compound)

	def __str__(self):
		res  ="scaffold: "+str(self.scaffold_name)
		res +="\nn: "+str(self.n)
		res +="\nR: "+str(self.R)	
		res +="\nmass_peak (w0): "+str(self.mass_peak)
		res +="\nmass_peak_min:"+str(self.mass_peak_min)
		res +="\nmass_peak_max: "+str(self.mass_peak_max)
		res +="\nscaffold_weight_rel2config:" +str(self.scaffold_weight_rel2config)	
		res +="\nscaffold_probability_rel2config: "+str(self.scaffold_probability_rel2config)
		res +="\nmin_possible_weight: "+str(self.min_possible_weight)
		res +="\nmax_possible_weight: "+str(self.max_possible_weight)
		res +="\nwmin: "+str(self.wmin)
		res +="\nwmax: "+str(self.wmax)
		res +="\nconfiguration: "+str(self.configuration)
		res +="\nnumber of sidechains (K): "+str(self.K)
		res +="\nnumber of compounds: "+str(self.number_of_compounds)
		res +="\nW,P,sidechain_smiles won't be displayed..."
		return res

	def prepare_matrix_CL(self,offset=0,ranks=0,precision=1,roundup=False):
		'''
		Creates the C and L matrix required for DP algorithms
		Args:
			int offset   : the maxtrix will use w <= wmax+offset+ 0.5*roundup
			int ranks    : if set>0 the value will override default R parameter of CSCCP (used for IDP)
			int precision: equal to 10^dec. Indicates how many digits are wanted for integer representation
			int roundup  : indicates wmax + offset should be rounded up when rescaled to integer 
		Returns:
			tuple array C,L  : the results
		'''
		l=max([len(x) for x in self.W])
		wmax=self.rescale_2int(self.wmax,precision=precision,roundup=roundup)

		if not ranks:
			C =zeros((self.n,wmax+offset+1,self.R),float)   # [0 .. n-1] x [0 .. wmax] x [0 .. R-1] 
			L =zeros((self.n,wmax+offset+1,self.R,2),int)
		else:
			C =zeros((self.n,wmax+offset+1,ranks),float)   # [0 .. n-1] x [0 .. wmax] x [0 .. R-1] 
			L =zeros((self.n,wmax+offset+1,ranks,2),int)
		return C,L

	def prepare_W(self,precision=0):
		'''
		Creates the W matrix for DP algorithms by rescaling to integer 
		Args:
			int precision: equal to 10^dec. Indicates how many digits are wanted for integer representation
			int roundup  : indicates wmax + offset should be rounded up when rescaled to integer 
		Returns:
			array W    : the result
		'''
		l=max([len(x) for x in self.W])
		if not precision:
			W    =array([(x+[0]*l)[:l] for x in self.W])
		else:
			W    =array([([int(y*precision) for y in x]+[0]*l)[:l] for x in self.W])
		return W

	def prepare_P(self):
		'''
		Creates the P matrix for DP algorithms 
		Returns:
			array P    : the result
		'''
		l=max([len(x) for x in self.W])
		return array([(x+[0]*l)[:l] for x in self.P])

	def prepare_K(self):
		'''
		Creates the K array for DP algorithms 
		Returns:
			array K    : the result
		'''
		return array(self.K)

	def rescale_2int(self,w,precision,roundup=False):
		'''
		Rescales a quatity w to integer
		Args:
			float w      : the quantity to rescale
			int precision: equal to 10^dec. Indicates how many digits are wanted for integer representation
			int roundup  : indicates wmax + offset should be rounded up when rescaled to integer 
		Returns:
			int res  : the results
		'''
		w =int(w*precision+0.5*int(roundup))
		return w

	def trivial(self):
		return len(self.K)==0

	def trivial_solution(self):
		weight=self.scaffold_weight_rel2config
		# cleanup the molecule be removing all atoms with null atomic number
		mol=pybel.Molecule(self.structure)
		todelete=[]
		for atom in mol:
			if not atom.atomicnum:
				todelete.append(atom)
		for td in todelete:
			mol.OBMol.DeleteAtom(td.OBAtom)	
		if self.mass_peak_min<weight and weight<self.mass_peak_max:
			return [dict(probability=1,weight=weight,mol=mol)]
		else: 
			return []

	def generate_output_molecules(self,csccp_result,options={}):
		'''
		Generates output molecules of the csccp
		Args:
			CSCCPResult csccp_result: the result returned by running a CSCCP solver
		Returns:
			dict: a dictionary
		'''

		# subroutine to connect the scaffold and side chain together
		def connect_sidechain_to_scaffold(scaffold,sidechain_smile,s):
			"""
			Connects the scaffold and side chain at desired position s. The scaffold is modified by reference
			Args:
				pybel.molecule scaffold: the structure of the scaffold
				string sidechain_smile : the smile of the sidechain 
				int s                  : position in the scaffold where the sidechain should be connected
			"""

			#index of substituted positions start at 0 but molecule index starts at 1. Therefore pos=s+1
			pos=s+1 
			smile=sidechain_smile
			# cleanup the sidechain smile
			smile=smile.replace("[H]","")
			smile=smile.replace("()","")
			r = pybel.readstring("smi",smile)
			#--

			# identify the atom to connect. His index is stored in toconnect
			toconnect=1
			for atom in r:
				if atom.atomicnum==0:
					break
				toconnect+=1	
			# removing the indicator of connection [*]
			smile=smile.replace("[*]","")
			r = pybel.readstring("smi",smile)
			#--

			# connecting the molecules
			# now we have r and s and know the molecule atom idx that need to be connected
			builder = OBBuilder()
			# Disable charge perception if desired
			scaffold.OBMol.SetAutomaticPartialCharge(False)
			# remember the size of molecule
			n = scaffold.OBMol.NumAtoms()
			# add sidechain to scaffold
			scaffold.OBMol += r.OBMol
			# connect the fragments
			builder.Connect(scaffold.OBMol, pos, n+toconnect)
		# -- end of connect_sidechain_to_scaffold

		count=0
		final_result=[]
		for res in csccp_result:
			prob=res[0]
			weight=res[1]
			comp=res[2]
			# we extract the smiles of the selected sidechains in the "sidechains"
			pos=0
			sidechains_smi=[]
			for scIDx in list(comp):
				sidechains_smi.append(self.sidechain_smiles[pos][scIDx])
				pos+=1
			#--
			# connecting all the sidechains to the scaffold (SOME ERRORS WITH THE INDICES !!!)
			pos=0
			mol=pybel.Molecule(self.structure)
			for sidechain_smi in sidechains_smi:
				connect_sidechain_to_scaffold(mol,sidechain_smi,self.configuration[pos])
				pos+=1
			#--
			# cleanup the molecule be removing all atoms with null atomic number
			todelete=[]
			for atom in mol:
				if not atom.atomicnum:
					todelete.append(atom)
			for td in todelete:
				mol.OBMol.DeleteAtom(td.OBAtom)			
			#--

			# generate final result and reajusting weight and probability
			final_result.append(dict(probability=prob,weight=weight,mol=mol))
			count+=1
		return final_result

# * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * 
# * * * * * * * * * * * * * *   CSCCPResult   * * * * * * * * * * * * * * * * * * * *

class CSCCPResult:
	'''This class represents the parameters to be given to a CSCCP solver.'''

	def __init__(self,data):
		if not str(data.dtype)=="[('p', '<f8'), ('w', '<f8'), ('comp', 'O')]":
			print "CSCCPResult.__init__(): Something is wrong with the dtype of the argument... terminating"
			exit()
		self.data=data

	def __iter__(self):
		for x in self.data:
			yield x

	def __len__(self):
		return len(data)

	def __str__(self):
		return str(self.data)

# * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * 
# * * * * * * * * * * * * * *   Scaffold   * * * * * * * * * * * * * * * * * * * *

class ScaffoldNullPopularity(Exception):
	def __init__(self,message):
		Exception.__init__(self,message)

class ScaffoldCorrupted(Exception):	
	def __init__(self,message):
		Exception.__init__(self,message)

class Scaffold:

	'''This class represents a scaffold and all the related data. It can be loaded using a cIdx file using the load_data() method.'''
	
	def __init__(self):
		#------------ data directly loaded from NPSD without further processing ------------------------------
		# name: corresponds to the name of the corresponding sdf file (without extension)
		self.name           = ""
		self.popularity     =0
		# structure; pybel.Molecule: pybel object representing a molecule (it is wrapper of open-babel.OBMol) 
		self.structure      =0 
		# scaffoldLength: number of positions on the scaffold (hydrogen and non hydrogen positions all together) ???????????
		self.scaffoldLength =0		
		#
		# min_possible_weight, max_possible_weight: min and max weight the scaffold for all possible configurations and combinations  
		self.min_possible_weight       =0
		self.max_possible_weight       =0
		# scaffold_weight_plus_hydrogens : "all hydrogen sidechains_weight" 
		# weight of the scaffold after adding hydrogens at all position of the scaffold. 
		# this is the wight that is obtained when opening the scaffold .sdf file.
		# indeed, for open-babel, the hydrogens are implicit.
		self.scaffold_weight_plus_hydrogens =0
		#
		self.configurations =[]
		self.number_of_compounds=[]
		self.sidechainlist  =[]
		#------------ data obtained by processing -----------------------------------------------------------
		# scaffold_weight_rel2configs : "config_hydrogene sidechains_weights" 
		# Contains for each configuration the weight of the scaffold excluding only the non hydrogen side chains
		self.scaffold_weight_rel2configs       =[]
		# scaffold_probability_rel2configs : "config_hydrogene sidechains_scaffold_probabilities" 
		# Contains for each configuration the cumulative probability of all hydrogen side chains  (excluding only the non hydrogen side chains)
		self.scaffold_probability_rel2configs =[]
		
	
	def __str__(self):
		res="name: "+str(self.name)+"\npopularity: "+str(self.popularity)+"\nscaffold_weight_plus_hydrogens: "+str(self.scaffold_weight_plus_hydrogens)
		res+="\nmin_possible_weight: "+str(self.min_possible_weight)+"\nmax_possible_weight: "+str(self.max_possible_weight)
		res+="\nnumber of configurations: "+str(len(self.configurations))
		res+="\nmax number of compounds: "
		for idx,num in self.number_of_compounds:
			res+="%s:%s " % (idx,num)
		res+="\n"
		return res
	
	def load_data_mr(self,text):
		'''
		This method loads the scaffold data from a mcIdx file <scaffold name>.mcIdx.
		Args:
			string scaffoldName: the name of the scaffold
		Raises:
			IOError
			ScaffoldNullPopularity
			ScaffoldCorrupted
		'''		
		#------------ data directly loaded from NPSD without further processing ------------------------------
		separator="@$@$@"
	
		cIdxf,sdff=text.split(separator)
		# Processing the cIdx file 
		Idx=0
		configNumber=0
		sidechainBlock=0
		positionsNumber=0
		position=-1
		
		for line in cIdxf.splitlines():
			if Idx==0:
				a=line.strip().split()
				self.name       = a[0]
				self.popularity = int(a[1])
				if not self.popularity:
					break
				Idx+=1
				continue
				
			if Idx==1:
				a=line.strip().split()
				self.max_possible_weight = float(a[0])
				self.min_possible_weight = float(a[1])
				Idx+=1
				continue
				
			if Idx==2:
				configNumber=float(line.strip().split()[0])
				sidechainBlock=3+configNumber
				Idx+=1
				continue
				
			if Idx>=3 and Idx<sidechainBlock:
				a=line.strip().split()
				nonHydrogenPositions=int(a[0])
				self.configurations.append([int(x) for x in a[1:nonHydrogenPositions+1]])
				Idx+=1
				continue
				
			if Idx>=3 and Idx==sidechainBlock:
				self.scaffoldLength=int(line.split()[0])
				for i in range(0,self.scaffoldLength):
					self.sidechainlist.append([]) 
				Idx+=1	
				continue
				
			if Idx>=3 and Idx>sidechainBlock:
				a=line.strip().split('\t')
				if len(a)==2: #new position
					position         =int(a[0])
					sidechainsNumber =int(a[1])
				else:
					self.sidechainlist[position].append([float(a[4]),float(a[3]),a[0]])
				Idx+=1
				continue
				
			Idx+=1
		if [] in self.sidechainlist: 
			raise ScaffoldCorrupted("Warning: some positions are missing: %d/%d. The file is probably corrupted." %(position+1,self.scaffoldLength))
			
		if not self.popularity:
			raise ScaffoldNullPopularity("Warning: the popularity of this scaffold is 0.")

		# Loading the scaffold molecule (.sdf file) thanks to pybel
		self.structure  = pybel.readstring("sdf",sdff)
		self.scaffold_weight_plus_hydrogens = self.structure.molwt

		#------------ data obtained by processing -----------------------------------------------------------
		c=0
		for conf in self.configurations:
			# we add hydrogens to be able to count them. Nevertheless it doesn't affect the molecular weight which always includes the hydrogens
			self.structure.OBMol.AddHydrogens()

			# \we compute the total number of hydrogens
			total_number_of_hydrogenes     =0
			for atom in [a for a in self.structure if not a.OBAtom.GetAtomicNum()==1]:
				if atom.idx-1 in conf:
					number_of_hydrogenes       =len([x for x in OBAtomAtomIter(atom.OBAtom) if x.GetAtomicNum()==1])
					total_number_of_hydrogenes +=number_of_hydrogenes
			#

			# \Dealing with the mass 
			# we remove the hydrogens because we don't need them anymore 
			self.structure.OBMol.DeleteHydrogens()
			# we compute the total mass of hydrogens at non hydrogen side chain positions    
			total_hydrogene_mass = total_number_of_hydrogenes * 1.00794
			# we compute the mass of the scaffold without the non hydrogen positions
			hs_weight            =self.scaffold_weight_plus_hydrogens - total_hydrogene_mass
			self.scaffold_weight_rel2configs.append(hs_weight)
			# 


			# \Dealing with the probability 
			# the condition a.atomicnum>0 is to ensure that we don't count the "atoms" R# that are just used to number the positions on the scaffold
			hs_probabilities=[1]
			for atom in [ a for a in self.structure if not a.idx-1 in conf and a.atomicnum>0] :
				idx =atom.idx-1
				# \now we find the hydrogen sidechain in the list and save its probability
				for sidechain in self.sidechainlist[idx]:
					smile       =sidechain[2]
					probability =sidechain[1]
					original_smile = smile
					# if not H, that SMILE cannot be an hydrogen sidechain
					if 'H' in smile:
						for char in smiles_char:
							#removing all reserved characters from the SMILE. Only atom names remain
							smile =smile.replace(char,'')
						# \now we check that only H are present in the string, otherwise it is not a hydrogen side chain 
						ishydrogen = True
						for char in smile:
							if not char=='H':
								ishydrogen=False
						#
						if ishydrogen:
							hs_probabilities.append(probability)
				#
			#//end for atom in [ a for a in self.structure if not a.idx-1 in conf and a.atomicnum>0]
			hs_cumulative_probability = reduce(mul,hs_probabilities)
			self.scaffold_probability_rel2configs.append(hs_cumulative_probability) 
			c+=1
		#//end for conf in self.configurations

		self.number_of_compounds=self.count_compounds()
		
	def load_data(self,scaffoldName,options={}):
		'''
		This method loads the scaffold data from a cIdx file <scaffold name>.cIdx and scaffold name>.sdf
		Args:
			string scaffoldName: the name of the scaffold
		Raises:
			IOError
			ScaffoldNullPopularity
			ScaffoldCorrupted
		Example:
			load_data("s000000001")
		'''
		
		verbose =options['verbose'] if options.has_key('verbose') else False
		debug   =options['debug'] if options.has_key('debug') else False
		#------------ data directly loaded from NPSD without further processing ------------------------------
		
		cIdxf=open(scaffoldName+".cIdx")
		
		# Processing the cIdx file 
		Idx=0
		configNumber=0
		sidechainBlock=0
		positionsNumber=0
		position=-1
		
		for line in cIdxf:	
			if Idx==0:
				a=line.strip().split()
				self.name       = a[0]
				self.popularity = int(a[1])
				if not self.popularity:
					break
				Idx+=1
				continue
				
			if Idx==1:
				a=line.strip().split()
				self.max_possible_weight = float(a[0])
				self.min_possible_weight = float(a[1])
				Idx+=1
				continue
				
			if Idx==2:
				configNumber=float(line.strip().split()[0])
				sidechainBlock=3+configNumber
				Idx+=1
				continue
				
			if Idx>=3 and Idx<sidechainBlock:
				a=line.strip().split()
				nonHydrogenPositions=int(a[0])
				self.configurations.append([int(x) for x in a[1:nonHydrogenPositions+1]])
				Idx+=1
				continue
				
			if Idx>=3 and Idx==sidechainBlock:
				self.scaffoldLength=int(line.split()[0])
				for i in range(0,self.scaffoldLength):
					self.sidechainlist.append([]) 
				Idx+=1	
				continue
				
			if Idx>=3 and Idx>sidechainBlock:
				a=line.strip().split('\t')
				if len(a)==2: #new position
					position         =int(a[0])
					sidechainsNumber =int(a[1])
				else:
					self.sidechainlist[position].append([float(a[4]),float(a[3]),a[0]])
				Idx+=1
				continue
				
			Idx+=1
		if [] in self.sidechainlist:
			raise ScaffoldCorrupted("Warning: some positions are missing: %d/%d. The file is probably corrupted." %(position+1,self.scaffoldLength))
			
		if not self.popularity:
			cIdxf.close()
			raise ScaffoldNullPopularity("Warning: the popularity of this scaffold is 0.")

		# Loading the scaffold molecule (.sdf file) thanks to pybel
		self.structure  = pybel.readfile("sdf",scaffoldName+".sdf").next()
		self.scaffold_weight_plus_hydrogens = self.structure.molwt

		#------------ data obtained by processing -----------------------------------------------------------
		c=0
		for conf in self.configurations:
			# we add hydrogens to be able to count them. Nevertheless it doesn't affect the molecular weight which always includes the hydrogens
			self.structure.OBMol.AddHydrogens()

			# \we compute the total number of hydrogens
			total_number_of_hydrogenes     =0
			for atom in [a for a in self.structure if not a.OBAtom.GetAtomicNum()==1]:
				if atom.idx-1 in conf:
					number_of_hydrogenes       =len([x for x in OBAtomAtomIter(atom.OBAtom) if x.GetAtomicNum()==1])
					total_number_of_hydrogenes +=number_of_hydrogenes
			#

			# \Dealing with the mass 
			# we remove the hydrogens because we don't need them anymore 
			self.structure.OBMol.DeleteHydrogens()
			# we compute the total mass of hydrogens at non hydrogen side chain positions    
			total_hydrogene_mass = total_number_of_hydrogenes * 1.00794
			# we compute the mass of the scaffold without the non hydrogen positions
			hs_weight            =self.scaffold_weight_plus_hydrogens - total_hydrogene_mass
			self.scaffold_weight_rel2configs.append(hs_weight)
			# 


			# \Dealing with the probability 
			# the condition a.atomicnum>0 is to ensure that we don't count the "atoms" R# that are just used to number the positions on the scaffold
			hs_probabilities=[1]
			for atom in [ a for a in self.structure if not a.idx-1 in conf and a.atomicnum>0] :
				idx =atom.idx-1
				# \now we find the hydrogen sidechain in the list and save its probability
				for sidechain in self.sidechainlist[idx]:
					smile       =sidechain[2]
					probability =sidechain[1]
					original_smile = smile
					# if not H, that SMILE cannot be an hydrogen sidechain
					if 'H' in smile:
						for char in smiles_char:
							#removing all reserved characters from the SMILE. Only atom names remain
							smile =smile.replace(char,'')
						# \now we check that only H are present in the string, otherwise it is not a hydrogen side chain 
						ishydrogen = True
						for char in smile:
							if not char=='H':
								ishydrogen=False
						#
						if ishydrogen:
							hs_probabilities.append(probability)
				#
			#//end for atom in [ a for a in self.structure if not a.idx-1 in conf and a.atomicnum>0]
			hs_cumulative_probability = reduce(mul,hs_probabilities)
			self.scaffold_probability_rel2configs.append(hs_cumulative_probability) 
			c+=1
		#//end for conf in self.configurations
		cIdxf.close()	

		self.number_of_compounds=self.count_compounds()
	
	def generate_csccp(self,R,mass_peak,ppm,configIdx):
		'''
		This method generates a CSCCP object that will be the input of all the CSCCP solver.
		Args:
			int R          : the csccp solver will output the top-R best solutions
			float mass_peak: the mass of the mass-spectrum peak that we want to identify
			int ppm        : resolution of the mass-spectrometer in ppm 
			int configIdx  : the index of the configuration that should be generated (@see scaffold.configurations)
		Returns:
			CSCCP csccp  : the result 
		raises:
			BadCSCCP
		'''

		csccp=CSCCP()

		csccp.structure=self.structure
		csccp.scaffold_name=self.name
		# we compute the number n of substituted positions 
		csccp.n =len(self.configurations[configIdx])
		csccp.R =R

		csccp.ppm=ppm
		csccp.mass_peak     =mass_peak
		csccp.mass_peak_max =(1.0 + float(ppm)/10000) * mass_peak
		csccp.mass_peak_min =(1.0 - float(ppm)/10000) * mass_peak 
		
		csccp.configIdx=configIdx
		csccp.scaffold_weight_rel2config      =self.scaffold_weight_rel2configs[configIdx] 
		csccp.scaffold_probability_rel2config =self.scaffold_probability_rel2configs[configIdx] 
		csccp.min_possible_weight     =self.min_possible_weight
		csccp.max_possible_weight     =self.max_possible_weight

		csccp.wmin          =csccp.mass_peak_min - csccp.scaffold_weight_rel2config
		csccp.wmax          =csccp.mass_peak_max - csccp.scaffold_weight_rel2config


		if csccp.wmin<0 or csccp.wmax<0:
			raise BadCSCCP("This CSCCP is not valid wmin<0 or wmax<0")

		def only_non_hydrogen_list(lst):
			return [sc for sc in lst if not is_hydrogen_sidechain(sc[2])]

		sidechains =[only_non_hydrogen_list(self.sidechainlist[position]) for position in self.configurations[configIdx]]
		csccp.W =[[x[0] for x in slist]for slist in sidechains]
		csccp.P =[[x[1] for x in slist]for slist in sidechains]
		csccp.K =[len(slist) for slist in sidechains]

		csccp.sidechain_smiles =[[x[2] for x in slist]for slist in sidechains]

		csccp.configuration    =self.configurations[configIdx]
		# we compute the number of possible compounds
		csccp.number_of_compounds=reduce(mul,csccp.K) if csccp.K else 1
		return csccp 

	def count_compounds(self):
		'''
		Computes the number of compounds that can be generated from the csccp
		Returns:
			int res: the result  
		'''
		def only_non_hydrogen_list(lst):
			return [sc for sc in lst if not is_hydrogen_sidechain(sc[2])]

		res=[]
		configIdx=0
		for configuration in self.configurations:
			sidechains =[only_non_hydrogen_list(self.sidechainlist[position]) for position in self.configurations[configIdx]]
			K =[len(slist) for slist in sidechains]
			if len(K):
				res.append((configIdx,reduce(mul,K)))
			else:
				res.append((configIdx,1))
			configIdx+=1
		return res
		



