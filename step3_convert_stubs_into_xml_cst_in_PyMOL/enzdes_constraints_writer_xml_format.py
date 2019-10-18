#Script to be run in PyMOL session where pdbs of networked residues are loaded.
#Each object represents different stereo isomer and conformation of LYS-LG2 conjugate
#with different configurations of networks represented in different states of the object.

#Script extracts 3D coordinates of anchor atoms of ligand (hardcoded in lig_ats dictionary)
#and protein side chains (hardcoded in get_upstream_anchors function).
#computes 6 dof between upstream and downstream anchors and writes them in the form of xml constraint file.

#(2017-07-12)
#To change which order constraints are written into individual files ( nuc-bridge-shuttle-support or support-shuttle-bridge-nuc)
#need to change (uncomment) two matching lines
#.....
##ress=sorted(ress, key=lambda x: int(x[1]))                              # for ["sup","sht","brd","nuc"] order of constraints
#ress=sorted(ress, key=lambda x: int(x[1]), reverse=True)                 # reversing for ["nuc","brd","sht","sup"] order of constraints
#.....
##cst_file_template=[("TYR","SER","THR"),("TYR",),("ASN","GLN"),("LYS",)] # for ["sup","sht","brd","nuc"] order of constraints
#cst_file_template=[("LYS",),("ASN","GLN"),("TYR",),("TYR","SER","THR")]  # reversing for ["nuc","brd","sht","sup"] order of constraints

#(2017-07-14)
#looks like can cut down on number of constraints by selecting only unique ones. Earlier was operating
#under assumption that number of constraint blocks for each interaction needs to be identical therefore was writing hundreds of identical
#blocks for LYS-LG2 to equalize their number with number of blocks for other constraints. Which is deleterious for matcher performance.


from pymol import cmd
import sys
from itertools import product

def get_6dof_cst(D3,D2,D1,U1,U2,U3):
	cur_cst={
			"distanceAB":[0.0, 0.2,80.0,  0,0],
			"angle_A":   [0.0, 5.0,10.0,360,0],
			"angle_B":   [0.0, 5.0,10.0,360,0],
			"torsion_A": [0.0,10.0,10.0,360,0],
			"torsion_AB":[0.0,10.0,10.0,360,0],
			"torsion_B": [0.0,10.0,10.0,360,0]
			}
	if U1.split("/")[-1] == "NZ" and U2.split("/")[-1] == "CE" and U3.split("/")[-1] == "2HZ" : #means processing constraint for LYS
		cur_cst={
				"distanceAB":[0.0, 0.1,80.0,  1,0], #the only difference is periodicity will be set at 1 indicating covalent interaction in the cst file
				"angle_A":   [0.0, 5.0,10.0,360,0],
				"angle_B":   [0.0, 5.0,10.0,360,0],
				"torsion_A": [0.0,10.0,10.0,360,0],
				"torsion_AB":[0.0,10.0,10.0,360,0],
				"torsion_B": [0.0,10.0,10.0,360,0]
				}
	cur_cst["distanceAB"][0]=round(cmd.get_distance(D1,U1),1)
	cur_cst["angle_A"][0]=round(cmd.get_angle(D2,D1,U1),1)
	cur_cst["angle_B"][0]=round(cmd.get_angle(D1,U1,U2),1)
	cur_cst["torsion_A"][0]=round(cmd.get_dihedral(D3,D2,D1,U1),1)
	cur_cst["torsion_AB"][0]=round(cmd.get_dihedral(D2,D1,U1,U2),1)
	cur_cst["torsion_B"][0]=round(cmd.get_dihedral(D1,U1,U2,U3),1)
	return cur_cst

def res3to1_res1to3(resn):
	dic3to1={
			"ALA":"A","CYS":"C","ASP":"D","GLU":"E","PHE":"F",
			"GLY":"G","HIS":"H","ILE":"I","LYS":"K","LEU":"L",
			"MET":"M","ASN":"N","PRO":"P","GLN":"Q","ARG":"R",
			"SER":"S","THR":"T","VAL":"V","TRP":"W","TYR":"Y"
			}
	dic1to3={
			"A":"ALA","C":"CYS","D":"ASP","E":"GLU","F":"PHE",
			"G":"GLY","H":"HIS","I":"ILE","K":"LYS","L":"LEU",
			"M":"MET","N":"ASN","P":"PRO","Q":"GLN","R":"ARG",
			"S":"SER","T":"THR","V":"VAL","W":"TRP","Y":"TYR"
			}
	if len(resn) == 3 and resn in dic3to1.keys():
		return dic3to1[resn]
	elif len(resn) == 1 and resn in dic1to3.keys():
		return dic1to3[resn]
	else:
		print "Something went wrong. Cannot convert name: {}".format(resn)

def get_upstream_anchors(resn,names_only=True,obj_name=None,idx=None): #can return names of atoms only, or selection string for these atoms
	if names_only:
		if resn == "LYS":
			U1='NZ'
			U2='CE'
			U3='2HZ'
		elif resn == "TYR":
			U1='OH'
			U2='CZ'
			U3='HH'
		elif resn == "SER":
			U1='OG'
			U2='CB'
			U3='HG'
		elif resn == "THR":
			U1='OG1'
			U2='CB'
			U3='HG1'
		elif resn == "ASN":
			U1='OD1'
			U2='CG'
			U3='1HD2'
		elif resn == "GLN":
			U1='OE1'
			U2='CD'
			U3='1HE2'
		else:
			sys.exit("Restype is not expected to be in RA theozyme")
	else:
		if resn == "LYS":
			U1='{}///{}`{}/NZ'.format(obj_name,resn,idx)
			U2='{}///{}`{}/CE'.format(obj_name,resn,idx)
			U3='{}///{}`{}/2HZ'.format(obj_name,resn,idx)
		elif resn == "TYR":
			U1='{}///{}`{}/OH'.format(obj_name,resn,idx)
			U2='{}///{}`{}/CZ'.format(obj_name,resn,idx)
			U3='{}///{}`{}/HH'.format(obj_name,resn,idx)
		elif resn == "SER":
			U1='{}///{}`{}/OG'.format(obj_name,resn,idx)
			U2='{}///{}`{}/CB'.format(obj_name,resn,idx)
			U3='{}///{}`{}/HG'.format(obj_name,resn,idx)
		elif resn == "THR":
			U1='{}///{}`{}/OG1'.format(obj_name,resn,idx)
			U2='{}///{}`{}/CB'.format(obj_name,resn,idx)
			U3='{}///{}`{}/HG1'.format(obj_name,resn,idx)
		elif resn == "ASN":
			U1='{}///{}`{}/OD1'.format(obj_name,resn,idx)
			U2='{}///{}`{}/CG'.format(obj_name,resn,idx)
			U3='{}///{}`{}/1HD2'.format(obj_name,resn,idx)
		elif resn == "GLN":
			U1='{}///{}`{}/OE1'.format(obj_name,resn,idx)
			U2='{}///{}`{}/CD'.format(obj_name,resn,idx)
			U3='{}///{}`{}/1HE2'.format(obj_name,resn,idx)
		else:
			sys.exit("Restype is not expected to be in RA theozyme")
	return (U1,U2,U3)

####################################( Main )###########################################################################	
obj_list=cmd.get_names(enabled_only=1) #selecting visible objects only
# print obj_list

lign="MTK"
#atom names in HBnetGenerated pdb MTK stub file are differen from atom names in params/pdb files of LG2 ligand used for matching or backbone building 
#tuple[0] is in MTK name, tuple[1] is in LG2
lig_ats={
		"D1":("C11","C12"),
		"D2":("C10","C11"),
		"D3":("C12","C13"),
		}
#need to decide how to hold all data for constraints before writing them into cst file.
#Using list of lists to keep order of all elements synchronized and order of cst blocks
#[ nuc=[cst1,...cst#networks],brd=[cst1,...cst#networks],sht=[cst1,...],sup=[cst1,...] ]

theozyme_initialized=False
count_skipped_networks=0

for obj_name in obj_list:
	for state in range(1,cmd.count_states(obj_name)+1):
		cmd.set("state",state)
		stored.lig=[]
		cmd.iterate("{} and chain X".format(obj_name),"stored.lig.append((resn,resi))",1) #for some reason ligand (== chain X) has different resi => have to read it each iteration
		lig=list(set(stored.lig)) #iterate runs over all atoms in the selection, applying "set" collapses into one pair (resn,resi)
		lign=lig[0][0]
		ligi=lig[0][1]
		D1='{}///{}`{}/{}'.format(obj_name,lign,ligi,lig_ats["D1"][0])
		D2='{}///{}`{}/{}'.format(obj_name,lign,ligi,lig_ats["D2"][0])
		D3='{}///{}`{}/{}'.format(obj_name,lign,ligi,lig_ats["D3"][0])

		stored.ress=[]
		cmd.iterate("{} and chain A".format(obj_name),"stored.ress.append((resn,resi))",1)
		ress=list(set(stored.ress)) #iterate runs over all atoms in the selection, applying "set" collapses into one pair (resn,resi)
		#ress=sorted(ress, key=lambda x: int(x[1])) # for ["sup","sht","brd","nuc"] order of constraints
		ress=sorted(ress, key=lambda x: int(x[1]), reverse=True) # reversing for ["nuc","brd","sht","sup"] order of constraints
		
		#resns=map(int,resns)
		#looks like order of residues in the theozyme (as appears in the base_stubs.pdb) is the following, but potentially can change as well as restypes
		# 1 -- support     (Y) can be S or T
		# 2 -- shuttle     (Y)
		# 3 -- bridge      (N) can be Q
		# 4 -- nucleophile (K)
		# initializing theozyme=["nuc"[],"brd"[],"sht"[],"sup"[]] if was not done previously
		if not theozyme_initialized:
			theozyme=[[] for i in range(len(ress))]
			theozyme_initialized=True
			
		for i in range(len(ress)):
			idx=ress[i][1]
			resn=ress[i][0]
			U1,U2,U3=get_upstream_anchors(resn,0,obj_name,idx)
			
			atoms={
					"Dres":D1.split("/")[-2].split("`")[0],
					"D1":  D1.split("/")[-1],
					"D2":  D2.split("/")[-1],
					"D3":  D3.split("/")[-1],
					"Ures":U1.split("/")[-2].split("`")[0],
					"U1":  U1.split("/")[-1],
					"U2":  U2.split("/")[-1],
					"U3":  U3.split("/")[-1],
					}
			
			try:
				#print(D3,D2,D1,U1,U2,U3)
				cur_cst_block=get_6dof_cst(D3,D2,D1,U1,U2,U3) #at this point of the runtime D1-3 and U1-3 are selection strings
			except CmdException:
				print "Something went wrong. Skipping current network #",state," in object ",obj_name #Some networks in current set of stub files have missing residues (ask Gerard ???). Skipping these networks
				count_skipped_networks+=1
				break
			
			#renaming D atoms to match names in LG2 (as opposed to MTK) representation of the ligand before writing constraint file
			atoms["Dres"]="LG2"
			atoms["D1"]=lig_ats["D1"][1]
			atoms["D2"]=lig_ats["D2"][1]
			atoms["D3"]=lig_ats["D3"][1]
			
			if 	(atoms,cur_cst_block) not in theozyme[i]:
# 				print "Appending to theozyme{i}: ".format(i=i), (atoms,cur_cst_block)
				theozyme[i].append((atoms,cur_cst_block)) #accumulates constraints definitions for each constraint and for for each state if they have not been seen yet

print "Number of skipped networks: ",count_skipped_networks

#cst_file_template=[("TYR","SER","THR"),("TYR",),("ASN","GLN"),("LYS",)] #each tuple in the list is a constraint, and each element in the tuple is what restype can be in that constraint
cst_file_template=[("LYS",),("TYR",),("TYR","SER","THR"),("ASN","GLN")]
#cst_file_template=[("TYR",),("LYS",),("ASN","GLN"),("TYR","SER","THR")]
#cst_file_template=[("TYR","SER","THR"),("TYR",),("ASN","GLN")]
#cst_file_template=[("TYR",),("TYR",),("ASN",),("LYS",)]

out_name_prefix="RA2_{ntwrkid}".format(ntwrkid=obj_list[0].split("_")[4]) #assumes obj names are of the form grp_ntwrk_base_stubs_0001
for template in product(*cst_file_template): #generate all individual combinations of residues in template
	out_name="{prefix}_{res_combination}.cst.txt".format(prefix=out_name_prefix,res_combination="".join(res3to1_res1to3(res3) for res3 in template))
	with open(out_name,"w") as out:
		for idx,cst in enumerate(theozyme[:]):
			Ures=template[idx]
			U1,U2,U3=get_upstream_anchors(Ures,1)
			out.write('<MatcherConstraint>\n')
			out.write('\t<DownstreamResidue   atom1="{U1}" atom2="{U2}" atom3="{U3}" name="{resn_XXX}" />\n'.format(U1=U1,U2=U2,U3=U3,resn_XXX=Ures)) #seems like BW switched the logic of what is called Up and Down, so have to adjust
			out.write('\t<UpstreamResidue atom1="{D1}" atom2="{D2}" atom3="{D3}" name="{resn_XXX}" />\n'.format(D1=cst[0][0]["D1"],D2=cst[0][0]["D2"],D3=cst[0][0]["D3"],resn_XXX=cst[0][0]["Dres"]))
			for rcrd in cst:
				out.write('\t<Combination>\n')
				out.write('\t\t<DistanceAB x0="{}" xtol="{}" k="{}" periodicity="{}" noSamples="{}" />\n'.format(*rcrd[1]["distanceAB"]))
				out.write('\t\t<AngleA    x0="{}" xtol="{}" k="{}" periodicity="{}" noSamples="{}" />\n'.format(*rcrd[1]["angle_A"]))
				out.write('\t\t<AngleB    x0="{}" xtol="{}" k="{}" periodicity="{}" noSamples="{}" />\n'.format(*rcrd[1]["angle_B"]))
				out.write('\t\t<TorsionA  x0="{}" xtol="{}" k="{}" periodicity="{}" noSamples="{}" />\n'.format(*rcrd[1]["torsion_A"]))
				out.write('\t\t<TorsionAB x0="{}" xtol="{}" k="{}" periodicity="{}" noSamples="{}" />\n'.format(*rcrd[1]["torsion_AB"]))
				out.write('\t\t<TorsionB  x0="{}" xtol="{}" k="{}" periodicity="{}" noSamples="{}" />\n'.format(*rcrd[1]["torsion_B"]))
				out.write('\t</Combination>\n')	
			out.write('</MatcherConstraint>\n')
