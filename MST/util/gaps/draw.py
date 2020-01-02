#! /usr/bin/env python

from time import asctime, localtime, time
from os import environ
from sys import stderr, exit, argv
from string import Template
import re

'''
This file is used for generating an x3d file using a crystal.
Usage:
./draw.py input [output]
If output is not specified, the output will be saved as imput.x3d

A cross platform viewer for x3d files is view3dscene. The binaries can be found
here:
http://vrmlengine.sourceforge.net/view3dscene.php
'''
def main(infile,outfile=''):
	if outfile == '':
		outfile = infile + '.x3d'
	crystal = parse(infile)
	x3dOut(crystal,outfile)

'''
Output an x3d file.
The format is xml encoded x3d. The template strings used can be found below
'''
def x3dOut(crystal,outfile):
	a,b,c,atomset = crystal
	elements = []
	d=dict(a=str(a),b=str(b),c=str(c),
		ab=str(a+b),ac=str(a+c),bc=str(b+c),
		abc=str(a+b+c))
	box = box_template.substitute(d)
	# Get a list of the elements that will be in use.
	for atom in atomset:
		if atom.sym not in elements:
			elements.append(atom.sym)
	# Define the color for each one of those elements.
	defs = ''
	for element in elements:
		d = dict(sym=element,color=colorstr(colors[element]))
		defs += defs_template.substitute(d)
	atoms = ''
	# Insert a sphere for each atom.
	for atom in atomset:
		d = dict(x=atom.x,y=atom.y,z=atom.z,sym=atom.sym,rad=atom.rad)
		atoms += atom_template.substitute(d)
	# Get the name of the creator.
	if 'NAME' in environ.keys():
		name = environ['NAME']
	elif 'HOSTNAME' in environ.keys():
		name = environ['USER'] + '@' + environ['HOSTNAME']
	elif 'LOGNAME' in environ.keys():
		name = environ['LOGNAME']
	else:
		name = environ['USER']
	timestamp=asctime(localtime(time()))	
	d = dict(timestamp=timestamp,box=box,defs=defs,atoms=atoms,name=name)
	out = main_template.substitute(d)
	fp = open(outfile,'w')
	fp.write(out)
	fp.close()

'''
Reads in the atomic structure and identifies the atoms and vectors as well as 
the redefined radii
'''
def parse(filename):
	fp = open(filename)
	direct = False
	atoms = []
	# Regex for matching vectors.
	vec = re.compile("""
		^[ \t]* 							# whitespace
		(-?[0-9]*[.]?[0-9]+[E]?[+-]?[0-9]*)	# floating point number 
		[ \t]*[,; \t][ \t]*					# delimiter
		(-?[0-9]*[.]?[0-9]+[E]?[+-]?[0-9]*)	# floating point number 
		[ \t]*[,; \t][ \t]*					# delimiter
		(-?[0-9]*[.]?[0-9]+[E]?[+-]?[0-9]*)	# floating point number 
		""",re.VERBOSE)
	# Regex for amtching atoms.
	atom = re.compile("""
		^[ \t]*(?:[0-9]*[ \t]*[;, \t])?[ \t]* # whitespace and a possible int
		(-?[0-9A-Z][0-9a-z]*[.]?[0-9]*[ED]?[+-]?[0-9]*) 
			# atomic radius, number or element name
		[ \t]*[,; \t][ \t]* 				# delimiter
		(-?[0-9]*[.]?[0-9]+[E]?[+-]?[0-9]*)	# floating point number 
		[ \t]*[,; \t][ \t]*					# delimiter
		(-?[0-9]*[.]?[0-9]+[E]?[+-]?[0-9]*)	# floating point number 
		[ \t]*[,; \t][ \t]*					# delimiter
		(-?[0-9]*[.]?[0-9]+[E]?[+-]?[0-9]*)	# floating point number 
		""",re.VERBOSE);
	# Regex for matching atomic radious definitions.
	redef = re.compile('''
		^[ \t]*\#[ \t]*						# Spacing and comment symbol
		([0-9A-Z][0-9a-z]?)					# Atomic symbol or number
		[ \t]*=[ \t]*						# Spacing and '='
		(-?[0-9]*[.]?[0-9]+[E]?[+-]?[0-9]*)	# floating point number
		''',re.VERBOSE)
	
	# Match the first bounding vector
	for line in fp:
		m = vec.match(line)
		if m:
			a = Vector(m.group(1),m.group(2),m.group(3))
			break
	# Match the second bounding vector
	for line in fp:
		m = vec.match(line)
		if m:
			b = Vector(m.group(1),m.group(2),m.group(3))
			break
	# Match the third bounding vector
	for line in fp:
		m = vec.match(line)
		if m:
			c = Vector(m.group(1),m.group(2),m.group(3))
			break
	# Match atoms, atom radii and whether the coordinates are direct
	for line in fp:
		# Check if the coordinates are direct
		if re.match("[ \t]directCoordinates",line):
			direct = True
		else:
			# Match atoms
			m = atom.match(line)
			if m:
				if direct:
					Ca,Cb,Cc = map(float,m.group(2,3,4))
					x = Ca*a.x + Cb*b.x + Cc*c.x
					y = Ca*a.y + Cb*b.y + Cc*c.y
					z = Ca*a.z + Cb*b.z + Cc*c.z
					atoms.append(Atom(x,y,z,m.group(1)))
				else:
					t,x,y,z = m.groups()
					atoms.append(Atom(x,y,z,t))
			else:
				# It's not an atom so try matching it to a radius definition
				m = redef.match(line)
				if m:
					t = m.group(1)
					if t.isdigit():
						t = symbols[int(t)]
					Atom.redefs[m.group(1)] = float(m.group(2))
	fp.close()	
	return (a,b,c,atoms)


class Vector:
	def __init__(self,x=0.0,y=0.0,z=0.0):
		self.x = float(x)
		self.y = float(y)
		self.z = float(z)
	# Override the str operation
	def __str__(self):
		return str(self.x)+' '+str(self.y)+' '+str(self.z)
	# Override +
	def __add__(self,other):
		return Vector(self.x+other.x,self.y+other.y,self.z+other.z)

class Atom:
	redefs = {}
	def __init__(self,x=0.0,y=0.0,z=0.0,t=''):
		self.x = float(x)
		self.y = float(y)
		self.z = float(z)
		if '.' in t:
			self.sym = 'Va'
			self.rad = float(t)
			return
		elif t.isalpha():
			self.sym = t
		elif t.isdigit():
			self.sym = symbols[abs(int(t))]
		else:
			print >> stderr, 'Invalid type ', t
			exit(2)

		if self.sym in Atom.redefs.keys():
			self.rad = Atom.redefs[self.sym]
		else:
			self.rad = radii[symbols.index(self.sym)]

# The template which all others are inserted in
main_template = Template('''<?xml version='1.0' encoding='UTF-8'?>

<X3D version='3.0' profile='Interchange'
xmlns:xsd='http://www.w3.org/2001/XMLSchema-instance'
xsd:noNamespaceSchemaLocation='http://www.web3d.org/specifications/x3d-3.0.xsd'
>
	<head>
		<meta content='$name' name='creator' />
		<meta content='$timestamp' name='created'/>
	</head>
	<Scene>
		<Background skyColor='0 0 0'/>
$box
		<!-- Defining the colors for the elements -->$defs

		<!-- The atoms themselves -->$atoms
	</Scene>
</X3D>''')

# Creates the bounding box for the crystal
box_template = Template('''
		<!-- The box made up of the three bounding vectors -->
		<Shape>
			<Appearance>
				<Material emissiveColor='1 1 1'/>
			</Appearance>
			<IndexedLineSet coordIndex='0 1 5 3 0 2 4 1 4 7 5 7 6 3 6 2'>
				<Coordinate point='0 0 0 
					$a $b $c 
					$ab $ac $bc 
					$abc'/>
			</IndexedLineSet>
		</Shape>
''')

# Defines an element color
defs_template = Template('''
		<Appearance DEF='$sym'>
			<Material ambientIntensity='1' shininess='0.234375'
				diffuseColor = '$color' 
				specularColor='0.85 0.85 0.85'/>
		</Appearance>''')

# Defines an atom by radius
atom_template = Template('''
		<Transform translation='$x $y $z'>
			<Shape>
				<Appearance USE='$sym' />
				<Sphere radius='$rad' />
			</Shape>
		</Transform>''')

'''
Transform an rgb color into a string of floats
'''
def colorstr(rgb):
	r = rgb/0x10000
	g = (rgb/0x100) % 0x100
	b = rgb % 0x100

	s = str(float(r)/float(0xff))
	s += " "
	s += str(float(g)/float(0xff))
	s += " "
	s += str(float(b)/float(0xff))
	return s

'''
Colors for each atom type. The colors are the same as those found at:
http://www.cgl.ucsf.edu/chimera/docs/UsersGuide/colortables.html#byelement
and
http://jmol.sourceforge.net/jscolors/
'''
colors =	{ "Va":0xffffff,"H" :0xffffff,"He":0xd9ffff,"Li":0xcc80ff,
"Be":0xc2ff00,"B" :0xffb5b5,"C" :0x909090,"N" :0x3050f8,"O" :0xff0d0d,
"F" :0x90e050,"Ne":0xb3e3f5,"Na":0xab5cf2,"Mg":0x8aff00,"Al":0xbfa6a6,
"Si":0xf0c8a0,"P" :0xff8000,"S" :0xffff30,"Cl":0x1ff01f,"Ar":0x80d1e3,
"K" :0x8f40d4,"Ca":0x3dff00,"Sc":0xe6e6e6,"Ti":0xbfc2c7,"V" :0xa6a6ab,
"Cr":0x8a99c7,"Mn":0x9c7ac7,"Fe":0xe06633,"Co":0xf090a0,"Ni":0x50d050,
"Cu":0xc88033,"Zn":0x7d80b0,"Ga":0xc28f8f,"Ge":0x668f8f,"As":0xbd80e3,
"Se":0xffa100,"Br":0xa62929,"Kr":0x5cb8d1,"Rb":0x702eb0,"Sr":0x00ff00,
"Y" :0x94ffff,"Zr":0x94e0e0,"Nb":0x73c2c9,"Mo":0x54b5b5,"Tc":0x3b9e9e,
"Ru":0x248f8f,"Rh":0x0a7d8c,"Pd":0x006985,"Ag":0xc0c0c0,"Cd":0xffd98f,
"In":0xa67573,"Sn":0x668080,"Sb":0x9e63b5,"Te":0xd47a00,"I" :0x940094,
"Xe":0x429eb0,"Cs":0x57178f,"Ba":0x00c900,"La":0x70d4ff,"Ce":0xffffc7,
"Pr":0xd9ffc7,"Nd":0xc7ffc7,"Pm":0xa3ffc7,"Sm":0x8fffc7,"Eu":0x61ffc7,
"Gd":0x45ffc7,"Tb":0x30ffc7,"Dy":0x1fffc7,"Ho":0x00ff9c,"Er":0x00e675,
"Tm":0x00d452,"Yb":0x00bf38,"Lu":0x00ab24,"Hf":0x4dc2ff,"Ta":0x4da6ff,
"W" :0x2194d6,"Re":0x267dab,"Os":0x266696,"Ir":0x175487,"Pt":0xd0d0e0,
"Au":0xffd123,"Hg":0xb8b8d0,"Tl":0xa6544d,"Pb":0x575961,"Bi":0x9e4fb5,
"Po":0xab5c00,"At":0x754f45,"Rn":0x428296,"Fr":0x420066,"Ra":0x007d00,
"Ac":0x70abfa,"Th":0x00baff,"Pa":0x00a1ff,"U" :0x008fff,"Np":0x0080ff,
"Pu":0x006bff,"Am":0x545cf2}

# Atomic radii. For citation see radii.c in the gaps source code.
radii	= [ 1.00, 0.25, 0.00, 1.45, 1.05, 0.85, 0.70, 0.65, 0.60, 0.50, 0.00, 
1.80, 1.50, 1.25, 1.10, 1.00, 1.00, 1.00, 0.00, 2.20, 1.80, 1.60, 1.40, 1.35, 
1.40, 1.40, 1.40, 1.35, 1.35, 1.35, 1.35, 1.30, 1.25, 1.15, 1.15, 1.15, 0.00, 
2.35, 2.00, 1.80, 1.55, 1.45, 1.45, 1.35, 1.30, 1.35, 1.40, 1.60, 1.55, 1.55, 
1.45, 1.45, 1.40, 1.40, 0.00, 2.60, 2.15, 1.95, 1.85, 1.85, 1.85, 1.85, 1.85, 
1.85, 1.80, 1.75, 1.75, 1.75, 1.75, 1.75, 1.75, 1.75, 1.55, 1.45, 1.35, 1.35, 
1.30, 1.35, 1.35, 1.35, 1.50, 1.90, 1.80, 1.60, 1.90, 0.00, 0.00, 0.00, 2.15, 
1.95, 1.80, 1.80, 1.75, 1.75, 1.75, 1.75 ]

# A list of atomic symbols.
symbols = [ "Va", "H" , "He", "Li", "Be", "B" , "C" , "N" , "O" , "F" , "Ne", 
"Na", "Mg", "Al", "Si", "P" , "S" , "Cl", "Ar", "K" , "Ca", "Sc", "Ti", "V" , 
"Cr", "Mn", "Fe", "Co", "Ni", "Cu", "Zn", "Ga", "Ge", "As", "Se", "Br", "Kr", 
"Rb", "Sr", "Y" , "Zr", "Nb", "Mo", "Tc", "Ru", "Rh", "Pd", "Ag", "Cd", "In", 
"Sn", "Sb", "Te", "I" , "Xe", "Cs", "Ba", "La", "Ce", "Pr", "Nd", "Pm", "Sm", 
"Eu", "Gd", "Tb", "Dy", "Ho", "Er", "Tm", "Yb", "Lu", "Hf", "Ta", "W" , "Re", 
"Os", "Ir", "Pt", "Au", "Hg", "Tl", "Pb", "Bi", "Po", "At", "Rn", "Fr", "Ra", 
"Ac", "Th", "Pa", "U" , "Np", "Pu", "Am"]

if __name__ == "__main__":
	c = len(argv)
	if c == 1:
		print >> stderr, "Not enough arguments"
	elif c == 2:
		main(argv[1])
	elif c >= 3:
		main(argv[1],argv[2])
