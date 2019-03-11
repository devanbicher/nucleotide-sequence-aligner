import sys
import getopt
import string
import os
import re
import glob

oldSeqs = []; ##this is a list consisting the of original sequences and their names
newSeqs = []; ##this is the new list containing the new alignments(backwards)

x = [];
y = [];
m = [];

########################################################################################################
########################################################################################################
###usage statement
########################################################################################################
########################################################################################################
def usage():
	print("This is Devan Bicher's Sequence Alignment Program for CSE 308, Project2");
	print("");
	print("The formatting for this is as follows:");
	print("\tpython Aligner.py [inputSequences.fasta] [outputfile]");
	print("where: inputSequences.fasta is the 2 sequences to be aligned in fasta format");
	print("       and outputfile is the file that the alignment will be written to");
	print("NOTE: [outputfile] is optional, either way the alignment will print to the screen");
	print("NOTE: Your Terminal window width may need to be adjusted to properly display the alignment*");
	print("          *The display shows the names + 100 bases at a time.");
	print("ONE MORE NOTE: A simple alignment runs first and it simply tests if the sequences are the same");
	print("                   or if one sequence is completely contained within the other.");
	print("");
	print("");

########################################################################################################
########################################################################################################
###File Reading in stuff
########################################################################################################
########################################################################################################
	
################################################
###gobble in a file
################################################
def gobble(fileName):
	#print("Gobbling: "+ fileName);
	
	inputFile = open(fileName, 'rt');
	inputLine = inputFile.readline();
	
	outputList = "".split();  #### list declaration
	
	while inputLine:
		outputList.append(inputLine);
		inputLine = inputFile.readline();

	inputFile.close();
	
	return outputList;

################################################
###parse the fasta file
################################################
def parseFasta(fileName):
	lineList = gobble(fileName);
	seqs = [];
	names = [];

	##first get the mutant sequences.
	tempLine = "";
	for line in lineList:
		if(line[0] == ">"):
			if(tempLine != "" ):
				seqs.append(tempLine);
			tempLine = "";
			names.append(line[1:-1]);
			continue;
		tempLine = tempLine + line[0:-1];
	seqs.append(tempLine);
	
	oldSeqs.append(seqs);
	oldSeqs.append(names);

########################################################################################################
########################################################################################################
###Aligning Methods
########################################################################################################
########################################################################################################	

################################################
###simple alignment test
################################################
def simpleAlign():
	## This simplest case alignment checks if the two sequences are equal or if one is completely contained in the other
	if oldSeqs[0][0] == oldSeqs[0][1]:
		newSeqs.append(oldSeqs[0][0]);
		newSeqs.append(oldSeqs[0][1]);
		return True;
		
	elif oldSeqs[0][0] in oldSeqs[0][1]:
		newSeqs.append(substring(oldSeqs[0][0],oldSeqs[0][1]));
		newSeqs.append(oldSeqs[0][1]);
		return True;
	
	elif oldSeqs[0][1] in oldSeqs[0][0]:
		newSeqs.append(oldSeqs[0][0]);
		newSeqs.append(substring(oldSeqs[0][1],oldSeqs[0][0]));
		return True;
	
	return False;

################################################
###simple Substring alignment
################################################	
def substring(sub,string):
	newSub = '';
	for i in range(string.index(sub)):
		newSub += '-';
	
	newSub += sub;
	
	for i in range(len(string)-len(newSub)):
		newSub += '-';

	return newSub;
	
################################################
###simple Substring alignment
################################################	
def match(a,b):
	if a == b:
		return 1; #match
	else: 
		return -4; #mismatch

################################################
###simple Substring alignment
################################################
def setupTables():
	
	## z = x,y,m --> z[i][j] 
	for i in range(len(oldSeqs[0][1])+1):
		x.append([]);
		y.append([]);
		m.append([]);
		for j in range(len(oldSeqs[0][0])+1):
			if i == 0 and j == 0:  ## top right corner *[0][0] = 0 for all graphs
				x[i].append(0);
				y[i].append(0);
				m[i].append(0);
			elif i == 0: ##top row of the graphs 
				y[i].append(-999999);
				m[i].append(-999999);
				if j == 1:
					x[i].append(-10.5);
				else:
					x[i].append(x[i][j-1] - .5);
			elif j == 0: ##left column of the graphs
				x[i].append(-999999);
				m[i].append(-999999);
				if i == 1:
					y[i].append(-10.5);
				else:
					y[i].append(y[i-1][j] - .5);
			else:
				##do something with x
				x[i].append(max(-10 -.5 + m[i][j-1],-10 + x[i][j-1],-10 -.5 +y[i][j-1]));
				##do something with y
				y[i].append(max(-10 -.5 + m[i-1][j],-10 -.5 + x[i-1][j],-.5 +y[i-1][j]));
				##do something with m
				m[i].append(match(oldSeqs[0][1][i-1],oldSeqs[0][0][j-1]) + max(m[i-1][j-1],x[i][j],y[i][j]));

################################################
###add bases to the newSeqs alignment sequences
################################################
def extendNewSeqs(z,i,j):
	if z == 'm':
		newSeqs[0] += oldSeqs[0][0][j-1];
		newSeqs[1] += oldSeqs[0][1][i-1];
	elif z == 'x':
		newSeqs[0] += oldSeqs[0][0][j-1];
		newSeqs[1] += '-';
	elif z == 'y':
		newSeqs[0] += '-';
		newSeqs[1] += oldSeqs[0][1][i-1];

################################################
###BackTracking
################################################
def backtracking():
	i = len(oldSeqs[0][1]);
	j = len(oldSeqs[0][0]);
	
	backList = [x[i][j], y[i][j], m[i][j]];
	
	while i != 0 and j != 0:
		
		node = max(backList);
		if node == backList[0]: 	##node on the X graph
			extendNewSeqs('x',i,j);
			j += -1;	##we are on the X graph so decrease the j index by one
			##Create a list of the possible backtrack possibilities
			backList = [];	
			backList.append(x[i][j] -.5);
			backList.append(y[i][j] -10 -.5);
			backList.append(m[i][j] -10 -.5);
			
		elif node == backList[1]: 	##node on the Y graph
			extendNewSeqs('y',i,j);
			i += -1;	##we are on the Y graph so decrease the i index by one
			##Create a list of the possible backtrack possibilities
			backList = [];
			backList.append(x[i][j] -10 -.5);
			backList.append(y[i][j] -.5);
			backList.append(m[i][j] -10 -.5);
			
		else:  						##node on the M graph
			extendNewSeqs('m',i,j);
			mscore = match(oldSeqs[0][0][j-2], oldSeqs[0][1][i-2]);
			backList = [];
			backList.append(x[i][j]);
			backList.append(y[i][j]);
			backList.append(mscore + m[i-1][j-1]);
			##this decides how to increment our i and j indexes
			biggest = max(backList);
			if biggest == backList[2]:
				i += -1;
				j += -1;

################################################
###fix the length of the names to be the same
################################################
def fixNames(name1,name2):
	## this simply makes the names of the sequences the same length so it looks pretty
	if len(name1) > len(name2):
		newName1 = name1;
		newName2 = '';
		for i in range(len(name1) - len(name2)):
			newName2 += ' ';
		newName2 += name2;	
		
	else:
		newName2 = name2;
		newName1 = '';
		for i in range(len(name2) - len(name1)):
			newName1 += ' ';
		newName1 += name1;
	
	names = [];
	names.append(newName1);
	names.append(newName2);
	
	return names;
				
########################################################################################################
########################################################################################################
###main method
########################################################################################################
########################################################################################################
def main():
	usage();
	if len(sys.argv) == 1:
		print("But you did not enter any sequence file to be read in, please run again with the usage described above");
		raise SystemExit, 5
	elif not os.path.exists(sys.argv[1]):
		print("The sequence file that you eneterd does not exist or is named incorrectly, please run again with a proper file");
		raise SystemExit, 5
	
	parseFasta(sys.argv[1]);	

	##create the output filename if it exists	
	fileName = '';
	if(len(sys.argv) == 3):
		fileName = sys.argv[2];
		finalFile = open(fileName, 'w');
	
	##simple alignment test
	print("Running Simple Alignment test...");
	if simpleAlign() == True:
		print("Simple Alignment Test successful!"); 
		## IF the simple alignment worked than we don't need to do the affine gap scoring
	else:
		newSeqs.append("");
		newSeqs.append("");
		print("The simple Alignment test was unsuccessful, performing affine gap scoring alignment...");
		print("Setting up X, Y, and M graphs...");
		##Setup the X,Y,M tables
		setupTables();
		print("Graphs Created, performing alignment with backtracking...");
		##now perform the backtracking
		backtracking();
			
	###we need to reverse the sequences since they were assembled in reverse order
	finalSeq1 = newSeqs[0][::-1];
	finalSeq2 = newSeqs[1][::-1];
	
	##fix the lengths of the names to be the same
	if len(oldSeqs[1][0]) != len(oldSeqs[1][1]):
		names = fixNames(oldSeqs[1][0],oldSeqs[1][1]);
	else:
		names = [];
		name1 = oldSeqs[1][0];
		name2 = oldSeqs[1][1];
		names.append(name1);
		names.append(name2);
	
	##final Printing stuff
	print("\nBacktracking done.  Printing Alignments...\n\n");
	pos = 0;
	while pos < len(finalSeq1):
		fpos = pos + 100;
		if fpos > len(finalSeq1):
			fpos = len(finalSeq1);
		print(names[0] + ":\t" +finalSeq1[pos:fpos]); ##print the 1t sequence
		print(names[1] + ":\t" +finalSeq2[pos:fpos]); ##print the 2nd sequence over the first
		if len(fileName) > 0:
			finalFile.write(names[0] + ":\t" +finalSeq1[pos:fpos] + "\n"); ##print the 1t sequence to the file
			finalFile.write(names[1] + ":\t" +finalSeq2[pos:fpos] + "\n"); ##print the 2nd sequence over the first to the file
		pos = fpos;
		
	if len(fileName) > 0:
		finalFile.close();	

if __name__ == "__main__":
	main()
