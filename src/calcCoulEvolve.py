#!/usr/bin/env python
"""
calcCoulEvolve.py

Calculates the total coulomb energy (kcal/mol) of an evolved protein sequence given
its 1-letter abbreviated sequence and corresponding pdb structure.  User inputs:
    PDB_FILE: pdb file containing structure
    DIELEC_CONST: dielectric constant of the system
    IONIC_STR: ionic strength of the system in molar
    TEMP_K: temperature of the system in Kelvin

- Read "ABD..." sequence string
- Record numbers of charged residues in the sequence (make a list) and the value
  of their charges (make a list)
- Get coordinates for CBs of all charged residues
- Calculate r and get energy as usual
"""

__author__ = "Abraham J. Rickett"
__date__ = "160719"


# *************************************************************************** #
#                                User inputs                                  #
# *************************************************************************** #
from sys import argv

#EVO_SEQ = argv[1]
PDB_FILE = "1stn.pdb"
DIELEC_CONST = 20
IONIC_STR = 0.1
TEMP_K = 300


# *************************************************************************** #
#                           Function defintions                               #
# *************************************************************************** #

from math import sqrt,exp

def readPDB_evo(pdb_file):
    """
    Takes a pdb file and reads in the coordinates of each beta carbon.
    **Alpha carbons are used instead for glycine residues.**
    """

    # Initalize list to hold coordinates
    all_coord = []

    # Open pdb_file and read each line into pdb (a list of lines)
    f = open(pdb_file,'r')
    pdb = f.readlines()
    f.close()

    # Go through each line in the pdb file
    for line in pdb:

        # Only look at a line if it starts with "ATOM". Record the type of
        # atom and the amino acid the atom is in.
        if line[0:4] == "ATOM":
            amino_acid = line[17:20]
            atom = line[13:16]

            # If the amino acid is not a glycine, append the xyz coordinates of
            # its CB to all_coord
            if amino_acid != "GLY":
                if atom == "CB ":
                    all_coord.append([float(line[30:38]),
                                      float(line[38:46]),
                                      float(line[46:54])])

            # For glycine residues, find the alpha carbon ...
            else:
                if atom == "CA ":

                    # ... and put the xyz coordinates in all_coord
                    all_coord.append([float(line[30:38]),
                                      float(line[38:46]),
                                      float(line[46:54])])

    # Return coordinates
    return all_coord

def readPDB_anc(pdb_file):
    """
    Takes a pdb file and reads in the coordinates of all CB atoms within
    titratable residues. This list of coordinates is used to calculate the
    energy of the ancestral structure.
    """
    # Initialize lists/dictionaries
    anc_coord = [] # coordinates of just the CBs for the charged residues in the
                   # ancestral structure
    anc_charge = [] # signs of charged residues in ancestral structure

    titr_res = ["ASP","GLU","TYR","ARG","HIS","LYS"]
    titr_charge = {"ASP":-1.0,"GLU":-1.0,"TYR":-1.0,"ARG":1.0,"HIS":1.0,"LYS":1.0}

    # Open pdb_file and read each line into pdb (a list of lines)
    f = open(pdb_file,'r')
    pdb = f.readlines()
    f.close()

    # Go through each line in the pdb file
    for line in pdb:

        # Record each amino acid and atom name
        if line[0:4] == "ATOM":
            amino_acid = line[17:20]
            atom = line[13:16]

            # If the amino acid is titratable, append its CB coordinates to
            # anc_coord and its charge to anc_charge
            if amino_acid in titr_res:
                if atom == "CB ":
                    anc_charge.append(titr_charge[amino_acid])
                    anc_coord.append([float(line[30:38]),
                                      float(line[38:46]),
                                      float(line[46:54])])

    # Return the coordinates
    print anc_charge
    return anc_coord, anc_charge

def readSEQ(evo_seq):
    """
    Reads a string containing an evolved sequence and records which residues are
    charged.
    """

    # Initialize lists to hold charge signs and the index of charged residues.
    # Start count at -1 so first iteration is 0.
    evo_charge = []
    charged_res = []
    titr_charge2 = {"D":-1.0,"E":-1.0,"H":1.0,"R":1.0,"K":1.0,"Y":-1.0}
    titr_aa = titr_charge2.keys()
    count = -1

    # Cycle through each letter of the evolved sequence and identify the charged
    # residues
    for a in evo_seq:
        count += 1
        if a in titr_aa:
            evo_charge.append(titr_charge2[a])
            charged_res.append(count)

    # Return the charge lists
    return evo_charge, charged_res

def calcCoulomb(q1,q2,r12,dielec_const,ionic_str,temp_k):
    """
    Calculates the Coulomb energy of a charge-charge interaction,incorporating
    the Debye-Huckel theory.
    """

    return 332*q1*q2/(r12*dielec_const)*exp(-50.29*sqrt(ionic_str/(dielec_const*temp_k)))


def calculateEnergy(anc_coord,all_coord,anc_charge,evo_charge,charged_res,dielec_const,ionic_str,temp_k):
    """
    Calculates the energy of a structure given the coordinates of each
    charged atom, their fractional charge, and the dielectric constant of the
    system.
        E = 332*sum_i[sum_j[q_i*q_j/(r_ij * dielec_const)]] (kcal/mol)
    """

    # Initialize variables
    num_evo_groups = len(evo_charge)
    evo_energy = 0.

    # Calculate energy of interaction of every ij interaction (making sure not
    # to double count; note we start j at i + 1).
    for i in range(num_evo_groups):
        for j in range(i+1,num_evo_groups):

            # Calculate distance between atom i and atom j
            evo_r = (all_coord[charged_res[i]][0] - all_coord[charged_res[j]][0])**2
            evo_r = evo_r + (all_coord[charged_res[i]][1] - all_coord[charged_res[j]][1])**2
            evo_r = evo_r + (all_coord[charged_res[i]][2] - all_coord[charged_res[j]][2])**2
            evo_r = sqrt(evo_r)
            #print (i,j,evo_r)
            # Add the energy of this interaction to the total energy
            evo_energy = evo_energy + calcCoulomb(evo_charge[i],evo_charge[j],
                                                  evo_r,dielec_const,ionic_str,
                                                  temp_k)

    # Initialize variables
    num_anc_groups = len(anc_charge)
    anc_energy = 0.

    # Calculate energy of interaction of every ij interaction (making sure not
    # to double count; note we start j at i + 1).
    for k in range(num_anc_groups):
        for l in range(k+1,num_anc_groups):

            # Calculate distance between atom i and atom j
            anc_r = (anc_coord[k][0] - anc_coord[l][0])**2
            anc_r = anc_r + (anc_coord[k][1] - anc_coord[l][1])**2
            anc_r = anc_r + (anc_coord[k][2] - anc_coord[l][2])**2
            anc_r = sqrt(anc_r)
            #print (k,l,anc_r)
            # Add the energy of this interaction to the total energy
            anc_energy = anc_energy + calcCoulomb(anc_charge[k],anc_charge[l],
                                                  anc_r,dielec_const,
                                                  ionic_str,temp_k)

    # Return energy
    return evo_energy, anc_energy

class AncCoul:

    def __init__(self,PDB_FILE):

         self.all_coord = readPDB_evo(PDB_FILE)
         self.anc_coord, self.anc_charge = readPDB_anc(PDB_FILE)

    def calcAncEnergy(self,some_seq):
        evo_charge, charged_res = readSEQ(some_seq)
        evo_energy, anc_energy = calculateEnergy(self.anc_coord,
                                                 self.all_coord,
                                                 self.anc_charge,
                                                 evo_charge,charged_res,
                                                 DIELEC_CONST,IONIC_STR,TEMP_K)

        # Find the difference between them
        deltaE = anc_energy - evo_energy
        deltaE_p = deltaE/anc_energy
        return deltaE_p


X = AncCoul(PDB_FILE)


# *************************************************************************** #
#                               Main code                                     #
# *************************************************************************** #

# Read in coordinates and charge lists
"""
all_coord = readPDB_evo(PDB_FILE)
anc_coord, anc_charge = readPDB_anc(PDB_FILE)
evo_charge, charged_res = readSEQ(EVO_SEQ)

# Determine electrostatic energies for ancestral and evolved sequences
evo_energy, anc_energy = calculateEnergy(anc_coord,all_coord,anc_charge,evo_charge,charged_res,DIELEC_CONST,IONIC_STR,TEMP_K)

# Find the difference between them
deltaE = anc_energy - evo_energy

return deltaE
"""
