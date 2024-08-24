# SDOCK2.0
Fast Fourier transform (FFT) based global protein-protein docking program
                               SDOCK2.0

   Global Protein-Protein Docking Using Stepwise Force-Field Potentials   
                Written by ZHANG Changsheng 
  Copyright (c) 2024 Molecular Design Lab at Peking University, Beijing, China  
 
    contact: Author: changshengzhang@pku.edu.cn,
             Head of the lab: professor LAI Luhua :lhlai@pku.edu.cn  
 
====================================NOTE========================================
	SDOCK2.0 is free for academic usage and please reference to:
(1) ZHANG Changsheng, LAI Luhua. SDOCK: A Global Protein-Protein Docking Program
 Using Stepwise Force-Field Potentials. Journal of Computational Chemistry 2011,
 32 (12), 2598-2612.
(2) ZHANG Changsheng, LAI Luhua. De novo design of cyclic peptide binder based 
on fragment docking and assembling. 2024
o==============================================================================o

==> PRINCIPLE and ALGORITHMS

SDOCK is a global protein-protein docking program which uses Fast Fourier 
Transform method. The aim of protein-protein docking is to predict the detail 
structure of protein complex from the structure of the unbound components.

The score function of the first version, SDOCK1.0, contains van der Waals 
attractive potential, geometric collision, screened electrostatic potential and 
Lazaridis_Karplus desolvation energy, and in SDOCK2.0, cation-pi energy, 
backbone hydrogen bonding energy, and explicit water exclusion terms. We 
generate stepwise potentials from the corresponding continuous forms to adopt 
the flexibility implicitly.  All the two-particle potential functions are 
exactly or approximately expressed as the form of the product of one atom's 
property and the interacting atom's field.  And atom properties and the fields 
are mapped to 3D grids.  the docking scoring function is expressed as multiple 
correlations of the grids which can be efficiently calculated by Fast Fourier 
Transform method.

The van der Waals and electrostatic parameters are Modified from GROMOS 96 force
field parameters. Solvation parameters are mainly from the parameter table in 
the article of Lazaridis and Karplus. cation-pi and hydrogen bonding energy are 
short range interactions. Cation-pi interaction is considered as the 
electrostatic interaction between a positive charge of amino group or guanidine 
group and the induced negative charge of the aromatic atom.  The water exclusion
 energy was then expressed as the energy summation of all the possible excluded 
waters. we developed a fast method to generate the solvent map, which describes 
the energy of water molecules in the first-shell target surface. The water 
exclusion energy term can be selected to include or exclude the scoring function. 
The optimization of the atomic parameters and the weights of different potential
 terms is based on a docking test set NLC2.0 that contain 486 cases with small 
or moderate conformational change upon binding.

Which protein is fixed and which is rotate and the size of the grids are 
optimized by the program. The rotational space is quasi-uniformly sampled using 
a deterministic layered Sukharev grid sequence suggested by Lindemann et al. 
The program reads the quaternions and then translates them to rotational 
matrices.

A given number of best docking conformations will be saved and then clustered. 
The distance of any two solutions is evaluated by the root mean square deviation
 of interface Ca atoms of the moving protein.

To learn more about the algorithm, please read the reference papers.

==> INSTALLATION

This software can be used on LINUX platform. 
(1) Download fftw from http://www.fftw.org/, for example fftw-3.3.8. 
    Enter to directory /fftw-3.3.8. Read INSTALL and Install the fftw-3.3.8 
    software for fast Fourier Transform calculation.
(2) Enter to directory /src, open Makefile. Change the FFTW_DIR variable to the
    correct fftw library route.
(3) Run Make for compiling the source codes, and get 4 executable files, 
    preprocess, watmap, sdock, and build.
(4) cp src/preprocess bin/; cp src/watmap bin/; cp src/sdock bin/; 
    cp src/build bin/;
(5) Enter to directory /example. Run ./1AY7dock.sh and ./2REXdock.sh for testing.

==> PROTEIN-PROTEIN DOCKING PROCEDURE

I.  Prepare: The receptor and ligand structure files (PDB format). Change N 
    terminal 'N' and 'CA' to 'NT' and 'CAT'; change C terminal 'C' and 'O' to  
    'CT' and 'OT'.
II. Preprocess: add force field parameters to protein atoms and assign surface 
    atoms.
III. Build water map, if choose the explicit water mode.
IV. Docking: docking the two protein together, save the rotational and 
    translational information of the best docking results, and cluster the 
    docking results.
V. Build Models: Get the 3D complex model for your interest docking results.

  The executable binary file 'preprocess', 'watmap', 'sdock' and 'build' can do 
  the task II, III, IV, and V respectively. See the following manual for detail
  usage. Please have a try with the two example cases under the example 
  directory.  Just run the commands in file '1AY7dock.sh' and '2REXdock.sh' to 
  quickly learn about this software.

==> MANUAL

PREPROCESS: preprocess the origin protein structure file (PDB format) to the 
SDOCK input file

 Usage: preprocess pdb_file 
    [-c chain_id (all chains in the pdb_file)]  
    [-o output_file (preprocessed.pdb)]  
    [-d dynamic_or_not (0)]  
    [-u uncommon_amino_acid (no)] 
    [-s small_molecule (no)]  
    [-a atom_parameter_file (ATM)]

 Format of uncommon_amino_acid or small_molecule: 
    chain_id1+name1:chain_id2+name2:chain_id3+name3:...  '@' means any chain,
    '_' means blank
 
 Examples:
 preprocess 1SBB.pdb
 preprocess 1SBB.pdb -u no -d 1    #consider dynamic using GNM model
 preprocess 1FKM.pdb -c A -o 2G77L.pdb -u AMSE
 preprocess 1EXB.pdb -o 1EXBR.pdb -s  @NDP -a myATM
 preprocess 1ZMC.pdb -c AB  -o 2F5ZR.pdb -s  AFAD:BFAD -a ../data/ATM  
           # A and B chains

WATMAP: build the water map on the protein surface
 Usage: watmap preprocessed_protein_structure watermap_for_protein

 Examples:
 watmap 1AY7R.pdb 1AY7Rwat.pdb 

SDOCK: protein-protein docking and save results.
 Usage: sdock preprocessed_protein_structure1 preprocessed_protein_structure2 
    [watermap_for_protein1] [watermap_for_protein2] 
    [-o sdock_result_file_name (sdock_record)]  
    [-p progress with percentage (0.1)]  
    [-r rotation_sampling_file (so3layer.qua)] 
    [-c weight_of_collision_term (0.17)]  
    [-e weight_of_electrostatic_term (1.4)]
    [-s weight_of_desolvation_term (0.67/0.60)]  
    [-i weight_of_induction_term (4.0)]  
    [-H weight_of_backboneO_HB_term] (7.0)
    [-h weight_of_backboneN_HB_term (0.5/1.0)]  
    [-w weight_of_explicit_water_term (0.0/1.0)]  
    [-n cluster_number (1000)]  
    [-d cluster_radius (3.0)]

 Examples:
 sdock 1AY7R.pdb 1AY7L.pdb                       
      #no watmap, explicit water term is not included
 sdock 1AY7R.pdb 1AY7L.pdb 1AY7Rwat.pdb 1AY7Lwat.pdb
      #watmaps are inputed, explicit water mode
 sdock 1AY7R.pdb 1AY7L.pdb -o 1AY7_dock_result -p 0
      #result file name and no progress is displayed
 sdock 1AY7R.pdb 1AY7L.pdb -o 1AY7_dock_result -H 8.0
      #a larger weight for backbone oxygen hydrogen bond term
 sdock 1AY7R.pdb 1AY7L.pdb 1AY7Rwat.pdb 1AY7Lwat.pdb -s 0.63
      #a larger desolvation weight in explicit water mode
 sdock 1AY7R.pdb 1AY7L.pdb -n 100 -d 2.0    
      #only output the top 100 solutions with a smaller cluster radius

BUILD: build models according to the docking records.
 Usage: preprocessed_protein_structure1 preprocessed_protein_structure2 
    [-o  sdock_result_file_name (sdock_record)]
    [-r rotation_sampling_file (so3layer.qua)] 
    [-c cluster_number (1)] 
    [-m cluster_member_number (1)]
    [-n four_letter_model_name (AAAA)]
    [-d model_file_route (./)]
 
 Examples:
 build 1AY7R.pdb 1AY7L.pdb
 build 1AY7R.pdb 1AY7L.pdb -n 1ay7
      #The first 4 letter of file name of the built complex model is 1ay7
 build 1AY7R.pdb 1AY7L.pdb -r ../data/so3layer.qua
 build 1AY7R.pdb 1AY7L.pdb -o 1AY7_dock_resultf -c 0 -m 1 -n 1AY7
      #The first structure of all clusters
 build 1AY7R.pdb 1AY7L.pdb -o 1AY7_dock_resultf -c 0 -m 0 -n 1AY7 -d model/
      #build the structure of all results to the direction model/

==================================NOTE===================================================

Any questions and problems, please contact the author: changshengzhang@pku.edu.cn.
Thanks for your support!

===================================END===================================================

