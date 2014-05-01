#-------------------------------------------------------------------------------
# Name:        pymol_calc_accessibility.py
#
# Purpose:     Calculates per-residue accessiblity from a pymol object for all
#              atoms and side-chain atoms only and saves to a csv file.
#
# Author:      Tet Woo Lee <tw.lee(at)auckland.ac.nz>
#
# Copyright:   (c) Tet Woo Lee 2011-2014
#
# Licence:     GNU General Public License version 3
#
#              This program is free software: you can redistribute it and/or modify
#              it under the terms of the GNU General Public License as published by
#              the Free Software Foundation, either version 3 of the License, or
#              (at your option) any later version.
#
#              This program is distributed in the hope that it will be useful,
#              but WITHOUT ANY WARRANTY; without even the implied warranty of
#              MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#              GNU General Public License for more details.
#
#              You should have received a copy of the GNU General Public License
#              along with this program.  If not, see <http://www.gnu.org/licenses/>.
#-------------------------------------------------------------------------------

# Tested with PyMOL(TM) Molecular Graphics System, Version 1.6.9.0.
import csv
from collections import OrderedDict
from pymol import stored

"""
Determines the surface area of each residue in a selection in square Angstroms and writes this data to csv file(s).

The surface area for a residue is calculated for all atoms and for side-chain atoms only; side-chain atoms are defined
as those that are not named N, C, or O for glycine, and those that are not named N, C, O or CA for other residues.

One output file is writen for each model (pymol object) in the selection, with the following columns:
* model : pymol object name
* resi : residue identifier (number)
* resn : residue name
* sc_acc : area of side-chain atoms
* all_acc : area of all atoms

Details:

Uses ``Get Area`` in PyMol with ``load_b`` set so this command will alter b-factors of the selection,
see documentation for this command for more information.

By default, calculates the solvent accessible surface area (change with the ``dot_solvent`` parameter,
0 will calculate the molecular surface area).

``dot_density`` controls the accuracy of the measurement,
defaults to ``3``.

``selection`` can be used to restrict surface area calculate to certain selection. Must be a single object for
the ``get_area`` command. By

For simplicity, all solvent, hydrogen and het atoms as well as alternate locations are deleted prior to the calculation.

Note:

Script may remove some atoms, change options and alter b-factors of selection, see above.

"""

def calc_residue_accessibility(model,resi,resn,name,b):
    """
    Processes surface area for an atom, adding appropriate position in the ``stored.acc_data`` dict.
    """
    all_acc = 0.0 # stores all atom acc area
    sc_acc = 0.0 # sc atom acc area

    initial_data = (model,resi,resn,all_acc,sc_acc) # default data for initializing if no data stored in dict
    res_data = stored.acc_data.get(resi,initial_data) # lookup dict for current data or use default
    model,resi,resn,all_acc,sc_acc = res_data # unpack looked-up data

    main_chain_names_GLY = set(("N","C","O")) # name of atoms to include as main chain for glycine
    main_chain_names_etc = set(("N","C","O","CA")) # name of atoms to include as main chain for other residues

    if resn=="GLY": # decide if atom is main chain or not
      main_chain = name in main_chain_names_GLY
    else:
      main_chain = name in main_chain_names_etc

    # add to areas
    all_acc += b
    if not main_chain: sc_acc += b

    print "Processing %s %s-%s atom %3s: exposed area = %5.1f, main-chain = %5r, current sc total area %5.1f, current total area %5.1f"%\
        (model,resn,resi,name,b,main_chain,sc_acc,all_acc)

    # store in dict
    altered_data = (model,resi,resn,all_acc,sc_acc)
    stored.acc_data[resi] = (altered_data)


def calc_accessibility(selection,dot_solvent = 1, dot_density = 3):

    print "Accessibility calculation for selection %s..."%selection

    print "Removing solvent, hydrogens, het atoms and alternate locations..."
    cmd.remove("solvent")
    cmd.remove("hydrogens")
    cmd.remove("het")
    cmd.remove("alt b-z")

    print "Setting options..."
    cmd.set("dot_solvent",dot_solvent)
    cmd.set("dot_density",dot_density)

    print "Calculating surface area..."
    cmd.get_area(selection=selection, load_b=1) # calc surface area: use get_area with load_b to store area as b_factor

    print "Processing each atom..."
    stored.acc_data = OrderedDict() # used to store accesibility data
    cmd.iterate(selection,'calc_residue_accessibility(model,resi,resn,name,b)')
        # iterate through each atom to process the area for the item

    print "Writing data..."
    outfile = None
    writer = None
    for cur_acc_data in stored.acc_data.values():
        model = cur_acc_data[0]
        if writer is None: # no writer for this model, make a new one and writer header
            filename = '%s_accessibility.csv'%model
            print "Opening file %s..."%filename
            outfile = open(filename,'wb')
            writer = csv.writer(outfile)
            writer.writerow(("model","resi","resn","all_acc","sc_acc")) # write header
        # write data for this residue
        writer.writerow(cur_acc_data)

    print "Closing files..."
    if outfile is not None: outfile.close()
    print "Done."

def calc_accessibility_all(dot_solvent = 1, dot_density = 3):
  for object_name in cmd.get_object_list(): # perform calculation for all objects
    calc_accessibility(object_name, dot_solvent, dot_density)