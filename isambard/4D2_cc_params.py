#Wood, Christopher W., et al. "ISAMBARD: an open-source computational environment for biomolecular analysis, modelling and design." Bioinformatics 33.19 (2017): 3043-3050.
#
#ISAMBARD analysis of 4D2 helical pairs
#
import isambard
fourD2ac = isambard.ampal.convert_pdb_to_ampal("4D2_finalac.pdb") reference_axis = isambard.analyse_protein.reference_axis_from_chains(fourD2ac) crangles_a = isambard.analyse_protein.crick_angles(fourD2a,reference_axis) crangles_c = isambard.analyse_protein.crick_angles(fourD2c,reference_axis)

#INTERFACE ANGLE
import numpy
phica_a_list = [crangles_a[x] for x in range(0,len(crangles_a),7) if crangles_a[x] is not None]
phica_c_list = [crangles_c[x] for x in range(0,len(crangles_c),7) if crangles_c[x] is not None]
phica_a = numpy.mean(phica_a_list)
phica_c = numpy.mean(phica_c_list)
i_angle = (phica_a+phica_c)/2

#RADIUS
radius_a_list = isambard.analyse_protein.polymer_to_reference_axis_distances(fourD2a, reference_axis)
radius_c_list = isambard.analyse_protein.polymer_to_reference_axis_distances(fourD2c, reference_axis)
radius_a = numpy.mean(radius_a_list)
radius_c = numpy.mean(radius_c_list)
radius = numpy.mean(radius_a_list+radius_c_list)

#PITCH
alpha_a_list = isambard.analyse_protein.alpha_angles(fourD2a, reference_axis) 
alpha_c_list = isambard.analyse_protein.alpha_angles(fourD2c, reference_axis) 
pitch_a_list = [(2* numpy.pi * radius) / numpy.tan(numpy.deg2rad(x)) for x in alpha_a_list if x is not None]
pitch_c_list = [(2* numpy.pi * radius) / numpy.tan(numpy.deg2rad(x)) for x in alpha_c_list if x is not None]
pitch = numpy.mean(pitch_a_list + pitch_c_list)
