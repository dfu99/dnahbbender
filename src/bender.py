# -*- coding: utf-8 -*-
"""
Created on Wed Mar 28 21:43:45 2018

@author: Dan

CadnanoManualEdit

changelog:
    12/26/2018:
        whitespace removed from json.dump by adding separator options
"""
import json
import numpy as np
from matplotlib import pyplot as plt
import argparse

AXIAL_RISE = 0.332
INTERHELICAL_DISTANCE = 2.6

def strandnum2idx(data, strandnum):
    """
    Returns the list index of the strand by its num
    """
    for idx, strand in enumerate(data['vstrands']):
        if strand['num'] == strandnum:
            return idx

def insertbase(data, strandidx, loc):
    """
    Edits a strand's 'loop' key to add insertions to the strand
    :param data: json data of DNA bundle
    :param strandidx: The strand index
    :param loc: The index of the base
    :return: json-formated dict with modified data
    """
    data['vstrands'][strandidx]['loop'][loc] = 1  # default 1 insertion
    return data


def deletebase(data, strandidx, loc):
    """
    Edits a strand's 'skip' key to add deletions to the strand
    :param data: json data of DNA bundle
    :param strandidx: The strand index of the base
    :param loc: The index of the base
    :return: json-formated dict with modified data
    """
    data['vstrands'][strandidx]['skip'][loc] = -1
    return data


def jsondatafromfile(filename):
    """
    Reads json data
    :param filename: string input filename, must be in same path as script
    :return: json data in json-formated dict
    """
    json_data = open(filename).read()
    data = json.loads(json_data)
    return data


def datatojson(write_data, filename):
    """
    Writes json-formated dict to .json file
    :param write_data: json-formated dict
    :param filename: string output filename
    :return: None, generates file
    """
    with open(filename, 'w') as outfile:
        json.dump(write_data, outfile, separators=(',', ':'))
    outfile.close()
    print("Written data to", filename)


def distro(n, i, x):
    """
    selects evenly distributed points throughout a list
    :param n: int, nt length of bend section
    :param i: int, index of first base starting the section
    :param x: int, number of elements to distribute
    :return: List of indices to be edited
    """
    x = abs(x)
    if x == 0:
        return []
    elif x == 1:
        return [int(i+n/2)]
    elif x == 2:
        return [int(i+n/3), int(i+2*n/3)]
    elif x >= 3:
        dist = [i, i+n-1]
        for y in range(int(x-2)):
            dist.append(int(i+(y+1)*(n/(x-1))))
        return dist
    else:
        return []

def calc_gradients(mystrands, align_axis=1, theta=0):
    """
    :param hb_locs: x-y coordinates of helical centers
    :param align_axis: x or y normal plane
    :param theta: angle of curvature
    """
    # Calculate the distance of each helix to the normal plane
    # The normal plane is the mean between the min and max of the selected alignment axis
    hb_locs = [strand['pos'] for strand in mystrands.values()]
    hb_locs = np.array(hb_locs)
    d = hb_locs[:, align_axis]
    normal_plane = np.mean(d)

    for key in mystrands.keys():
        mystrands[key]['grad'] = round(np.radians(theta) * (mystrands[key]['pos'][align_axis] - normal_plane) / AXIAL_RISE, 0)
    
    return mystrands

def create_strands(template_data):
    """
    Convert the strands to a dictionary with col, row as keys
    """
    mystrands = {}
    for i in range(len(template_data['vstrands'])):
        this_strand = template_data['vstrands'][i]
        mystrands[(this_strand['col'], this_strand['row'])] = {'strandnum': this_strand['num']}
    return mystrands

def cadnano2_to_xy(mystrands, template_data, lattice_type="hc"):
    # Set the root to start from
    first_strand = template_data['vstrands'][0]
    first_key = (first_strand['col'], first_strand['row'])
    mystrands[first_key]['pos'] = np.array([0, 0])

    # Set the positions of the other strands
    # until all are defined
    num_defined = 1
    curr_key = first_key
    new_keys = []
    while num_defined < len(template_data['vstrands']):
        if lattice_type == "hc":
            mod_key = (curr_key[0] + curr_key[1]) % 2
            offsets = {0: [[0, -5.2], [0, 2.6], [2.25, -1.3], [-2.25, -1.3]], 1: [[0, -2.6], [0, 5.2], [2.25, 1.3], [-2.25, 1.3]]}
        elif lattice_type == "sq":
            mod_key = 0
            offsets = {0: [[0, -2.6], [0, 2.6], [2.6, 0], [-2.6, 0]]}
        else:
            raise ValueError("Invalid lattice type. Use 'hc' for honeycomb or 'sq' for square.")
        directions = [(0, 1), (0, -1), (1, 0), (-1, 0)]

        # Enumerate each adjacent helix
        for i, direction in enumerate(directions):
            new_key = (curr_key[0] + direction[0], curr_key[1] + direction[1])
            # If the new key is in the dictionary and the position is not defined
            if new_key in mystrands and 'pos' not in mystrands[new_key]:
                mystrands[new_key]['pos'] = mystrands[curr_key]['pos'] + np.array(offsets[mod_key][i])
                num_defined += 1
                new_keys.append(new_key) # build stack of new keys
        curr_key = new_keys.pop(0) # pop the next key from the stack
    return mystrands

def plot_lattice(mystrands):
    # Sample data structure
    data = mystrands
    # Create a new figure and axis
    fig, ax = plt.subplots()

    # Iterate through the dictionary
    for item in data.values():
        x, y = item['pos']
        strandnum = item['strandnum']
        
        # Plot the point
        ax.plot(x, y, 'o')
        
        # Add the text label
        ax.text(x, y, strandnum, fontsize=9, ha='right', va='bottom')

    # Set labels and title
    ax.set_xlabel('X coordinate')
    ax.set_ylabel('Y coordinate')
    ax.set_title('Strand Positions')

    # Adjust the plot layout
    plt.tight_layout()

    # Show the plot
    plt.show()

def main():
    # Create parser
    parser = argparse.ArgumentParser(description="Adding insertions and deletions into cadnano2 json file to bend a DNA bundle")

    # Add arguments
    parser.add_argument("filename", help="Input filename of cadnano2 json file")
    parser.add_argument("output_filename", help="Output filename of cadnano2 json file")
    parser.add_argument("lattice_type", help="Type of lattice, honeycomb or square")
    parser.add_argument("angle", help="Desired bend angle")
    parser.add_argument("bend_length", help="Length of bundle to add insertions/deletions")
    parser.add_argument("bend_start", help="Column ID of the first base in the bend")

    # Parse arguments
    args = parser.parse_args()

    # Load the template data
    template_data = jsondatafromfile(args.filename)

    # Create the dictionary of strand positions
    mystrands = create_strands(template_data)
    # Convert the cadnano2 coordinates to x-y coordinates
    mystrands = cadnano2_to_xy(mystrands, template_data, args.lattice_type)
    # Get the edits
    mystrands = calc_gradients(mystrands, theta=int(args.angle))

    bendstarts = [int(args.bend_start)]
    for bst in bendstarts:
        for strand in mystrands.values():
            strandidx = strandnum2idx(template_data, strand['strandnum'])
            for pos in distro(int(args.bend_length), bst, abs(strand['grad'])):
                if strand['grad'] < 0: # deletions
                    template_data = deletebase(template_data, strandidx, pos)
                else: # insertions
                    template_data = insertbase(template_data, strandidx, pos)

    datatojson(template_data, args.output_filename)

if __name__ == "__main__":
    main()
