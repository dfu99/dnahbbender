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


def insertbase(data, strandnum, loc):
    """
    Edits a strand's 'loop' key to add insertions to the strand
    # inputs data - json data of DNA bundle
    #       strandnum - the strand number
    #       loc - the base pair position to be edited
    :param data: json data of DNA bundle
    :param strandnum: The strand index of the base
    :param loc: The index of the base
    :return: json-formated dict with modified data
    """
    data['vstrands'][strandnum]['loop'][loc] = 1  # default 1 insertion
    return data


def deletebase(data, strandnum, loc):
    """
    Edits a strand's 'skip' key to add deletions to the strand
    :param data: json data of DNA bundle
    :param strandnum: The strand index of the base
    :param loc: The index of the base
    :return: json-formated dict with modified data
    """
    data['vstrands'][strandnum]['skip'][loc] = -1
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
    # N is the nt length of the section
    # i is the starting base id
    # x is the number of edits
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


if __name__ == "__main__":
    # Choose and import the tempalte
    template_name = "templates/4hb-512.json"
    template_data = jsondatafromfile(template_name)

    # Copy data to new location
    # (Isn't this useless because dict is mutable?)
    new_data = template_data

    # Set where bends will start
    bendstarts = [211]

    # Set the number of edits to be placed in each strand
    # Keys are strand indices. Values are number of edits
    # Edits: (-) are deletions, (+) are insertions
    edits = {0: 12, 1: 12, 2: -12, 3: -12}
    N = 90
    for bst in bendstarts:
        for strand in edits.keys():
            for pos in distro(N, bst, abs(edits[strand])):
                if edits[strand] < 0:  # deletions
                    new_data = deletebase(new_data, strand, pos)
                else:  # insertions
                    new_data = insertbase(new_data, strand, pos)
        
    datatojson(new_data, "output/4hb-512-180-N90.json")
