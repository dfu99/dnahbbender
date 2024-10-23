# -*- coding: utf-8 -*-
"""
Created on Wed Mar 28 16:18:00 2018

@author: Dan

ShiftCorrect.py moves inserts or deletions off crossover positions

Shifts to the closest non-crossover position

changelog:
    12/26/2018:
        whitespace removed from json.dump by adding separator options
"""

import json


def getlength(data):
    """
    Gets the length of the json data
    :param data: json data
    :return: int, length
    """
    return len(data['vstrands'][0]['scaf'])


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


def delbase(data, strandnum, loc):
    """
    Edits a strand's 'skip' key to add deletions to the strand
    :param data: json data of DNA bundle
    :param strandnum: The strand index of the base
    :param loc: The index of the base
    :return: json-formated dict with modified data
    """
    data['vstrands'][strandnum]['skip'][loc] = -1
    return data


def resetbase(data, strandnum, loc):
    """
    Edits a strand's 'skip' and 'loop' keys to 0
    :param data: json data of DNA bundle
    :param strandnum: The strand index of the base
    :param loc: The index of the base
    :return: json-formated dict with modified data
    """
    data['vstrands'][strandnum]['skip'][loc] = 0
    data['vstrands'][strandnum]['loop'][loc] = 0
    return data


def getcrossoverplane(data, pos, numstrands):
    """
    Checks if crossovers are placed in a cross-section of the bundles
    :param data: json data
    :param pos: base index position
    :return: True or False
    """
    for strandnum in range(numstrands):
        connectingstrands = []
        connectingstrands.append(data['vstrands'][strandnum]['scaf'][pos][0])
        connectingstrands.append(data['vstrands'][strandnum]['scaf'][pos][2])
        connectingstrands.append(data['vstrands'][strandnum]['stap'][pos][0])
        connectingstrands.append(data['vstrands'][strandnum]['stap'][pos][2])
        if any(x not in [-1, data['vstrands'][strandnum]['num']] for x in connectingstrands):
            print("Crossover plane @ ", pos)
            return True
    return False


def getshiftdirection(data, pos, numstrands):
    """
    Checks whether the base should be moved left or right
    :param data: json data
    :param pos: current base index position
    :return: Ture or False
    """
    for strandnum in range(numstrands):
        connectingstrands = []
        connectingstrands.append(data['vstrands'][strandnum]['scaf'][pos+1][0])
        connectingstrands.append(data['vstrands'][strandnum]['scaf'][pos+1][2])
        connectingstrands.append(data['vstrands'][strandnum]['stap'][pos+1][0])
        connectingstrands.append(data['vstrands'][strandnum]['stap'][pos+1][2])
        if any(x not in [-1, data['vstrands'][strandnum]['num']] for x in connectingstrands):
            return True  # Crossover is to the right, shift to the left
    return False  # Crossover is to the left, shift to the right


if __name__ == "__main__":
    
    # 0 deletions - 3 insertions
    # 1,5 deletions - 2,4 insertions

    filename = "output/4hb-512-180-N60"
    json_data = open(filename+".json").read()
    data = json.loads(json_data)
    
    # the workspace size
    jsonlength = getlength(data)

    # the number of strands
    number_of_strands = len(data['vstrands'])
    print("Number of strands =", number_of_strands)

    # look at each base pair left to right
    for pos in range(jsonlength): 
        # look at each helix
        for i, row in enumerate(data['vstrands'], start=0):
            # check if there is an edit
            if row['skip'][pos] == -1 or row['loop'][pos] == 1:
                # check if we're on a crossover plane
                if getcrossoverplane(data, pos, number_of_strands):
                    # check which direction to shift
                    if getshiftdirection(data, pos, number_of_strands):
                        # shift left
                        # if i in [0,1,5]:
                        if row['skip'][pos] == -1:
                            resetbase(data, i, pos)
                            delbase(data, i, pos-1)
                        else:
                            resetbase(data, i, pos)
                            insertbase(data, i, pos-1)
                    else:
                        # shift right
                        # if i in [0,1,5]:
                        if row['skip'][pos] == -1:
                            resetbase(data, i, pos)
                            delbase(data, i, pos+1)
                        else:
                            resetbase(data, i, pos)
                            insertbase(data, i, pos+1)
                else:
                    pass
            else:
                pass
            
    with open(filename+'-shiftcorrected.json', 'w') as outfile:
        json.dump(data, outfile, separators=(',', ':'))
    outfile.close()
