"""
Created on Wed Mar 28 16:18:00 2018

@author: Dan

Moves inserts or deletions off crossover positions
Shifts to the closest non-crossover position
"""

from .editor import insertbase, deletebase, resetbase

def getlength(data):
    """
    Gets the length of the json data
    :param data: json data
    :return: int, length
    """
    return len(data['vstrands'][0]['scaf'])

def getcrossoverplane(data, pos, numstrands, verbose=False):
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
            if verbose:
                print("Shifting edits away from a crossover plane @ ", pos)
            return True
    return False

def getshiftdirection(data, pos, numstrands):
    """
    Checks whether the base should be moved left or right
    :param data: json data
    :param pos: current base index position
    :return: True or False
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

def shiftcorrect(data, verbose=False):
    """
    Shifts insertions and deletions off crossover planes
    :param data: cadnano2 format json data
    :return: json-formated dict with modified data
    """
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
                if getcrossoverplane(data, pos, number_of_strands, verbose):
                    # check which direction to shift
                    if getshiftdirection(data, pos, number_of_strands):
                        # shift left
                        if row['skip'][pos] == -1:
                            resetbase(data, i, pos)
                            deletebase(data, i, pos-1)
                        else:
                            resetbase(data, i, pos)
                            insertbase(data, i, pos-1)
                    else:
                        # shift right
                        if row['skip'][pos] == -1:
                            resetbase(data, i, pos)
                            deletebase(data, i, pos+1)
                        else:
                            resetbase(data, i, pos)
                            insertbase(data, i, pos+1)
    return data
