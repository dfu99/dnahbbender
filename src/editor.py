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

def resetbase(data, strandidx, loc):
    """
    Edits a strand's 'skip' and 'loop' keys to 0
    :param data: json data of DNA bundle
    :param strandnum: The strand index of the base
    :param loc: The index of the base
    :return: json-formated dict with modified data
    """
    data['vstrands'][strandidx]['skip'][loc] = 0
    data['vstrands'][strandidx]['loop'][loc] = 0
    return data