'''
Bag handling functions to allow for updating and sorting of bags for the
various representations
'''


def bag_updater(bag, bag_sizes):
    """
    Checks if new bag created is larger than previous bags and updates the
    bag size if it is larger
    Parameters:
    -----------
    bag : dict
        dictionary of all the bags for the current molecule
    bag_sizes : dict
        dictionary of the largest bag sizes in the dataset
    """
    # grab the keys for the bag
    bag_keys = list(bag.keys())
    for i in range(len(bag_keys)):
        key = bag_keys[i]
        # check if the new bag has a larger size than what is in bag_sizes
        #   update the bag_sizes if it is larger, pass if it is not
        if key in bag_sizes:
            if bag[key] > bag_sizes[key]:
                bag_sizes[key] = bag[key]
            else:
                pass
        # if the bag is not in bag_sizes, add it to bag_sizes
        else:
            bag_sizes[key] = bag[key]


def bag_organizer(bag_set, bag_sizes):
    """
    Sorts bags by magnitude, pads, and concactenates into one feature list
    Parameters:
    -----------
    bag_set : dict
        dictionary filled with all of the current molecules information
    bag_sizes : dict
        dictionary of the largest bag sizes in the dataset
    Returns:
    --------
    feat_list : list
        sorted and padded feature list of the current molecule
    """
    feat_list = []
    bag_keys = list(bag_set.keys())
    for i in range(len(bag_keys)):
        # grab the size of the largest bag and length of current molecule bag
        size = bag_sizes[bag_keys[i]] + 1
        baglen = len(bag_set[bag_keys[i]])
        if baglen > (size - 1):
            raise Exception(
                '{}-bag size is too small. Increase size to {}.'.format(bag_keys[i], baglen))
        pad = size - baglen
        # sort the bag by magnitude and pad with zeros to make all same length
        bag_set[bag_keys[i]] = sorted(bag_set[bag_keys[i]], reverse=True)
        bag_set[bag_keys[i]].extend([0.] * pad)
        feat_list.append(bag_set[bag_keys[i]])

    return feat_list
