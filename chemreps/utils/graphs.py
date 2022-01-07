"""
Graph based functions for angles and torsions
"""


def gen_graph(connect, n_atom):
    '''
    Generates connectivity graph

    Parameters
    ---------
    connect: list
        list of bond connectivity 
    n_atom : int
        number of atoms

    Returns
    -------
    graph: dict
        graph of entire system
    '''
    graph = {}
    for atom in range(1, n_atom+1):
        aconn = []
        for conn in connect:
            if atom in conn:
                c = list(conn)
                c.remove(float(atom))
                aconn.append(int(c[0]))
        graph.update({atom: aconn})

    return graph


def dfs_connections(graph, start, length, connections, connection=[]):
    '''
    Searches for all connected atoms in an angle or torsion using a depth first search
    Adapted from:
    https://stackoverflow.com/a/42232789

    Parameters
    ---------
    graph: dict
        graph of entire system
    start: int
        position of start of the torsion
    length: int
        number of atoms searched over (eg 3 for angle, 4 for torsion)
    connections: list
        empty list of all connections to fill in depth first search
    connection: list
        empty list for each individual connection

    '''
    # start or append connection list
    connection = connection + [start]
    # check for desired length
    if len(connection) == length:
        # ensure no looping over the same atom
        if len(list(set(connection))) == length:
            # check if reversed torsion is in list
            if connection[::-1] not in connections:
                connections.append(connection)
    else:
        # loop over bonds
        for bond in graph[start]:
            dfs_connections(graph, bond, length, connections, connection)
