import numpy as np


# -------------- edgelist methods --------------------------------------------------------------
def nodelist_like(shape, nhood):
    nEdge = nhood.shape[0]
    nodes = np.arange(np.prod(shape), dtype=np.uint64).reshape(shape)
    node1 = np.tile(nodes, (nEdge, 1, 1, 1))
    node2 = np.full(node1.shape, -1, dtype=np.uint64)

    for e in range(nEdge):
        node2[e, \
        max(0, -nhood[e, 0]):min(shape[0], shape[0] - nhood[e, 0]), \
        max(0, -nhood[e, 1]):min(shape[1], shape[1] - nhood[e, 1]), \
        max(0, -nhood[e, 2]):min(shape[2], shape[2] - nhood[e, 2])] = \
            nodes[max(0, nhood[e, 0]):min(shape[0], shape[0] + nhood[e, 0]), \
            max(0, nhood[e, 1]):min(shape[1], shape[1] + nhood[e, 1]), \
            max(0, nhood[e, 2]):min(shape[2], shape[2] + nhood[e, 2])]

    return (node1, node2)


def affgraph_to_edgelist(aff, nhood):
    num_vert = aff.shape[1] * aff.shape[2] * aff.shape[3]
    node1, node2 = nodelist_like(aff.shape[1:], nhood)
    node1, node2, aff = node1.ravel(), node2.ravel(), aff.ravel()
    # discard illegal vertices
    valid = np.logical_and.reduce((node1 > 0, node2 > 0, node1 < num_vert,
                                   node2 < num_vert, node1 != node2))
    return node1[valid], node2[valid], aff[valid]


def mknhood3d(radius=1):
    ceilrad = np.ceil(radius)
    x = np.arange(-ceilrad, ceilrad + 1, 1)
    y = np.arange(-ceilrad, ceilrad + 1, 1)
    z = np.arange(-ceilrad, ceilrad + 1, 1)
    [i, j, k] = np.meshgrid(z, y, z)

    idxkeep = (i ** 2 + j ** 2 + k ** 2) <= radius ** 2
    i = i[idxkeep].ravel();
    j = j[idxkeep].ravel();
    k = k[idxkeep].ravel();
    zeroIdx = np.ceil(len(i) / 2).astype(np.int32);

    nhood = np.vstack((k[:zeroIdx], i[:zeroIdx], j[:zeroIdx])).T.astype(np.int32)
    return np.ascontiguousarray(np.flipud(nhood))
