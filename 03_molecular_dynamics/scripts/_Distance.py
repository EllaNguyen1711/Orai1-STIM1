import os
import numpy as np
import MDAnalysis as mda
import pickle as pkl

def dis_Tyr208_His57(psf, trajs, start, end, pairs = None):
    """
    Calculates the distance between the centroids of Tyr208 and His57 
    aromatic rings across 6 fixed subunit pairs.
    """
    u = mda.Universe(psf, trajs)

    if pairs is None:
        pairs = [
                ("PROA", "PROQ"), ("PROB", "PROG"), ("PROC", "PROI"),
                ("PROD", "PROK"), ("PROE", "PROM"), ("PROF", "PROO"),
            ]
    else:
        pairs = pkl.load(open(pairs, 'rb'))

    # Pre-select the ring atom groups
    tyr_rings = []
    his_rings = []

    for seg_tyr, seg_his in pairs:
        # Tyrosine ring: 6 Carbons
        tyr = u.select_atoms(f"segid {seg_tyr} and resid 208 and name CG CD1 CD2 CE1 CE2 CZ")
        # Histidine ring: 3 Carbons, 2 Nitrogens
        his = u.select_atoms(f"segid {seg_his} and resid 57 and name CG ND1 CD2 CE1 NE2")

        if len(tyr) != 6 or len(his) != 5:
            print(f"Warning: Expected 6 Tyr/5 His atoms, got {len(tyr)}/{len(his)} for {seg_tyr}-{seg_his}")
            # Depending on forcefield, names might slightly differ (e.g., HSC vs HIS)
        
        tyr_rings.append(tyr)
        his_rings.append(his)

    n_frames = end - start
    # We will store the minimum distance found among the 6 pairs for each frame
    min_distances = np.empty(n_frames, dtype=np.float32)

    for i, ts in enumerate(u.trajectory[start:end]):
        dvals = []

        for tyr, his in zip(tyr_rings, his_rings):
            # Calculate Centroids (geometric center of the ring atoms)
            # .center_of_geometry() is faster than center_of_mass() for this purpose
            tyr_centroid = tyr.center_of_geometry()
            his_centroid = his.center_of_geometry()
            
            # Euclidean distance between centroids
            d = np.linalg.norm(tyr_centroid - his_centroid)
            dvals.append(d)

        min_distances[i] = np.min(dvals)

    return min_distances


def dis_Tyr208_Phe53(psf, trajs, start, end, pairs = None):
    """
    6 fixed cross-subunit Tyr208–Phe53 pairs
    Returns minimum COM distance among the 6 pairs per frame.
    Tyr208 aromatic ring: CG CD1 CD2 CE1 CE2 CZ
    Phe53  aromatic ring: CG CD1 CD2 CE1 CE2 CZ
    """
    u = mda.Universe(psf, trajs)
    if pairs is None:
        pairs = [
            ("PROA", "PROQ"),
            ("PROB", "PROG"),
            ("PROC", "PROI"),
            ("PROD", "PROK"),
            ("PROE", "PROM"),
            ("PROF", "PROO"),
            ]
    else:
        pairs = pkl.load(open(pairs, 'rb'))
    
    tyr_ring_atoms = "name CG CD1 CD2 CE1 CE2 CZ"
    phe_ring_atoms = "name CG CD1 CD2 CE1 CE2 CZ"

    tyr_atoms = []
    phe_atoms = []
    for seg_tyr, seg_phe in pairs:
        tyr = u.select_atoms(f"segid {seg_tyr} and resid 208 and ({tyr_ring_atoms})")
        phe = u.select_atoms(f"segid {seg_phe} and resid 53  and ({phe_ring_atoms})")
        if len(tyr) != 6 or len(phe) != 6:
            raise ValueError(
                f"Expected 6 ring atoms each, got Tyr={len(tyr)} Phe={len(phe)} "
                f"for segids {seg_tyr}-{seg_phe}"
            )
        tyr_atoms.append(tyr)
        phe_atoms.append(phe)

    n_frames = end - start
    dist = np.empty(n_frames, dtype=np.float32)

    for i, ts in enumerate(u.trajectory[start:end]):
        dvals = []
        for tyr, phe in zip(tyr_atoms, phe_atoms):
            com_tyr = tyr.center_of_mass()
            com_phe = phe.center_of_mass()
            d = np.linalg.norm(com_tyr - com_phe)
            dvals.append(d)
        dist[i] = np.min(dvals)
        # dist[i] = np.mean(dvals)
    return dist
