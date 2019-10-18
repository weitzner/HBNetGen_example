#Make sure Tyr resid is set appropriately in line
#tyr_coords = np.array([mdl[2].coords for mdl in models])
#to cluster on Tyr

from collections import namedtuple
import attr
import logging
import numpy as np
import os


logging.basicConfig(level=logging.INFO)


@attr.s
class AtomRecord(object):
    record_name = attr.ib()
    serial = attr.ib()
    name = attr.ib()
    altLoc = attr.ib()
    resName = attr.ib()
    chainID = attr.ib()
    resSeq = attr.ib()
    iCode = attr.ib()
    x = attr.ib()
    y = attr.ib()
    z = attr.ib()
    occupancy = attr.ib()
    tempFactor = attr.ib()
    element = attr.ib()
    charge = attr.ib()

    @classmethod
    def from_str(cls, record):
        record_name = record[:6].strip()
        serial = int(record[6:11].strip())
        name = record[12:16].strip()
        altLoc = record[16].strip()
        resName = record[17:20].strip()
        chainID = record[21].strip()
        resSeq = int(record[22:26].strip())
        iCode = record[26].strip()
        x = float(record[30:38].strip())
        y = float(record[38:46].strip())
        z = float(record[46:54].strip())
        occupancy, tempFactor, element, charge = [None] * 4
        """
        try:
            occupancy = float(record[54:60].strip())
            tempFactor = float(record[60:66].strip())
            element = record[76:78].strip()
            charge = record[78:80].strip()
        except (IndexError, ValueError):
            pass
        """
        return cls(
            record_name,
            serial,
            name,
            altLoc,
            resName,
            chainID,
            resSeq,
            iCode,
            x,
            y,
            z,
            occupancy,
            tempFactor,
            element,
            charge,
        )

    def __str__(self):
        return "{record_name:6}{serial:5d} {name:4}{altLoc:1}{resName:3} {chainID:1}{resSeq:>4}{iCode:1}   {x:8.3f}{y:8.3f}{z:8.3f}{occupancy:6.2}{tempFactor:6.2}          {element:2}{charge:2}".format(
            record_name=self.record_name,
            serial=self.serial,
            name=self.name,
            altLoc=self.altLoc,
            resName=self.resName,
            chainID=self.chainID,
            resSeq=self.resSeq,
            iCode=self.iCode,
            x=self.x,
            y=self.y,
            z=self.z,
            occupancy=self.occupancy if self.occupancy is not None else "",
            tempFactor=self.tempFactor if self.tempFactor is not None else "",
            element=self.element if self.element is not None else "",
            charge=self.charge if self.charge is not None else "",
        )


@attr.s
class Residue(object):
    atom_records = attr.ib()
    coords = attr.ib()

    @classmethod
    def from_records(cls, atom_records):
        coords = np.array([[ar.x, ar.y, ar.z] for ar in atom_records])
        return cls(atom_records, coords)

    @property
    def name(self):
        names = set(ar.resName for ar in self.atom_records)
        assert(len(names) == 1)
        return names.pop()

    def __str__(self):
        return "\n".join([str(ar) for ar in self.atom_records])


def read_in_stubs_file(fname):
    models = []
    model = []
    curr_res = []
    resSeq = 0

    with open(fname, "r") as f:
        for l in f:
            if not l:
                continue

            if l.startswith("MODEL"):
                # clear previous model data
                model = []
                curr_res = []
                resSeq = 0
                continue

            if l.startswith("ENDMDL"):
                # finalize current model
                model.append(Residue.from_records(curr_res))
                models.append(model)
                continue

            record = AtomRecord.from_str(l)
            if record.resSeq == resSeq:
                curr_res.append(record)
            else:
                if curr_res:
                    model.append(Residue.from_records(curr_res))
                curr_res = [record]
                resSeq = record.resSeq
    return models

def process_stubs(fname, n_clusters=100):

    base_name, ext = os.path.splitext(fname)
    if not os.path.exists(base_name):
        os.makedirs(base_name)

    logging.info("Reading in file...")
    # filter out silly broken models from the list as soon as they are returned
    models = [mdl for mdl in read_in_stubs_file(fname) if len(mdl) == 5]
    logging.info("File processed!")

    logging.info("Extracting TYR coords...")
    tyr_coords = np.array([mdl[2].coords for mdl in models])
    logging.info("TYR coords exrtracted!")

    # flatten array so each residue is one coordinate in a higher-dimensional space
    logging.info("Resizing coord array for clustering...")
    tyr_coords = tyr_coords.reshape(
        tyr_coords.shape[:-2] + (tyr_coords.shape[-1] * tyr_coords.shape[-2],)
    )
    logging.info("Done!")

    # cluster on the TYR coordinates
    logging.info("clustering...")
    from sklearn.cluster import DBSCAN
    clustering = DBSCAN(eps=0.25, metric="euclidean", n_jobs=-1).fit(tyr_coords)
    logging.info("Done!")

    from itertools import groupby

    # identify top N clusters
    top_n_clusters = sorted(
        [(key, len(list(group))) for key, group in groupby(sorted(clustering.labels_))],
        key=lambda x: x[-1],
        reverse=True,
    )[:n_clusters]
    # "noisy" samples are given the label -1, so if -1 is included in the top N,
    # we should remove it
    top_n_clusters = [cluster for cluster in top_n_clusters if cluster[0] >= -1]

    # collect groups of model IDs for the top N clusters
    # ensure they are consistent with the previous step
    ntwrks = []
    for clusID, n_occ in top_n_clusters:
        model_indices_for_network = [
            i for i, elem in enumerate(clustering.labels_) if elem == clusID
        ]
        assert len(model_indices_for_network) == n_occ
        ntwrks.append(model_indices_for_network)

    # write each network cluster to a separate PDB-formatted file
    for i, ntwrk in enumerate(ntwrks):
        with open(os.path.join(base_name, "grp_ntwrk_{}_{:04d}".format(base_name, i) + ext), "w") as f:
            f.write("MODEL\n")
            f.write(
                "\nENDMDL\nMODEL\n".join(
                    ["\n".join([str(rsd) for rsd in models[mdlNo]]) for mdlNo in ntwrk]
                )
            )
            f.write("\nENDMDL")

[process_stubs(fn) for fn in os.listdir(".") if fn.endswith(".pdb")]
