import os
import sys
import numpy as np
from . import db_management as dbm
from .plot_NIST import load_NIST_spectra
from . import get_file
from . import jdx_Reader as jdx

def NIST_Smile_List(expect="all"):
    kwargs = {"db_name": get_file("Molecule_DB.db"),
              "user": "azariven",
              "dir": "",
              "DEBUG": False, "REMOVE": False, "BACKUP": False, "OVERWRITE": False}
    "../../../BiosigSeas2/BioSig_SEAS-master/input/molecule_info"
    cross_db = dbm.database(**kwargs)
    cross_db.access_db()

    cmd = 'SELECT ID.SMILES, ID.inchikey, Spectra.CAS, ID.Formula, ID.IUPAC_chemical_name  \
            FROM ID,Spectra WHERE ID.SMILES=Spectra.Smiles AND Spectra.Is_Gas="Y"'

    result = cross_db.c.execute(cmd)
    data = np.array(result.fetchall()).T

    smiles = data[0]
    inchikeys = data[1]
    CAS = data[2]
    formula = data[3]
    name = data[4]

    if expect == "all":
        return smiles, inchikeys, CAS, formula, name

def test_nist():
    file = "c2h2_test.jdx"
    data = jdx.JdxFile(file)
    info = NIST_Smile_List()
    smiles = info[0]
    print(smiles, len(smiles))

def nist_spectrum(molecule_smile):
    NIST_data = NIST_Smile_List()
    NIST_Smiles = NIST_data[0]
    CAS = NIST_data[2]

    if molecule_smile == "CSC":
        cas = "C75183"
    else:
        cas = CAS[list(NIST_Smiles).index(molecule_smile)]

    nu, absorb = np.load(get_file("NIST_Spectra_Smile_Calibrated/%s.npy" % cas))
    return nu, absorb

def check_NIST_db():
    Smiles = "OC(=O)C(Cl)"
    Absorption_Boost = 3
    nu1, absorb1 = load_NIST_spectra(Smiles, ["wn", "A"], is_smile=True, NIST_Spectra="NIST_Spectra")

    abfix = []
    for i in absorb1:
        if i > 2:
            abfix.append(2)
        else:
            abfix.append(i)
    absorb1 = np.array(abfix)
    absorb_baseline = 0.001  # baseline_als(absorb1, 10**6, 0.001, niter=10)
    absorb2 = absorb1 - absorb_baseline
    absorb2_calibrated = absorb2 / max(absorb2) * Absorption_Boost

    return nu1, absorb1, absorb2_calibrated

if __name__ == "__main__":
    nu1, absorb1, absorb2_calibrated = check_NIST_db()
    print("Baseline corrected and calibrated spectra available.")


