"""Download the public per-star catalogues listed in data/literature_census.csv from VizieR.

Saves every table of each catalogue to ~/data/old_disks/vizier/<catalogue with / -> _>/<table>.ecsv and writes
data/vizier_manifest.csv (catalogue, table, rows, columns, sha256, server). Catalogues the server does not return
are listed with rows = 0.
Server: the CfA mirror; CDS Strasbourg timed out from this machine on 2026-11-06.
"""
import hashlib

import pandas as pd
from astroquery.vizier import Vizier

from olddisks.paths import DATA

SERVER = "vizier.cfa.harvard.edu"
CATALOGUES = [
    "J/MNRAS/456/4315", "J/MNRAS/439/1015", "J/ApJ/779/102", "J/ApJ/834/9", "J/A+A/634/A10", "J/A+A/618/A122",
    "J/A+A/635/A152", "J/MNRAS/466/2006", "J/A+A/675/A49", "J/AJ/137/3100", "J/A+A/459/423", "J/ApJ/923/77",
    "J/AJ/131/895", "J/MNRAS/411/1013", "J/A+A/698/A53", "J/MNRAS/448/2717", "J/MNRAS/495/3022", "J/ApJ/836/202",
    "J/ApJS/191/352", "J/ApJ/921/32", "J/ApJ/920/92", "J/MNRAS/516/2348", "J/ApJS/268/19", "J/ApJ/711/361",
    "J/ApJ/752/45", "J/AJ/159/46",
    "J/ApJ/798/77", "J/AJ/131/332",  # Ho+2015 M31 dE CaT [Fe/H]; Geha+2006 NGC 205 velocities
    # reported missing on the mirror; retried here
    "J/ApJ/750/33", "J/ApJ/903/10", "J/A+A/689/A88", "J/MNRAS/457/1299", "J/A+A/672/A131", "J/ApJ/681/L13",
    "J/ApJ/839/20", "J/ApJ/758/124",
]


def main():
    viz = Vizier(row_limit=-1, columns=["**"], vizier_server=SERVER)
    rows = []
    for cat in CATALOGUES:
        tables = viz.get_catalogs(cat)
        if len(tables) == 0:
            rows.append(dict(catalogue=cat, table="", rows=0, columns="", sha256="", server=SERVER))
            continue
        for name, t in zip(tables.keys(), tables):
            path = DATA / "vizier" / cat.replace("/", "_") / (name.replace("/", "_") + ".ecsv")
            path.parent.mkdir(parents=True, exist_ok=True)
            t.write(path, overwrite=True)
            rows.append(dict(catalogue=cat, table=name, rows=len(t), columns=" ".join(t.colnames),
                             sha256=hashlib.sha256(path.read_bytes()).hexdigest(), server=SERVER))
        print(cat, [(k, len(x)) for k, x in zip(tables.keys(), tables)])
    pd.DataFrame(rows).to_csv("data/vizier_manifest.csv", index=False)


if __name__ == "__main__":
    main()
