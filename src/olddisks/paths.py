"""Locations of external data (see data/sources.yaml)."""
from pathlib import Path

DATA = Path.home() / "data/old_disks"
LVDB = DATA / "lvdb/comb_all.csv"  # Local Volume Database v1.1.2
SURVEY_MEMBERS = DATA / "survey_members.parquet"
