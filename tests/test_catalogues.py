"""Harmonised catalogues must reproduce the member counts published with each table."""
import numpy as np
import pytest

from olddisks.catalogues import load


@pytest.fixture(scope="module")
def stars():
    return load()


@pytest.mark.parametrize("source, galaxy, n_members, n_members_feh", [
    ("Tolstoy+2023", "Sculptor", 1604, 1339),   # A&A 675, A49: 1604 members, 1339 with CaT [Fe/H]
    ("Kacharov+2017", "Phoenix", 194, 194),     # MNRAS 466, 2006: 194 members
    ("Kirby+2014", None, 792, None),            # MNRAS 439, 1015: Mmb = Y for 792 of 862 stars
])
def test_published_member_counts(stars, source, galaxy, n_members, n_members_feh):
    s = stars[(stars.source == source) & stars.member_paper.eq(True)]
    if galaxy:
        s = s[s.galaxy == galaxy]
    assert len(s) == n_members
    if n_members_feh is not None:
        assert np.isfinite(s.feh).sum() == n_members_feh


def test_kirby2020_m31_members(stars):
    s = stars[(stars.source == "Kirby+2020") & stars.member_paper.eq(True)]
    assert len(s) == 256  # AJ 159, 46: 256 members in And I, III, V, VII, X
    assert set(s.galaxy) == {"Andromeda I", "Andromeda III", "Andromeda V", "Andromeda VII", "Andromeda X"}


def test_one_row_per_star_in_walker2023(stars):
    s = stars[stars.source == "Walker+2023"]
    assert not s.duplicated(["instrument", "star_id"]).any()


def test_galaxy_names_are_lvdb_names(stars):
    import pandas as pd
    from olddisks.paths import LVDB
    assert set(stars.galaxy) <= set(pd.read_csv(LVDB).name)


def test_members_centred_on_assigned_galaxy(stars):
    """Catches tables assigned to the wrong galaxy: members' median position within 2 R_half of the LVDB centre."""
    import pandas as pd
    from olddisks.geometry import elliptical_radius
    from olddisks.paths import LVDB
    lvdb = pd.read_csv(LVDB).drop_duplicates("name").set_index("name")
    bad = []
    for (src, gal), s in stars[stars.member_paper.eq(True)].groupby(["source", "galaxy"]):
        g = lvdb.loc[gal]
        r = elliptical_radius(s.ra.median(), s.dec.median(), g.ra, g.dec, 0.0, 0.0) * 60 / g.rhalf
        if r > 2:
            bad.append((src, gal, round(float(r), 1)))
    assert bad == []
