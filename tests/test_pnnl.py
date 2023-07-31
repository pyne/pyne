import pytest
from pyne.pnnl.materials_compendium import Root, Datum, Element, Isotope, Contact, Mol


@pytest.fixture
def sample_data():
    """
    Sample data for testing.
    """
    data = {
        "siteVersion": "0.1.1",
        "data": [
            {
                "Comment": [
                    "The above density is estimated to be accurate to 3 significant digits. Uncertainties are not addressed.",
                    "MatWeb also lists reactor grade boron carbide with 2.65 density, density varies between 2.5 to 2.65g/cc.",
                    "The MCNP primer lists density as 2.51.",
                ],
                "Density": 2.52,
                "Elements": [
                    {
                        "WeightFraction_whole": 0.7826710045852137,
                        "NonIsotopic": True,
                        "Element": "B",
                        "WeightFraction": 0.782671,
                        "AtomicMass": 10.8135,
                        "ZAID": "5000",
                        "AtomFraction": 0.8,
                        "AtomDensity_whole": 0.10984098113584465,
                        "AtomFraction_whole": 0.8,
                        "id": "B",
                        "Isotopes": [
                            {
                                "WeightPercent": 0.184267,
                                "Isotope": "B10",
                                "WeightFraction_whole": 0.14422067312891074,
                                "IsotopicWeightFraction_whole": 0.18426730041614653,
                                "WeightFraction": 0.144221,
                                "Abundance": 0.199,
                                "IsotopicAtomDensity": 0.021858,
                                "AtomicNumber_whole": 5,
                                "ZAID": "5010",
                                "AtomFraction": 0.1592,
                                "AtomicNumber": 5,
                                "IsotopicWeightFraction": 0.184267,
                                "RelativeAtomicMass": 10.012937,
                                "RelativeAtomicMass_whole": 10.01293695,
                                "IsotopicAtomFraction": 0.199,
                                "Abundance_whole": 0.199,
                                "IsotopicAtomFraction_whole": 0.199,
                                "AtomFraction_whole": 0.1592,
                                "IsotopicAtomDensity_whole": 0.021858355246033086,
                            },
                            {
                                "WeightPercent": 0.815504,
                                "Isotope": "B11",
                                "WeightFraction_whole": 0.638271413770117,
                                "IsotopicWeightFraction_whole": 0.8155041007407409,
                                "WeightFraction": 0.638271,
                                "Abundance": 0.801,
                                "IsotopicAtomDensity": 0.087983,
                                "AtomicNumber_whole": 5,
                                "ZAID": "5011",
                                "AtomFraction": 0.6408,
                                "AtomicNumber": 5,
                                "IsotopicWeightFraction": 0.815504,
                                "RelativeAtomicMass": 11.009305,
                                "RelativeAtomicMass_whole": 11.00930536,
                                "IsotopicAtomFraction": 0.801,
                                "Abundance_whole": 0.801,
                                "IsotopicAtomFraction_whole": 0.801,
                                "AtomFraction_whole": 0.6408,
                                "IsotopicAtomDensity_whole": 0.08798262588981155,
                            },
                        ],
                        "AtomDensity": 0.109841,
                        "AtomicMass_whole": 10.8135,
                    },
                ],
                "Source": "PNNL",
                "References": [
                    "Density and weight fractions from http://physics.nist.gov/cgi-bin/Star/compos.pl?matno=121.",
                    "Formula also from http://www.matweb.com/search/datasheet.aspx?MatGUID=45fd34d496fe48e3ab513bcbc4079430",
                    "Also listed at  LA-UR-09-0380, Criticality Calculations with MCNP5: A Primer by R. Brewer, LANL, Jan 2009.",
                ],
                "Contact": {
                    "Phone": "352-484-4040",
                    "Name": "Rebecca Detwiler",
                    "Email": "Rebecca.Detwiler@pnnl.gov",
                },
                "MaterialAtomDensity": 0.137301,
                "Mols": [
                    {"Mols": 4, "Element": "B", "Isotope": "B4"},
                    {"Mols": 1, "Element": "C", "Isotope": "C"},
                ],
                "Formula": "B4C",
                "MatNum": 42,
                "MaterialWeight": "55.264600",
                "Name": "Boron Carbide",
            },
        ],
    }
    return data


def test_root_from_dict(sample_data):
    root = Root.from_dict(sample_data)
    assert root.siteVersion == "0.1.1"
    assert len(root.data) == 1


def test_datum_from_dict(sample_data):
    datum_dict = sample_data["data"][0]
    datum = Datum.from_dict(datum_dict)
    assert datum.Name == "Boron Carbide"
    assert datum.Density == 2.52
    assert len(datum.Elements) == 1


def test_element_from_dict(sample_data):
    element_dict = sample_data["data"][0]["Elements"][0]
    element = Element.from_dict(element_dict)
    assert element.Element == "B"
    assert element.NonIsotopic is True
    assert len(element.Isotopes) == 2
