import pytest
from pyne.pnnl.materials_compendium import Root, Datum, Element, Isotope, Contact, Mol
from pyne.pnnl.utils import (
    ContactInfo,
    MolsInfo,
    IsotopeInfo,
    ElementInfo,
    Material,
)


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
    """
    Test the 'from_dict' method of the 'Root' class to ensure it correctly converts
    a sample data dictionary into a 'Root' object and checks basic attributes.
    """
    root = Root.from_dict(sample_data)
    assert root.siteVersion == "0.1.1"
    assert len(root.data) == 1


def test_datum_from_dict(sample_data):
    """
    Test the 'from_dict' method of the 'Datum' class to ensure it correctly converts
    a sample data dictionary into a 'Datum' object and checks basic attributes.
    """
    datum_dict = sample_data["data"][0]
    datum = Datum.from_dict(datum_dict)
    assert datum.Name == "Boron Carbide"
    assert datum.Density == 2.52
    assert len(datum.Elements) == 1


def test_element_from_dict(sample_data):
    """
    Test the 'from_dict' method of the 'Element' class to ensure it correctly converts
    a sample data dictionary into an 'Element' object and checks basic attributes.
    """
    element_dict = sample_data["data"][0]["Elements"][0]
    element = Element.from_dict(element_dict)
    assert element.Element == "B"
    assert element.NonIsotopic is True
    assert len(element.Isotopes) == 2


def test_isotope_from_dict(sample_data):
    """
    Test the 'from_dict' method of the 'Isotope' class to ensure it correctly converts
    a sample data dictionary into an 'Isotope' object and checks basic attributes.
    """
    isotope_dict1 = sample_data["data"][0]["Elements"][0]["Isotopes"][0]
    isotope_dict2 = sample_data["data"][0]["Elements"][0]["Isotopes"][1]
    isotope1 = Isotope.from_dict(isotope_dict1)
    isotope2 = Isotope.from_dict(isotope_dict2)
    assert isotope1.Isotope == "B10"
    assert isotope1.WeightPercent == 0.184267
    assert isotope2.Isotope == "B11"
    assert isotope2.WeightPercent == 0.815504


def test_contact_from_dict(sample_data):
    """
    Test the 'from_dict' method of the 'Contact' class to ensure it correctly converts
    a sample data dictionary into a 'Contact' object and checks basic attributes.
    """
    contact_dict = sample_data["data"][0]["Contact"]
    contact = Contact.from_dict(contact_dict)
    assert contact.Name == "Rebecca Detwiler"
    assert contact.Phone == "352-484-4040"


def test_mol_from_dict(sample_data):
    """
    Test the 'from_dict' method of the 'Mol' class to ensure it correctly converts
    a sample data dictionary into a 'Mol' object and checks basic attributes.
    """
    mol_dict1 = sample_data["data"][0]["Mols"][0]
    mol_dict2 = sample_data["data"][0]["Mols"][1]
    mol1 = Mol.from_dict(mol_dict1)
    mol2 = Mol.from_dict(mol_dict2)
    assert mol1.Element == "B"
    assert mol1.Mols == 4
    assert mol1.Isotope == "B4"
    assert mol2.Element == "C"
    assert mol2.Mols == 1
    assert mol2.Isotope == "C"


def test_dataclass_equality():
    """
    Test the equality of two dataclass objects of the 'Element' class.
    """
    element1 = Element(
        0.782671,
        True,
        "B",
        0.782671,
        10.8135,
        "5000",
        0.8,
        0.109841,
        0.8,
        "B",
        [],
        0.109841,
        10.8135,
        "0.0",
    )
    element2 = Element(
        0.782671,
        True,
        "B",
        0.782671,
        10.8135,
        "5000",
        0.8,
        0.109841,
        0.8,
        "B",
        [],
        0.109841,
        10.8135,
        "0.0",
    )
    assert element1 == element2


def test_dataclass_not_equality():
    """
    Test the inequality of two dataclass objects of the 'Element' class.
    """
    element1 = Element(
        0.782671,
        True,
        "B",
        0.782671,
        10.8135,
        "5000",
        0.8,
        0.109841,
        0.8,
        "B",
        [],
        0.109841,
        10.8135,
        "0.0",
    )
    element2 = Element(
        0.782671,
        False,
        "B",
        0.782671,
        10.8135,
        "5000",
        0.8,
        0.109841,
        0.8,
        "B",
        [],
        0.109841,
        10.8135,
        "0.0",
    )
    assert element1 != element2


def test_contact_info():
    """
    Test the ContactInfo class to ensure it returns the correct contact information.
    """
    contact_data = Contact(
        Name="Ahnaf Tahmid Chowdhury", Phone="123-456-7890", Email="tahmid@example.com"
    )
    contact_info = ContactInfo(contact_data)

    assert contact_info.get_name() == "Ahnaf Tahmid Chowdhury"
    assert contact_info.get_phone() == "123-456-7890"
    assert contact_info.get_email() == "tahmid@example.com"


def test_mols_info():
    """
    Test the MolsInfo class to ensure it returns the correct molecular information.
    """
    mol_data = Mol(Mols=42, Isotope="C-14", Element="Carbon")
    mols_info = MolsInfo(mol_data)

    assert mols_info.get_mols() == 42
    assert mols_info.get_isotope() == "C-14"
    assert mols_info.get_element() == "Carbon"


def test_isotope_info():
    """
    Test the IsotopeInfo class to ensure it returns the correct isotope information.
    """
    isotope_data = Isotope(
        WeightPercent=0.011,
        Isotope="C-12",
        WeightFraction_whole=12.0,
        IsotopicWeightFraction_whole=13.0,
        WeightFraction=0.012,
        Abundance=98.9,
        IsotopicAtomDensity=0.013,
        AtomicNumber_whole=6,
        ZAID="12000",
        AtomFraction=0.014,
        AtomicNumber=6,
        IsotopicWeightFraction=0.015,
        RelativeAtomicMass=12.01,
        RelativeAtomicMass_whole=12.0,
        IsotopicAtomFraction=0.016,
        Abundance_whole=99.0,
        IsotopicAtomFraction_whole=16.0,
        AtomFraction_whole=14.0,
        IsotopicAtomDensity_whole=13.0,
    )
    isotope_info = IsotopeInfo(isotope_data)

    assert isotope_info.get_all() == (
        "Isotope: C-12 \n"
        " ZAID: 12000 \n"
        " Weight Percent: 0.011 \n"
        " Weight Fraction (Whole): 12.0 \n"
        " Isotopic Weight Fraction (Whole): 13.0 \n"
        " Weight Fraction: 0.012 \n"
        " Abundance: 98.9 \n"
        " Isotopic Atom Density: 0.013 \n"
        " Atomic Number (Whole): 6 \n"
        " Atom Fraction: 0.014 \n"
        " Atomic Number: 6 \n"
        " Isotopic Weight Fraction: 0.015 \n"
        " Relative Atomic Mass: 12.01 \n"
        " Relative Atomic Mass (Whole): 12.0 \n"
        " Isotopic Atom Fraction: 0.016 \n"
        " Abundance (Whole): 99.0 \n"
        " Isotopic Atom Fraction (Whole): 16.0 \n"
        " Atom Fraction (Whole): 14.0 \n"
        " Isotopic Atom Density (Whole): 13.0 \n"
    )


def test_element_info():
    """
    Test the ElementInfo class to ensure it returns the correct element information.
    """
    isotope_data = Isotope(
        WeightPercent=0.011,
        Isotope="C-12",
        WeightFraction_whole=12.0,
        IsotopicWeightFraction_whole=13.0,
        WeightFraction=0.012,
        Abundance=98.9,
        IsotopicAtomDensity=0.013,
        AtomicNumber_whole=6,
        ZAID="12000",
        AtomFraction=0.014,
        AtomicNumber=6,
        IsotopicWeightFraction=0.015,
        RelativeAtomicMass=12.01,
        RelativeAtomicMass_whole=12.0,
        IsotopicAtomFraction=0.016,
        Abundance_whole=99.0,
        IsotopicAtomFraction_whole=16.0,
        AtomFraction_whole=14.0,
        IsotopicAtomDensity_whole=13.0,
    )
    element_data = Element(
        WeightFraction_whole=12.0,
        NonIsotopic=False,
        Element="C",
        WeightFraction=0.012,
        AtomicMass=12.01,
        ZAID="12000",
        AtomFraction=0.014,
        AtomDensity_whole=1.0,
        AtomFraction_whole=14.0,
        id="42",
        Isotopes=[isotope_data],
        AtomDensity=1.5,
        AtomicMass_whole=12,
        Abundances="98.9",
    )
    element_info = ElementInfo(element_data)

    assert element_info.get_all() == (
        "Element: C \n"
        " Id: 42 \n"
        " ZAID: 12000 \n"
        " Atomic Mass: 12.01 \n"
        " Atom Density: 1.5 \n"
        " Atomic Mass (Whole): 12 \n"
        " Atom Fraction: 0.014 \n"
        " Weight Fraction: 0.012 \n"
        " Atom Fraction (Whole): 14.0 \n"
        " Weight Fraction (Whole): 12.0 \n"
        " Non Isotopic: False \n"
        " Isotopes: C-12 \n"
        " Abundances: 98.9"
    )


def test_material():
    """
    Test the Material class to ensure it returns the correct material information.
    """
    datum = Datum(
        Name="Sample Material",
        Formula="H2O",
        Comment=["Comment 1", "Comment 2"],
        Density=1.0,
        Acronym=["ABC", "DEF"],
        Elements=[
            Element(
                WeightFraction_whole=12.0,
                NonIsotopic=False,
                Element="C",
                WeightFraction=0.012,
                AtomicMass=12.01,
                ZAID="12000",
                AtomFraction=0.014,
                AtomDensity_whole=1.0,
                AtomFraction_whole=14.0,
                id="42",
                Isotopes=[
                    Isotope(
                        WeightPercent=0.011,
                        Isotope="C-12",
                        WeightFraction_whole=12.0,
                        IsotopicWeightFraction_whole=13.0,
                        WeightFraction=0.012,
                        Abundance=98.9,
                        IsotopicAtomDensity=0.013,
                        AtomicNumber_whole=6,
                        ZAID="12000",
                        AtomFraction=0.014,
                        AtomicNumber=6,
                        IsotopicWeightFraction=0.015,
                        RelativeAtomicMass=12.01,
                        RelativeAtomicMass_whole=12.0,
                        IsotopicAtomFraction=0.016,
                        Abundance_whole=99.0,
                        IsotopicAtomFraction_whole=16.0,
                        AtomFraction_whole=14.0,
                        IsotopicAtomDensity_whole=13.0,
                    )
                ],
                AtomDensity=1.5,
                AtomicMass_whole=12,
                Abundances="98.9",
            )
        ],
        Source="Some source information",
        References=["Ref 1", "Ref 2"],
        Contact={
            "Name": "Ahnaf Tahmid Chowdhury",
            "Phone": "123-456-7890",
            "Email": "tahmid@example.com",
        },
        MaterialAtomDensity=2.5,
        Mols=[Mol(Mols=42, Isotope="C12", Element="Carbon")],
        MatNum=12,
        MaterialWeight="10",
        Verification_Notes=["Some verification notes"],
    )
    material = Material(datum)

    assert material.get_name() == "Sample Material"
    assert material.get_formula() == "H2O"
