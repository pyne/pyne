"""
This Python script parses the JSON data from the "Compendium of Material
Composition Data for Radiation Transport Modeling" provided by the Pacific
Northwest National Laboratory (PNNL). The JSON file contains material properties,
including material density, atom density, elemental and isotopic compositions
provided as weight, atom fraction, and atom density. The data in the JSON file
is intended for radiation transport modeling and can be used with various radiation
transport codes.

Usage:
1. Ensure you have downloaded the JSON file "Materials Compendium" from PNNL.
2. Place the JSON file in the same directory as this script.
3. Run this script using a Python interpreter.

Disclaimer:
- The material composition data in the JSON file is based on Revision 2 of the Compendium
    and is provided with references. Users should be aware of variabilities in composition
    or densities for certain materials, and ranges are shown where possible in the references.
- Users may need to supply application-specific impurities for certain materials if they are
    not found in the known references.
- It is essential to ensure that the chosen radiation transport code uses appropriate
    simulation parameters (e.g., reaction cross sections) to ensure simulation accuracy.

Note:
- This script is specifically designed to parse the JSON file for radiation transport
    modeling but can potentially be adapted for other purposes.


For more information about the Compendium, visit https://compendium.cwmd.pnnl.gov/
"""

import os
from typing import List
from typing import Any
from dataclasses import dataclass
import json


# Data class for representing contact information
@dataclass
class Contact:
    Phone: str
    Name: str
    Email: str

    # Static method to create a Contact object from a dictionary
    @staticmethod
    def from_dict(obj: Any) -> 'Contact':
        _Phone = str(obj.get("Phone", ""))
        _Name = str(obj.get("Name", ""))
        _Email = str(obj.get("Email", ""))
        return Contact(_Phone, _Name, _Email)


# Data class for representing information about an isotope
@dataclass
class Isotope:
    # List all attributes of the Isotope class here
    WeightPercent: float
    Isotope: str
    WeightFraction_whole: float
    IsotopicWeightFraction_whole: float
    WeightFraction: float
    Abundance: float
    IsotopicAtomDensity: float
    AtomicNumber_whole: int
    ZAID: str
    AtomFraction: float
    AtomicNumber: int
    IsotopicWeightFraction: float
    RelativeAtomicMass: float
    RelativeAtomicMass_whole: float
    IsotopicAtomFraction: float
    Abundance_whole: float
    IsotopicAtomFraction_whole: float
    AtomFraction_whole: float
    IsotopicAtomDensity_whole: float

    @staticmethod
    def from_dict(obj: Any) -> 'Isotope':
        _WeightPercent = float(obj.get("WeightPercent"))
        _Isotope = str(obj.get("Isotope"))
        _WeightFraction_whole = float(obj.get("WeightFraction_whole"))
        _IsotopicWeightFraction_whole = float(
            obj.get("IsotopicWeightFraction_whole"))
        _WeightFraction = float(obj.get("WeightFraction"))
        _Abundance = float(obj.get("Abundance"))
        _IsotopicAtomDensity = float(obj.get("IsotopicAtomDensity"))
        _AtomicNumber_whole = int(obj.get("AtomicNumber_whole"))
        _ZAID = str(obj.get("ZAID"))
        _AtomFraction = float(obj.get("AtomFraction"))
        _AtomicNumber = int(obj.get("AtomicNumber"))
        _IsotopicWeightFraction = float(obj.get("IsotopicWeightFraction"))
        _RelativeAtomicMass = float(obj.get("RelativeAtomicMass"))
        _RelativeAtomicMass_whole = float(obj.get("RelativeAtomicMass_whole"))
        _IsotopicAtomFraction = float(obj.get("IsotopicAtomFraction"))
        _Abundance_whole = float(obj.get("Abundance_whole"))
        _IsotopicAtomFraction_whole = float(
            obj.get("IsotopicAtomFraction_whole"))
        _AtomFraction_whole = float(obj.get("AtomFraction_whole"))
        _IsotopicAtomDensity_whole = float(
            obj.get("IsotopicAtomDensity_whole"))
        return Isotope(_WeightPercent, _Isotope, _WeightFraction_whole,
                       _IsotopicWeightFraction_whole, _WeightFraction, _Abundance,
                       _IsotopicAtomDensity, _AtomicNumber_whole, _ZAID, _AtomFraction,
                       _AtomicNumber, _IsotopicWeightFraction, _RelativeAtomicMass,
                       _RelativeAtomicMass_whole, _IsotopicAtomFraction, _Abundance_whole,
                       _IsotopicAtomFraction_whole, _AtomFraction_whole,
                       _IsotopicAtomDensity_whole)


# Data class for representing an element
@dataclass
class Element:
    # List all attributes of the Element class here
    WeightFraction_whole: float
    NonIsotopic: bool
    Element: str
    WeightFraction: float
    AtomicMass: float
    ZAID: str
    AtomFraction: float
    AtomDensity_whole: float
    AtomFraction_whole: float
    id: str
    Isotopes: List[Isotope]
    AtomDensity: float
    AtomicMass_whole: float
    Abundances: str

    @staticmethod
    def from_dict(obj: Any) -> 'Element':
        _WeightFraction_whole = float(obj.get("WeightFraction_whole"))
        _NonIsotopic = bool(obj.get("NonIsotopic"))
        _Element = str(obj.get("Element"))
        _WeightFraction = float(obj.get("WeightFraction"))
        _AtomicMass = float(obj.get("AtomicMass"))
        _ZAID = str(obj.get("ZAID"))
        _AtomFraction = float(obj.get("AtomFraction"))
        _AtomDensity_whole = float(obj.get("AtomDensity_whole"))
        _AtomFraction_whole = float(obj.get("AtomFraction_whole"))
        _id = str(obj.get("id"))
        _Isotopes = [Isotope.from_dict(y) for y in obj.get("Isotopes")]
        _AtomDensity = float(obj.get("AtomDensity"))
        _AtomicMass_whole = float(obj.get("AtomicMass_whole"))
        _Abundances = str(obj.get("Abundances"))
        return Element(_WeightFraction_whole, _NonIsotopic, _Element, _WeightFraction,
                       _AtomicMass, _ZAID, _AtomFraction, _AtomDensity_whole,
                       _AtomFraction_whole, _id, _Isotopes, _AtomDensity,
                       _AtomicMass_whole, _Abundances)


# Data class for representing a molecul
@dataclass
class Mol:
    # List all attributes of the Mol class here
    Mols: int
    Element: str
    Isotope: str

    @staticmethod
    def from_dict(obj: Any) -> 'Mol':
        _Mols = int(obj.get("Mols"))
        _Element = str(obj.get("Element"))
        _Isotope = str(obj.get("Isotope"))
        return Mol(_Mols, _Element, _Isotope)


# Data class for representing data associated with a material
@dataclass
class Datum:
    # List all attributes of the Datum class here
    Comment: List[str]
    Density: float
    Acronym: object
    Elements: List[Element]
    Source: str
    References: List[str]
    Contact: Contact
    MaterialAtomDensity: float
    Mols: List[Mol]
    MatNum: int
    MaterialWeight: str
    Name: str
    Verification_Notes: List[str]
    Formula: str

    @staticmethod
    def from_dict(obj: Any) -> 'Datum':
        _Comment = obj.get("Comment")
        _Density = float(obj.get("Density"))
        _Acronym = obj.get("Acronym")
        _Elements = [Element.from_dict(y) for y in obj.get("Elements")]
        _Source = str(obj.get("Source"))
        _References = obj.get("References")
        _Contact = obj.get("Contact")
        _MaterialAtomDensity = float(obj.get("MaterialAtomDensity"))
        _Mols = [Mol.from_dict(y) for y in obj.get("Mols")]
        _MatNum = int(obj.get("MatNum"))
        _MaterialWeight = str(obj.get("MaterialWeight"))
        _Name = str(obj.get("Name"))
        _Verification_Notes = obj.get("Verification Notes")
        _Formula = str(obj.get("Formula"))
        return Datum(_Comment, _Density, _Acronym, _Elements, _Source, _References,
                     _Contact, _MaterialAtomDensity, _Mols, _MatNum, _MaterialWeight,
                     _Name, _Verification_Notes, _Formula)


# Data class for representing the root of the JSON data
@dataclass
class Root:
    # List all attributes of the Datum class here
    siteVersion: str
    data: List[Datum]

    @staticmethod
    def from_dict(obj: Any) -> 'Root':
        _siteVersion = str(obj.get("siteVersion"))
        _data = [Datum.from_dict(y) for y in obj.get("data")]
        return Root(_siteVersion, _data)


# Get the current directory of the script
current_directory = os.path.dirname(os.path.abspath(__file__))

# Construct the path to the JSON file
json_file_path = os.path.join(current_directory, "MaterialsCompendium.json")

# Read the JSON data from the file
with open(json_file_path, "r") as file:
    jsonstring = file.read()

# Convert the JSON data into a Python dictionary
json_data = json.loads(jsonstring)

# Convert the Python dictionary into structured objects using the Root class
MaterialsCompendium = Root.from_dict(json_data)


"""
Test

for datum in MaterialsCompendium.data:
    if datum.Density == 5.5:
        print(f"Material Name: {datum.Name}")
        print("Elements with density 5.5:")
        for element in datum.Elements:
            print(element.Element)
        print("-------------------------")
"""
