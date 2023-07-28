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
    """
    Represents a contact with phone, name, and email attributes.
    """
    Phone: str
    Name: str
    Email: str

    # Static method to create a Contact object from a dictionary
    @staticmethod
    def from_dict(obj: Any) -> 'Contact':
        """
        Create a Contact object from a dictionary representation.

        Parameters:
        obj (dict): A dictionary containing keys "Phone", "Name", and "Email".

        Returns:
        Contact: A new Contact object with attributes initialized using the values from the dictionary.
        """
        _Phone = str(obj.get("Phone", ""))
        _Name = str(obj.get("Name", ""))
        _Email = str(obj.get("Email", ""))
        return Contact(_Phone, _Name, _Email)


# Data class for representing information about an isotope
@dataclass
class Isotope:
    """
    Represents an isotope with various attributes related to its properties.

    Attributes:
        WeightPercent (float): Weight percent of the isotope.
        Isotope (str): The name or identifier of the isotope.
        WeightFraction_whole (float): Weight fraction of the isotope as a whole number.
        IsotopicWeightFraction_whole (float): Isotopic weight fraction of the
        isotope as a whole number.
        WeightFraction (float): Weight fraction of the isotope.
        Abundance (float): Abundance of the isotope.
        IsotopicAtomDensity (float): Isotopic atom density of the isotope.
        AtomicNumber_whole (int): Atomic number of the isotope as a whole number.
        ZAID (str): ZAID (unique identifier) of the isotope.
        AtomFraction (float): Atom fraction of the isotope.
        AtomicNumber (int): Atomic number of the isotope.
        IsotopicWeightFraction (float): Isotopic weight fraction of the isotope.
        RelativeAtomicMass (float): Relative atomic mass of the isotope.
        RelativeAtomicMass_whole (float): Relative atomic mass of the isotope as a whole number.
        IsotopicAtomFraction (float): Isotopic atom fraction of the isotope.
        Abundance_whole (float): Abundance of the isotope as a whole number.
        IsotopicAtomFraction_whole (float): Isotopic atom fraction of the isotope as a whole number.
        AtomFraction_whole (float): Atom fraction of the isotope as a whole number.
        IsotopicAtomDensity_whole (float): Isotopic atom density of the isotope as a whole number.

    Methods:
        from_dict(obj: Any) -> 'Isotope':
            Create an Isotope object from a dictionary representation.
    """
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
        """
        Create an Isotope object from a dictionary representation.

        Parameters:
            obj (dict): A dictionary containing attributes of the Isotope.

        Returns:
            Isotope: A new Isotope object with attributes initialized using the values 
            from the dictionary.
        """
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
    """
    Represents an element with various attributes related to its properties.

    Attributes:
        WeightFraction_whole (float): Weight fraction of the element as a whole number.
        NonIsotopic (bool): Whether the element is non-isotopic.
        Element (str): The name or identifier of the element.
        WeightFraction (float): Weight fraction of the element.
        AtomicMass (float): Atomic mass of the element.
        ZAID (str): ZAID (unique identifier) of the element.
        AtomFraction (float): Atom fraction of the element.
        AtomDensity_whole (float): Atom density of the element as a whole number.
        AtomFraction_whole (float): Atom fraction of the element as a whole number.
        id (str): Identifier of the element.
        Isotopes (List[Isotope]): List of Isotope objects representing the isotopes of the element.
        AtomDensity (float): Atom density of the element.
        AtomicMass_whole (float): Atomic mass of the element as a whole number.
        Abundances (str): Abundances of the element.

    Methods:
        from_dict(obj: Any) -> 'Element':
            Create an Element object from a dictionary representation.
    """
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
        """
        Create an Element object from a dictionary representation.

        Parameters:
            obj (dict): A dictionary containing attributes of the Element.

        Returns:
            Element: A new Element object with attributes initialized using the values from the dictionary.
        """
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
    """
    Represents a molecular entity with various attributes related to its properties.

    Attributes:
        Mols (int): The number of moles of the molecular entity.
        Element (str): The name or identifier of the element in the molecular entity.
        Isotope (str): The name or identifier of the isotope in the molecular entity.

    Methods:
        from_dict(obj: Any) -> 'Mol':
            Create a Mol object from a dictionary representation.
    """
    Mols: int
    Element: str
    Isotope: str

    @staticmethod
    def from_dict(obj: Any) -> 'Mol':
        """
        Create a Mol object from a dictionary representation.

        Parameters:
            obj (dict): A dictionary containing attributes of the Mol.

        Returns:
            Mol: A new Mol object with attributes initialized using the values from the dictionary.
        """
        _Mols = int(obj.get("Mols"))
        _Element = str(obj.get("Element"))
        _Isotope = str(obj.get("Isotope"))
        return Mol(_Mols, _Element, _Isotope)


# Data class for representing data associated with a material
@dataclass
class Datum:
    """
    Represents a data entry with various attributes related to its properties.

    Attributes:
        Comment (List[str]): List of comments associated with the data entry.
        Density (float): The density value of the data entry.
        Acronym (object): An acronym or abbreviation associated with the data entry.
        Elements (List[Element]): List of Element objects representing the elements in the data entry.
        Source (str): The source of the data entry.
        References (List[str]): List of references associated with the data entry.
        Contact (Contact): A Contact object representing the contact information related to the data entry.
        MaterialAtomDensity (float): The material atom density value of the data entry.
        Mols (List[Mol]): List of Mol objects representing the molecular entities in the data entry.
        MatNum (int): The material number associated with the data entry.
        MaterialWeight (str): The material weight information of the data entry.
        Name (str): The name of the data entry.
        Verification_Notes (List[str]): List of verification notes related to the data entry.
        Formula (str): The formula associated with the data entry.

    Methods:
        from_dict(obj: Any) -> 'Datum':
            Create a Datum object from a dictionary representation.
    """
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
        """
        Create a Datum object from a dictionary representation.

        Parameters:
            obj (dict): A dictionary containing attributes of the Datum.

        Returns:
            Datum: A new Datum object with attributes initialized using the values from the dictionary.
        """
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
    """
    Represents a root object that contains site version and a list of data entries.

    Attributes:
        siteVersion (str): The site version associated with the root object.
        data (List[Datum]): List of Datum objects representing the data entries.

    Methods:
        from_dict(obj: Any) -> 'Root':
            Create a Root object from a dictionary representation.
    """
    siteVersion: str
    data: List[Datum]

    @staticmethod
    def from_dict(obj: Any) -> 'Root':
        """
        Create a Root object from a dictionary representation.

        Parameters:
            obj (dict): A dictionary containing attributes of the Root.

        Returns:
            Root: A new Root object with attributes initialized using the values from the dictionary.
        """
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
MaterialsCompendium = Root.from_dict(json_data).data
