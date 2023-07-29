"""
Separate classes and functions to perform specific tasks
to retrieve specific information from the Materials Compendium
"""
import difflib
from .materials_compendium import (
    MaterialsCompendium,
    Datum,
    Isotope,
    Element,
    Contact,
    Mol,
)


class ContactInfo:
    def __init__(self, contact_data: Contact):
        """
        Initialize a ContactInfo object with contact data.

        Parameters:
            contact_data (Contact): An instance of the Contact class containing name, phone, and email information.

        Example:
            contact_data = Contact(name="Ahnaf Tahmid Chowdhury", phone="123-456-7890", email="tahmid@example.com")
            contact_info = ContactInfo(contact_data)
        """
        self.name = contact_data.Name
        self.phone = contact_data.Phone
        self.email = contact_data.Email

    def __str__(self):
        """
        Return a string representation of the ContactInfo object.

        Returns:
            str: A formatted string containing contact information.

        Example:
            contact_data = Contact(name="Ahnaf Tahmid Chowdhury", phone="123-456-7890", email="tahmid@example.com")
            contact_info = ContactInfo(contact_data)
            print(contact_info)
            # Output: "Name: Ahnaf Tahmid Chowdhury, Phone: 123-456-7890, Email: tahmid@example.com"
        """
        return f"Name: {self.name}, Phone: {self.phone}, Email: {self.email}"

    def get_name(self):
        """
        Get the name associated with the contact.

        Returns:
            str: The name of the contact.

        Example:
            print(contact_info.get_name())
            # Output: "Ahnaf Tahmid Chowdhury"
        """
        return self.name

    def get_phone(self):
        """
        Get the phone number associated with the contact.

        Returns:
            str: The phone number of the contact.

        Example:
            print(contact_info.get_phone())
            # Output: "123-456-7890"
        """
        return self.phone

    def get_email(self):
        """
        Get the email address associated with the contact.

        Returns:
            str: The email address of the contact.

        Example:
            print(contact_info.get_email())
            # Output: "tahmid@example.com"
        """
        return self.email


class MolsInfo:
    def __init__(self, mol_data: Mol):
        """
        Initialize a MolsInfo object with molecular data.

        Parameters:
            mol_data (Mol): An instance of the Mol class containing mols, isotope, and element information.

        Example:
            mol_data = Mol(mols=42, isotope="C-14", element="Carbon")
            mols_info = MolsInfo(mol_data)
        """
        self.mols = mol_data.Mols
        self.isotope = mol_data.Isotope
        self.element = mol_data.Element

    def __str__(self):
        """
        Return a string representation of the MolsInfo object.

        Returns:
            str: A formatted string containing molecular information.

        Example:
            print(mols_info)
            # Output: "Mols: 42, Isotopes: C-14, Elements: Carbon"
        """
        return f"Mols: {self.mols}, Isotopes: {self.isotope}, Elements: {self.element}"

    def get_mols(self):
        """
        Get the mols associated with the molecular data.

        Returns:
            int: The quantity of mols.

        Example:
            print(mols_info.get_mols())
            # Output: 42
        """
        return self.mols

    def get_isotope(self):
        """
        Get the isotope associated with the molecular data.

        Returns:
            str: The isotope information.

        Example:
            print(mols_info.get_isotope())
            # Output: "C-14"
        """
        return self.isotope

    def get_element(self):
        """
        Get the element associated with the molecular data.

        Returns:
            str: The element information.

        Example:
            print(mols_info.get_element())
            # Output: "Carbon"
        """
        return self.element


class IsotopeInfo:
    def __init__(self, isotope_data: Isotope):
        """
        Initialize an IsotopeInfo object with isotope data.

        Parameters:
            isotope_data (Isotope): An instance of the Isotope class containing various isotope information.

        Example:
            isotope_data = Isotope(
                WeightPercent=0.011, Isotope="C-12", WeightFraction_whole=12.0, IsotopicWeightFraction_whole=13.0,
                WeightFraction=0.012, Abundance=98.9, IsotopicAtomDensity=0.013, AtomicNumber_whole=6,
                ZAID="12000", AtomFraction=0.014, AtomicNumber=6, IsotopicWeightFraction=0.015,
                RelativeAtomicMass=12.01, RelativeAtomicMass_whole=12.0, IsotopicAtomFraction=0.016,
                Abundance_whole=99.0, IsotopicAtomFraction_whole=16.0, AtomFraction_whole=14.0,
                IsotopicAtomDensity_whole=13.0
            )
            isotope_info = IsotopeInfo(isotope_data)
        """
        self.weight_percent = isotope_data.WeightPercent
        self.isotope = isotope_data.Isotope
        self.weight_fraction_whole = isotope_data.WeightFraction_whole
        self.isotopic_weight_fraction_whole = isotope_data.IsotopicWeightFraction_whole
        self.weight_fraction = isotope_data.WeightFraction
        self.abundance = isotope_data.Abundance
        self.isotopic_atom_density = isotope_data.IsotopicAtomDensity
        self.atomic_number_whole = isotope_data.AtomicNumber_whole
        self.zaid = isotope_data.ZAID
        self.atom_fraction = isotope_data.AtomFraction
        self.atomic_number = isotope_data.AtomicNumber
        self.isotopic_weight_fraction = isotope_data.IsotopicWeightFraction
        self.relative_atomic_mass = isotope_data.RelativeAtomicMass
        self.relative_atomic_mass_whole = isotope_data.RelativeAtomicMass_whole
        self.isotopic_atom_fraction = isotope_data.IsotopicAtomFraction
        self.abundance_whole = isotope_data.Abundance_whole
        self.isotopic_atom_fraction_whole = isotope_data.IsotopicAtomFraction_whole
        self.atom_fraction_whole = isotope_data.AtomFraction_whole
        self.isotopic_atom_density_whole = isotope_data.IsotopicAtomDensity_whole

    def __str__(self):
        """
        Return a string representation of the IsotopeInfo object.

        Returns:
            str: A formatted string containing isotope information.
        """
        return f"Isotope: {self.isotope}, Weight Percent: {self.weight_percent}, ZAID: {self.zaid}"

    def get_all(self):
        """
        Get all available isotope information in a formatted string.

        Returns:
            str: A formatted string containing all available isotope information.
        """
        return (
            f"Isotope: {self.isotope} \n"
            + f" ZAID: {self.zaid} \n"
            + f" Weight Percent: {self.weight_percent} \n"
            + f" Weight Fraction (Whole): {self.weight_fraction_whole} \n"
            + f" Isotopic Weight Fraction (Whole): {self.isotopic_weight_fraction_whole} \n"
            + f" Weight Fraction: {self.weight_fraction} \n"
            + f" Abundance: {self.abundance} \n"
            + f" Isotopic Atom Density: {self.isotopic_atom_density} \n"
            + f" Atomic Number (Whole): {self.atomic_number_whole} \n"
            + f" Atom Fraction: {self.atom_fraction} \n"
            + f" Atomic Number: {self.atomic_number} \n"
            + f" Isotopic Weight Fraction: {self.isotopic_weight_fraction} \n"
            + f" Relative Atomic Mass: {self.relative_atomic_mass} \n"
            + f" Relative Atomic Mass (Whole): {self.relative_atomic_mass_whole} \n"
            + f" Isotopic Atom Fraction: {self.isotopic_atom_fraction} \n"
            + f" Abundance (Whole): {self.abundance_whole} \n"
            + f" Isotopic Atom Fraction (Whole): {self.isotopic_atom_fraction_whole} \n"
            + f" Atom Fraction (Whole): {self.atom_fraction_whole} \n"
            + f" Isotopic Atom Density (Whole): {self.isotopic_atom_density_whole} \n"
        )


class ElementInfo:
    def __init__(self, element_data: Element):
        """
        Initialize an ElementInfo object with element data.

        Parameters:
            element_data (Element): An instance of the Element class containing various element information.

        Example:
            element_data = Element(
                WeightFraction_whole=12.0, NonIsotopic=False, Element="C",
                WeightFraction=0.012, AtomicMass=12.01, ZAID="12000", AtomFraction=0.014,
                AtomDensity_whole=1.0, AtomFraction_whole=14.0, id="42",
                Isotopes=[
                    Isotope(
                WeightPercent=0.011, Isotope="C-12", WeightFraction_whole=12.0, IsotopicWeightFraction_whole=13.0,
                WeightFraction=0.012, Abundance=98.9, IsotopicAtomDensity=0.013, AtomicNumber_whole=6,
                ZAID="12000", AtomFraction=0.014, AtomicNumber=6, IsotopicWeightFraction=0.015,
                RelativeAtomicMass=12.01, RelativeAtomicMass_whole=12.0, IsotopicAtomFraction=0.016,
                Abundance_whole=99.0, IsotopicAtomFraction_whole=16.0, AtomFraction_whole=14.0,
                IsotopicAtomDensity_whole=13.0
                ),],
                AtomDensity=1.5, AtomicMass_whole=12, Abundances="98.9"
            )
            element_info = ElementInfo(element_data)
        """
        self.weight_fraction_whole = element_data.WeightFraction_whole
        self.non_isotopic = element_data.NonIsotopic
        self.element = element_data.Element
        self.weight_fraction = element_data.WeightFraction
        self.atomic_mass = element_data.AtomicMass
        self.zaid = element_data.ZAID
        self.atom_fraction = element_data.AtomFraction
        self.atom_density_whole = element_data.AtomDensity_whole
        self.atom_fraction_whole = element_data.AtomFraction_whole
        self.id = element_data.id
        self.isotopes = [
            IsotopeInfo(isotope_data) for isotope_data in element_data.Isotopes
        ]
        self.atom_density = element_data.AtomDensity
        self.atomic_mass_whole = element_data.AtomicMass_whole
        self.abundances = element_data.Abundances

    def __str__(self):
        """
        Return a string representation of the ElementInfo object.

        Returns:
            str: A formatted string containing element information.

        Example:
            print(element_info)
            # Output: "Element: Carbon \n ZAID: 12000 \n Isotopes: C-12, C-13, ..."
        """
        return f"Element: {self.element} \n ZAID: {self.zaid} \n Isotopes:{', '.join([isotope.isotope for isotope in self.isotopes])}"

    def get_all(self):
        """
        Get all available element information in a formatted string.

        Returns:
            str: A formatted string containing all available element information.

        Example:
            print(element_info.get_all())
            # Output: "Element: Carbon
            #  Id: 42
            #  ZAID: 12000
            #  Atomic Mass: 12.01
            #  Atom Density: 1.5
            #  Atomic Mass (Whole): 12
            #  Atom Fraction: 0.014
            #  Weight Fraction: 0.012
            #  Atom Fraction (Whole): 14
            #  Weight Fraction (Whole): 12
            #  Non Isotopic: False
            #  Isotopes: C-12, C-13, ...
            #  Abundances: 98.9"
        """
        return (
            f"Element: {self.element} \n"
            + f" Id: {self.id} \n"
            + f" ZAID: {self.zaid} \n"
            + f" Atomic Mass: {self.atomic_mass} \n"
            + f" Atom Density: {self.atom_density} \n"
            + f" Atomic Mass (Whole): {self.atomic_mass_whole} \n"
            + f" Atom Fraction: {self.atom_fraction} \n"
            + f" Weight Fraction: {self.weight_fraction} \n"
            + f" Atom Fraction (Whole): {self.atom_fraction_whole} \n"
            + f" Weight Fraction (Whole): {self.weight_fraction_whole} \n"
            + f" Non Isotopic: {self.non_isotopic} \n"
            + f" Isotopes: {', '.join([isotope.isotope for isotope in self.isotopes])} \n"
            + f" Abundances:{self.abundances}"
        )

    # Methods for accessing specific attributes

    def get_weight_fraction(self):
        """
        Get the weight fraction associated with the element.

        Returns:
            float: The weight fraction.

        Example:
            print(element_info.get_weight_fraction())
            # Output: 0.012
        """
        return self.weight_fraction

    def get_weight_fraction_whole(self):
        """
        Get the whole number weight fraction associated with the element.

        Returns:
            float: The whole number weight fraction.

        Example:
            print(element_info.get_weight_fraction_whole())
            # Output: 12.0
        """
        return self.weight_fraction_whole

    def get_non_isotopic(self):
        """
        Get the non-isotopic information associated with the element.

        Returns:
            bool: The non-isotopic information.

        Example:
            print(element_info.get_non_isotopic())
            # Output: False
        """
        return self.non_isotopic

    def get_element(self):
        """
        Get the name of the element.

        Returns:
            str: The name of the element.

        Example:
            # Output: "C"
        """
        return self.element

    def get_atomic_mass(self):
        """
        Get the atomic mass associated with the element.

        Returns:
            float: The atomic mass.

        Example:
            # Output: 12.01
        """
        return self.atomic_mass

    def get_zaid(self):
        """
        Get the ZAID (unique identifier) associated with the element.

        Returns:
            str: The ZAID.

        Example:
            # Output: "12000"
        """
        return self.zaid

    def get_atom_fraction(self):
        """
        Get the atom fraction associated with the element.

        Returns:
            float: The atom fraction.

        Example:
            print(element_info.get_atom_fraction())
            # Output: 0.014
        """
        return self.atom_fraction

    def get_atom_density_whole(self):
        """
        Get the whole number atom density associated with the element.

        Returns:
            float: The whole number atom density.

        Example:
            print(element_info.get_atom_density_whole())
            # Output: 1.0
        """
        return self.atom_density_whole

    def get_atom_fraction_whole(self):
        """
        Get the whole number atom fraction associated with the element.

        Returns:
            float: The whole number atom fraction.

        Example:
            print(element_info.get_atom_fraction_whole())
            # Output: 14
        """
        return self.atom_fraction_whole

    def get_id(self):
        """
        Get the identifier associated with the element.

        Returns:
            str: The element identifier.

        Example:
            print(element_info.get_id())
            # Output: "42"
        """
        return self.id

    def get_isotopes(self):
        """
        Get information about isotopes associated with the element.

        Returns:
            str: A formatted string containing information about isotopes.

        Example:
            print(element_info.get_isotopes())
            # Output: "Index: 0 \n Isotope: C-12, Weight Percent: 0.011, ZAID: 12000 \n
            #          Index: 1 \n Isotope: C-13, Weight Percent: ..., ZAID: ..."
        """
        return "\n".join(
            [
                f"Index: {index} \n {isotope}"
                for index, isotope in enumerate(self.isotopes)
            ]
        )

    def get_atom_density(self):
        """
        Get the atom density associated with the element.

        Returns:
            float: The atom density.

        Example:
            print(element_info.get_atom_density())
            # Output: 1.5
        """
        return self.atom_density

    def get_atomic_mass_whole(self):
        """
        Get the whole number atomic mass associated with the element.

        Returns:
            float: The whole number atomic mass.

        Example:
            print(element_info.get_atomic_mass_whole())
            # Output: 12.0
        """
        return self.atomic_mass_whole

    def get_abundances(self):
        """
        Get the abundances associated with the element.

        Returns:
            str: A formatted string containing the abundances.

        Example:
            print(element_info.get_abundances())
            # Output: "98.9"
        """
        return self.abundances


class Material:
    def __init__(self, datum: Datum):
        """
        Initialize a Material object with material data.

        Parameters:
            datum (Datum): An instance of the Datum class containing various material information.

        Example:
            datum = Datum(Name="Sample Material", Formula="H2O", ...)
            material = Material(datum) #or
            material = Material(MaterialsCompendium[0])
        """
        self.comment = datum.Comment
        self.density = datum.Density
        self.acronym = datum.Acronym
        self.elements = [ElementInfo(element_data) for element_data in datum.Elements]
        self.source = datum.Source
        self.references = datum.References
        self.contact = ContactInfo(Contact.from_dict(datum.Contact))
        self.material_atom_density = datum.MaterialAtomDensity
        self.mols = [MolsInfo(mol_data) for mol_data in datum.Mols]
        self.mat_num = datum.MatNum
        self.material_weight = datum.MaterialWeight
        self.name = datum.Name
        self.verification_notes = datum.Verification_Notes
        self.formula = datum.Formula

    def __str__(self):
        """
        Return a string representation of the Material object.

        Returns:
            str: A formatted string containing material information.

        Example:
            print(material)
            # Output: "Material Name: Sample Material \n Formula: H2O \n Density: 1.0"
        """
        return (
            f"Material Name: {self.name} \n"
            + f"Formula: {self.formula} \n"
            + f"Density: {self.density}"
        )

    def get_all(self):
        """
        Get all available material information in a formatted string.

        Returns:
            str: A formatted string containing all available material information.

        Example:
            print(material.get_all())
            # Output: "Material Name: Sample Material
            #  Acronym: ...
            #  Formula: H2O
            #  Density: 1.0
            #  Material Atom Density: ...
            #  Elements: ...
            #  Mat Num: ...
            #  Material Weight: ...
            #  Comments:
            #  ...
            #  Source: ...
            #  Verification Notes: ...
            #  References: ...
            #  Contact:
            #   Name: ...
            #   Phone: ...
            #   Email: ..."
        """
        element_data = []
        for element in self.elements:
            element_data.append(element.element)
        elements_str = ", ".join(element_data)
        comments = "\n ".join(self.comment) if self.comment else ""
        acronym = ", ".join(self.acronym) if self.acronym else ""
        verification_notes = (
            "\n ".join(self.verification_notes) if self.verification_notes else ""
        )
        contact = f"Contact:\n Name: {self.contact.name} \n Phone: {self.contact.phone} \n Email: {self.contact.email}"
        references = "\n ".join(self.references)
        return (
            f"Material Name: {self.name} \n"
            + f"Acronym: {acronym} \n"
            + f"Formula: {self.formula} \n"
            + f"Density: {self.density} \n"
            + f"Material Atom Density: {self.material_atom_density} \n"
            + f"Elements: {elements_str} \n"
            + f"Mat Num: {self.mat_num} \n"
            + f"Material Weight: {self.material_weight} \n"
            + f"Comments:\n {comments} \n"
            + f"Source: {self.source} \n"
            + f"Verification Notes: {verification_notes} \n"
            + f"References: {references} \n"
            + contact
        )

    # Methods for accessing specific attributes

    def get_comment(self):
        """
        Get the comments associated with the material.

        Returns:
            list: A list of comments.

        Example:
            datum = Datum(Comment=["Comment 1", "Comment 2"])
            material = Material(datum)
            print(material.get_comment())
            # Output: ["Comment 1", "Comment 2"]
        """
        return self.comment

    def get_density(self):
        """
        Get the density associated with the material.

        Returns:
            float: The density.

        Example:
            datum = Datum(Density=1.0)
            material = Material(datum)
            print(material.get_density())
            # Output: 1.0
        """
        return self.density

    def get_acronym(self):
        """
        Get the acronym associated with the material.

        Returns:
            object: The acronym.

        Example:
            print(material.get_acronym())
            # Output: ["ABC", "DEF"]
        """
        return self.acronym

    def get_elements(self):
        """
        Get information about elements associated with the material.

        Returns:
            str: A formatted string containing information about elements.

        Example:
            print(material.get_elements())
            # Output: "Index: 0 \n Element: Carbon, ZAID: 12000, ..."
        """
        return "\n".join(
            [
                f"Index: {index} \n {element}"
                for index, element in enumerate(self.elements)
            ]
        )

    def get_source(self):
        """
        Get the source of the material.

        Returns:
            str: The source of the material.

        Example:
            print(material.get_source())
            # Output: "Some source information"
        """
        return self.source

    def get_references(self):
        """
        Get the references associated with the material.

        Returns:
            list: A list of references.

        Example:
            datum = Datum(References=["Ref 1", "Ref 2"])
            material = Material(datum)
            print(material.get_references())
            # Output: ["Ref 1", "Ref 2"]
        """
        return self.references

    def get_contact(self):
        """
        Get the contact information associated with the material.

        Returns:
            str: A formatted string containing the contact information.

        Example:
            print(material.get_contact())
            # Output: "Contact:
            #          Name: Ahnaf Tahmid Chowdhury
            #          Phone: 1234567890
            #          Email: tahmid@example.com"
        """
        contact = f"Contact:\n Name: {self.contact.name} \n Phone: {self.contact.phone} \n Email: {self.contact.email}"
        return contact

    def get_material_atom_density(self):
        """
        Get the material atom density.

        Returns:
            float: The material atom density.

        Example:
            print(material.get_material_atom_density())
            # Output: 2.5
        """
        return self.material_atom_density

    def get_mols(self):
        """
        Get information about mols associated with the material.

        Returns:
            str: A formatted string containing information about mols.

        Example:
            print(material.get_mols())
            # Output: "Mols: 100, Isotopes: H1, C12, O16"
        """
        return "\n".join(
            [
                f"Mols: {mol.get_mols()}, Isotopes: {mol.get_isotope()}"
                for mol in self.mols
            ]
        )

    def get_mat_num(self):
        """
        Get the material number.

        Returns:
            int: The material number.

        Example:
            print(material.get_mat_num())
            # Output: 12
        """
        return self.mat_num

    def get_material_weight(self):
        """
        Get the material weight.

        Returns:
            str: The material weight.

        Example:
            print(material.get_material_weight())
            # Output: "10"
        """
        return self.material_weight

    def get_name(self):
        """
        Get the name of the material.

        Returns:
            str: The name of the material.

        Example:
            print(material.get_name())
            # Output: "Sample Material"
        """
        return self.name

    def get_verification_notes(self):
        """
        Get the verification notes associated with the material.

        Returns:
            list: A list of verification notes.

        Example:
            print(material.get_verification_notes())
            # Output: ["Some verification notes"]
        """
        return self.verification_notes

    def get_formula(self):
        """
        Get the formula of the material.

        Returns:
            str: The formula of the material.

        Example:
            print(material.get_formula())
            # Output: "H2O"
        """
        return self.formula

    @classmethod
    def from_name(cls, material_name):
        """
        Create a Material object by searching for a material using its name.

        Parameters:
            material_name (str): The name of the material to search for.

        Returns:
            Material or None: The Material object if found, or None if not found.

        Example:
            material = Material.from_name("Sample Material")
            if material:
                print(material.get_all())
            else:
                print("Material not found.")
        """
        matched_materials = []
        for datum in MaterialsCompendium:
            if datum.Name == material_name:
                return cls(datum)
            elif difflib.SequenceMatcher(None, material_name, datum.Name).ratio() > 0.6:
                # Consider a match if the similarity ratio is greater than 0.6 (adjust the threshold as needed)
                matched_materials.append(datum.Name)

        if matched_materials:
            suggestions = "\n".join(matched_materials)
            print(f"Material '{material_name}' not found. Did you mean:\n{suggestions}")
        else:
            print(f"Material '{material_name}' not found in the data.")
        return None

    @classmethod
    def from_formula(cls, material_formula):
        """
        Create a Material object by searching for a material using its formula.

        Parameters:
            material_formula (str): The formula of the material to search for.

        Returns:
            Material or None: The Material object if found, or None if not found.

        Example:
            material = Material.from_formula("H2O")
            if material:
                print(material.get_all())
            else:
                print("Material not found.")
        """
        matched_materials = []

        for datum in MaterialsCompendium:
            if datum.Formula is not None and datum.Formula == material_formula:
                return cls(datum)
            elif (
                datum.Formula is not None
                and difflib.SequenceMatcher(
                    None, material_formula, datum.Formula
                ).ratio()
                > 0.6
            ):
                # Consider a match if the similarity ratio is greater than 0.6 (adjust the threshold as needed)
                matched_materials.append(datum.Formula)

        if matched_materials:
            suggestions = "\n".join(matched_materials)
            print(
                f"Material with formula '{material_formula}' not found. Did you mean:\n{suggestions}"
            )
        else:
            print(f"Material with formula '{material_formula}' not found in the data.")
        return None

    @classmethod
    def from_acronym(cls, material_acronym):
        """
        Create a Material object by searching for a material using its acronym.

        Parameters:
            material_acronym (str): The acronym of the material to search for.

        Returns:
            Material or None: The Material object if found, or None if not found.

        Example:
            material = Material.from_acronym("ABC")
            if material:
                print(material.get_all())
            else:
                print("Material not found.")
        """
        matched_materials = []
        for datum in MaterialsCompendium:
            if datum.Acronym is not None and datum.Acronym == material_acronym:
                return cls(datum)
            elif (
                datum.Acronym is not None
                and difflib.SequenceMatcher(
                    None, material_acronym, datum.Acronym
                ).ratio()
                > 0.6
            ):
                # Consider a match if the similarity ratio is greater than 0.6 (adjust the threshold as needed)
                matched_materials.append(datum.Acronym)

        if matched_materials:
            suggestions = "\n".join(matched_materials)
            print(
                f"Material with acronym '{material_acronym}' not found. Did you mean:\n{suggestions}"
            )
        else:
            print(f"Material with acronym '{material_acronym}' not found in the data.")
        return None
