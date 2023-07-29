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
        self.name = contact_data.Name
        self.phone = contact_data.Phone
        self.email = contact_data.Email

    def __str__(self):
        return f"Name: {self.name}, Phone: {self.phone}, Email: {self.email}"

    def get_name(self):
        return self.name

    def get_phone(self):
        return self.phone

    def get_email(self):
        return self.email


class MolsInfo:
    def __init__(self, mol_data: Mol):
        self.mols = mol_data.Mols
        self.isotope = mol_data.Isotope
        self.element = mol_data.Element

    def __str__(self):
        return f"Mols: {self.mols}, Isotopes: {self.isotope}, Elements: {self.element}"

    def get_mols(self):
        return self.mols

    def get_isotope(self):
        return self.isotope

    def get_element(self):
        return self.element


class IsotopeInfo:
    def __init__(self, isotope_data: Isotope):
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
        return f"Isotope: {self.isotope}, Weight Percent: {self.weight_percent}, ZAID: {self.zaid}"

    def get_all(self):
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
        return f"Element: {self.element} \n ZAID: {self.zaid} \n Isotopes:{', '.join([isotope.isotope for isotope in self.isotopes])}"

    def get_all(self):
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

    def get_weight_fraction(self):
        return self.weight_fraction

    def get_weight_fraction_whole(self):
        return self.weight_fraction_whole

    def get_non_isotopic(self):
        return self.non_isotopic

    def get_element(self):
        return self.element

    def get_atomic_mass(self):
        return self.atomic_mass

    def get_zaid(self):
        return self.zaid

    def get_atom_fraction(self):
        return self.atom_fraction

    def get_atom_density_whole(self):
        return self.atom_density_whole

    def get_atom_fraction_whole(self):
        return self.atom_fraction_whole

    def get_id(self):
        return self.id

    def get_isotopes(self):
        return "\n".join(
            [
                f"Index: {index} \n {isotope}"
                for index, isotope in enumerate(self.isotopes)
            ]
        )

    def get_atom_density(self):
        return self.atom_density

    def get_atomic_mass_whole(self):
        return self.atomic_mass_whole

    def get_abundances(self):
        return self.abundances


class Material:
    def __init__(self, datum: Datum):
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
        return (
            f"Material Name: {self.name} \n"
            + f"Formula: {self.formula} \n"
            + f"Density: {self.density}"
        )

    def get_all(self):
        element_data = []
        for element in self.elements:
            element_data.append(element.element)
        elements_str = ", ".join(element_data)
        comments = "\n ".join(self.comment)
        acronym = ", ".join(self.acronym)
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
            + f"Verification Notes: {self.verification_notes} \n"
            + f"References: {references} \n"
            + contact
        )

    def get_comment(self):
        return self.comment

    def get_density(self):
        return self.density

    def get_acronym(self):
        return self.acronym

    def get_elements(self):
        return "\n".join(
            [
                f"Index: {index} \n {element}"
                for index, element in enumerate(self.elements)
            ]
        )

    def get_source(self):
        return self.source

    def get_references(self):
        return self.references

    def get_contact(self):
        contact = f"Contact:\n Name: {self.contact.name} \n Phone: {self.contact.phone} \n Email: {self.contact.email}"
        return contact

    def get_material_atom_density(self):
        return self.material_atom_density

    def get_mols(self):
        return self.mols

    def get_mat_num(self):
        return self.mat_num

    def get_material_weight(self):
        return self.material_weight

    def get_name(self):
        return self.name

    def get_verification_notes(self):
        return self.verification_notes

    def get_formula(self):
        return self.formula

    @classmethod
    def from_name(cls, material_name):
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
