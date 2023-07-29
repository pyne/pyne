"""
Separate classes and functions to perform specific tasks
to retrieve specific information from the Materials Compendium
"""
from dataclasses import dataclass
import difflib
from .materials_compendium import MaterialsCompendium, Datum

@dataclass
class Material:
    def __init__(self, datum: Datum):
        self.comment = datum.Comment
        self.density = datum.Density
        self.acronym = datum.Acronym
        self.elements = datum.Elements
        self.source = datum.Source
        self.references = datum.References
        self.contact = datum.Contact
        self.material_atom_density = datum.MaterialAtomDensity
        self.mols = datum.Mols
        self.mat_num = datum.MatNum
        self.material_weight = datum.MaterialWeight
        self.name = datum.Name
        self.verification_notes = datum.Verification_Notes
        self.formula = datum.Formula

    def __str__(self):
        return f"Material Name: {self.name} \nFormula: {self.formula} \nDensity: {self.density} \nSource: {self.source}"

    def get_comment(self):
        return self.comment
    
    def get_density(self):
        return self.density
    
    def get_acronym(self):
        return self.acronym
    
    def get_elements(self):
        return self.elements
    
    def get_source(self):
        return self.source
    
    def get_references(self):
        return self.references
    
    def get_contact(self):
        return self.contact
    
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
            print(f"Material '{material_name}' not found. Did you mean:\n{suggestions}?")
        else:
            print(f"Material '{material_name}' not found in the data.")
        return None

    @classmethod
    def from_formula(cls, material_formula):
        matched_materials = []
        for datum in MaterialsCompendium.data:
            if datum.Formula == material_formula:
                return cls(datum)
            elif difflib.SequenceMatcher(None, material_formula, datum.Formula).ratio() > 0.6:
                # Consider a match if the similarity ratio is greater than 0.6 (adjust the threshold as needed)
                matched_materials.append(datum.Name)

        if matched_materials:
            suggestions = "\n".join(matched_materials)
            print(f"Material with formula '{material_formula}' not found. Did you mean:\n{suggestions}?")
        else:
            print(f"Material with formula '{material_formula}' not found in the data.")
        return None