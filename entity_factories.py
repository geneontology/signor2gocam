import csv
import os
from naming_conventions import NamingConvention
from entity_models import SignorEntity, SignorProtein, SignorMicroRNA, SignorComplex, SignorProteinFamily, SignorSmallMolecule
from download import SignorDownloader


class SignorGroupingFactory:
    NAME_FIELD = None
    GROUPING_CLASS = None

    def __init__(self, filename):
        self.grouping = {}
        with open(filename, "r") as f:
            data = list(csv.DictReader(f, delimiter=";"))

            for line in data:
                entities = []

                for entity in line['LIST OF ENTITIES'].split(", "):
                    entities.append(entity.strip())

                args = {
                    "signor_id": line['SIGNOR ID'],
                    "name": line[self.NAME_FIELD],
                    "entities": entities
                }
                sig_grouping = self.GROUPING_CLASS(**args)
                self.grouping[sig_grouping.id] = sig_grouping


class SignorComplexFactory(SignorGroupingFactory):
    FILENAME = "resources/SIGNOR_complexes.csv"

    def __init__(self, filename=None):
        if filename is None:
            filename = self.FILENAME
        if not os.path.exists(filename):
            filename = SignorDownloader.download_complexes()

        self.NAME_FIELD = "COMPLEX NAME"
        self.GROUPING_CLASS = SignorComplex
        SignorGroupingFactory.__init__(self, filename)
        self.complexes = self.grouping


class SignorProteinFamilyFactory(SignorGroupingFactory):
    FILENAME = "resources/SIGNOR_PF.csv"

    def __init__(self, filename=None):
        if filename is None:
            filename = self.FILENAME
        if not os.path.exists(filename):
            filename = SignorDownloader.download_families()

        self.NAME_FIELD = "PROT. FAMILY NAME"
        self.GROUPING_CLASS = SignorProteinFamily
        SignorGroupingFactory.__init__(self, filename)
        self.families = self.grouping


class SignorEntityFactory:
    complex_factory = SignorComplexFactory()
    family_factory = SignorProteinFamilyFactory()
    entity_type_map = {
        'complex': SignorComplex,
        'protein': SignorProtein,
        'mirna': SignorMicroRNA,
        'smallmolecule': SignorSmallMolecule
    }

    @classmethod
    def determine_entity(cls, entity_id: str, entity_name: str, entity_type: str) -> SignorEntity:
        if NamingConvention.is_complex(entity_id):
            return SignorEntityFactory.complex_from_id(entity_id)
        else:
            return cls.entity_type_map[entity_type](entity_id, entity_name)

    @classmethod
    def complex_from_id(cls, entity_id: str):
        possible_complexes = cls.complex_factory.complexes
        if entity_id in possible_complexes:
            return possible_complexes[entity_id]


def main():
    # TODO: parameterize SIGNOR_complexes.csv path
    complex_list = SignorComplexFactory("SIGNOR_complexes.csv")
    complex = complex_list.complexes["SIGNOR-C87"]
    complex_cc = "GO:0032991"
    # create individual for go term
    for entity in complex.entities:
        relation = "complex_cc has_part " + entity
        # create individual for entity
        # state relation
        print(relation)
    pathway_line_mech_or_reg_activity = "pc.whatevs"
    # state "pathway_line_mech_or_reg_activity ENABLED_BY complex_cc"

    # TODO: parameterize SIGNOR_PF.csv path
    pf_list = SignorProteinFamilyFactory("SIGNOR_PF.csv")
    pf = pf_list.families["SIGNOR-PF11"]
    print("PF list for " + pf.name + ":")
    for entity in pf.entities:
        print(entity)


if __name__ == "__main__":
    main()