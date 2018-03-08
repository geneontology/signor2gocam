import csv

complexes = []

class SignorComplex():
    def __init__(self, signor_id, name, entities):
        self.id = signor_id
        self.name = name
        self.entities = entities

class SignorComplexFactory():
    def __init__(self, filename):
        self.complexes = []
        with open(filename, "r") as f:
            data = list(csv.DictReader(f, delimiter=";"))


            for line in data:
                entities = []

                for entity in line['LIST OF ENTITIES'].split(", "):
                    entities.append(entity.strip())
                # print(line)

                sig_complex = SignorComplex(line['SIGNOR ID'], line['COMPLEX NAME'], entities)
                complexes.append(sig_complex)