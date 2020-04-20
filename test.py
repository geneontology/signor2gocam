import unittest
import yaml
from pathway_connections import MechanismToGoMappingSet

M_FILE = "metadata/signor_mechanism_go_mapping.yaml"

class TestSignor2Gocam(unittest.TestCase):

    def test_mechanism_map_loading(self):
        with open(M_FILE) as mf:
            mappings = yaml.load(mf)
        self.assertGreater(len(mappings), 1)

    def test_mechanism_mapping_set(self):
        mapping_set = MechanismToGoMappingSet(M_FILE)
        go_term = mapping_set.go_id_by_mechanism("catalytic activity")
        self.assertEqual(go_term, "GO:0003824")

if __name__ == '__main__':
    unittest.main()