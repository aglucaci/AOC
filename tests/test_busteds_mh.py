import unittest

from scripts.busteds_mh import count_busteds_mh_branches


class CountBustedsMhBranchesTests(unittest.TestCase):
    def test_current_hyphy_labeled_mapping(self):
        data = {
            "tested": {
                "0": {
                    "branch1": "background",
                    "branch2": "background",
                    "branch3": "test",
                }
            },
            "background": 1,
        }

        self.assertEqual(count_busteds_mh_branches(data), (1, 2))

    def test_labeled_mapping_across_partitions(self):
        data = {
            "tested": {
                "0": {"branch1": "test", "branch2": "background"},
                "1": {"branch3": "test"},
            }
        }

        self.assertEqual(count_busteds_mh_branches(data), (2, 1))

    def test_legacy_separate_branch_lists(self):
        data = {
            "tested": {"0": ["branch1", "branch3"]},
            "background": {"0": ["branch2"]},
        }

        self.assertEqual(count_busteds_mh_branches(data), (2, 1))


if __name__ == "__main__":
    unittest.main()
