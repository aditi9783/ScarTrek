import filecmp
from scartrek import find_scars
import unittest


class TestScarTrekIntegration(unittest.TestCase):
    @unittest.skip("Skipping in remote build until we get test data committed")
    def test_assert_values(self):
        basedir = "./test1"
        find_scars.main(["-i", basedir,
                         "--geneseq", "../../reference/H37Rv_genes.txt",
                         "--protseq", "../../reference/H37Rv_proteins_from_genbank.txt"])

        samples = ["sample1", "sample2", "sample3"]
        files_to_compare = ["frameshifts2", "indels", "genewise.mutations2", "stopcodon_causing_mut2"]
        scar_results = {}
        for sample in samples:
            for filename in files_to_compare:
                sample_basedir = basedir + "/" + sample + "/mapped"

                expected_file = sample_basedir + "/oldfiles/" + sample + "." + filename
                actual_file = sample_basedir + "/" + sample + "." + filename

                scar_results["{} {}".format(sample, filename)] = filecmp.cmp(expected_file, actual_file)

        failed_tests = [sample for sample, result in scar_results.items() if not result]
        if failed_tests:
            print("Failed ")
            print([sample for sample in failed_tests])
        else:
            print("All passed!")
