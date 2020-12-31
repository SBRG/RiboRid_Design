
import hashlib
import unittest
import os.path

from riborid import rRNA  # ,design_oligos, experiment,


class TestRrna(unittest.TestCase):
    """ Test initialization of rRNA class """

    def setUP(self):
        print('Setting Up rrna tests...')
        self.rrna_test = rRNA.RRNA(rrna_fa='tests/example_data/Saureus_TCH1516_23S.fa',
                                   name='ex', rtype='23S', outdir='test_out',
                                   pre=150, oligos_df=None)
        # checksum of properly created files
        self.consesus_hex = 'b3c97db85f8210009be98751bf076f99'
        self.rrna_hex = '6e5c23c08f91ceab5967a22fcd34feb9'

    def tearDown(self):
        print('Cleaning Up')
        pass

    def test_rrna_init(self):
        """ test if rrna object is created properly with rrna fasta input """

        print('Testing rrna init')
        self.assertEqual(self.rrna_test.name, 'ex')
        self.assertEqual(self.rrna_test.rtype, '23S')
        self.assertEqual(self.rrna_test.outdir, 'rrd')
        self.assertEqual(self.rrna_test.pre, 150)
        self.assertIsNone(self.rrna_test.oligos_df)
        self.assertTrue(os.path.exists(self.rrna_test.rrna_fa))

        # check if the files are created properly
        self.assertEqual(hashlib.md5(open(self.rrna_test.rrna_fa, 'rb').read()).hexdigest(),
                         self.rrna_hex)
        self.assertEqual(hashlib.md5(open(self.rrna_test.consensus, 'rb').read()).hexdigest(),
                         self.consesus_hex)

    # TODO: generate a test oligo_df and write this test.
    def test_rrna_oligodf(self):
        """ test initializing the rrna object with exisiting oligo_df """

        print('Testing rrna init with oligodf')
        pass


class TestExperiment(unittest.TestCase):
    """ Test initialization of experiment class """
    pass


class TestDesignOligos(unittest.TestCase):
    """ Test design oligos """
    pass


if __name__ == '__main__':
    TestRrna.main()
    TestExperiment.main()
    TestDesignOligos.main()

# TODO: Need to write tests for rRNA.py
