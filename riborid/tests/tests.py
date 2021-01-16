
import hashlib
import unittest
import os.path

from riborid import rRNA  # ,design_oligos, experiment,


class TestRrna(unittest.TestCase):
    """ Test initialization of rRNA class """

    def setUP(self):
        print('Setting Up rrna tests...')
        self.rrna_test = rRNA.RRNA(rrna_fa='tests/test_data/Saureus_TCH1516_23S.fa',
                                   name='ex', rtype='23S', outdir='test_out',
                                   pre=150, oligos_df=None)
        # checksum of properly created files
        self.consesus_hex = 'b3c97db85f8210009be98751bf076f99'
        self.rrna_hex = '6e5c23c08f91ceab5967a22fcd34feb9'

    def tearDown(self):
        print('Cleaning Up')
        pass

    def test_init_name(self):
        self.assertEqual(self.rrna_test.name, 'ex')

    def test_init_rtype(self):
        self.assertEqual(self.rrna_test.rtype, '23S')

    def test_init_outdir(self):
        self.assertEqual(self.rrna_test.outdir, 'rrd')

    def test_init_pre(self):
        self.assertEqual(self.rrna_test.pre, 150)

    def test_init_na_oligosdf(self):
        self.assertIsNone(self.rrna_test.oligos_df)

    def test_init_outpath(self):
        self.assertTrue(os.path.exists(self.rrna_test.rrna_fa))

    # check if the files are created properly
    def test_init_rrnafa(self):
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
