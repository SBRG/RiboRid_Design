
import hashlib
import os.path
import pytest
from riborid import rRNA  # ,design_oligos, experiment,
import shutil


@pytest.fixture(scope='class')
def create_rrna(request):
    request.cls.test_rrna = rRNA.RRNA(rrna_fa='test_data/Saureus_TCH1516_23S.fa',
                     name='ex', rtype='23S', outdir='test_out',
                     pre=150, oligos_df=None)

    yield

    shutil.rmtree('test_out')


@pytest.mark.usefixtures("create_rrna")
class TestRRNA:

    def test_init_name(self):
        assert self.test_rrna.name == 'ex'

    def test_init_rtype(self):
        assert self.test_rrna.rtype == '23S'

    def test_init_outdir(self):
        assert self.test_rrna.outdir == 'test_out'

    def test_init_pre(self):
        assert self.test_rrna.pre == 150

    def test_init_na_oligosdf(self):
        assert self.test_rrna.oligos_df is None

    def test_init_outpath(self):
        assert os.path.exists(self.test_rrna.rrna_fa)

    # check if the files are created properly
    def test_init_rrnafa(self):
        # checksum of properly created files
        consesus_hex = '61228d5c3eeb6342805a4a4d1eed2aea'
        assert hashlib.md5(open(self.test_rrna.consensus, 'rb').read()).hexdigest() == consesus_hex

# # TODO: generate a test oligo_df and write this test.
# def test_rrna_oligodf(self):
#     """ test initializing the rrna object with exisiting oligo_df """
#
#     print('Testing rrna init with oligodf')
#     pass
#
#
# class TestExperiment():
#     """ Test initialization of experiment class """
#     pass
#
#
# class TestDesignOligos():
#     """ Test design oligos """
#     pass


# if __name__ == '__main__':
#     TestRrna.main()
#     TestExperiment.main()
#     TestDesignOligos.main()

# TODO: Need to write tests for rRNA.py
