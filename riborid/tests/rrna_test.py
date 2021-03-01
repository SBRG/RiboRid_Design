
import hashlib
import os
import pytest
from riborid import rRNA, experiment  # ,design_oligos, experiment,
import shutil
import pandas as pd

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


@pytest.fixture(scope='class')
def create_experiment(request):

    os.mkdir('exp_test')
    request.cls.test_exp = experiment.Experiment(max_gap=50, max_shift=10, oligo_len=32, mt_thresh=65,
                                                 mt_err=3, na=100, mg=4, oligoc=150)
    yield
    shutil.rmtree('exp_test')


@pytest.mark.usefixtures("create_experiment", "create_rrna")
class TestExp:

    def test_init_mgap(self):
        assert self.test_exp.max_gap, 50

    def test_init_mshift(self):
        assert self.test_exp.max_shift, 10

    def test_init_oligolen(self):
        assert self.test_exp.oligo_len, 32

    def test_init_mtthresh(self):
        assert self.test_exp.mt_thresh, 65

    def test_init_mterr(self):
        assert self.test_exp.mt_err, 3

    def test_init_na(self):
        assert self.test_exp.na, 100

    def test_init_mg(self):
        assert self.test_exp.mg, 4

    def test_init_oligoc(self):
        assert self.test_exp.oligoc, 150

    def test_log_exp(self):
        log_file = 'exp_test/test_log.txt'
        self.test_exp.log_exp(log_file)
        assert os.path.isfile(log_file)

        exp_hex = '1bd04aae1fc5df33cf57d140c0ad1c05'
        assert hashlib.md5(open(log_file, 'rb').read()).hexdigest() == exp_hex

    def test_exp_gapfill(self):
        df = self.test_exp.gapfill(self.test_rrna)
        df.melting_temp = pd.Series.round(df.melting_temp)
        test_df = pd.read_csv('test_data/Saureus_oligosdf_none.csv', index_col=0)
        test_df.melting_temp = pd.Series.round(test_df.melting_temp)
        assert df.equals(test_df)
