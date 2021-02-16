import hashlib
import os
import pytest
from riborid import experiment
import shutil

@pytest.fixture(scope='class')
def create_experiment(request):
    request.cls.test_exp = experiment.Experiment(max_gap=50, max_shift=10, oligo_len=32, mt_thresh=65,
                 mt_err=3, na=100, mg=4, oligoc=150)

@pytest.mark.usefixtures("create_experiment")
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