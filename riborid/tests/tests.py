#location of files needed for tests
DATA_DIR = 'example_data'
import sys
sys.path.append('../')
import hashlib

from riborid import design_oligos, experiment, rRNA


def test_rrna_wgbk():
    "test if rrna object is created properly with gbk input"

    rrna_test = rRNA.RRNA('ex', '23S', 'tests/example_data/Saureus_TCH1516.gb', 'genbank', outdir='test_out')
    rrna_checksum(rrna_test)

def test_rrna_wfa():
    " test if rrna object is created properly with rrna fasta input "
    rrna_test = rRNA.RRNA('ex', '23S', 'tests/example_data/Saureus_TCH1516_23S.fa', 'fasta', outdir='test_out')
    rrna_checksum(rrna_test)

def rrna_checksum(rrna_test):
    assert hashlib.md5(open(rrna_test.rrna_fa, 'rb').read()).hexdigest() ==
    assert hashlib.md5(open(rrna_test.consensus, 'rb').read()).hexdigest() ==


# TODO: Need to write tests for rRNA.py



