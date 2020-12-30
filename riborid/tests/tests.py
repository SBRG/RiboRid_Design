#location of files needed for tests
DATA_DIR = 'example_data'
import sys
sys.path.append('../')
import hashlib

from riborid import design_oligos, experiment, rRNA


def test_rrna_wgbk():
    "test if rrna object is created properly with gbk input"

    rrna_test = rRNA.RRNA('ex', '23S', 'tests/example_data/Saureus_TCH1516.gb',
                          'genbank', outdir='test_out')
    rrna_checksum(rrna_test)

def test_rrna_wfa():
    " test if rrna object is created properly with rrna fasta input "
    rrna_test = rRNA.RRNA('ex', '23S', 'tests/example_data/Saureus_TCH1516_23S.fa',
                          'fasta', outdir='test_out')
    rrna_checksum(rrna_test)

def rrna_checksum(rrna_test):
    consesus_hex = 'b3c97db85f8210009be98751bf076f99'
    rrna_hex = '6e5c23c08f91ceab5967a22fcd34feb9'
    assert hashlib.md5(open(rrna_test.rrna_fa, 'rb').read()).hexdigest() == rrna_hex
    assert hashlib.md5(open(rrna_test.consensus, 'rb').read()).hexdigest() == consesus_hex


# TODO: Need to write tests for rRNA.py



