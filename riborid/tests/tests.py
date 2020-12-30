#location of files needed for tests
DATA_DIR = 'example_data'
import sys
sys.path.append('../')

from riborid import design_oligos, experiment, rRNA


def test_loadgb():
    # test if genbank file is loaded properly
    gb_dir = 'example_data/Saureus_TCH1516.gb'


def test_consensus():
    # test if correct consensus sequence is generated from genbank
    pass


def test_faout():
    # test if proper output fasta file is generated
    pass
# TODO: Need to write tests for rRNA.py



