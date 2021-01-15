#!/usr/bin/env nextflow

""" Processing pipeline:
1. Download the genbank files from SRA. Info in ftpfilepaths.txt

1a. If fasta file, run prokka
2. Push the genbank files through oligo_design.
3. Save the rRNA fasta file. And the OligoDF.
4. Delete everything else.
"""

/// *******************************
// * Step 1: Download genbank files *
// *******************************

params.ftpfilepath = '/home/saugat/Desktop/gb_test/testpaths.txt'
params.stage = '/home/saugat/Desktop/gb_test/stage'

Channel
	.fromPath(params.ftpfilepath)
	.splitText()
	.map{ it.trim() }
	.set{ ftplinks_ch }

process download_gb {
	
	publishDir "${params.stage}"

	input:
	val download_path from ftplinks_ch
	
	output:
	file '*gbff.gz' into result

	shell:

	'''
	OUTPUT=$(basename !{download_path})
	curl !{download_path} -o rando_text.gbff.gz
	'''
}
//curl !{download_path} -o "/home/saugat/Desktop/gb_test/stage/$OUTPUT"	
//result.subscribe {println it}

