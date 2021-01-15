#!/usr/bin/env nextflow

""" Processing pipeline:
1. Download the genbank files from SRA. Info in ftpfilepaths.txt

1a. If fasta file, run prokka
2. Push the genbank files through oligo_design.
3. Save the rRNA fasta file. And the OligoDF.
4. Delete everything else.
"""

/// *******************************
// * Step 1: Parse metadata file *
// *******************************

params.ftpfilepath = '/home/saugat/Desktop/gb_test/testpaths.txt'
params.stage = '/home/saugat/Desktop/gb_test/stage'

// Ensure file exists
File ftplist = new File(params.ftpfilepath)
assert(ftplist.exists())

allLines = ftplist.readLines()

/*proces inandout {
	input:
	file download_path from allLines

	output:
	set 
}
*/
process download_gb {
	
	input:
	val download_path from allLines
	
	output:
	stdout into result

	shell:

	'''
	OUTPUT=$(basename !{download_path})
	echo  "$OUTPUT" 
	echo !{download_path} 
	'''
}
	
result.subscribe {println it}

