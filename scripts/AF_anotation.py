#!/usr/bin/env python
"""
Introduction: Allele frequency counter
------------

Input a path to a bed or vcf file, this script will annotate the with the following values:

AF_MAX       (maximum MAF from 1000genomes)
UK10K_cohort (alternative allele frequency in UK10K_COHORT release version 20130116)
ESP_ALL      (MAF from ESP project 6500 exomes)
rsID_137     (the rs ID from the dbSNP version 137)


Source files:
-------------

/nfs/users/nfs_r/rr11/codes/resources/MAF/RELEASE-21-04-2019/EPS/final/ESP6500SI_21012013_MAF.bed.gz
/nfs/users/nfs_r/rr11/codes/resources/MAF/RELEASE-21-04-2019/UK10K_COHORT/final/UK10K_COHORT_20130116_AF.bed.gz
/nfs/users/nfs_r/rr11/codes/resources/MAF/RELEASE-21-04-2019/1KG/final/annots-rsIDs-AFs.2012-07-19.tab.gz
/nfs/users/nfs_r/rr11/codes/resources/MAF/RELEASE-21-04-2019/dbSNP/final/annots-rsIDs-dbSNPv137.2012-09-13.tab.gz


How the source files were generated?
------------------------------------

/nfs/users/nfs_r/rr11/codes/resources/MAF/RELEASE-21-04-2019/EPS/README
/nfs/users/nfs_r/rr11/codes/resources/MAF/RELEASE-21-04-2019/UK10K_COHORT/README
/nfs/users/nfs_r/rr11/codes/resources/MAF/RELEASE-21-04-2019/1KG/README
/nfs/users/nfs_r/rr11/codes/resources/MAF/RELEASE-21-04-2019/dbSNP/README


Input format:
-------------

*BED format

This could be the output of the DNG pipeline (version 0.7.0) or any other tabular format as long as it include
chr --> first column
pos --> second column
ref --> fourth column
alt --> fifth column

To change the order of these columns, use option -c (or --columnIndexes ) to indicate the order of chr,pos,ref,alt
(1,2,3,4)


*VCF format

The script accept both *.vcf or *.vcf.gz formats (version 4.0 or 4.1). It will try to automatically detect if the
input file is a bed or vcf by reading the line of the VCF file (which should start with ##fileformat=VCFv4.0
or ##fileformat=VCFv4.1). If not, then it will consider as a tab separated files.


Output format:
--------------

If the input file is a bed file (or any tab-based format): the script will add four columns. However, in
case of VCF file, the script will add the values to the INFO columns with the following tag names:

	AF_MAX
	UK10K_cohort
	ESP_ALL
	rsID_137

Also, in the VCF file case, four lines will be added to the header section:
	##INFO=<ID=1KG_MAF,Number=1,Type=String,Description="AF_MAX from 1000g,/nfs/users/nfs_s/sa9/chd/resources/MAF/RELEASE-21-01-2013/1KG/README">
	##INFO=<ID=UK10K_cohort_MAF,Number=1,Type=String,Description="UK10K_cohort_MAF,/nfs/users/nfs_s/sa9/chd/resources/MAF/RELEASE-21-01-2013/UK10K_COHORT/README">
	##INFO=<ID=UK10K_cohort_MAF,Number=1,Type=String,Description="ESP_ALL_MAF,/nfs/users/nfs_s/sa9/chd/resources/MAF/RELEASE-21-01-2013/UK10K_COHORT/README">
	##INFO=<ID=rsID_137,Number=1,Type=String,Description="rsID_137,/nfs/users/nfs_s/sa9/chd/resources/MAF/RELEASE-21-01-2013/1KG/README">


Indel matching:
---------------

There are two different matching strategy as follows:

A) An exact match which work as follows:
    - First , the script tried to match by (chrom, pos, ref, alt) for both SNVs and INDELs (the only option for SNVs)
    - For indels, if not exact I match found , the script tries to match by (chrom, pos, slice, direction)
     where  slice is DNA nuecoltide difference between ref and alt and direction is ins or del.

B) 'Lenient' match for indels where
    - first script tries to match by the above exact matching, then it looks for matching for indels  within +/- 10 bp
    to match slice and direction.
	- if there are multiple indels with the same slice and direction in +/- 10 bp of the source indel then
		- chose the closer
		- if still there are multiple indels at the same distance with same slice and direction
			- get the one with lower Allele frequency

Usage:
------

	python maf.py -i file.bed -c 0,1,3,4 -m exact

	[-i] option is the path to the input file. Could be bed , vcf, bed.gz or vcf.gz (if empty, script will read stdin)
	[-c] option is the indexes of 4 important columns used by the script (chr, pos, ref and alt). This is a zero-based index.
	[-t] option is the matching strategy for indels (see above for details.)


"""

import pyVCFtoolbox as tlbx
import sys
import os
import argparse


try:
    import pysam
except ImportError:
    try:
        sys.path.append("/nfs/users/nfs_r/rr11/.local/lib/python2.6/site-packages/pysam-0.6-py2.6-linux-x86_64.egg/pysam")
        import pysam
    except ImportError:
        print >> sys.stderr, "pysam module is not installed. Please try to install it and start again."
        print >> sys.stderr, "http://code.google.com/p/pysam/"
        sys.exit(1)

MAF_BED_dict = {}

def checkMAFFiles(MAFBedPath):
    """
    This function check if the Master MAF table file exist then
    check if all the paths in the master table are also exist.
    """
    try:
        f = tlbx.openFile(MAFBedPath)
        for line in f:
            if not line.startswith("#"):
                line = line.strip().split('\t')
                try:
                    j = tlbx.openFile(line[1])
                    j.close()
                except IOError:
                    print >> sys.stderr, "Couldn't open the [%s] MAF at [%s] " % line
                    sys.exit(1)
    except Exception as e:
        print >> sys.stderr, e
        sys.exit(1)

    return True


def loadAFs(MAF_BED_path):
    """
    For all the paths in the Master MAF table file, load the their tabix index file
    to randomly access overlapping variants in the source file.
    """
    global MAF_BED_dict
    f = open(MAF_BED_path,'r')
    for line in f:
        if not line.startswith("#") and line != '':
            line = line.strip().split("\t")
            if len(line) != 4:
                #TODO : I should have validated BED file before starting the pipeline. This is a belated checkpoint
                print >> sys.stderr, line
                print >> sys.stderr, "Please check the format of your file for BED paths %s. Error len(line) != 4" % MAF_BED_path
                sys.exit(1)
            else:
                if os.path.isfile(line[1]):
                    if not MAF_BED_dict.has_key(line[0]):
                        try:
                            MAF_BED_dict[line[0]] = {'tabix':pysam.Tabixfile(line[1]),
                                                     'indices':line[2],
                                                     'vcfHeader':line[3]
                                                     }
                        except Exception as e:
                            print >> sys.stderr, e
                            print >> sys.stderr, "Unexpected file format [%s]" % line[1]
                            sys.exit(1)
                else:
                    print >> sys.stderr, "Can't find the [%s] MAF file [%s] in the MAF table [%s]" % (line[0], line[1] , MAF_BED_path)
                    sys.exit(1)
    return MAF_BED_dict

def findMatch(sourceVarInfo, targetVariants):
    """
    For every variant in a list of variants found in the target file:
        Find an exact match I if the source variant is SNV or
        if the source variant is an indel then try to find an exact
        match then exact match II.
    """
    foundExactMatch = None
    # if source var type is snv then get exact match only
    if sourceVarInfo['type'] == 'snv':
        for targetVarInfo, impValue in targetVariants:
            if targetVarInfo['type'] == 'snv':
                if sourceVarInfo['chr'] == targetVarInfo['chr'] and \
                                sourceVarInfo['pos'] == targetVarInfo['pos'] and \
                                sourceVarInfo['ref'] == targetVarInfo['ref'] and \
                                sourceVarInfo['alt'] == targetVarInfo['alt']:
                    return impValue
        return foundExactMatch # return None here. Couldn't find any exact match for the source snv

    # if source var type is indel then first,
    # try to find the first exact match and return its value (MAF value or rs Id)
    if sourceVarInfo['type'] == 'indel':
        for targetVarInfo,impValue in targetVariants:
            if targetVarInfo['type'] == 'indel':
                if sourceVarInfo['chr'] == targetVarInfo['chr'] and \
                                sourceVarInfo['pos'] == targetVarInfo['pos']:
                    if sourceVarInfo['ref'] == targetVarInfo['ref'] and \
                                    sourceVarInfo['alt'] == targetVarInfo['alt']:
                        return impValue
                    elif sourceVarInfo['slice'] == targetVarInfo['slice'] and \
                                    sourceVarInfo['direction'] == targetVarInfo['direction']:
                        return impValue
        # Return None here , couldn't find any exact I or II
        return foundExactMatch

def parseMAF(line, indices):
    """
    This function extracts the chr,pos,ref and alt from four different pre-processed MAF files.
    The impValue are usually the MAF score but it could the rsId from the dbSNP file.
    :param line: # line from  a MAF bed file
    :param indices:  the chrom,pos,ref,alt abd af column indices
    """
    chrom    = line[indices[0]]
    pos      = line[indices[1]]
    ref      = line[indices[2]]
    alt      = line[indices[3]]
    impValue = line[indices[4]] # the impValue here is the AF_MAX
    return [chrom,pos,ref,alt,impValue]

def getMAF(sourceVarInfo, mafType):
    """
    Given the start and end of a source variants, this function will return all
    overlapping variants in the target file (MAF or dbSNP).
    """
    targetVariants = []
    lines = []
    mafTabixFile = MAF_BED_dict[mafType]['tabix']
    chrom = sourceVarInfo['chr']
    start = sourceVarInfo['start']
    end   = sourceVarInfo['end']

    try:
        lines = mafTabixFile.fetch(chrom,start,end)
    except:
        pass
    for line in lines:
        line = line.strip().split("\t")
        indices = [int(x) for x in MAF_BED_dict[mafType]['indices'].split(",")]
        chrom,pos,ref,alt,impValue = parseMAF(line, indices)
        targetVarInfo = tlbx.extractVarInfo(chrom,pos,ref,alt)
        targetVariants.append([targetVarInfo,impValue])
    return targetVariants

def printVCFHeader():
    """
    Get the path to the file with additional VCF headers and print them.
    """
    for key in sorted(MAF_BED_dict):
        print MAF_BED_dict[key]['vcfHeader'].strip()



def parseInputFile( f_path, columnIndexes):
    f = tlbx.openFile(f_path)

    isValidColumns = False
    lineCounter = 0
    for line in f:
        lineCounter += 1
        originalLine = line
        if lineCounter == 1:
            tlbx.isHeaderExists(line)
            f_type  = tlbx.isBEDorVCF(line)
        if not line.startswith("#"):
            line  = line.strip().split("\t")
            if not isValidColumns:
                # executed once
                isValidColumns = tlbx.isValidColumns(columnIndexes,line)
                chrIdx,posIdx,refIdx,altIdx = [int(x) for x in columnIndexes.split(',')]
            chrom = tlbx.removeChrPrefix(line[chrIdx]).strip(' ') #if any
            pos = line[posIdx]
            ref = line[refIdx]
            # in the case of multi-allelic variants at a single locus, deal with each alt separately
            alts = tlbx.getAlts(line[altIdx])
            all_AFLists = []
            for alt in alts:
                # get the var info (start, end, type, slice, direction ,etc)
                sourceVarInfo = tlbx.extractVarInfo(chrom,pos,ref,alt)

                #get MAF values and dbSNP rsId
                current_alt_AFList = []
                for mafType in sorted(MAF_BED_dict):
                    mafReturned = None
                    targetVariants = getMAF(sourceVarInfo, mafType)
                    if targetVariants: # not empty
                        mafReturned = findMatch(sourceVarInfo, targetVariants)
                        pass
                    if mafReturned:
                        current_alt_AFList.append(mafReturned)
                    else:
                        current_alt_AFList.append('.')
                all_AFLists.append(current_alt_AFList)
            # print results
            if f_type == 'bed':
                newline = originalLine.strip("\n") + "\t"
                newline += '\t'.join(all_AFLists[0])

                print newline
            elif f_type == 'vcf':
                newInfo = ""
                newColumnNames = [x for x in sorted(MAF_BED_dict)]
                for i,tag in enumerate(newColumnNames):
                    newInfo += ";%s=%s" % (tag,  ",".join([x[i] for x in all_AFLists]))
                line[7] += newInfo
                print "\t".join(line)
        else:
            if f_type == 'bed':
                newColumnNames = [x for x in sorted(MAF_BED_dict)]
                newline = originalLine.strip() + "\t" + "\t".join(newColumnNames)
                print newline
            elif f_type == 'vcf':
                if not line.startswith("#CHROM"):
                    print line.strip()
                else:
                    printVCFHeader()
                    print line.strip()
    f.close()


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('-i', '--inputFile', help='A path to the VCF or BED file or stdin if empty', default=None)
    parser.add_argument('-c', '--columnIndexes', help='A 0-based index for chrom, pos, ref and alt columns'
                                                    ' (comma separated)', default='0,1,3,4')
    parser.add_argument('-m', '--mafBedPath', help='A path to tab-separated text file with MAF name and path',
                        required=True)

    args           = parser.parse_args()
    f_path         = args.inputFile
    columnIndexes  = args.columnIndexes
    mafBedPath     = args.mafBedPath

    if args.mafBedPath and (args.inputFile or not sys.stdin.isatty()):
        checkMAFFiles(mafBedPath)   # check if the master MAF table exist and all of its content or exit
        loadAFs(mafBedPath)         # load the tabix index files for the MAF files
        parseInputFile( f_path , columnIndexes) # annotate the input file
    else:
        parser.print_help()
        sys.exit(0)