#!/usr/bin/python2
# 
# 
# 
# 
# Kim Brugger (04 Jan 2018), contact: kim@brugger.dk

import sys
import subprocess
import os
import argparse
import re 
import gzip
import json

import pprint
pp = pprint.PrettyPrinter(indent=4)

sys.path.append("/software/dev/lib//python2.7/site-packages/")
import pysam


def open_file( filename ):

    fh = None
    if ".gz" in filename:
        fh = gzip.open(filename, 'rb')
    else:
        fh = open(filename, 'rb')

    return fh


def get_next_json_line_as_dict(json_fh):
    """
    Return the next line in json_fh as a dict
    """
    line = json_fh.readline()
    line = line.strip().rstrip(",")
    return json.loads(line)


def merge_files( vcf_file, json_file, out_file):
    """ merge the content of the two files into a new vcf

    Args:
        vcf_file (str): original vcf file
        json_file (str): json nirvana annotation file

    Returns:
        None
    
    Raises:
        No exceptions are caught by the function
    """

    out_fh = open(out_file, "w")

    json_fh = open_file( json_file )
    json_fh.readline()  # Skip header line

    vcf_in = pysam.VariantFile( vcf_file )

    # add CSQ to the info field in header
    vcf_in.header.add_line( '##INFO=<ID=CSQ,Number=.,Type=String,Description="Consequence type as predicted by VEP. Format: Allele|Gene|HGNC|RefSeq|Feature|Consequence|cDNA_position|Protein_position|Amino_acids|Existing_variation|SIFT|PolyPhen|HGVSc">')

    # print the vcf header
    out_fh.write(str( vcf_in.header ))

    vcf_position_counter = 0
    json_position_counter = 0

    mnv_adjustment = 0  # Nirvana merges SNVs into MNVs in an extra position record. We need to count these to keep vcf aligned with json 

    # Loop through each of the vcf entries.
    for rec in vcf_in.fetch():

        vcf_position_counter += 1

        # Get the equivalent record from the json
        position_info = get_next_json_line_as_dict(json_fh)
        json_position_counter += 1

        vcf_chrom = rec.chrom
        json_chrom = position_info["chromosome"]

        vcf_pos = rec.pos
        json_pos = position_info["position"]

        # 2.0.10 has MNVs represented as both as individual SNVs and as MNVs
        # These are in separate json records
        # Here we skip those MNV records
        # Ideally we would include these, but this will require a major redesign of this code
        while not ((vcf_chrom == json_chrom) and (vcf_pos == json_pos)):
            mnv_adjustment += 1
            position_info = get_next_json_line_as_dict(json_fh)
            json_position_counter += 1
            json_chrom = position_info["chromosome"]
            json_pos = position_info["position"]

        # Make sure vcf record and json record positions are aligned 
        assert (vcf_chrom == json_chrom) and (vcf_pos == json_pos), "vcf and json out of sync at vcf line %d! Aborting." % vcf_position_counter

        CSQs = []

        for variant in position_info[ 'variants' ]:

            # No transcript annotation so include the unannotated vcf line
            if 'transcripts' not in variant.keys() or 'refSeq' not in variant [ 'transcripts' ]:
                continue

            for transcript in variant [ 'transcripts' ][ 'refSeq' ]:
                CSQ = []

                # Some kind of stupid representation of in/dels
                if (variant[ 'altAllele'] == "-"):
                    variant['altAllele'] = position_info[ 'refAllele' ].replace( variant['refAllele'], "", 1 )
                elif (variant[ 'refAllele'] == "-"):
                    variant['altAllele'] = position_info[ 'refAllele' ] + variant['altAllele']
                    variant['refAllele'] = position_info[ 'refAllele' ]

                # Add CSQ values from json to vcf CSQ string 
                CSQ.append( variant[ 'altAllele'])
                CSQ.append( transcript[ 'geneId'])
                CSQ.append( transcript[ 'hgnc'])
                CSQ.append( transcript[ 'transcript'])
                CSQ.append( "")
                CSQ.append( "&".join(transcript[ 'consequence']))

                CSQ.append( transcript.get('cdsPos', ""))
                CSQ.append( transcript.get( 'proteinPos', ""))
                CSQ.append( transcript.get( 'aminoAcids', ""))
                CSQ.append( "")

                sift_pred = transcript.get( 'siftPrediction', "")
                sift_score = transcript.get( 'siftScore', "")
                if sift_pred or sift_score:
                    sift = "%s(%s)" % (sift_pred, sift_score)
                else:
                    sift = ""
                CSQ.append( sift)

                poly_pred = transcript.get( 'polyPhenPrediction', "")
                poly_score = transcript.get( 'polyPhenScore', "")
                if poly_pred or poly_score:
                    poly = "%s(%s)" % (poly_pred, poly_score)
                else:
                    poly = ""
                CSQ.append( poly)

                CSQ.append( transcript.get( 'hgvsc', ""))

                CSQs.append( "|".join(CSQ) )

        rec.info['CSQ' ] = ",".join(CSQs)

        out_fh.write(str( rec )) 

    # Check that there are no more variants in json
    next_line = json_fh.readline()
    assert next_line.strip() == """],"genes":[""", "Error - unprocessed variants positions in json\n%s" % next_line

    # A check that both files contained tme same number of records (taking into account extra records in json for MNVs)
    print "Annotations: ", json_position_counter - mnv_adjustment
    print "VCF Records:", vcf_position_counter
    print "MNVs:", mnv_adjustment
    out_fh.close()

if __name__ == "__main__":

    parser = argparse.ArgumentParser(description='tool for creating a VEP like output from a vcf+nirvan json file ')

    parser.add_argument('-j', '--json-file',  help="nirvana json file")
    parser.add_argument('-v', '--vcf-file',    help="vcf file")
    parser.add_argument('-o', '--out-file',    help="out filename, otherwise vcf-file.nirvana.vcf")

    args = parser.parse_args()

    merge_files( vcf_file = args.vcf_file, json_file = args.json_file, out_file = args.out_file)
