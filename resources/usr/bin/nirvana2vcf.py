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



def readin_json( filename ):

    json_data = open_file( filename )

    annotations = json.load(json_data)

    return annotations


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

    annotations = readin_json( json_file )
    
    vcf_in = pysam.VariantFile( vcf_file )

    print vcf_file
    print json_file

    # add CSQ to the info field
    vcf_in.header.add_line( '##INFO=<ID=CSQ,Number=.,Type=String,Description="Consequence type as predicted by VEP. Format: Allele|Gene|HGNC|RefSeq|Feature|Consequence|cDNA_position|Protein_position|Amino_acids|Existing_variation|SIFT|PolyPhen|HGVSc">')

#    vcf_in.header.remove_line("CSQT")

#    vcf_in.header.add_line( '##source_nirvana= Nirvana X.XXXX')

    # print the vcf header
    #sys.stdout.write(str( vcf_in.header ))
    out_fh.write(str( vcf_in.header ))

    variant_counter = 0

    # Loop through each of the vcrf entries.
    for rec_index, rec in enumerate(vcf_in.fetch()):
#        print rec

        variant_info =  annotations[ "positions" ][ variant_counter ]

        vcf_chrom = rec.chrom
        json_chrom = variant_info["chromosome"]

        vcf_pos = rec.pos
        json_pos = variant_info["position"]

        assert (vcf_chrom == json_chrom) and (vcf_pos == json_pos), "vcf and json out of sync at vcf line %d! Aborting." % variant_counter
            #print "Mismatch:"
            #print "rec_index", rec_index
            #print "variant_counter", variant_counter
            #print
            #print "VCF CHR", vcf_chrom
            #print "JSON CHR", json_chrom
            #print
            #print "VCF POS", vcf_pos
            #print "JSON POS", json_pos
            #print

 #       my ($Allele,$ENS_gene, $HGNC,$RefSeq,$feature,$effects,$CDS_position,$Protein_position,$Amino_acid,$Existing_variation,$SIFT,$PolyPhen,$HGVSc,$Distance) = @$CSQ;

        CSQs = []
        
        for variant in variant_info[ 'variants' ]:
            #pp.pprint( variant)
            
            if 'transcripts' not in variant or 'refSeq' not in variant [ 'transcripts' ]:
                #sys.stdout.write(str( rec ))
                out_fh.write(str( rec ))
                continue

           
#            pp.pprint( variant )
            for transcript in  variant [ 'transcripts' ][ 'refSeq' ]:
                CSQ = []

                # Some kind of stupid representation of in/dels
                if (variant[ 'altAllele'] == "-"):
                    variant['altAllele'] = variant_info[ 'refAllele' ].replace( variant['refAllele'], "", 1 )
                elif (variant[ 'refAllele'] == "-"):
                    variant['altAllele'] = variant_info[ 'refAllele' ] + variant['altAllele']
                    variant['refAllele'] = variant_info[ 'refAllele' ]
                    
                

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

        
        #print rec.info['CSQ' ]

        #sys.stdout.write(str( rec ))
        out_fh.write(str( rec ))  # Does an output file stop truncated records? 
        # sys.stdout.flush()  # Does flush stop truncated records? Nope

        variant_counter += 1


    print "Annotations: ", len(annotations[ "positions" ])
    print "VCF Records:", variant_counter
    out_fh.close()


if __name__ == "__main__":

    parser = argparse.ArgumentParser(description='tool for creating a VEP like output from a vcf+nirvan json file ')

    parser.add_argument('-j', '--json-file',  help="nirvana json file")
    parser.add_argument('-v', '--vcf-file',    help="vcf file")
    parser.add_argument('-o', '--out-file',    help="out filename, otherwise vcf-file.nirvana.vcf")

    args = parser.parse_args()

#    args.json_file = '/data/projects/matt/variant_comparison/nirvana_annotation/X004856.json.gz'
#    args.vcf_file  = '/data/projects/matt/variant_comparison/nirvana_annotation/X004856.vcf.gz'
#    args.json_file = 'test.json'


    merge_files( vcf_file = args.vcf_file, json_file = args.json_file, out_file = args.out_file)

