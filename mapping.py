import argparse
import subprocess


def main(args):
    subprocess.call("bwa mem ~/Desktop/project/GCF_000001215.4_Release_6_plus_ISO1_MT_genomic.fna " + args.map_seq + " > " + args.map_seq + ".sam", shell=True)
    subprocess.call("samtools sort " + args.map_seq + ".sam > " + args.map_seq + ".bam", shell=True)
    subprocess.call("samtools index " + args.map_seq + ".bam", shell=True)
    print "Done"


parser = argparse.ArgumentParser(description='')
parser.add_argument('map_seq', metavar = 'filename.FASTA', help='Sequence file to be mapped')
args = parser.parse_args()
main(args)
