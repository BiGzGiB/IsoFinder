import argparse
import subprocess


def main(args):
    subprocess.call("bwa mem " + args.reference_seq + " " + args.map_seq + " > " + args.map_seq + ".sam", shell=True)
    subprocess.call("samtools sort " + args.map_seq + ".sam > " + args.map_seq + ".bam", shell=True)
    subprocess.call("samtools index " + args.map_seq + ".bam", shell=True)
    print "Done"


parser = argparse.ArgumentParser(description='')
parser.add_argument('reference_seq', metavar = 'filename.FASTA', help='Reference sequence file to be mapped onto')
parser.add_argument('map_seq', metavar = 'filename.FASTA', help='Sequence file to be mapped on')
args = parser.parse_args()
main(args)
