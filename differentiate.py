import os
import subprocess
import argparse

def main(seq):

    if not os.path.exists("./Final_sequences.fasta.nin"):
        subprocess.call("makeblastdb -dbtype nucl -in Final_sequences.fasta", shell=True) # Contradicting with test.py!!
        # Make blast database with the list of transcripts
    subprocess.call("blastn -query Final_sequences.fasta -db Final_sequences.fasta -outfmt 6 > b_Final_sequences.txt", shell=True)
  # Run sequence on blast to database made above, and redirect output into separate file

    original_mRNA_transcripts = seq
    transcript_infos = dict()
    exon_counter = 1
    # Blast between the genome walked transcript and original transcript?

    subprocess.call("blastn -query " + original_mRNA_transcripts + " -db Final_sequences.fasta -outfmt 6 > " + original_mRNA_transcripts + "_exon_coordinates.txt", shell=True)

    with open(original_mRNA_transcripts + "_exon_coordinates.txt", 'r') as ins:
        # 1. Open the sequences.FASTA file and split each sequences to separate files
        # Sort the file content first by alphabetically, then numerically on s.start(which is local_alignment[8])
        for line in ins:
            local_alignment = line.rstrip().split("\t")

            if (local_alignment[0] == local_alignment[1]):
                if (local_alignment[0] not in transcript_infos):
                    transcript_infos[local_alignment[0]] = dict()

                transcript_infos[str(local_alignment[0])][exon_counter] = list((local_alignment[8], local_alignment[9]))
                exon_counter = 1 + exon_counter

            elif (local_alignment[0] != local_alignment[1]):
                exon_counter = 1

        print transcript_infos



    group_dictionary = dict()
    groups = list()
    group_counter = 1
    with open("b_Final_sequences.txt", 'r') as ins:
        for line in ins:
            local_alignment = line.rstrip().split("\t") # Open the result file and split the tabs and make each line a vector
            # Redirect into a new file? Make the split file

            if (local_alignment[0] == local_alignment[1]):
                # print "self alignment"
                continue
            elif (local_alignment[0] != local_alignment[1]) and (float(local_alignment[2]) >= 95.000):
                #print ("found alignment between " + local_alignment[0] + ' and ' + local_alignment[1])


                matching_exons = list()
                tscript_exon_coord = transcript_infos[local_alignment[0]]


                for exons in transcript_infos[local_alignment[0]]:
                    intron_match = (int(local_alignment[7]) - int(tscript_exon_coord[exons][1])) + (int(tscript_exon_coord[exons][0]) - int(local_alignment[6]))
                    if (int(local_alignment[6]) > int(tscript_exon_coord[exons][0])) and (int(local_alignment[7]) < int(tscript_exon_coord[exons][1])):
                        continue
                        # Ignore this alignment because its shorter than the exon. Useless
                    elif (int(local_alignment[6]) <= int(tscript_exon_coord[exons][0])) and (int(local_alignment[7]) >= int(tscript_exon_coord[exons][1])):
                        # Matching exon and possibly some introns
                        #print (local_alignment[6] + "  " + local_alignment[7] + "\n" + transcript_infos[local_alignment[0]][exons][0] + "  " + transcript_infos[local_alignment[0]][exons][1])
                        # The above is just to compare between exon coordinates and alignment coordinate
                        matching_exons.append(str(exons))
                        print ("Exon" + str(exons) + " has " + str(intron_match) + " extra intronic matches")
                if (len(matching_exons) == 0):
                    print "no matching exons"
                elif (len(matching_exons) >= 1):
                    print ("Transcript " + local_alignment[1] + " matches " + str(matching_exons) + " of " + local_alignment[0] + "'s exon(s)")
                    q_transc = list()
                    s_transc = list()

                    for group_num in groups :
                        # Checks if the query and subject transcript is already grouped
                        if (local_alignment[0] in group_dictionary[group_num]):
                            q_transc.append(group_num)
                                #print (local_alignment[0] + " is in Group " + str(group_num))
                        if (local_alignment[1] in group_dictionary[group_num]):
                            s_transc.append(group_num)
                                #print (local_alignment[1] + " is in Group " + str(group_num))


                    if (len(q_transc) == 0) and (len(s_transc) == 0):
                        # Group the two transcripts into a new group if either one are not already grouped
                        print ("Both transcripts are not within any group")
                            # Do they have good intronic matches?
                        if ((len(matching_exons) == 1) and (intron_match >= 30)) or (len(matching_exons) > 1):
                            group_dictionary[group_counter] = [local_alignment[0], local_alignment[1]]
                            groups.append(group_counter)
                            print (local_alignment[0] + " and " + local_alignment[1] + " grouped to Group" + str(group_counter))
                            group_counter = group_counter + 1
                        elif (len(matching_exons) == 0):
                            print "Not grouped due to no matching exons or insufficient intron matches"



                                    # You need to solve the issue where one transcript is already added to a group so whether to add the non-grouped transcript to the original group or to make a new group
                                    # I've made it so that if transcript A is grouped while transcript B isn't, then transcript B is grouped into all the groups that has transcript A
                                    # Two problems here. 1) Order of the blast input file, 2)
                    elif (len(q_transc) == 0) and (len(s_transc) != 0):
                        print (local_alignment[1] + ' is in Group(s) ' + str(s_transc))
                        print ('Adding ' + local_alignment[0] + ' into group(s)' + str(s_transc))
                        for num in s_transc :
                            group_dictionary[num].append(str(local_alignment[0]))


                    elif (len(q_transc) != 0) and (len(s_transc) == 0):
                        print (local_alignment[0] + ' is in Group(s) ' + str(q_transc))
                        print ('Adding ' + local_alignment[1] + ' into group(s)' + str(q_transc))
                        for num in q_transc :
                            group_dictionary[num].append(str(local_alignment[1]))


    print group_dictionary


parser = argparse.ArgumentParser(description='')
parser.add_argument('seq', metavar = 'filename.FASTA', help='Original sequence file before running Genome Walker')
args = parser.parse_args()

main(args.seq)

# After blast, take each start and end coordinates from each transcripts and store them as exon coordinates

# dict[transcript_name][exon_num][start_coordinate, end_coordinate]


#
# After I've grouped all the potential transcripts together, I want to run MAFFT on each other to see the sequence alignment and determine whether or not they are isoforms of transcript from the same gene
#
