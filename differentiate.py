import os
import subprocess

def main():

    if not os.path.exists("./Final_sequences.fasta.nin"):
        subprocess.call("makeblastdb -dbtype nucl -in Final_sequences.fasta", shell=True) # Contradicting with test.py!!
        # Make blast database with the list of transcripts


    subprocess.call("blastn -query Final_sequences.fasta -db Final_sequences.fasta -outfmt 6 > Blast_result_Final_sequences.txt", shell=True)
  # Run sequence on blast to database made above, and redirect output into separate file





    original_mRNA_transcripts = "dnc_sequences.FASTA.txt" # sequences is the argument to pass
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

            if (local_alignment[0] != local_alignment[1]):
                exon_counter = 1

        print transcript_infos


    group_dictionary = dict()
    groups = list()
    group_counter = 1
    with open("Blast_result_Final_sequences.txt", 'r') as ins:
        for line in ins:
            local_alignment = line.rstrip().split("\t") # Open the result file and split the tabs and make each line a vector
            # Redirect into a new file? Make the split file

            if (local_alignment[0] == local_alignment[1]):
                print "self alignment"
            elif (local_alignment[0] != local_alignment[1]) and (local_alignment[2] == '100.000'):
                print ("found alignment between " + local_alignment[0] + ' and ' + local_alignment[1])

                q_transc = list()
                s_transc = list()

                for group_num in groups :
                    # Checks if the query and subject transcript is already grouped
                    if (local_alignment[0] in group_dictionary[group_num]):
                        q_transc.append(group_num)
                        print (local_alignment[0] + " is in Group " + str(group_num))
                    if (local_alignment[1] in group_dictionary[group_num]):
                        s_transc.append(group_num)
                        print (local_alignment[1] + " is in Group " + str(group_num))


                if (len(q_transc) == 0) and (len(s_transc) == 0):
                    # Group the two transcripts into a new group if either one are not already grouped
                    print ("Both transcripts are not within any group")
                    group_dictionary[group_counter] = [local_alignment[0], local_alignment[1]]
                    groups.append(group_counter)
                    group_counter = group_counter + 1
                    print group_dictionary


                    # You need to solve the issue where one transcript is already added to a group so whether to add the non-grouped transcript to the original group or to make a new group
                    # I've made it so that if transcript A is grouped while transcript B isn't, then transcript B is grouped into all the groups that has transcript A
                    # Two problems here. 1) Order of the blast input file, 2) the group will get super long if the % identity parameter is too lenient.
                elif (len(q_transc) == 0) and (len(s_transc) != 0):
                    print (local_alignment[1] + ' is in Group(s) ' + str(s_transc))
                    print ('Adding ' + local_alignment[0] + ' into group(s)' + str(s_transc))
                    for num in s_transc :
                        group_dictionary[num].append(str(local_alignment[0]))
                    print group_dictionary

                elif (len(q_transc) != 0) and (len(s_transc) == 0):
                    print (local_alignment[0] + ' is in Group(s) ' + str(q_transc))
                    print ('Adding ' + local_alignment[1] + ' into group(s)' + str(q_transc))
                    for num in q_transc :
                        group_dictionary[num].append(str(local_alignment[1]))
                    print group_dictionary





# After blast, take each start and end coordinates from each transcripts and store them as exon coordinates

# dict[transcript_name][exon_num][start_coordinate, end_coordinate]


#
# After I've grouped all the potential transcripts together, I want to run MAFFT on each other to see the sequence alignment and determine whether or not they are isoforms of transcript from the same gene
#
main()
