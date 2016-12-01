import os
import subprocess
import argparse
from Bio import SeqIO
from Bio import AlignIO
from Bio.Align import AlignInfo
from StringIO import StringIO
import re
from collections import OrderedDict
import collections


def main(seq):

#    if not os.path.exists("./Final_sequences.fasta.nin"):
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
        # Exon numbers aren't neccessarily in the correct order.
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


    #is_there_any_match = False
    group_dictionary = dict() #
    groups = list()
    group_counter = 1
    hundred_percent_coordinates = dict() # To store all exon coordinates
    with open("b_Final_sequences.txt", 'r') as ins:
        for line in ins:
            local_alignment = line.rstrip().split("\t") # Open the result file and split the tabs and make each line a vector
            # Redirect into a new file? Make the split file

            if (local_alignment[0] == local_alignment[1]) and (float(local_alignment[2]) == 100.000) and (local_alignment[6] == local_alignment[8]) and (local_alignment[7] == local_alignment[9]):
                if (local_alignment[0] not in hundred_percent_coordinates):
                    hundred_percent_coordinates[local_alignment[0]] = dict()
                    hundred_percent_coordinates[local_alignment[0]]['start'] = list()
                    hundred_percent_coordinates[local_alignment[0]]['end'] = list()
                hundred_percent_coordinates[local_alignment[0]]['start'].append(local_alignment[6])
                hundred_percent_coordinates[local_alignment[0]]['end'].append(local_alignment[7])
                hundred_percent_coordinates[local_alignment[0]]['start'].sort(key=int)
                hundred_percent_coordinates[local_alignment[0]]['end'].sort(key=int)

            # Get all the alignment infomation, sort the alignment coordinates, then calculate the 500Ns coordinate.


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
                            group_dictionary[group_counter] = dict()
                            group_dictionary[group_counter][local_alignment[0]] = dict()
                            group_dictionary[group_counter][local_alignment[1]] = dict()
                            groups.append(group_counter)
                            print (local_alignment[0] + " and " + local_alignment[1] + " grouped to Group" + str(group_counter))
                            group_counter = group_counter + 1
                            #is_there_any_match = True
                        elif (len(matching_exons) == 0):
                            print "Not grouped due to no matching exons or insufficient intron matches"



                                    # You need to solve the issue where one transcript is already added to a group so whether to add the non-grouped transcript to the original group or to make a new group
                                    # I've made it so that if transcript A is grouped while transcript B isn't, then transcript B is grouped into all the groups that has transcript A
                                    # Two problems here. 1) Order of the blast input file, 2)
                    elif (len(q_transc) == 0) and (len(s_transc) != 0):
                        print (local_alignment[1] + ' is in Group(s) ' + str(s_transc))
                        print ('Adding ' + local_alignment[0] + ' into group(s)' + str(s_transc))
                        for num in s_transc :
                            group_dictionary[num][local_alignment[0]] = dict()


                    elif (len(q_transc) != 0) and (len(s_transc) == 0):
                        print (local_alignment[0] + ' is in Group(s) ' + str(q_transc))
                        print ('Adding ' + local_alignment[1] + ' into group(s)' + str(q_transc))
                        for num in q_transc :
                            group_dictionary[num][local_alignment[1]] = dict()



    if (len(group_dictionary) == 0):
        print "No isoforms. All transcripts are expressed from separate genes"

    elif (len(group_dictionary) >= 1):
        print group_dictionary
        print groups
        print "\n \n \n \n Starting gene model construction \n \n \n \n \n \n \n \n \n \n "



        #
        # Each transcript have been grouped. Not time to stitch the information up to build a gene model for each group
        #

        # Save each transcript sequence in individual variables
        for group_num in group_dictionary:
            for record in SeqIO.parse("Final_sequences.fasta", "fasta"):
                if (record.name in group_dictionary[group_num]):
                    group_dictionary[group_num][record.name] = record.seq
        print group_dictionary



        # Name genes by group numbers
        # Initialise gene model by first taking the biggest transcript
        if not os.path.exists("./gene_models/"):
            subprocess.call("mkdir gene_models", shell=True) # Contradicting with test.py!!


        gene_model_dict = dict()
        largest_transcripts = dict()
        if (len(gene_model_dict) == 0):
            for gene_num in group_dictionary:
                if (gene_num not in largest_transcripts):
                    largest_transcripts[gene_num] = str()
                for transcript in group_dictionary[gene_num]:
                    if (largest_transcripts[gene_num] == ''):
                        print " YES IT IS EMPTY "
                        largest_transcripts[gene_num] = transcript
                        print ("Largest transcript for group " + str(gene_num) + " is " + str(largest_transcripts[gene_num]))
                        print len(group_dictionary[gene_num][str(largest_transcripts[gene_num])])

                    elif (len(group_dictionary[gene_num][str(largest_transcripts[gene_num])]) >= 1) and (len(group_dictionary[gene_num][transcript]) > len(group_dictionary[gene_num][str(largest_transcripts[gene_num])])):
                        largest_transcripts[gene_num] = transcript
                        print ("Largest transcript for group " + str(gene_num) + " is " + str(largest_transcripts[gene_num]))

                # Save the largest transcript sequence of each gene as the gene_model
                with open('gene_models/gene_model' + str(gene_num) + '.fasta.txt', 'w') as ins:
                    ins.write(">gene_model" + str(gene_num) + "\n")
                    ins.write(str(group_dictionary[gene_num][str(largest_transcripts[gene_num])])+"\n")
                    ins.close()

        for gene_num in group_dictionary:


#
#   Start individual gene_model stitching process
#
        #    count = 0
            Current_gene_model = ("gene_model" + str(gene_num))
            for strand in group_dictionary[gene_num]:
                First_chunk = True # To detect if any
                if (strand == largest_transcripts[gene_num]):
                    print "THIS IS THE LARGEST STRAND!!"
                    continue
                else :
                    subprocess.call("blastn -query gene_models/" + Current_gene_model + ".fasta.txt -subject gene_models/" + Current_gene_model + ".fasta.txt -perc_identity 100 -outfmt 6 > n_coordinates.fasta.txt", shell=True)
                    with open("n_coordinates.fasta.txt", 'r') as ins:
                        hundred_percent_coordinates[Current_gene_model] = dict()
                        hundred_percent_coordinates[Current_gene_model]['start'] = list()
                        hundred_percent_coordinates[Current_gene_model]['end'] = list()
                        for line in ins:
                            local_alignment = line.rstrip().split("\t") # Open the result file and split the tabs
                            if (local_alignment[0] == local_alignment[1]) and (float(local_alignment[2]) == 100.000) and (local_alignment[6] == local_alignment[8]) and (local_alignment[7] == local_alignment[9]):
                                hundred_percent_coordinates[local_alignment[0]]['start'].append(local_alignment[6])
                                hundred_percent_coordinates[local_alignment[0]]['end'].append(local_alignment[7])
                        hundred_percent_coordinates[Current_gene_model]['start'].sort(key=int)
                        hundred_percent_coordinates[Current_gene_model]['end'].sort(key=int)


                    for record in SeqIO.parse("gene_models/gene_model" + str(gene_num) + ".fasta.txt", "fasta"):
                        gene_model_dict[record.name] = record.seq


                    # Save the single transcript data into a separate file to run blast against gene_model
                    with open('tmp_compare_seq.fasta.txt', 'w') as ins:
                        ins.write(">" + strand + "\n")
                        ins.write(str(group_dictionary[gene_num][strand])+"\n")
                        ins.close()

                    print ("length of " + strand + " " + str(len(group_dictionary[gene_num][strand])))
                    subprocess.call("blastn -query gene_models/" + Current_gene_model + ".fasta.txt -subject tmp_compare_seq.fasta.txt -outfmt 6 -culling_limit 1 > tmp_blast.txt", shell=True)
#blastn -query gene_models/gene_model1.fasta.txt -subject tmp_compare_seq.fasta.txt -outfmt 7


                    with open("tmp_blast.txt", 'r') as ins:
                        chunk_counter = 1
                        matching_coordinates = dict()
                        current_transcript = str()
                        for line in ins:
                            local_alignment = line.rstrip().split("\t")
                            current_transcript = local_alignment[1]
                            if (float(local_alignment[2]) >= 95.000): # MANIPULATE THIS NUMBER!!!
                                matching_coordinates[int(local_alignment[6])] = list((local_alignment[6], local_alignment[7], local_alignment[8], local_alignment[9], chunk_counter)) #                        # Do I need [6] and [7] coordinates??
                                chunk_counter = chunk_counter + 1

                        #print hundred_percent_coordinates['var_L']
                        print matching_coordinates


                        matching_coordinates = collections.OrderedDict(sorted(matching_coordinates.items()))


                        print ("All Matching Coordinates \n" + str(matching_coordinates) + "\n")



                        # Identifying what parts to insert
                        insert_start_coordinates = list()
                        insert_end_coordinates = list()
                        base_insert_coordinate = int()
                        inserting_seq = list()
                        matched_chunk = str()
                        b_matched_already_added = dict() # To prevent blast result on the same exon causing a duplicate exon adding
                        b_matched_already_added['start'] = list()
                        b_matched_already_added['end'] = list()
                        c_matched_already_added = dict() # To prevent blast result on the same exon causing a duplicate exon adding
                        c_matched_already_added['start'] = list()
                        c_matched_already_added['end'] = list()
                        new_seq_info = True
                        chunks_covered = 0
                        Done = 0


                        for matching_self in range(0,len(hundred_percent_coordinates[current_transcript]['start'])):
                            for chunks in matching_coordinates:
                                if (int(matching_coordinates[chunks][2]) >= int(hundred_percent_coordinates[current_transcript]['start'][matching_self])) and (int(matching_coordinates[chunks][3]) <= int(hundred_percent_coordinates[current_transcript]['end'][matching_self])):
                                    print (str(int(hundred_percent_coordinates[current_transcript]['start'][matching_self])) + " <= " + str(int(matching_coordinates[chunks][2])) + "-" + str(int(matching_coordinates[chunks][3])) + " <= " + str(int(hundred_percent_coordinates[current_transcript]['end'][matching_self])))
                                    print "YES! Hit matched"
                                    if (First_chunk == True):
                                        print "This is the first match"
                                        First_chunk = False

                                        # Was there any new sequence information before this?
                                        print insert_start_coordinates
                                        print insert_end_coordinates

                                    if (int(matching_coordinates[chunks][0]) <= Done) and (int(matching_coordinates[chunks][1]) <= Done):
                                        print (str(matching_coordinates[chunks][0]) + " <= " + str(Done))
                                        print (str(matching_coordinates[chunks][1]) + " <= " + str(Done))
                                        print "False alarm. Already added"
                                    new_seq_info = False
                                    matched_chunk = int(matching_coordinates[chunks][0])
                                    a = hundred_percent_coordinates[Current_gene_model]['start']


                                    # Find the insert coordinate in master transcript
                                    for index in range(0, len(a)): # Hard coded
                #                        print matched_chunk
                                        if (int(a[index]) <= matched_chunk) and (matched_chunk <= int(hundred_percent_coordinates[Current_gene_model]['end'][index])):

                                            #assert base_insert_coordinate <= int(a[index])
                                            base_insert_coordinate = int(a[index]) - 1 # -1 as well. DOUBLE CHECK THIS!!!
                                            print (str(a[index]) + " <= " + str(matched_chunk) + " <= " + str(hundred_percent_coordinates[Current_gene_model]['end'][index]))
                                            print ("Base insert coordinate: " + str(base_insert_coordinate))

                                            break


                                    # ADD LEFT
                                    inserting_seq.append(str(gene_model_dict[Current_gene_model])[Done:base_insert_coordinate]) # !! Check if I am adding a single base in between !!
                                    print ("ADDING " + str(len(str(gene_model_dict[Current_gene_model])[Done:base_insert_coordinate])) + " bases")



                                    print "Added Left"



                                    # ADD NEW
                                    for index in range(0, len(insert_start_coordinates)):
                                        insert_start = int(insert_start_coordinates[index]) - 1  # The -1 is needed to not miss the first base
                                        insert_end = int(insert_end_coordinates[index])

                                        inserting_seq.append(str(group_dictionary[gene_num][current_transcript])[insert_start:insert_end])
                                        inserting_seq.append("N"*500)
                                        print ("ADDING " + str(len(str(group_dictionary[gene_num][current_transcript])[insert_start:insert_end])) + " + 500N bases")
                                        print ("Length of gene_model: " + str(len("".join(inserting_seq))))
                                        print "Added New"

                                    # ADD MATCHED
                                    # Calculate the size of base exon
                                    print hundred_percent_coordinates[Current_gene_model]
                                    print ("matching coordinates[chunks]: " + str(matching_coordinates[chunks]))
                                    end_of_match = int(matching_coordinates[chunks][1])
                                    print ("BEFORE end_of_match: " + str(end_of_match))
                                    print (str(int(hundred_percent_coordinates[Current_gene_model]['start'][index])) + " <= " + str(end_of_match) + " <= " + str(int(hundred_percent_coordinates[Current_gene_model]['end'][index])))


                                    if (int(hundred_percent_coordinates[Current_gene_model]['start'][index]) <= end_of_match) and (end_of_match <= int(hundred_percent_coordinates[Current_gene_model]['end'][index])):
                                        end_of_match = int(hundred_percent_coordinates[Current_gene_model]['end'][index])
                                        print ("AFTER end_of_match: " + str(end_of_match))

                                    # Calculate the size of compare exon
                                    comp_ex_start_site = int()
                                    comp_ex_end_site = int()
                                    base_ex_start_site = int()
                                    base_ex_end_site = int()
                                    for coord in range(0, len(hundred_percent_coordinates[Current_gene_model]['start'])):
                                        base_ex_start_site = int(hundred_percent_coordinates[Current_gene_model]['start'][coord])
                                        base_ex_end_site = int(hundred_percent_coordinates[Current_gene_model]['end'][coord])

                                        #print "OVER HERE!!!"
                                        #print (str(base_ex_start_site) + ' <= ' + str(matching_coordinates[chunks][0])) + ' - ' + str((matching_coordinates[chunks][1]) + ' <= ' + str(base_ex_end_site))
                                        if (base_ex_start_site <= int(matching_coordinates[chunks][0])) and (int(matching_coordinates[chunks][1]) <= base_ex_end_site):
                                            print ("Length of Gene Model Exon: " + str(len(str(gene_model_dict[Current_gene_model])[base_ex_start_site:base_ex_end_site])))
                                            break
                                    for coord in range(0, len(hundred_percent_coordinates[strand]['start'])):
                                        comp_ex_start_site = int(hundred_percent_coordinates[strand]['start'][coord])
                                        comp_ex_end_site = int(hundred_percent_coordinates[strand]['end'][coord])
                                        #print "OVER HERE!!!22222"

                                        #print (str(comp_ex_start_site) + ' <= ' + str(matching_coordinates[chunks][2])) + ' - ' + str((matching_coordinates[chunks][3]) + ' <= ' + str(comp_ex_end_site))
                                        if (comp_ex_start_site <= int(matching_coordinates[chunks][2])) and (int(matching_coordinates[chunks][3]) <= comp_ex_end_site):
                                            print ("Length of " + str(strand) + " Exon: " + str(len(str(group_dictionary[gene_num][strand])[comp_ex_start_site:comp_ex_end_site])))
                                            break


                                    # Check to see if there is a difference in length between base[exon] and compare[exon]
                                    # Two options
                                    # IF the base matched is larger or same, then add that
                                    # ELIF the compare matched is larger, Run MAFFT and add consensus sequence
                                    print comp_ex_start_site
                                    print comp_ex_end_site
                                    if ((comp_ex_start_site in c_matched_already_added['start']) and (comp_ex_end_site in c_matched_already_added['end'])) or ((base_ex_start_site in c_matched_already_added['start']) and (base_ex_end_site in c_matched_already_added['end'])):
                                        print "Already added. Move on"
                                        print "Second Condition Met"
                                        print c_matched_already_added
                                        Done = end_of_match + 500

                                    elif (len(str(gene_model_dict[Current_gene_model])[base_ex_start_site:base_ex_end_site]) != len(str(group_dictionary[gene_num][strand])[comp_ex_start_site:comp_ex_end_site])):
                                    #    print (str(gene_model_dict[Current_gene_model])[base_insert_coordinate:end_of_match] + " \n " + str(group_dictionary[gene_num][strand])[comp_ex_start_site:comp_ex_end_site])
                                        c_matched_already_added['start'].append(comp_ex_start_site)
                                        c_matched_already_added['end'].append(comp_ex_end_site)
                                        c_matched_already_added['start'].append(base_ex_start_site)
                                        c_matched_already_added['end'].append(base_ex_end_site)
                                        print "First Condition Met"


                                        with open('gene_models/tmp_mafft.fasta.txt', 'w') as ins:
                                            ins.write(">Gene_model_exon" + "\n" + str(str(gene_model_dict[Current_gene_model])[base_insert_coordinate:end_of_match]) + "\n")
                                            ins.write(">" + str(strand) + "_exon \n" + str((str(group_dictionary[gene_num][strand])[comp_ex_start_site:comp_ex_end_site])) + "\n")
                                            ins.close()


                                        # MAFFT
                                        p1 = subprocess.Popen("mafft --quiet --auto gene_models/tmp_mafft.fasta.txt", shell=True, universal_newlines = True, stdout=subprocess.PIPE)
                                        stdout, stderr = p1.communicate()
                                        align = AlignIO.read(StringIO(stdout), "fasta")
#                                        print(align)
                                        summary_align = AlignInfo.SummaryInfo(align)
                                        consensus = summary_align.dumb_consensus(ambiguous='N')
                                        print "MAFFT CONSENSUS!!!"
#                                        print str(consensus).upper()
                                        inserting_seq.append(str(consensus).upper())
                                        inserting_seq.append("N"*500)
                                        print ("ADDING " + str(len(str(consensus).upper())) + " +500"  + " bases")
                                        print "Added Matched (MAFFT)"
                                        print c_matched_already_added
                                        Done = end_of_match + 500

#




                                    # Only add MATCHED if not already added.
                                    elif (end_of_match not in b_matched_already_added['end']) and (base_insert_coordinate not in b_matched_already_added['start']):
                                        inserting_seq.append(str(gene_model_dict[Current_gene_model])[base_insert_coordinate:end_of_match]) # THIS LINE OF CODE IS MESSING THIS UP!!
                                        inserting_seq.append("N"*500)
                                        b_matched_already_added['end'].append(end_of_match)
                                        b_matched_already_added['start'].append(base_insert_coordinate)
                                        print "Third Condition Met"
                                        print b_matched_already_added
                                        print ("ADDING " + str(len(str(gene_model_dict[Current_gene_model])[base_insert_coordinate:end_of_match])) + " +500"  + " bases")
                                        print "Added Matched"

                                        Done = end_of_match + 500 # Added 500 Ns
                                    else :
                                        print "Fourth Condition Met"
                                        print ("Length of gene_model: " + str(len("".join(inserting_seq))))
                                        print "Matched exon already added. Skipping."
                                    # GET DONE COORDINATE







                                    # WIPE EVERY LIST AND UPDATE
                                    insert_start_coordinates = list()
                                    insert_end_coordinates = list()
                                    chunks_covered = chunks_covered + 1

                                    print (" \n \n \n \n \nBase insert coordinates: "+ str(base_insert_coordinate))
                                    print ("End of match: "+ str(end_of_match))
                                    #print ("Matching coordinates: " + str(matching_coordinates))
                                    print ("Chunks covered: " + str(chunks_covered))
                                    #print ("Length of matching_coordinates: " + str(len(matching_coordinates)))
                                    print ("Length of gene_model: " + str(len("".join(inserting_seq))))
                                    print ("Done: "+ str(Done) + " \n \n \n \n")


                                if (chunks_covered == len(matching_coordinates)):
                                    print "No more new sequence to add. Adding remaining sequences"

                                    # ADD REMAINING SEQUENCE OF GENE_MODEL
                                    inserting_seq.append(str(gene_model_dict[Current_gene_model])[Done:])
                                    chunks_covered = 0
                                    print ("Length of gene_model: " + str(len("".join(inserting_seq))))
                                    print ("Length of gene_model: " + str(len(inserting_seq)))
                                    print "Finish \n \n \n \n \n \n \n \n \n \n \n \n \n \n \n"

                                    break


                            if (new_seq_info == True):
                                insert_start_coordinates.append(int(hundred_percent_coordinates[current_transcript]['start'][matching_self]))
                                insert_end_coordinates.append(int(hundred_percent_coordinates[current_transcript]['end'][matching_self]))
                                print "Found new sequence information!"
                #            else:
                #                if (len(inserting_seq) > 0): #??
                #                    continue

                        # WRITE TO A NEW FILE
                #        ite = ite + 1
                        with open('gene_models/' + Current_gene_model + '.fasta.txt', 'w') as ins:
                            ins.write(">" + Current_gene_model + "\n")
                #            ins.write(">gene_model1.iter" + str(ite) +"\n")
                            ins.write("".join(inserting_seq)+"\n")
                            ins.close()

                        # Reset after writing gene_model to a file
                        inserting_seq = list()
                #    count = count + 1
                #    if (count == 3):
                #        break
        #print group_dictionary




parser = argparse.ArgumentParser(description='')
parser.add_argument('seq', metavar = 'filename.FASTA', help='Original sequence file before running Genome Walker')
args = parser.parse_args()

main(args.seq)

# After blast, take each start and end coordinates from each transcripts and store them as exon coordinates

# dict[transcript_name][exon_num][start_coordinate, end_coordinate]




#
# After I've grouped all the potential transcripts together, I want to run MAFFT on each other to see the sequence alignment and determine whether or not they are isoforms of transcript from the same gene
#
