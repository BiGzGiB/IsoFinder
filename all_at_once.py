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


# Things to fix and review
# line 529-530 and 666-669 about making decisions on adding chunk when there are more than one possible insertional location
# Update 8/2/17: I've added the arguement and tested it. It works! "Line 121 for how many intronic match it should have should also be a optional arguement"
# Line 875. Am I missing one transcript on the list of transcripts that have not been grouped with any other transcripts?
# Problem identified with Drosophila na transcript. It should group but it isn't grouping.
def main(args):

    if not os.path.exists("./gene_models/"):
        subprocess.call("mkdir gene_models", shell=True)
#    if not os.path.exists("./Final_sequences.fasta.nin"):
    subprocess.call("makeblastdb -dbtype nucl -in " + args.Genome_Walked_seq, shell=True) # Contradicting with test.py!!
        # Make blast database with the list of transcripts
    subprocess.call("blastn -query " + args.Genome_Walked_seq + " -db " + args.Genome_Walked_seq + " -outfmt 6 > gene_models/tmp_blast.txt", shell=True)
  # Output is to use for grouping isoforms together.

    original_mRNA_transcripts = args.Original_seq
    exon_counter = 1
    # Blast between the genome walked transcript and original transcript?

    subprocess.call("blastn -query " + original_mRNA_transcripts + " -db " + args.Genome_Walked_seq + " -outfmt 6 > gene_models/exon_coordinates.txt", shell=True)
    # For finding exon coordinates and calculating how many intronic sequences there are



    transcript_infos = dict() # This contains how many exons each transcripts contain and its coordinates
    with open("gene_models/exon_coordinates.txt", 'r') as ins:
        # 1. Open the sequences.FASTA file and split each sequences to separate files
        # Exon numbers aren't neccessarily in the correct order.
        for line in ins:
            local_alignment = line.rstrip().split("\t")

            if (local_alignment[0] == local_alignment[1]):
                if (local_alignment[0] not in transcript_infos):
                    transcript_infos[local_alignment[0]] = dict()  # This contains how many exons each transcripts contain and its coordinates
                transcript_infos[str(local_alignment[0])][exon_counter] = list((local_alignment[8], local_alignment[9]))
                exon_counter = 1 + exon_counter

            elif (local_alignment[0] != local_alignment[1]): # Basically Else:
                exon_counter = 1

        print ("Tranascript_infos: \n" + str(transcript_infos))



    group_dictionary = dict() # This contains the number of genes, its corresponding transcripts, and transcript's corresponding sequence.
    groups = list()
    group_counter = 1
    hundred_percent_coordinates = dict() # To store all self_matching/chunk coordinates
    # Chunk = exonic sequence + intronic sequence
    with open("gene_models/tmp_blast.txt", 'r') as ins:
        for line in ins:
            local_alignment = line.rstrip().split("\t") # Open the result file and split the tabs and make each line a vector

            if (local_alignment[0] == local_alignment[1]) and (float(local_alignment[2]) == 100.000) and (local_alignment[6] == local_alignment[8]) and (local_alignment[7] == local_alignment[9]):
                if (local_alignment[0] not in hundred_percent_coordinates):
                    hundred_percent_coordinates[local_alignment[0]] = dict()
                hundred_percent_coordinates[local_alignment[0]][int(local_alignment[6])] = (local_alignment[6], local_alignment[7], local_alignment[0])
                # Saving chunk coordinates for each transcripts
                # This information is used in the ordering of chunks on gene_model construction. Its not used in the process of grouping transcript isoforms


            # Get the exon coordinates for each transcripts
            elif (local_alignment[0] != local_alignment[1]) and (float(local_alignment[2]) >= args.p) and (local_alignment[3] >= 50):
                # (local_alignment[3] >= 50??)
                matching_exons = list()

                # WHY IS THE local_alignment[0] NOT IN transcript_infos??? Because the original transcript file did not contain the transcript data, OR genome walked file did not contain the genome_walked data
                # POSSIBLE ERROR!!!!
                if (local_alignment[0] not in transcript_infos):
                    print (str(local_alignment[0]) + " is not in transcript_infos.")
                    print ("Please be aware that your original sequence input file does not contain the mRNA sequence information of transcript " + str(local_alignment[0]))
                    # Error message
                    continue
                else :
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
                        # No matching exons between Transcript A and B so read

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
                            if ((len(matching_exons) == 1) and (intron_match >= args.i)) or (len(matching_exons) > 1):
                                # Group the transcripts together if the transcript shares one exon with at least 30 intronic sequence match, or there is more than 1 exon matching.
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
                        #else :
                        #    assert False, "Both transcripts are grouped. Why are these two transcripts being analysed again?"
            else :
                print "Not good enough match for grouping"


    if (len(group_dictionary) == 0):
        print "No isoforms. All transcripts are expressed from separate genes"
        subprocess.call("rm gene_models/tmp_blast.txt", shell=True)
        subprocess.call("rm gene_models/exon_coordinates.txt", shell=True)


    elif (len(group_dictionary) >= 1): # Basically Else:
        print group_dictionary
        print groups
        print "\n \n \n \n Starting gene model construction \n \n \n \n"


        for transcripts in hundred_percent_coordinates:
            hundred_percent_coordinates[transcripts] = collections.OrderedDict(sorted(hundred_percent_coordinates[transcripts].items()))
        # Order the self_matching chunks. Chunks are separated by 500Ns


        #
        # Each transcript have been grouped. Now time to stitch the information up to build a gene model for each group
        #


        # Save each transcript sequence in individual variables
        for group_num in group_dictionary:
            for record in SeqIO.parse(args.Genome_Walked_seq, "fasta"):
                if (record.name in group_dictionary[group_num]):
                    group_dictionary[group_num][record.name] = record.seq
#        group_dictionary now has group numbers, transcript_IDs that are within the groups, and its sequences


        Gene_model = dict() # This is like the insert_seq in the differentiate.py version for each gene_model
        for group_num in group_dictionary:
            Gene_model[group_num] = dict()

        # Start building the gene models by finding out which chunks overlap and find out how many chunks each gene_model have
        assertation_list = list() # To save chunks to loop through and check if there are more than one occurrences
        chunk_num = 1

        with open("gene_models/tmp_blast.txt", 'r') as ins:
            for line in ins:
                chunk_counter = str(chunk_num)
                local_alignment = line.rstrip().split("\t")
                query = local_alignment[0]
                subject = local_alignment[1]
                q_matching_chunk = tuple()
                s_matching_chunk = tuple()
                # Check if the q and s are in the same group
                for groups in group_dictionary:
                    if (query in group_dictionary[groups]) and (subject in group_dictionary[groups]):
                        # Both q and s in the same group. Now check through all chunks if it has already been added

                        # Check if the alignment is acceptable and save the coordinates if it is.
                        if (query != subject) and (float(local_alignment[2]) >= args.p):
                            for self_alignment in hundred_percent_coordinates[query]:
                                if (int(hundred_percent_coordinates[query][self_alignment][0]) <= int(local_alignment[6])) and (int(local_alignment[7]) <= int(hundred_percent_coordinates[query][self_alignment][1])):
                                    q_matching_chunk = hundred_percent_coordinates[query][self_alignment]
                            for self_alignment in hundred_percent_coordinates[subject]:
                                if (int(hundred_percent_coordinates[subject][self_alignment][0]) <= int(local_alignment[8])) and (int(local_alignment[9]) <= int(hundred_percent_coordinates[subject][self_alignment][1])):
                                    s_matching_chunk = hundred_percent_coordinates[subject][self_alignment]

                            s_identity = s_matching_chunk[2] + "_" + s_matching_chunk[0]
                            q_identity = q_matching_chunk[2] + "_" + q_matching_chunk[0]
                            not_in_any = True
                            for chunks in Gene_model[groups]:
                                assertation_list.append(s_identity) # To assert
                                assertation_list.append(q_identity) # To assert
                                # e.g. grh_L_1

                                if (q_identity in Gene_model[groups][chunks]) and (s_identity not in Gene_model[groups][chunks]):
                                    Gene_model[groups][chunks][s_identity] = s_matching_chunk
                                    not_in_any = False

                                elif (q_identity not in Gene_model[groups][chunks]) and (s_identity in Gene_model[groups][chunks]):
                                    Gene_model[groups][chunks][q_identity] = q_matching_chunk
                                    not_in_any = False

                                elif (q_identity in Gene_model[groups][chunks]) and (s_identity in Gene_model[groups][chunks]):
                                    not_in_any = False

                                elif (q_identity not in Gene_model[groups][chunks]) and (s_identity not in Gene_model[groups][chunks]):
                                    continue

                            if (not_in_any == True): # Both q_identity and s_identity are not in any other chunk. Make a new chunk with these in it
                                Gene_model[groups][chunk_counter] = dict()
                                Gene_model[groups][chunk_counter][q_identity] = q_matching_chunk
                                Gene_model[groups][chunk_counter][s_identity] = s_matching_chunk
                                chunk_num = chunk_num + 1

                        else : # Break the loop if q and s are not in the same group or % identity is too low
                            continue


            # Add the remaining chunks of transcripts that are not shared between other transcripts within the same group/gene_model
            lonely_chunks = dict() # To save those chunks that are on its own, to find possible unique chunks
            for Gene in Gene_model:
                lonely_chunks[Gene] = list()
                for transcript in hundred_percent_coordinates:
                    if (transcript in group_dictionary[Gene]):
                        for self_alignment in hundred_percent_coordinates[transcript]:
                            q_matching_chunk = hundred_percent_coordinates[transcript][self_alignment]
                            q_identity = q_matching_chunk[2] + "_" + q_matching_chunk[0]
                            stand_alone_chunk = True
                            for chunks in Gene_model[Gene]:
                                if (q_identity in Gene_model[Gene][chunks]):
                                    stand_alone_chunk = False
                                    break
                            if (stand_alone_chunk == True):
                                chunk_counter = str(chunk_num)
                                Gene_model[Gene][chunk_counter] = dict()
                                Gene_model[Gene][chunk_counter][q_identity] = q_matching_chunk
                                chunk_num = chunk_num + 1
                                lonely_chunks[Gene].append(chunk_counter)

            # REMOVING DUPLICATED CHUNKS
            duplicates = dict()
            print "Finding duplicate chunks"
            print duplicates
            for Gene in Gene_model:
                print ("Gene: " + str(Gene))
                for chunks in sorted(Gene_model[Gene], key=int):
                    print ("Chunk" + str(chunks) + ": " + str(sorted(Gene_model[Gene][chunks], key=str)))
                    for chunk2 in sorted(Gene_model[Gene], key=int):
                        if (sorted(Gene_model[Gene][chunks].items()) == sorted(Gene_model[Gene][chunk2].items())) and (chunks != chunk2):
                            if (len(duplicates) == 0):
                                duplicates[chunks + "_" + chunk2] = ((chunks, chunk2))
                            else:
                                if (str(chunk2 + "_" + chunks) not in duplicates):
                                    duplicates[chunks + "_" + chunk2] = ((chunks, chunk2))

            # REMOVING DUPLICATED CHUNKS
            print "Deleting duplicated chunks"
            for Gene in Gene_model:
                for dupes in duplicates:
                    delete = duplicates[dupes][1]
                    if (delete in Gene_model[Gene]):
                        print ('deleted: ' + str(sorted(Gene_model[Gene][delete], key=str)))
                        del Gene_model[Gene][delete]

            duplicates = dict()
            print "Finding duplicate chunks"
            for Gene in Gene_model:
                print ("Gene: " + str(Gene))
                for chunks in sorted(Gene_model[Gene], key=int):
                    for single_chunk in Gene_model[Gene][chunks]:
                        for chunk2 in sorted(Gene_model[Gene], key=int):
                            if (single_chunk in Gene_model[Gene][chunk2]):
                                if (len(duplicates) == 0):
                                    if (chunks != chunk2):
                                        duplicates[chunks + "_" + chunk2] = [chunks, chunk2]
                                else:
                                    if (chunks != chunk2):
                                        if (str(chunk2 + "_" + chunks) not in duplicates):
                                            duplicates[chunks + "_" + chunk2] = [chunks, chunk2]
                                break
            # MERGE DUPLICATED CHUNKS

            print "-- Merging duplicated chunks --"
            print duplicates
            for Gene in Gene_model:
                print ("Gene: " + str(Gene))
                for dupes in duplicates.values():
                    if (dupes[0] in Gene_model[Gene]) and (dupes[1] in Gene_model[Gene]):
                        chunk_counter = str(chunk_num)
                        print "BEFORE"
                        print sorted(Gene_model[Gene], key=int)
                        if (chunk_counter not in Gene_model[Gene]):
                            Gene_model[Gene][chunk_counter] = dict()
                            print chunk_counter
                        else :
                            assert False, "Merging Groups not properly merging!"

                        for single_chunk in Gene_model[Gene][dupes[0]].values():
                            c_identity = single_chunk[2] + "_" + single_chunk[0]
                            Gene_model[Gene][chunk_counter][c_identity] = single_chunk

                        for single_chunk in Gene_model[Gene][dupes[1]].values():
                            c_identity = single_chunk[2] + "_" + single_chunk[0]
                            Gene_model[Gene][chunk_counter][c_identity] = single_chunk

                        print ("Chunk" + str(dupes[0]) + ": " + str(sorted(Gene_model[Gene][dupes[0]], key=str)))
                        print ("Chunk" + str(dupes[1]) + ": " + str(sorted(Gene_model[Gene][dupes[1]], key=str)))
                        del Gene_model[Gene][dupes[0]]
                        del Gene_model[Gene][dupes[1]]
                        print ("Merged chunks " + str(dupes[0]) + " and " + str(dupes[1]))
                        print ("Chunk" + str(chunk_counter) + ": " + str(sorted(Gene_model[Gene][chunk_counter], key=str)))
                        print "AFTER"
                        print sorted(Gene_model[Gene], key=int)
                        chunk_num = chunk_num + 1

            # ASSERTATION: To see if there are more than one chunk that contains the same exon(?)
            assertation_list = list(set(sorted(assertation_list)))
            for chunk_assert in assertation_list:
                assert_counter = 0
#                print chunk_assert
                for Gene in Gene_model:
                    for Gene_model_chunk in Gene_model[Gene]:
                        if (chunk_assert in Gene_model[Gene][Gene_model_chunk]):
                            assert_counter = assert_counter + 1
#                            print assert_counter
                            assert assert_counter == 1, (str(chunk_assert) + " is in more than one chunks")



            # Start sorting them in order
            print "-- Start sorting them in order -- \n"
            Order = dict()
            for Gene in Gene_model:
                Order[Gene] = dict()
                for transcript in group_dictionary[Gene]:
                    Order[Gene][transcript] = list()
                    for self_alignment in hundred_percent_coordinates[transcript]:
                        q_matching_chunk = hundred_percent_coordinates[transcript][self_alignment]
                        q_identity = q_matching_chunk[2] + "_" + q_matching_chunk[0]

                        for chunks in Gene_model[Gene]:
                            if (q_identity in Gene_model[Gene][chunks]):
                                if (len(Order[Gene][transcript]) == 0):
                                    Order[Gene][transcript].append(chunks)
                                elif (len(Order[Gene][transcript]) != 0) and (chunks not in Order[Gene][transcript]):
                                    Order[Gene][transcript].append(chunks)

                                else: # Continue. This particular chunk is already added into the transcript's order
                                    continue
                            else: # This is normal. Obviously, q_identity might not be in chunks
                                continue


            # List for the ambiguous chunks that should be added AFTER the most left and right ends
            # Two conditions of unique chunks
            # 1. That chunk is only found in one transcript, AND that chunks on its either side are ALL lonely chunks
            # 2. The chunk is positionally on end of the gene model, BUT it was added within the gene model, not at the end, because of the other relative chunk position
            ambiguous = dict() # appending will be done from other_uniq and left right chunks that went in between gene_model
            for Gene in Gene_model:
                Order[Gene]["relative_order"] = dict()
                Order[Gene]["left_uniq"] = list()
                Order[Gene]["right_uniq"] = list()
                Order[Gene]["other_uniq"] = list()
                ambiguous[Gene] = list()

                # To make the relative order of each chunk information
                for chunks in Gene_model[Gene]:
                    if (chunks not in Order[Gene]["relative_order"]):
                        Order[Gene]["relative_order"][chunks] = dict()
                        Order[Gene]["relative_order"][chunks]['left'] = list()
                        Order[Gene]["relative_order"][chunks]['right'] = list()
                        for transcript in Order[Gene]:
                            if (transcript != "relative_order") and (transcript != "right_uniq") and (transcript != "left_uniq") and (transcript != "other_uniq"):
                                for chunk_order in Order[Gene][transcript]:
                                    if (chunks in Order[Gene][transcript]) and (chunk_order in Order[Gene][transcript]):

                                        q = Order[Gene][transcript].index(chunks)
                                        s = Order[Gene][transcript].index(chunk_order)

                                        if (s < q):
                                            if (chunk_order not in Order[Gene]["relative_order"][chunks]['left']):
                                                Order[Gene]["relative_order"][chunks]['left'].append(chunk_order)
                                        elif (s > q):
                                            if (chunk_order not in Order[Gene]["relative_order"][chunks]['right']):
                                                Order[Gene]["relative_order"][chunks]['right'].append(chunk_order)


                # Save possible uniques and the very left and right chunks
                for chunks in sorted(Order[Gene]["relative_order"], key=int):
                    if (len(Order[Gene]["relative_order"][chunks]['left']) == 0):
                        Order[Gene]['left_uniq'].append(chunks)

                    elif (len(Order[Gene]["relative_order"][chunks]['right']) == 0):
                        Order[Gene]['right_uniq'].append(chunks)

                # Remove non-ambiguous most left or right lonely_chunks from the lonely_chunks list so that they don't get identified as unique chunks
                if (len(Order[Gene]['left_uniq']) == 1):
                    lonely_left = Order[Gene]['left_uniq'][0] # STR
                    if (lonely_left in lonely_chunks[Gene]):
                        lonely_left_idx = lonely_chunks[Gene].index(lonely_left)
                        del lonely_chunks[Gene][lonely_left_idx]
                if (len(Order[Gene]['right_uniq']) == 1):
                    lonely_right = Order[Gene]['right_uniq'][0] # STR
                    if (lonely_right in lonely_chunks[Gene]):
                        lonely_right_idx = lonely_chunks[Gene].index(lonely_right)
                        del lonely_chunks[Gene][lonely_right_idx]

                # Save those chunks where it appears only in one transcripts and has no relative position on one side, but is still ambiguous

                for chunks in Order[Gene]['left_uniq']:
                    if (chunks in lonely_chunks[Gene]):
                        if (chunks not in ambiguous[Gene]):
                            ambiguous[Gene].append(chunks)
                for chunks in Order[Gene]['right_uniq']:
                    if (chunks in lonely_chunks[Gene]):
                        if (chunks not in ambiguous[Gene]):
                            ambiguous[Gene].append(chunks)

                # Save other uniques (All uniques excluding the uniques in the left and right)
                for chunks in lonely_chunks[Gene]:
                    if (chunks not in Order[Gene]['left_uniq']) and (chunks not in Order[Gene]['right_uniq']): # ??? Get rid of this?
                        left_chunks_lonely = True
                        right_chunks_lonely = True
                        for left_chunks in Order[Gene]["relative_order"][chunks]['left']:
                            if (left_chunks not in lonely_chunks[Gene]):
                                left_chunks_lonely = False
                        for right_chunks in Order[Gene]["relative_order"][chunks]['right']:
                            if (right_chunks not in lonely_chunks[Gene]):
                                right_chunks_lonely = False
                        if (left_chunks_lonely == True) or (right_chunks_lonely == True):
                            Order[Gene]["other_uniq"].append(chunks)
                            if (chunks not in ambiguous[Gene]):
                                ambiguous[Gene].append(chunks)

                        # ASSERT
                        if (left_chunks_lonely == True) and (right_chunks_lonely == True):
                            assert False, "Why is this transcript even added into this group??"


                # Start sorting the chunks in order
                print "\n -- Final Ordering for each Gene_Model -- \n"
                print "\n  - Ordering Inner Gene Model Order - \n"
                for chunks in sorted(Order[Gene]["relative_order"], key=int):
                    if ("final_order" not in Order[Gene]):
                        Order[Gene]['final_order'] = list()

                    if (chunks not in Order[Gene]['left_uniq']) and (chunks not in Order[Gene]['right_uniq']) and (chunks not in Order[Gene]['other_uniq']):
                        # Ordering non-ambiguous inner chunks first
                        if (len(Order[Gene]['final_order']) == 0):
                            Order[Gene]['final_order'].append(chunks)
                        elif (len(Order[Gene]['final_order']) == 1):
                            element = Order[Gene]['final_order'][0]
                            if (element in Order[Gene]["relative_order"][chunks]['left']):
                                Order[Gene]['final_order'].insert(1, chunks)
                            elif (element in Order[Gene]["relative_order"][chunks]['right']):
                                Order[Gene]['final_order'].insert(0, chunks)
                        else: # if the length is now above 1,
                            left_counter = 0 # Index to be inserted from left of list
                            right_counter = 0 # Index to be inserted from right of list
                            left_chunk = Order[Gene]["relative_order"][chunks]['left']
                            right_chunk = Order[Gene]["relative_order"][chunks]['right']

                            for chunks_final in Order[Gene]['final_order']:

                                for left_chunks in left_chunk:
                                    if (chunks_final in left_chunk):
                                        if (left_counter < (Order[Gene]['final_order'].index(chunks_final)+1)):
                                            left_counter = int(Order[Gene]['final_order'].index(chunks_final)) + 1
                                for right_chunks in right_chunk:
                                    if (chunks_final in right_chunk):
                                        if (right_counter == 0):
                                            right_counter = int(Order[Gene]['final_order'].index(chunks_final))
                                        elif (right_counter > Order[Gene]['final_order'].index(chunks_final)):
                                            right_counter = int(Order[Gene]['final_order'].index(chunks_final))
                                            print ("Right:" + str(right_counter))
                            print left_counter
                            print right_counter

                            # DOUBLE CHECK
                            print ("Chunk: " + str(chunks))
                            print Order[Gene]['relative_order'][chunks]
                            if (left_counter==right_counter):
                                Order[Gene]['final_order'].insert(left_counter, chunks)
                                print Order[Gene]['final_order']
                            elif (left_counter == 0) and (right_counter == 1):
                                Order[Gene]['final_order'].insert(left_counter, chunks)
                                print Order[Gene]['final_order']
                            elif (left_counter == len(Order[Gene]['final_order'])) and (right_counter == 0):
                                Order[Gene]['final_order'].insert(left_counter, chunks)
                                print Order[Gene]['final_order']
#                            elif (left_counter == 0) and (right_counter > 1):

                            else:
                                print "UNKNOWN WHERE TO INSERT CHUNK"
                                print Order[Gene]['relative_order'][chunks]

                                # There's more than one possible position this chunk can be added in.
                                # Make a decision!!!

#                                assert False, "Unknown where to insert chunk"
                print "\n  - Ordering Left and Right Most Chunks - \n"
                print Order[Gene]['final_order']
                print "LONELY CHUNKS"
                print lonely_chunks[Gene]
                print "OTHER UNIQUES"
                print Order[Gene]["other_uniq"]


                # Now add the left and right most chunks
                for chunks in sorted(Order[Gene]["relative_order"], key=int):
                    left_chunks = Order[Gene]["relative_order"][chunks]['left']
                    right_chunks = Order[Gene]["relative_order"][chunks]['right']

                    # Adding left most chunk and uniques
                    if (chunks in Order[Gene]['left_uniq']):
                        insert_index = len(Order[Gene]['final_order'])
                        for right_chunk in Order[Gene]['final_order']:

                            if (right_chunk in right_chunks):
                                if (insert_index > int(Order[Gene]['final_order'].index(right_chunk))):
                                    insert_index = int(Order[Gene]['final_order'].index(right_chunk))
                        Order[Gene]['final_order'].insert(insert_index, chunks)
                        if (insert_index != 0):
                            print "AMBIGUOUS ADDED"
                            print chunks
                            if (chunks not in ambiguous[Gene]):
                                ambiguous[Gene].append(chunks)

                    # Adding right most chunk and uniques
                    elif (chunks in Order[Gene]['right_uniq']):
                        insert_index = 0
                        for left_chunk in Order[Gene]['final_order']:

                            if (left_chunk in left_chunks):
                                if (insert_index < int(Order[Gene]['final_order'].index(left_chunk)) + 1):
                                    insert_index = int(Order[Gene]['final_order'].index(left_chunk)) + 1
                        if (insert_index != len(Order[Gene]['final_order'])):
                            print "AMBIGUOUS ADDED"
                            print chunks
                            if (chunks not in ambiguous[Gene]):
                                ambiguous[Gene].append(chunks)
                        Order[Gene]['final_order'].insert(insert_index, chunks)
                        print Order[Gene]['final_order']

                print ("Order of chunks in Gene_model" + str(Gene) + ': '+ str(Order[Gene]['final_order']))

                print "\n  - Ordering Remaining Ambiguous Chunks - \n"
                # Now add the other uniques in relation to the left and right most chunks
                for unique_chunks in Order[Gene]["other_uniq"]:

                    insert_index = 0
                    left_chunks_lonely = True
                    right_chunks_lonely = True
                    for left_chunks in Order[Gene]["relative_order"][unique_chunks]['left']:
                        if (left_chunks not in lonely_chunks[Gene]):
                            left_chunks_lonely = False # Left chunks are not lonely!
                    for right_chunks in Order[Gene]["relative_order"][unique_chunks]['right']:
                        if (right_chunks not in lonely_chunks[Gene]):
                            right_chunks_lonely = False # Right chunks are not lonely!

                    print unique_chunks
                    print left_chunks_lonely
                    print right_chunks_lonely
                    if (right_chunks_lonely == True) and (left_chunks_lonely == False): #
                        for right_chunk in Order[Gene]["relative_order"][unique_chunks]['right']:
                            if right_chunk in Order[Gene]['final_order']:
                                insert_idx_compare = Order[Gene]['final_order'].index(right_chunk)# INT.
                                if (insert_index < insert_idx_compare):
                                    insert_index = insert_idx_compare
                    elif (right_chunks_lonely == False) and (left_chunks_lonely == True):
                        for left_chunk in Order[Gene]["relative_order"][unique_chunks]['left']:
                            insert_idx_compare = Order[Gene]['final_order'].index(left_chunk) + 1 # INT. + 1 to get the correct insertion point
                            if (insert_index < insert_idx_compare):
                                insert_index = insert_idx_compare

                    else :
                        assert False, "Why is this unique chunk even in here?"
                    print insert_index
                    Order[Gene]['final_order'].insert(insert_index, unique_chunks)

                print Order[Gene]['final_order']





                # Identify chunks that have not yet been added because of its ambiguous position
                missing_chunks = list()
                for chunk in Gene_model[Gene]:
                    if (chunk not in Order[Gene]['final_order']):
                        missing_chunks.append(chunk)





                # Duplicated and modified version of Inner Chunk Ordering
                for chunks in missing_chunks:
                    left_counter = 0 # Index to be inserted from left of list
                    right_counter = 0 # Index to be inserted from right of list
                    left_chunk = Order[Gene]["relative_order"][chunks]['left']
                    right_chunk = Order[Gene]["relative_order"][chunks]['right']

                    for chunks_final in Order[Gene]['final_order']:
                        for left_chunks in left_chunk:
                            if (chunks_final in left_chunk):
                                if (left_counter < (Order[Gene]['final_order'].index(chunks_final)+1)):
                                    left_counter = int(Order[Gene]['final_order'].index(chunks_final)) + 1
                        for right_chunks in right_chunk:
                            if (chunks_final in right_chunk):
                                if (right_counter == 0):
                                    right_counter = int(Order[Gene]['final_order'].index(chunks_final))
                                elif (right_counter > Order[Gene]['final_order'].index(chunks_final)):
                                    right_counter = int(Order[Gene]['final_order'].index(chunks_final))
                                    print ("Right:" + str(right_counter))
                    print left_counter
                    print right_counter

                    # DOUBLE CHECK
                    print ("Chunk: " + str(chunks))
                    print Order[Gene]['relative_order'][chunks]
                    if (left_counter==right_counter):
                        Order[Gene]['final_order'].insert(left_counter, chunks)
                        print Order[Gene]['final_order']
                    elif (left_counter == 0) and (right_counter == 1):
                        Order[Gene]['final_order'].insert(left_counter, chunks)
                        print Order[Gene]['final_order']
                    elif (left_counter == len(Order[Gene]['final_order'])) and (right_counter == 0):
                        Order[Gene]['final_order'].insert(left_counter, chunks)
                        print Order[Gene]['final_order']
#                            elif (left_counter == 0) and (right_counter > 1):

                    else:
                        print "UNKNOWN WHERE TO INSERT CHUNK"
                        print Order[Gene]['relative_order'][chunks]
#                                assert False, "Unknown where to insert chunk"
                        #


                if (len(missing_chunks) >= 1):
                    print missing_chunks
                    print Order[Gene]['final_order']
                    assert False, (chunk + " is not in the final order")



# list.insert(index, value).
# list.index(value). However, keep in mind that the return value will be the smallest index that the value matches to.


            # Start saving the ordered sequence information into text file
            if os.path.exists("./gene_models/gene_models.fasta.txt"):
                subprocess.call("rm gene_models/gene_models.fasta.txt", shell=True) # Remove because I need to append.

            if os.path.exists("gene_models/Summary_info_human.txt"):
                subprocess.call("rm gene_models/Summary_info_human.txt", shell=True) # Remove because I need to append.

            if os.path.exists("gene_models/Summary_info_computer.txt"):
                subprocess.call("rm gene_models/Summary_info_computer.txt", shell=True) # Remove because I need to append.

            inserting_seq = list() # For Gene Model
            inserting_separate_seq = list() # For Sequences for separate chunks
            for Gene in Gene_model:
                for ordered_chunks in Order[Gene]['final_order']:
                    first_chunk_idx = Order[Gene]['final_order'].index(Order[Gene]['final_order'][0]) # INT
                    first_chunk_id = Order[Gene]['final_order'][0] # STR

                    last_idx = Order[Gene]['final_order'].index(Order[Gene]['final_order'][-1]) # INT
                    last_id = Order[Gene]['final_order'][-1] # STR

                    current_chunk_idx = Order[Gene]['final_order'].index(ordered_chunks)
                    current_chunk_id = Order[Gene]['final_order'][current_chunk_idx]

                    next_chunk_idx = Order[Gene]['final_order'].index(ordered_chunks) + 1 # INT
                #    next_chunk_id = Order[Gene]['final_order'][next_chunk_idx] # STR NOT NEEDED
                    print ordered_chunks

                #    if (Order[Gene]['final_order'].index(next_chunk_id) < Order[Gene]['final_order'].index(last_id)):
                #        next_chunk_id = Order[Gene]['final_order'][next_chunk_idx] # STR

                #    elif (Order[Gene]['final_order'].index(next_chunk_id)== last_idx):
                #        print "Next chunk is the last chunk"


                    # Add the chunk
                    if (len(Gene_model[Gene][ordered_chunks]) == 1):
                        for indv_chunks in Gene_model[Gene][ordered_chunks]:
                            t_id = Gene_model[Gene][ordered_chunks][indv_chunks]
                            inserting_seq.append(str(str(group_dictionary[Gene][t_id[2]])[int(t_id[0])-1:int(t_id[1])]))
                            inserting_separate_seq.append(str(str(group_dictionary[Gene][t_id[2]])[int(t_id[0])-1:int(t_id[1])]))
                            print ("Length of chunk" + str(ordered_chunks) + ": " + str(len(str(str(group_dictionary[Gene][t_id[2]])[int(t_id[0])-1:int(t_id[1])]))))


                            # If CURRENT chunk is the last chunk in the order, don't add any Ns or Us
                            if (ordered_chunks == last_id):
                                print "LAST CHUNK"

                            # If CURRENT chunk is a unique chunk, add 500 Us
                            elif ((ordered_chunks in Order[Gene]["left_uniq"]) or (ordered_chunks in Order[Gene]["right_uniq"]) or (ordered_chunks in Order[Gene]['other_uniq'])) and (0 < current_chunk_idx) and (current_chunk_idx < last_idx):
                                inserting_seq.append("U"*500)
                                print "Us ADDED"

                            # Elif NEXT chunk is a unique chunk, add 500 Us
                            elif (next_chunk_idx < len(Order[Gene]['final_order'])):
                                if ((Order[Gene]['final_order'][next_chunk_idx] in Order[Gene]["left_uniq"]) or (Order[Gene]['final_order'][next_chunk_idx] in Order[Gene]["right_uniq"]) or (Order[Gene]['final_order'][next_chunk_idx] in Order[Gene]['other_uniq'])) and (0 < next_chunk_idx) and (next_chunk_idx < last_idx):
                                    inserting_seq.append("U"*500)
                                    print "Us ADDED"
                                # Why is there an else?
                                else:
                                    inserting_seq.append("N"*500)
                                    print "Ns ADDED"

                            # Else, just add 500 Ns
                            else:
                                inserting_seq.append("N"*500)
                                print "Ns ADDED"

                    # Run MAFFT on chunks that align to the same region
                    elif (len(Gene_model[Gene][ordered_chunks]) > 1):
                        for indv_chunks in Gene_model[Gene][ordered_chunks]:
                            t_id = Gene_model[Gene][ordered_chunks][indv_chunks]
                            with open('gene_models/tmp_mafft.fasta.txt', 'a') as ins:
                                ins.write(">" + str(indv_chunks) + '\n' + str(str(group_dictionary[Gene][t_id[2]])[int(t_id[0])-1:int(t_id[1])]) + "\n") # -1 on start to not leave out the single sequence
                                ins.close()
                        # MAFFT
                        p1 = subprocess.Popen("mafft --quiet --auto gene_models/tmp_mafft.fasta.txt", shell=True, universal_newlines = True, stdout=subprocess.PIPE)
                        stdout, stderr = p1.communicate()
                        align = AlignIO.read(StringIO(stdout), "fasta")
                        summary_align = AlignInfo.SummaryInfo(align)
                        consensus = summary_align.dumb_consensus(ambiguous='N')
                        inserting_seq.append(str(consensus).upper())
                        inserting_separate_seq.append(str(consensus).upper())
                        print ("Length of chunk" + str(ordered_chunks) + ": " + str(len(str(consensus).upper())))

                        # If CURRENT chunk is the last chunk in the order, don't add any Ns or Us
                        if (ordered_chunks == last_id):
                            print "LAST CHUNK"

                        elif ((ordered_chunks in Order[Gene]["left_uniq"]) or (ordered_chunks in Order[Gene]["right_uniq"]) or (ordered_chunks in Order[Gene]['other_uniq'])) and ((Order[Gene]['final_order'].index(ordered_chunks) != 0) and (Order[Gene]['final_order'].index(ordered_chunks) != len(Order[Gene]['final_order'])-1)):
                            inserting_seq.append("U"*500)
                            print "Us ADDED"
                        # If NEXT chunk is a unique chunk, add 500 Us
                        elif (next_chunk_idx < len(Order[Gene]['final_order'])):
                            if ((Order[Gene]['final_order'][next_chunk_idx] in Order[Gene]["left_uniq"]) or (Order[Gene]['final_order'][next_chunk_idx] in Order[Gene]["right_uniq"]) or (Order[Gene]['final_order'][next_chunk_idx] in Order[Gene]['other_uniq'])) and (0 < next_chunk_idx) and (next_chunk_idx < last_idx):
                                inserting_seq.append("U"*500)
                                print "Us ADDED"
                            # Why is there an else?
                            else:
                                inserting_seq.append("N"*500)
                                print "Ns ADDED"
                        # Else, just add 500 Ns
                        else:
                            inserting_seq.append("N"*500)
                            print "Ns ADDED"


                        with open('gene_models/tmp_mafft.fasta.txt', 'w') as ins:
                            ins.close()
                            # To reset the file for the next round of MAFFT

                print "\n Writing to files now \n "

                # Write all gene models into one file
                with open('gene_models/gene_models.fasta.txt', 'a') as ins:
                    ins.write(">gene_model" + str(Gene) + "\n")
                    ins.write("".join(inserting_seq)+"\n \n")
                    print ("Length: " + str(len("".join(inserting_seq))))
                    print ("No. Chunks: " + str(len(inserting_seq)) + " \n \n \n \n")

                    ins.close()

                # Write Human readable summary file for all gene model
                with open('gene_models/Summary_info_human.txt', 'a') as ins:
                    ins.write("-- Gene Model " + str(Gene) + " -- \n")
                    ins.write('Transcript IDs: \n')
                    for transcript_id in sorted(group_dictionary[Gene], key=str):
                        if sorted(group_dictionary[Gene], key=str).index(transcript_id) +1 == len(group_dictionary[Gene]):
                            ins.write(transcript_id)
                        else:
                            ins.write(transcript_id + ', ')
                    ins.write('\n')
                    ins.write('Length of Gene Model (including gaps): ' + str(len("".join(inserting_seq))) + '\n')
                    ins.write('Length of Gene Model (excluding gaps): ' + str(len("".join(inserting_separate_seq))) + '\n')
                    ins.write('No. Chunks: ' + str(len(Order[Gene]['final_order'])) + ' \n')
                    ins.write('Order of chunks: ' + str(Order[Gene]['final_order']) + '\n')
                    ins.write('Ambiguous chunks: ')
                    for chunks in sorted(ambiguous[Gene], key=int):
                        if (sorted(ambiguous[Gene], key=int).index(chunks) + 1 == len(ambiguous[Gene])):
                            ins.write(chunks)
                        else:
                            ins.write(chunks + ', ')
                    ins.write('\n \n \n \n')


                # Write computer readable summary file
                with open('gene_models/Summary_info_computer.txt', 'a') as ins:
                    ins.write("Gene_Model" + str(Gene) + "\t")
                    for transcript_id in group_dictionary[Gene]:
                        ins.write(transcript_id + '\t')
                    ins.write('\n')


                if (args.separate == True):
                    # For IF(Optional args) the user want information separated for each gene models
                    if not os.path.exists("gene_models/gene_model" + str(Gene) + "_folder/"):
                        subprocess.call("mkdir gene_models/gene_model" + str(Gene) + "_folder", shell=True)
                    if os.path.exists("gene_models/gene_model" + str(Gene) + "_folder/final_all_seq.fasta.txt"):
                        subprocess.call('rm gene_models/gene_model' + str(Gene) + '_folder/final_all_seq.fasta.txt', shell=True)
                    if os.path.exists("gene_models/gene_model" + str(Gene) + "_folder/all_in_one_seq.fasta.txt"):
                        subprocess.call('rm gene_models/gene_model' + str(Gene) + '_folder/all_in_one_seq.fasta.txt', shell=True)

                    with open('gene_models/gene_model' + str(Gene) + '_folder/final_all_seq.fasta.txt', 'w') as ins:
                        ins.write("")
                        for idx in range(0, len(inserting_separate_seq)):
                            ins.write(">chunk" + str(Order[Gene]['final_order'][idx]) + "_" + str(Gene) + "\n")
                            ins.write(inserting_separate_seq[idx] + "\n \n")
                        ins.close()

                    with open('gene_models/gene_model' + str(Gene) + '_folder/all_in_one.fasta.txt', 'w') as ins:
                        ins.write(">gene_model" + str(Gene) + "\n")
                        ins.write("".join(inserting_seq)+"\n \n")
                        ins.close()




                # RESET FOR THE NEXT GENE MODEL!
                inserting_seq = list()
                inserting_separate_seq = list()


            # Add additional info in Summary file with the transcripts that have not been grouped at the end of the file
            with open('gene_models/Summary_info_human.txt', 'a') as ins:
                ins.write("-- Only variant transcripts -- \n")
                for transcript in sorted(hundred_percent_coordinates, key=str):
                    grouped = False
                    for Gene in group_dictionary:
                        if (transcript in group_dictionary[Gene]):
                            grouped = True
                            break
                    if (grouped == False):
                        # Does this mean one transcript is missing from the list??
                        if (sorted(hundred_percent_coordinates, key=str).index(transcript) + 1 == len(hundred_percent_coordinates)):
                            ins.write(transcript + '\n\n\n\n\n')
                        else:
                            ins.write(transcript + ', ')
                ins.close()

            with open('gene_models/Summary_info_computer.txt', 'a') as ins:
                ins.write("lonely_variants\t")
                for transcript in sorted(hundred_percent_coordinates, key=str):
                    grouped = False
                    for Gene in group_dictionary:
                        if (transcript in group_dictionary[Gene]):
                            grouped = True
                            break
                    if (grouped == False):
                        ins.write(transcript + '\t')
                ins.close()

            # Remove all temporary files that are no longer needed
            subprocess.call("rm gene_models/exon_coordinates.txt", shell=True)
            subprocess.call("rm gene_models/tmp_mafft.fasta.txt", shell=True)
            if (args.blast == False):
                subprocess.call("rm gene_models/tmp_blast.txt", shell=True)
            print "\n\n\n -- Finish -- \n\n\n"



parser = argparse.ArgumentParser(description='')
parser.add_argument('Genome_Walked_seq', metavar = 'filename.FASTA', help='Genome Walked sequence file')
parser.add_argument('Original_seq', metavar = 'filename.FASTA', help='Original sequence file before running Genome Walker')
parser.add_argument('-p', nargs='?', metavar='FLOAT', default=95.000, type=float, help='Minimum percentage identity for the grouping to occur based on the BLAST result (default: %(default)s)')
parser.add_argument('-i', nargs='?', metavar='INT', default=30, type=int, help='Minimum number of intronic sequence match for transcripts to be grouped together if transcripts only have one matching exon (default: %(default)s)')
parser.add_argument('-s','--separate', action='store_true', help='Creates a new folder that contains information about individual gene models (default: %(default)s)')
parser.add_argument('-b','--blast', action='store_true', help='Does not delete the BLAST output (default: %(default)s)')
args = parser.parse_args()
# Line 121 for how many intronic match it should have should also be a optional arguement
main(args)
