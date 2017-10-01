import os, sys, subprocess, argparse, collections, logging
from Bio import SeqIO
from Bio import AlignIO
from Bio.Align import AlignInfo
from StringIO import StringIO


# Should learn to use list/dict comprehension to make the script more readable

def Obtain_exon_indices(BLAST_file):
    # Obtain the indices for exon sequences within each transcript.
    # BLAST sequence alignment output between mRNA transcript and Genome Walked transcript required.

    exon_indices = dict()
    exon_counter = 1
    try:
        ins = open(BLAST_file, 'r')
    except IOError:
        print("%s does not exist" % (BLAST_file))

    for line in ins:
        local_alignment = line.rstrip().split("\t")

        if (local_alignment[0] == local_alignment[1]):
            if (local_alignment[0] not in exon_indices):
                exon_indices[local_alignment[0]] = dict()
            exon_indices[str(local_alignment[0])][exon_counter] = list((int(local_alignment[8]), int(local_alignment[9]))) # Start and end index of the exon
            exon_counter = 1 + exon_counter

        elif (local_alignment[0] != local_alignment[1]):
            exon_counter = 1
    # Note: Each exons are assigned an integer, but they do not represent the position of the exon within the transcript (e.g. Exon_1 does not neccessarily mean that that exon is on the 5' or 3' most position)
    # print ("Tranascript_infos: \n" + str(exon_indices))
    return exon_indices

def Add_to_group(group_dictionary, transcript_0, transcript_1, group_counter):
    # This function is inside the Cluster_transcripts function. It is to group the two transcripts together that has been predicted to be transcribed from the same gene.
    q_transc = list()
    s_transc = list()

    for group_num in group_dictionary:
        # Checks if the query and subject transcript is already grouped
        if (transcript_0 in group_dictionary[group_num]):
            q_transc.append(group_num)
                #print (local_alignment[0] + " is in Group " + str(group_num))
        if (transcript_1 in group_dictionary[group_num]):
            s_transc.append(group_num)
                #print (local_alignment[1] + " is in Group " + str(group_num))


    if (len(q_transc) == 0) and (len(s_transc) == 0):
        # Group the two transcripts into a new group if either one are not already grouped
        #print ("Both transcripts are not within any group")

        group_dictionary[group_counter] = dict()
        group_dictionary[group_counter][transcript_0] = dict()
        group_dictionary[group_counter][transcript_1] = dict()
        print("%s and %s grouped to Group_%d" % (transcript_0, transcript_1, group_counter))
        group_counter = group_counter + 1

                           # elif (len(matching_exons) == 0):
                               # print "Not grouped due to no matching exon or insufficient intron matches"
    elif (len(q_transc) == 0) and (len(s_transc) != 0):
        print ('%s is in Group(s) %s' % (transcript_1, str(s_transc)))
        print ('Adding %s into group(s) %s' % (transcript_0, str(s_transc)))
        for num in s_transc :
            group_dictionary[num][transcript_0] = dict()


    elif (len(q_transc) != 0) and (len(s_transc) == 0):
        print ('%s is in Group(s) %s' % (transcript_0, str(q_transc)))
        print ('Adding %s into group(s) %s' % (transcript_1, str(q_transc)))
        for num in q_transc :
            group_dictionary[num][transcript_1] = dict()
    else : # (len(q_transc) != 0) and (len(s_transc) != 0)
        tmp_groups = list()
        tmp_transcripts = list()
        for group in q_transc:
            tmp_groups.append(group)
            for transcript in group_dictionary[group]:
                tmp_transcripts.append(transcript)
        for group in s_transc:
            tmp_groups.append(group)
            for transcript in group_dictionary[group]:
                tmp_transcripts.append(transcript)

        print("Merging groups %s" % (tmp_groups))
        logging.info("Merged groups %s" % (tmp_groups))
        del tmp_groups[0]

        for transcript in tmp_transcripts:
            group_dictionary[q_transc[0]][transcript] = dict()

        for group in tmp_groups:
            del group_dictionary[group]

    return group_dictionary, group_counter

def Group_transcripts(exon_indices, BLAST_file, min_perc_identity=95.000, min_intron_match=30):
    # Function to measure whether the two transcripts should be grouped together.
    # Transcripts which are grouped together are transcripts predicted to be transcribed from the same gene.

    group_dictionary = dict() # This contains the number of genes, its corresponding transcripts, and transcript's corresponding sequence.
    group_counter = 1
    try:
        ins = open(BLAST_file, 'r')
    except IOError:
        print("%s does not exist" % (BLAST_file))

    for line in ins:
        local_alignment = line.rstrip().split("\t")
        # Get the exon coordinates for each transcripts
        already_grouped = False
        for num in group_dictionary:
            if (local_alignment[0] in group_dictionary[num]) and (local_alignment[1] in group_dictionary[num]):
                already_grouped = True
            else:
                continue

        if (already_grouped == False) and (local_alignment[0] != local_alignment[1]) and (float(local_alignment[2]) >= min_perc_identity) and (local_alignment[3] >= 50):
            # local_alignment[3] = alignment length

            matching_exons = list()
            if (local_alignment[0] not in exon_indices):
                # Q: Why would local_alignment[0] not be in exon_indices?
                # A: Because the original transcript file does not contain the transcript data, OR genome walked file did not contain the genome_walked sequence data.
                logging.warning("Please be aware that your mRNA sequence file does not contain the mRNA sequence information of transcript %s" % (local_alignment[0]))
                # Error log message

            else :
                exons = exon_indices[local_alignment[0]]
                for exon in exons:
                    intron_match = (int(local_alignment[7]) - exons[exon][1]) + (exons[exon][0] - int(local_alignment[6]))
                    if (int(local_alignment[6]) > exons[exon][0]) and (int(local_alignment[7]) < exons[exon][1]): # Make this exon >= local and local >= exon?
                        # Ignore this alignment because its shorter than the exon. Useless
                        continue

                    elif (int(local_alignment[6]) <= exons[exon][0]) and (int(local_alignment[7]) >= exons[exon][1]): # Make this exon >= local and local >= exon?
                        # Matching exon and possibly some introns
                        # print (local_alignment[6] + "  " + local_alignment[7] + "\n" + transcript_infos[local_alignment[0]][exon][0] + "  " + transcript_infos[local_alignment[0]][exon][1])
                        # The above is just to compare between exon coordinates and alignment coordinate
                        matching_exons.append(str(exon))
                        print ("Exon%s has %s extra intronic matches" % (str(exon), str(intron_match)))
                    else :
                        # This is in the case of where only one end of the contig matches both intron and exons while the other end only matches exons
                        #
                        # exon = =; intron = -
                        # transc_1_contig_1 ----=======----
                        # transc_2_contig_1 ----=====
                        #                   or
                        # transc_1_contig_1 ----=======----
                        # transc_2_contig_1        ====----
                        # Will make an optional argument whether to include these matches for clustering later.
                        pass

                # Do they have good intronic matches?
                # Group the transcripts together if the transcript shares one exon with at least the user-defined (Default=30) number of intronic sequence match, or there is more than 1 exon matching.
                if (len(matching_exons) == 0):
                    print "no matching exons"
                    # No matching exons between Transcript A and B

                elif (len(matching_exons) == 1) and (intron_match >= min_intron_match):
                    logging.info("Condition (len(matching_exons) == 1) and (intron_match >= min_intron_match) met between %s and %s" % (local_alignment[0], local_alignment[1]))
                    logging.info("intron_match: %d" % (intron_match))
                    group_dictionary, group_counter = Add_to_group(group_dictionary, local_alignment[0], local_alignment[1], group_counter)

                elif (len(matching_exons) > 1):
                #elif ((len(matching_exons) == 1) and (intron_match >= min_intron_match)) or (len(matching_exons) > 1):
                    print ("Transcript %s matches %s of %s's exon(s)" % (local_alignment[1], str(matching_exons), local_alignment[0]))
                    group_dictionary, group_counter = Add_to_group(group_dictionary, local_alignment[0], local_alignment[1], group_counter)

                else: # len(matching_exons) == 1:
                    print ("Not good enough alignment for grouping between %s and %s" % (local_alignment[0], local_alignment[1]))

        elif already_grouped == True:
            #print ("%s and %s already clustered" % (local_alignment[0], local_alignment[1]))
            continue
        else :
            print ("Not good enough alignment for grouping between %s and %s" % (local_alignment[0], local_alignment[1]))

    return group_dictionary

def Obtain_contig_indices(BLAST_file):
    # To distinguish each contigs and sequence gaps
    # BLAST_file == BLAST between GW_transcript sequence to GW_transcript sequence
    hundred_percent_coordinates = dict() # Will contain sequence information of contigs
    with open(BLAST_file, 'r') as ins:
        for line in ins:
            local_alignment = line.rstrip().split("\t")

            if (local_alignment[0] == local_alignment[1]) and (float(local_alignment[2]) == 100.000) and (local_alignment[6] == local_alignment[8]) and (local_alignment[7] == local_alignment[9]):
                if (local_alignment[0] not in hundred_percent_coordinates):
                    hundred_percent_coordinates[local_alignment[0]] = dict()
                hundred_percent_coordinates[local_alignment[0]][int(local_alignment[6])] = (local_alignment[6], local_alignment[7], local_alignment[0])

    # Because MAFFT in Genome_Walker outputs a sequence based on consensus, some transcripts may have a few N sequences that BLAST may not align.
    # The below block of code is to treat block of sequences separated by a few N sequences to be one contig rather than several.
    for transcript in hundred_percent_coordinates:
        hundred_percent_coordinates[transcript] = collections.OrderedDict(sorted(hundred_percent_coordinates[transcript].items()))
        temp_dict = collections.OrderedDict()
        skip_next_iter = False
        for i in range(0, len(hundred_percent_coordinates[transcript].items())):
            if skip_next_iter == True:
                skip_next_iter = False
                continue
            if (i != (len(hundred_percent_coordinates[transcript].items()) - 1)):
                current_contig = hundred_percent_coordinates[transcript].items()[i]
                next_contig = hundred_percent_coordinates[transcript].items()[i+1]

                if (int(next_contig[1][0]) - int(current_contig[1][1])) >= 501:
                    temp_dict[current_contig[0]] = current_contig[1]

                else :
                    temp_dict[current_contig[0]] = [current_contig[1][0], next_contig[1][1], transcript]
                    skip_next_iter = True
            elif (i == (len(hundred_percent_coordinates[transcript].items()) - 1)):
                current_contig = hundred_percent_coordinates[transcript].items()[i]
                temp_dict[current_contig[0]] = current_contig[1]

        hundred_percent_coordinates[transcript] = temp_dict
    return hundred_percent_coordinates

def Save_transcript_sequences(GW_fasta, group_dictionary):
    # Store transcript sequences to a dictionary.
    # Store only those transcript sequences which have been grouped. This information will be used to create the gene model later.
    for group_num in group_dictionary:
        try:
            records = SeqIO.parse(GW_fasta, "fasta")
        except IOError:
            print "%s file does not exist" % (GW_fasta)

        for record in records:
            if (record.name in group_dictionary[group_num]):
                group_dictionary[group_num][record.name] = record.seq
    return group_dictionary

def Merge_matching_contigs(BLAST_file, group_dictionary, hundred_percent_coordinates, min_perc_identity=95.000):
    # Contigs between transcripts that closely matched are merged together.

    contig_counter = 1
    Gene_model = dict()
    for group_num in group_dictionary:
        Gene_model[group_num] = dict()

    # Start building the gene models by finding out which contigs overlap and find out how many contigs each gene_model have
    with open(BLAST_file, 'r') as ins:
        for line in ins:
            local_alignment = line.rstrip().split("\t")
            query = local_alignment[0]
            subject = local_alignment[1]
            q_values = tuple()
            s_values = tuple()
            # Check if the q and s are in the same group
            for group in group_dictionary:
                if (query in group_dictionary[group]) and (subject in group_dictionary[group]):
                    # Both q and s in the same group. Now check through all contigs if it has already been added

                    # Check if the alignment is acceptable and save the index if it is.
                    if (query != subject) and (float(local_alignment[2]) >= min_perc_identity):
                        for contig in hundred_percent_coordinates[query]:
                            if (int(hundred_percent_coordinates[query][contig][0]) <= int(local_alignment[6])) and (int(local_alignment[7]) <= int(hundred_percent_coordinates[query][contig][1])):
                                q_values = hundred_percent_coordinates[query][contig]

                        for contig in hundred_percent_coordinates[subject]:
                            if (int(hundred_percent_coordinates[subject][contig][0]) <= int(local_alignment[8])) and (int(local_alignment[9]) <= int(hundred_percent_coordinates[subject][contig][1])):
                                s_values = hundred_percent_coordinates[subject][contig]

                        s_id = "%s_%s" % (s_values[2], s_values[0])
                        q_id = "%s_%s" % (q_values[2], q_values[0])


                        q_already_in = list()
                        s_already_in = list()
                        for contig in Gene_model[group]:

                            if (q_id in Gene_model[group][contig]) and (s_id not in Gene_model[group][contig]):
                                #Gene_model[group][contig][s_id] = s_values
                                q_already_in.append(contig)

                            elif (q_id not in Gene_model[group][contig]) and (s_id in Gene_model[group][contig]):
                                #Gene_model[group][contig][q_id] = q_values
                                s_already_in.append(contig)

                            elif (q_id in Gene_model[group][contig]) and (s_id in Gene_model[group][contig]):
                                q_already_in.append(contig)
                                s_already_in.append(contig)
                                break

                            else: # (q_id not in Gene_model[group][contig]) and (s_id not in Gene_model[group][contig])
                                continue

                        if (len(q_already_in) == 0) and (len(s_already_in) == 0): # Both q_id and s_id are not in any other contig. Make a new contig with these in it
                            Gene_model[group][contig_counter] = dict()
                            Gene_model[group][contig_counter][q_id] = q_values
                            Gene_model[group][contig_counter][s_id] = s_values
                            contig_counter = contig_counter + 1


                        elif (len(q_already_in) >= 1) and (len(s_already_in) == 0):
                            for contig in q_already_in:
                                Gene_model[group][contig][s_id] = s_values
                            if len(q_already_in) > 1:
                                logging.debug("Bug1 %s is in more than one contig, so %s has been added into more than one contigs too" % (q_id, s_id))
                        elif (len(q_already_in) == 0) and (len(s_already_in) >= 1):
                            for contig in s_already_in:
                                Gene_model[group][contig][q_id] = q_values
                            if len(s_already_in) > 1:
                                logging.debug("Bug2 %s is in more than one contig, so %s has been added into more than one contigs too" % (s_id, q_id))

                        elif (len(q_already_in) >= 1) and (len(s_already_in) >= 1):
                            tmp_list = q_already_in + s_already_in
                            if len(tmp_list) > 2:
                                logging.debug("tmp_list bigger than 2: %s" % (tmp_list))
                            elif len(tmp_list) == 1:
                                logging.debug("tmp_list is 1: %s" % (tmp_list))

                            if tmp_list[0] != tmp_list[1]:
                                print("Merging contigs %s" % (tmp_list))
                                for contig in tmp_list:
                                    for contig_id in Gene_model[group][contig]:
                                        Gene_model[group][tmp_list[0]][contig_id] = Gene_model[group][contig][contig_id]
                                        logging.info("%s = %s" % (str(contig_id), str(Gene_model[group][contig][contig_id])))
                                logging.info("Merged contigs %s" % (tmp_list))
                                del tmp_list[0]
                                for contig in tmp_list:
                                    del Gene_model[group][contig]




                        else:
                            assert False


                    else : # Break the loop if q and s are not in the same group or % id is too low
                        continue
                elif ((query not in group_dictionary[group]) and (subject in group_dictionary[group])) or ((query in group_dictionary[group]) and (subject not in group_dictionary[group])):
                    logging.info("%s and %s is not in the same group" % (query, subject))
    return Gene_model, contig_counter

def Add_remaining_contigs(contig_counter, Gene_model, hundred_percent_coordinates, group_dictionary):
    # Add the remaining contigs of transcripts that are not shared between other transcripts within the same group/gene_model
    for Gene in Gene_model:
        for transcript in hundred_percent_coordinates:
            if (transcript in group_dictionary[Gene]):
                for contig in hundred_percent_coordinates[transcript]:
                    q_values = hundred_percent_coordinates[transcript][contig]
                    q_id = q_values[2] + "_" + q_values[0]
                    stand_alone_contig = True
                    for contig2 in Gene_model[Gene]:
                        if (q_id in Gene_model[Gene][contig2]):
                            stand_alone_contig = False
                            break
                            # This contig is already merged with a matching contig with the other transcript
                    if (stand_alone_contig == True):
                        #logging.info("Adding into contig%s (%s = %s)" % (contig_counter, q_id, q_values))
                        Gene_model[Gene][contig_counter] = dict()
                        Gene_model[Gene][contig_counter][q_id] = q_values
                        contig_counter = contig_counter + 1




    return Gene_model

def Check_duplication(Gene_model):
    # To check if a contig is found more than once.
    tmp_list = list([contig2 for group in Gene_model for contig in Gene_model[group] for contig2 in Gene_model[group][contig]])
    duplicate_dict = collections.Counter(tmp_list)
    for contig in duplicate_dict:
        if duplicate_dict[contig] > 1:
            logging.debug("Contig %s is part of %d contigs" % (contig, duplicate_dict[contig]))

def _check_correct_grouping(group_dictionary):
    # Not in use of main(). Ignore.
    for group in group_dictionary:
        print group
        for transcript in group_dictionary[group]:
            print transcript

def Save_transcript_contig_order(group_dictionary, Gene_model, hundred_percent_coordinates):
    # The contig order for each transcripts are stored in a dictionary.
    # This information is later used to order the contigs to create the gene_model.
    Order = dict()
    for group in group_dictionary:
        Order[group] = dict()

        for transcript in group_dictionary[group]:
            Order[group][transcript] = list()
            previous_contig = str()
            #print "\n"
            #print ("%s: %d" % (transcript, len(hundred_percent_coordinates[transcript])))

            for contig in hundred_percent_coordinates[transcript]:
                #print hundred_percent_coordinates[transcript][contig]
                for contig2 in Gene_model[group]:
                    q_values = hundred_percent_coordinates[transcript][contig]
                    q_id = "%s_%s" % (q_values[2], q_values[0])

                    if q_id in Gene_model[group][contig2]:
                        #print ("%s's contig%s is in Group%s's contig%s" % (transcript, contig, group, contig2))

                        if contig2 not in Order[group][transcript]:
                            Order[group][transcript].append(contig2)

                        elif (contig2 in Order[group][transcript]) and (previous_contig != contig2):
                            # Error message
                            logging.warning("contig%s" % (contig2))
                        previous_contig = contig2


    '''
    for group in Order:
        for transcript in Order[group]:
            print ("%s = %s" % (transcript, Order[group][transcript]))
    '''
    return Order

def Contig_relatives(Gene_model, Order):
    # Stores the information about a position of a contig relative to all other contigs.
    relative_order = dict()
    for group in Gene_model:
        relative_order[group] = dict()

        for contig in Gene_model[group]:
            relative_order[group][contig] = dict()
            relative_order[group][contig]["left"] = set()
            relative_order[group][contig]["right"] = set()

            for group2 in Order:
                for transcript in Order[group2]:
                    if contig in Order[group2][transcript]:
                        for left_element in range(0, Order[group2][transcript].index(contig), 1):
                            relative_order[group][contig]["left"].add(Order[group2][transcript][left_element])

                        for right_element in range(len(Order[group2][transcript])-1, Order[group2][transcript].index(contig), -1):
                            relative_order[group][contig]["right"].add(Order[group2][transcript][right_element])
    '''for group in relative_order:
        for contig in relative_order[group]:
            print contig
            print ("left: %s" % (relative_order[group][contig]["left"]))
            print ("right: %s" % (relative_order[group][contig]["right"]))'''
    return relative_order

def Finalise_Gene_Model(Order, relative_order):
    # Constructs the Gene_model for each group.
    Final_order = dict()

    for group in Order:
        Final_order[group] = list()
        for transcript in Order[group]:
            if len(Final_order[group]) == 0:
                Final_order[group].extend(Order[group][transcript])
            else:
                for element in Order[group][transcript]:
                    if element not in Final_order[group]:

                        #forward
                        s_L = -1
                        index = 0
                        for left_element in Final_order[group]:

                            '''if type(left_element) == list:
                                for element2 in left_element:
                                    if element2 in relative_order[group][element]["left"]:
                                        s_L = Final_order.index(left_element)'''

                            if left_element in relative_order[group][element]["left"]:
                                s_L = Final_order[group].index(left_element)

                        #backward
                        s_R = len(Final_order[group])
                        for right_element_idx in range(len(Final_order[group])-1, -1, -1):

                            '''if type(Final_order[right_element_index]) == list:
                                for element2 in Final_order[right_element_index]:
                                    if element2 in relative_order[group][element]["left"]:
                                        s_L = Final_order.index(Final_order[right_element_index])'''

                            if Final_order[group][right_element_idx] in relative_order[group][element]["right"]:
                                s_R = right_element_idx

                        if s_L > s_R:
                            #print s_L, s_R
                            assert False, "s_L > s_R"

                        #print relative_order[1][element]["left"]
                        #print relative_order[1][element]["right"]
                        #print s_L, s_R


                        if (s_L == -1) and (s_R == 0):
                            Final_order[group].insert(0, element)

                        elif (s_L == len(Final_order[group])-1) and (s_R == len(Final_order[group])):
                            Final_order[group].insert(s_R, element)

                        else :
                            if (s_R - s_L) == 1:
                                Final_order[group].insert(s_R, element)

                            elif (s_R - s_L) > 1:
                                Final_order[group].insert(s_R, element)


    return Final_order

    '''
        for group in Gene_model:
            print ("Group: %s" % (group))
            for contig in Gene_model[group]:
                print ("Contig%s = %s" % (contig, Gene_model[group][contig]))
        for transcript in hundred_percent_coordinates:
            print hundred_percent_coordinates[transcript]'''

def Group_into_clusters(Final_order, Order):
    # For the contigs to be considered as a cluster, three conditions must be met.
    # 1. The contigs are unique (i.e. the contig is only found in one transcript)
    # 2. The contigs must be next to one another within the Final_order
    # 3. The contigs are found from the same transcript

    tmp_list = list([contig for group in Order for transcript in Order[group] for contig in Order[group][transcript]])
    duplicate_dict = collections.Counter(tmp_list) # A dictionary to check which contigs are unique

    contig_clusters = dict()
    previous_contig = int()
    cluster_counter = 1

    for group in Final_order:
        contig_clusters[group] = dict()
        previous_contig = int()
        cluster_counter = 1

        for contig in Final_order[group]:
            if (previous_contig == 0):
                pass

            elif (duplicate_dict[previous_contig] == 1) and (duplicate_dict[contig] == 1):
                same_transcript = False
                for transcript in Order[group]:
                    if (contig in Order[group][transcript]) and (previous_contig in Order[group][transcript]):
                        same_transcript = True
                        break
                if same_transcript == True:
                    if cluster_counter not in contig_clusters[group]:
                        contig_clusters[group][cluster_counter] = list()
                        contig_clusters[group][cluster_counter].append(previous_contig)
                        contig_clusters[group][cluster_counter].append(contig)

                    elif cluster_counter in contig_clusters[group]:
                        contig_clusters[group][cluster_counter].append(contig)

                elif same_transcript == False:
                    if cluster_counter in contig_clusters[group]:
                        cluster_counter = cluster_counter + 1


            elif (duplicate_dict[previous_contig] == 1) and (duplicate_dict[contig] >= 2):
                if cluster_counter in contig_clusters[group]:
                    cluster_counter = cluster_counter + 1


            elif (duplicate_dict[previous_contig] >= 2) and (duplicate_dict[contig] == 1):
                pass
            elif (duplicate_dict[previous_contig] >= 2) and (duplicate_dict[contig] >= 2):
                pass

            previous_contig = contig
    return contig_clusters

def Detect_ambiguous_contigs(Final_order, relative_order, contig_clusters):
    # Checks which contigs are ambiguously positioned.
    # These contigs are flagged as ambiguous so that the user is notified that there are >1 possible positions for these contigs with the current information given.
    ambiguous_contigs = dict()

    for group in Final_order:
        ambiguous_contigs[group] = list()
        iter = 0
        for element in Final_order[group]:
            iter = iter + 1
            if iter == 1:
                end_counter = 0
                for contig in relative_order[group]:
                    if len(relative_order[group][contig]["left"]) == 0:
                        end_counter = end_counter + 1
                        if end_counter == 2:
                            ambiguous_contigs[group].append(element)
                            break

            elif iter == len(Final_order[group]):
                end_counter = 0
                for contig in relative_order[group]:
                    if len(relative_order[group][contig]["right"]) == 0:
                        end_counter = end_counter + 1
                        if end_counter == 2:
                            ambiguous_contigs[group].append(element)
                            break
            else :
                s_L = -1
                s_R = len(Final_order[group])
                for element2 in Final_order[group]:

                    if element2 in relative_order[group][element]["left"]:
                        s_L = Final_order[group].index(element)
                        continue
                    if element2 in relative_order[group][element]["right"]:
                        s_R = Final_order[group].index(element)
                        continue

                if (s_R - s_L) == 2:
                    continue

                elif (s_R - s_L) > 2:
                    if s_L == -1:
                        for index in range(0, s_R):
                            if Final_order[group][index] not in ambiguous_contigs[group]:
                                ambiguous_contigs[group].append(Final_order[group][index])
                    ambiguous_contigs[group].append(element)
                    print ("Contig %d Added" % (element))
                    print ("s_L = %d, s_R = %d" % (s_L, s_R))
                    print relative_order[group][element]

        for cluster in contig_clusters[group]:
            for contig in contig_clusters[group][cluster]:
                if contig in ambiguous_contigs[group]:
                    for contig2 in contig_clusters[group][cluster]:
                        if contig2 not in ambiguous_contigs[group]:
                            ambiguous_contigs[group].append(contig2)
                    break

    return ambiguous_contigs

def Generate_sequences(Final_order, Gene_model, group_dictionary, ambiguous_contigs):
    # Contig sequences saved with appropriate sequences gaps to create output gene_model FASTA file.
    # Sequences gaps are represented by either 500 "Z" or "N" sequences.
    # Contigs with Z sequences on either side (5' and 3') are indicative that this contig is ambiguously placed.

    sequences = dict()
    for group in Final_order:
        sequences[group] = list()


        for current_contig in Final_order[group]:
            if current_contig == Final_order[group][-1]:
                sequences = Added_contig(Gene_model, group_dictionary, current_contig, group, sequences)
                continue

            else:
                next_contig_idx = Final_order[group].index(current_contig) + 1
                next_contig = Final_order[group][next_contig_idx]

                sequences = Added_contig(Gene_model, group_dictionary, current_contig, group, sequences)
                if (current_contig or next_contig) in ambiguous_contigs[group]:
                    sequences[group].append("Z"*500)
                    print "Zs ADDED"
                else:
                    sequences[group].append("N"*500)
                    print "Ns ADDED"
    return sequences

def Added_contig(Gene_model, group_dictionary, current_contig, group, sequences):
    # This function is nested inside the Generate_sequences function.
    # Contigs which match between different transcripts are merged and the consensus sequence is obtained in this function.

    # If the current contig is found in more than 1 transcripts, gain the consensus then add the consensus sequence.
    if len(Gene_model[group][current_contig]) > 1:
        # MAFFT
        for indv_contig in Gene_model[group][current_contig]:
            t_id = Gene_model[group][current_contig][indv_contig]
            with open('gene_models/tmp_mafft.fasta.txt', 'a') as ins:
                ins.write(">%s\n%s\n" % (indv_contig, group_dictionary[group][t_id[2]][int(t_id[0])-1:int(t_id[1])])) # -1 on start to not leave out the single sequence
                ins.close()
        p1 = subprocess.Popen("mafft --quiet --auto gene_models/tmp_mafft.fasta.txt", shell=True, universal_newlines = True, stdout=subprocess.PIPE)
        stdout, stderr = p1.communicate()
        align = AlignIO.read(StringIO(stdout), "fasta")
        summary_align = AlignInfo.SummaryInfo(align)
        consensus = summary_align.dumb_consensus(ambiguous='N')
        sequences[group].append(str(consensus).upper())


        with open('gene_models/tmp_mafft.fasta.txt', 'w') as ins:
            ins.close()
            # To reset the file for the next round of MAFFT

    # If the current contig is found in only one transcript.
    elif len(Gene_model[group][current_contig]) == 1:
        # Just Add
        for indv_contig in Gene_model[group][current_contig]:
            t_id = Gene_model[group][current_contig][indv_contig]
            group_dictionary[group][t_id[2]]
            sequences[group].append(str(group_dictionary[group][t_id[2]][int(t_id[0])-1:int(t_id[1])]).upper())


    return sequences

def Write_to_files(sequences, group_dictionary, Final_order, ambiguous_contigs, hundred_percent_coordinates):
    # Creates three output files
    # gene_models.fasta.txt contains all the FASTA sequences of the gene models.
    # Summary_info_human.txt is a human friendly file with basic information about each gene model
    # Summary_info_computer.txt is a tab deliminated file with the gene_model number on the first column and the transcript IDs in the consecutive columns used for making the gene_models.
    for group in sequences:

        with open('gene_models/gene_models.fasta.txt', 'a') as ins:
            ins.write(">gene_model%d \n" % (group))
            ins.write("".join(sequences[group])+"\n \n")
            print ("Length: " + str(len("".join(sequences[group]))))
            print ("No. Contigs: " + str(len(sequences[group])) + " \n \n \n \n")

            ins.close()

        # Write Human readable summary file for all gene model
        with open('gene_models/Summary_info_human.txt', 'a') as ins:
            ins.write("-- Gene Model %d -- \n" % (group))
            ins.write('Transcript IDs: \n')
            for transcript_id in sorted(group_dictionary[group], key=str):
                if sorted(group_dictionary[group], key=str).index(transcript_id) + 1 == len(group_dictionary[group]):
                    ins.write(transcript_id)
                else:
                    ins.write(transcript_id + ', ')
            ins.write('\n')
            ins.write('Length of Gene Model (including gaps): ' + str(len("".join(sequences[group]))) + '\n')
            #ins.write('Length of Gene Model (excluding gaps): ' + str(len("".join(inserting_separate_seq))) + '\n')
            ins.write('No. Contigs: %s \n' % (str(len(Final_order[group]))))
            ins.write('Order of contigs: %s \n' % (str(Final_order[group])))
            ins.write('Ambiguous contigs: %s \n' % (str(ambiguous_contigs[group])))
            #for contigs in sorted(ambiguous[Gene], key=int):
            #    if (sorted(ambiguous[Gene], key=int).index(contigs) + 1 == len(ambiguous[Gene])):
            #        ins.write(str(contigs))
            #    else:
            #        ins.write(str(contigs) + ', ')
            ins.write('\n \n \n \n')
            ins.close()


        # Write computer readable summary file
        with open('gene_models/Summary_info_computer.txt', 'a') as ins:
            ins.write("Gene_Model" + str(group) + "\t")
            for transcript_id in group_dictionary[group]:
                ins.write(transcript_id + '\t')
            ins.write('\n')
            ins.close()

    with open('gene_models/Summary_info_human.txt', 'a') as ins:
        ins.write("-- Only variant transcripts -- \n")
        ins.close()

    with open('gene_models/Summary_info_computer.txt', 'a') as ins:
        ins.write("lonely_transcripts\t")
        ins.close()

    for transcript in sorted(hundred_percent_coordinates, key=str):
        grouped = False
        for group in group_dictionary:
            if (transcript in group_dictionary[group]):
                grouped = True
                break
        if (grouped == False):
            with open('gene_models/Summary_info_human.txt', 'a') as ins:
                if (sorted(hundred_percent_coordinates, key=str).index(transcript) + 1 == len(hundred_percent_coordinates)):
                    ins.write(transcript + '\n\n\n\n\n')
                else:
                    ins.write(transcript + ', ')
                ins.close()
            with open('gene_models/Summary_info_computer.txt', 'a') as ins:
                ins.write(transcript + '\t')
                ins.close()

def Create_separate_gene_model_folder(sequences, Final_order, path):
    # For IF(Optional args) the user want information for each individual gene model.
    for group in sequences:

        if not os.path.exists("%sgene_model%d_folder/" % (path, group)):
            subprocess.call("mkdir %sgene_model%d_folder" % (path, group), shell=True)
        if os.path.exists("%sgene_model%d_folder/final_all_seq.fasta.txt" % (path, group)):
            subprocess.call("rm %sgene_model%d_folder/final_all_seq.fasta.txt" % (path, group), shell=True)
        if os.path.exists("%sgene_model%d_folder/all_in_one_seq.fasta.txt"):
            subprocess.call('rm %sgene_model%d_folder/all_in_one_seq.fasta.txt' % (path, group), shell=True)

        with open("%sgene_model%d_folder/final_all_seq.fasta.txt" % (path, group), 'w') as ins:
            contig_idx = 0
            for seq_idx in range(0, len(sequences[group]), 2):
                ins.write(">contig%s_%s\n" % (str(Final_order[group][contig_idx]), str(group)))
                ins.write(sequences[group][seq_idx] + "\n \n")
                contig_idx = contig_idx + 1
            ins.close()

        with open("%sgene_model%d_folder/all_in_one.fasta.txt" % (path, group), 'w') as ins:
            ins.write(">gene_model%s\n" % (str(group)))
            ins.write("".join(sequences[group])+"\n \n")
            ins.close()


def main(args):
    logging.basicConfig(filename='Error.log',level=logging.DEBUG)
    with open('Error.log', 'a') as ins:
        ins.write("\n\n")
        ins.close()

    # Create directory if not already exist
    if not os.path.exists("./gene_models/"):
        subprocess.call("mkdir gene_models", shell=True)

    # Run BLAST between the mRNA and Genome Walked transcripts files to obtain the exon indices of all the transcripts
    subprocess.call("makeblastdb -dbtype nucl -in %s" % (args.Genome_Walked_seq), shell=True) # Contradicting with test.py!!
    subprocess.call("blastn -query %s -db %s -outfmt 6 > gene_models/exon_indices.txt" % (args.Original_mRNA_seq, args.Genome_Walked_seq), shell=True)

    exon_indices = Obtain_exon_indices("gene_models/exon_indices.txt")

    subprocess.call("blastn -query %s -db %s -outfmt 6 > gene_models/tmp_blast.txt" % (args.Genome_Walked_seq, args.Genome_Walked_seq), shell=True)

    group_dictionary = Group_transcripts(exon_indices, "gene_models/tmp_blast.txt", args.p, args.i)

    # If any transcripts were clustered, initiate the second phase. If not, terminate the process.
    if (len(group_dictionary) == 0):
        print "\nNo isoforms. All transcripts are expressed from separate genes. \n"
        subprocess.call("rm gene_models/tmp_blast.txt", shell=True)
        subprocess.call("rm gene_models/exon_coordinates.txt", shell=True)
        print "Terminating process"
        sys.exit() # stops the execution of the script
    elif (len(group_dictionary) <= -1):
        assert False, "The value of group_dictionary is < 0. Please contact the developer to debug"
    else : #  group_dictionary is >= 1
        print " \n \nTranscript variants found. \n Starting gene model construction \n \n"


    # Start of the second phase
    hundred_percent_coordinates = Obtain_contig_indices("gene_models/tmp_blast.txt")

    group_dictionary = Save_transcript_sequences(args.Genome_Walked_seq, group_dictionary)

    Gene_model, contig_counter = Merge_matching_contigs("gene_models/tmp_blast.txt", group_dictionary, hundred_percent_coordinates, args.p)

    Gene_model = Add_remaining_contigs(contig_counter, Gene_model, hundred_percent_coordinates, group_dictionary)
    # The lonely_contigs dictionary stores contigs which are not shared between other transcripts within the same group/gene_model. This dictionary is used in a later function, but will be removed in the future.

    Check_duplication(Gene_model)

    Order = Save_transcript_contig_order(group_dictionary, Gene_model, hundred_percent_coordinates)

    relative_order = Contig_relatives(Gene_model, Order)

    Final_order = Finalise_Gene_Model(Order, relative_order)

    for group in Final_order:
        print ("Group: %s" % (group))
        print Final_order[group]
        for transcript in Order[group]:
            print ("%s = %s" % (transcript, Order[group][transcript]))


    contig_clusters = Group_into_clusters(Final_order, Order)
    print "\n \n Clusters"
    for group in contig_clusters:
        print ("Group: %s" % (group))
        for cluster in contig_clusters[group]:
            print ("Cluster%d: %s" % (cluster, contig_clusters[group][cluster]))


    ambiguous_contigs = Detect_ambiguous_contigs(Final_order, relative_order, contig_clusters)
    print "\n\n\n"

    sequences = Generate_sequences(Final_order, Gene_model, group_dictionary, ambiguous_contigs)

    Write_to_files(sequences, group_dictionary, Final_order, ambiguous_contigs, hundred_percent_coordinates)

    if (args.separate == True):
        Create_separate_gene_model_folder(sequences, Final_order, "gene_models/")

    subprocess.call("rm gene_models/exon_indices.txt", shell=True)
    subprocess.call("rm gene_models/tmp_mafft.fasta.txt", shell=True)
    if (args.blast == False):
        subprocess.call("rm gene_models/tmp_blast.txt", shell=True)
    print "\n\n\n -- Finish -- \n\n\n"
    print "Please read the generated summary files"

if __name__ == "__main__":
    # Positional and Optional Arguments
    parser = argparse.ArgumentParser(description='')
    parser.add_argument('Genome_Walked_seq', metavar = 'filename.FASTA', help='Genome Walked sequence file')
    parser.add_argument('Original_mRNA_seq', metavar = 'filename.FASTA', help='Original sequence file before running Genome Walker')
    parser.add_argument('-p', nargs='?', metavar='FLOAT', default=95.000, type=float, help='Minimum percentage id for the grouping to occur based on the BLAST result (default: %(default)s)')
    parser.add_argument('-i', nargs='?', metavar='INT', default=30, type=int, help='Minimum number of intronic sequence match for transcripts to be grouped together if transcripts only have one matching exon (default: %(default)s)')
    parser.add_argument('-s','--separate', action='store_true', help='Creates a new folder that contains information about individual gene models (default: %(default)s)')
    parser.add_argument('-b','--blast', action='store_true', help='Does not delete the BLAST output (default: %(default)s)')
    args = parser.parse_args()

    main(args)
