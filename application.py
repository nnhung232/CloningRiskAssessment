# -*- coding: utf-8 -*-
"""
Created on Fri Nov 18 08:46:27 2022

@author: PTHAJD

"""
from distutils.log import debug
from fileinput import filename
import re
import regex
from flask import Flask,render_template,request
from Bio import SeqIO


application = Flask(__name__)

@application.route("/")
def hello():
    return render_template('form.html')

@application.route('/form')
def form():
   return render_template('form.html')

@application.route('/data/', methods = ['POST', 'GET'], strict_slashes=False)
def data():
    if request.method == 'GET':
        return "The URL /data is accessed directly. Try going to '/form' to submit form"
    if request.method == 'POST':
        print("***", request.form)
        ## get input from form
        form_data = request.form
        print("****", form_data)
        promoterSeq = form_data['Sequence']
        print("*****", promoterSeq)
        lowThresh = form_data['lowThreshold']
        highThresh = form_data['highThreshold']
        print("$$$$", form_data, promoterSeq, lowThresh, highThresh)
        web_output_dict = {}

            #promoterSeq = 'ATGGCGCAAGTTAGCAGAATCTGCAATGGTGTGCAGAACCCATCTCTTATCTCCAATCTCTCGAAATCCAGTCAACGCAAATCTCCCTTATCGGTTTCTCTGAAGACGCAGCAGCATCCACGAGCTTATCCGATTTCGTCGTCGTGGGGATTGAAGAAGAGTGGGATGACGTTAATTGGCTCTGAGCTTCGTCCTCTTAAGGTCATGTCTTCTGTTTCCACGGCGTGCATGCTTCACGGTGCAAGCAGCCGGCCCGCAACCGCCCGCAAATCCTCTGGCCTTTCCGGAACCGTCCGCATTCCCGGCGACAAGTCGATCTCCCACCGGTCCTTCATGTTCGGCGGTCTCGCGAGCGGTGAAACGCGCATCACCGGCCTTCTGGAAGGCGAGGACGTCATCAATACGGGCAAGGCCATGCAGGCGATGGGCGCCCGCATCCGTAAGGAAGGCGACACCTGGATCATCGATGGCGTCGGCAATGGCGGCCTCCTGGCGCCTGAGGCGCCGCTCGATTTCGGCAATGCCGCCACGGGCTGCCGCCTGACGATGGGCCTCGTCGGGGTCTACGATTTCGACAGCACCTTCATCGGCGACGCCTCGCTCACAAAGCGCCCGATGGGCCGCGTGTTGAACCCGCTGCGCGAAATGGGCGTGCAGGTGAAATCGGAAGACGGTGACCGTCTTCCCGTTACCTTGCGCGGGCCGAAGACGCCGACGCCGATCACCTACCGCGTGCCGATGGCCTCCGCACAGGTGAAGTCCGCCGTGCTGCTCGCCGGCCTCAACACGCCCGGCATCACGACGGTCATCGAGCCGATCATGACGCGCGATCATACGGAAAAGATGCTGCAGGGCTTTGGCGCCAACCTTACCGTCGAGACGGATGCGGACGGCGTGCGCACCATCCGCCTGGAAGGCCGCGGCAAGCTCACCGGCCAAGTCATCGACGTGCCGGGCGACCCGTCCTCGACGGCCTTCCCGCTGGTTGCGGCCCTGCTTGTTCCGGGCTCCGACGTCACCATCCTCAACGTGCTGATGAACCCCACCCGCACCGGCCTCATCCTGACGCTGCAGGAAATGGGCGCCGACATCGAAGTCATCAACCCGCGCCTTGCCGGCGGCGAAGACGTGGCGGACCTGCGCGTTCGCTCCTCCACGCTGAAGGGCGTCACGGTGCCGGAAGACCGCGCGCCTTCGATGATCGACGAATATCCGATTCTCGCTGTCGCCGCCGCCTTCGCGGAAGGGGCGACCGTGATGAACGGTCTGGAAGAACTCCGCGTCAAGGAAAGCGACCGCCTCTCGGCCGTCGCCAATGGCCTCAAGCTCAATGGCGTGGATTGCGATGAGGGCGAGACGTCGCTCGTCGTGCGTGGCCGCCCTGACGGCAAGGGGCTCGGCAACGCCTCGGGCGCCGCCGTCGCCACCCATCTCGATCACCGCATCGCCATGAGCTTCCTCGTCATGGGCCTCGTGTCGGAAAACCCTGTCACGGTGGACGATGCCACGATGATCGCCACGAGCTTCCCGGAGTTCATGGACCTGATGGCCGGGCTGGGCGCGAAGATCGAACTCTCCGATACGAAGGCTGCCTGA'

        ## check for file upload
        f = request.files['file']
        print("****$$", f)
        
        ### user uploaded a fasta file
        if (f):
            f.save(f.filename)  
            for seq_record in SeqIO.parse(f.filename, "fasta"):
                promoterSeq = seq_record.seq
                promoterName = seq_record.id
        ## verify DNA string
                verifyStatus, promoterSeq = verify_dna_string(promoterSeq)
                if(not verifyStatus):
                    print("character other than ATGC found in input string")
                    return render_template('error.html')
        ## run the analysis code
                output_dict = run_prom_seq_code(promoterSeq,lowThresh,highThresh,promoterName)
                web_output_dict = {**web_output_dict, **output_dict}

                
        else:
            ## user pasted in a single fasta sequence 
            promoterSeq = promoterSeq
            print("*****$", promoterSeq)
            ## verify DNA string
            verifyStatus, promoterSeq = verify_dna_string(promoterSeq)
            print("@@@@@", verifyStatus, promoterSeq)
            if(not verifyStatus):
                print("character other than ATGC found in input string")
                return render_template('error.html')
        ## run the analysis code
            output_dict = run_prom_seq_code(promoterSeq,lowThresh,highThresh,"")
            web_output_dict = {**output_dict}

        return render_template('data.html',form_data = web_output_dict)


        ## generate ouput dictionary and print
#        return render_template('data.html',form_data = form_data)
#        return render_template('data.html',form_data = output_dict)

### combine two dicts into one
#        web_output_dict = {**form_data,**output_dict}

#        web_output_dict = {**output_dict}



def verify_dna_string(s1):
    checkStatus = False
    s2 = s1.upper()
#    print("upper case string:")
#    print(s2)
    if bool(re.search("[^ATGC]", str(s2))):
         #Character other than ATGC was found
        checkStatus = False
    else:
         #No character other then ATGC was found
        checkStatus = True

    return checkStatus, str(s2)
 
    
 


def run_prom_seq_code(seq1,t1,t2, seqname):
    
    outd = {}
    promoter_sequence = seq1
    dict_key = "User Sequence " + seqname
    outd[dict_key] = promoter_sequence
    dict_key = "Low Threshold " + seqname
    outd[dict_key] = t1
    dict_key = "High Threshold " + seqname
    outd[dict_key] = t2
#   promoter_sequence = 'TTTTTTTTTTTTTTTTTTTTTTTTTGCTCACTCACCACTCACTATCCCCCCCCCCCGCTCACTCACCACTCACTATCACTCAAAAAGAATTTTTTTTTTTTTTTTTTTTTTTTTTCTTATCTCTTGAAACAACAACGGACCGCTCAGGTAAGTACTCCTCTTCTCTCTCTCGGCAAAAAATTAAAAAACCCTAGGCGGATTCGTTTCTTCTCTCTCCCTCTCGACCGCGTTTTCTCTCTTCTGTTAGCGATCACTCGTTGTTTCGTTTCTTCACTTTTTTCGCCTTCGTGTGTCTCAATTCTCTCTCTCTTCTCGGAAAGAGTTGAGATCTCTCGTTGCTAGGGTTTCTCTCTCTCTTTCTCAAACCCCACCCAAAAAAAACCCTAGGGTTTTGGTGATTTGGGGTTTTGGATTTCGATTTTTCAATTCGATCTGCGATTTTCTCGGATCTGTGCCTTTTTTTCTCGTCGTTCGTTTTCGCTCTTCTCTCCTTTTCCTTCGATTTCTCGTCGCTCTTTCTCTCAATTCGTTTCCGTTTTTCCACTTCACTCTGCCTAGGGTTTGGTGAAATTCTGGGTTTCTCGATCGTGTTTACTCTCTTCTCCGATTTTCGATTCTCTGGGTTTTTCCTCAACTCTCTTCGATTTTCCTTCACTCTCTTCGTTCGTTTTCCTCTTTCGATTTATTCGGATCCGATCTCTTTCGTTTCTCTCCACTGACTAATTACTCTCTCTCCTACCTTTTCAGGTTTTACTCGAGGTACCACCATGGCGGCGTCATCTTCGTCCGTCGTGAGCTTCTCGGGCATCTCGTTGTGCAGTACTCACTCGATCTCCAACAAGACCTATCTATTCTCCGCCCACCCGCGCATTTCGGTGTCGTTCCCCAGTAAGCCCAATAGTTTGAAGTCCTTCAAGCAGCTCCAGCTGAAGAAGAACGGACTCTTTGAGAAGTTCTCTCGTACCTCCAGTCGGAGCTTCGTGGTGAGGTGCGACGCGTCGAAGGCCTTGGTACTGTACTCGACGCGGGACGGCCAGACCCACGCAATTGCTTCATACATCGCCTCCTGCATGAAGGAGAAGGCCGAATGCGACGTGATCGACCTCACCCACGGGGAGCACGTGAACCTCACCCAATACGATCAGGTGCTAATCGGTGCGAGTATTCGTTACGGCCACTTCAACGCCGTGCTTGACAAGTTCATCAAGAGAAACGTGGATCAGCTGAACAACATGCCAAGCGCGTTCTTCTGCGTAAACCTCACAGCAAGGAAGCCCGAGAAGCGTACTCCCCAGACAAACCCTTATGTCCGAAAATTCTTGCTTGCTACCCCCTGGCAGCCCGCGTTGTGCGGAGTGTTCGCAGGGGCCCTTCGGTACCCGCGATACCGGTGGATCGACAAGGTGATGATCCAGCTAATAATGCGGATGACTGGGGGAGAGACAGACACGAGCAAGGAGGTCGAGTACACGGATTGGGAGCAGGTTAAGAAGTTCGCGGAGGATTTTGCAAAGCTATCGTACAAGAAGGCCCTCTAGTAGGCCGGCCTAATGGATATATATATATTTGTGCGGCGGATCTGTGTGGTTTGAACACACCGCCGCCGTGCAAGCTAGCTAGTTTTTCTTTGTTTGTTTTGTGTATTTGTGTTGTAAAGAACTATATTTCTTCTGATTGGAATAAATAAAATTATTATTTCATACCTCGTCTTGTGTGGTTCTTGTGCCGCCACAAATACCACCGAGAGAATTGTAGTATCTGTTTGTGTTTCTTGTTGCTACTAGCTGGTTTTGTTTGTTGTGCCGCCGGCGAGGATATATCCTTTTTGGTTCCACTTCTGAGTGATTGTATCAAACACTTGTAAACTCTTCTGATTTGTAATAAAGTTATTGTTTTGTCATACCTTTCTTGTTTGGGTCTAGTATATAATTTTTGCGCCACGCTGCGGTTT'
#promoter_sequence = 'ATGGCGCAAGTTAGCAGAATCTGCAATGGTGTGCAGAACCCATCTCTTATCTCCAATCTCTCGAAATCCAGTCAACGCAAATCTCCCTTATCGGTTTCTCTGAAGACGCAGCAGCATCCACGAGCTTATCCGATTTCGTCGTCGTGGGGATTGAAGAAGAGTGGGATGACGTTAATTGGCTCTGAGCTTCGTCCTCTTAAGGTCATGTCTTCTGTTTCCACGGCGTGCATGCTTCACGGTGCAAGCAGCCGGCCCGCAACCGCCCGCAAATCCTCTGGCCTTTCCGGAACCGTCCGCATTCCCGGCGACAAGTCGATCTCCCACCGGTCCTTCATGTTCGGCGGTCTCGCGAGCGGTGAAACGCGCATCACCGGCCTTCTGGAAGGCGAGGACGTCATCAATACGGGCAAGGCCATGCAGGCGATGGGCGCCCGCATCCGTAAGGAAGGCGACACCTGGATCATCGATGGCGTCGGCAATGGCGGCCTCCTGGCGCCTGAGGCGCCGCTCGATTTCGGCAATGCCGCCACGGGCTGCCGCCTGACGATGGGCCTCGTCGGGGTCTACGATTTCGACAGCACCTTCATCGGCGACGCCTCGCTCACAAAGCGCCCGATGGGCCGCGTGTTGAACCCGCTGCGCGAAATGGGCGTGCAGGTGAAATCGGAAGACGGTGACCGTCTTCCCGTTACCTTGCGCGGGCCGAAGACGCCGACGCCGATCACCTACCGCGTGCCGATGGCCTCCGCACAGGTGAAGTCCGCCGTGCTGCTCGCCGGCCTCAACACGCCCGGCATCACGACGGTCATCGAGCCGATCATGACGCGCGATCATACGGAAAAGATGCTGCAGGGCTTTGGCGCCAACCTTACCGTCGAGACGGATGCGGACGGCGTGCGCACCATCCGCCTGGAAGGCCGCGGCAAGCTCACCGGCCAAGTCATCGACGTGCCGGGCGACCCGTCCTCGACGGCCTTCCCGCTGGTTGCGGCCCTGCTTGTTCCGGGCTCCGACGTCACCATCCTCAACGTGCTGATGAACCCCACCCGCACCGGCCTCATCCTGACGCTGCAGGAAATGGGCGCCGACATCGAAGTCATCAACCCGCGCCTTGCCGGCGGCGAAGACGTGGCGGACCTGCGCGTTCGCTCCTCCACGCTGAAGGGCGTCACGGTGCCGGAAGACCGCGCGCCTTCGATGATCGACGAATATCCGATTCTCGCTGTCGCCGCCGCCTTCGCGGAAGGGGCGACCGTGATGAACGGTCTGGAAGAACTCCGCGTCAAGGAAAGCGACCGCCTCTCGGCCGTCGCCAATGGCCTCAAGCTCAATGGCGTGGATTGCGATGAGGGCGAGACGTCGCTCGTCGTGCGTGGCCGCCCTGACGGCAAGGGGCTCGGCAACGCCTCGGGCGCCGCCGTCGCCACCCATCTCGATCACCGCATCGCCATGAGCTTCCTCGTCATGGGCCTCGTGTCGGAAAACCCTGTCACGGTGGACGATGCCACGATGATCGCCACGAGCTTCCCGGAGTTCATGGACCTGATGGCCGGGCTGGGCGCGAAGATCGAACTCTCCGATACGAAGGCTGCCTGA'
    pro_length = len(promoter_sequence) # define the length of a sequence and use a variable which will be used later
    print("Sequence length = %d" %(pro_length))# d = integer, f=float
    dict_key = "Sequence Length " + seqname
    outd[dict_key] = pro_length
    sequence_gccontent = getGCcontent(promoter_sequence)#variable for call
    print(sequence_gccontent)
    dict_key = "GC Content " + seqname
    outd[dict_key] = sequence_gccontent

#print process details function.
#total_checked = 0
#total_matches_found = 0
#for y in range(0, pro_length - 201): # 0, pro_length (1870), minus 20 bases is 1850
#    print("subsequence starting at %d" %(y))
#    promoter_substring = promoter_sequence[y:y+200:1] #Python substring search = string[begin: end: step], counting up.
##    print(promoter_substring)
#    subsequence_gccontent = getGCcontent(promoter_substring)
     
    homo_matches = gethomopolymer(promoter_sequence)
    dict_key = "Homopolymer Matches " + seqname
    outd[dict_key] = homo_matches
    gc_match = {}
    final_ranges = {}
    final_gc_content = {}
    gc_string_matches = getGCsubstring(promoter_sequence, pro_length, t1, t2)
    dict_key = "GC Matches Above and Below Threshold " + seqname
    outd[dict_key] = gc_string_matches
    print(final_ranges)
    print(final_gc_content)
    for key in final_ranges:
        stop_coordinate = key+final_ranges[key]
        range_substring = promoter_sequence[key-1:final_ranges[key]:1]
        print(key, final_ranges[key], final_gc_content[key], range_substring)
        corrected_subsequence_gccontent = getGCcontent(range_substring)
        print(corrected_subsequence_gccontent)
    print("getting seqrepeats")
    seqrepeatslist = getseqrepeats(promoter_sequence, pro_length)
    print("getting seqrepeats1")
    seqrepeatslist1, seqschecked1, seqstringrepeats = getseqrepeats1(promoter_sequence, pro_length)
    dict_key = "Repeated Strings " + seqname
    outd[dict_key] = seqstringrepeats
         
    
    return outd



def getGCcontent(seq):#first function ever created!
    c=0
    a=0
    g=0
    t=0
    for x in seq:
        if "C" in x:
            c+=1    
        elif "G" in x:
            g+=1
        elif "A" in x:
            a+=1    
        elif "T" in x:
            t+=1
    #print("Base count - C=%d, G=%d, A=%d, T=%d" %(c,g,a,t))
    gc_content=(g+c)*100/(a+t+g+c)
    #print("Overall GC_content= %f" %(gc_content))
    return gc_content

def gethomopolymer(seq):
    space = '<&nbsp>'
    output_string = ""
    total_matches_found = []

    promoter_substring = "TTTTTTTTTTTTTTTTTTTT|CCCCCCCCCCCCCCCCCCCC|GGGGGGGGGGGGGGGGGGGG|AAAAAAAAAAAAAAAAAAAA|ATATATATATATATATATATATATATATAT" #Python substring search = string[begin: end: step], counting up.
    #print(promoter_substring)
    y = ([[m.start(), m.end()] for m in regex.finditer(promoter_substring, seq, overlapped=1)])
    matches_found = len(y)
    if matches_found > 1:
        total_matches_found.append(y)
        #print(y)
        for coordinate_pair in y:
            substring = seq[coordinate_pair[0]:coordinate_pair[1]]
            match_coords = "[" + str(coordinate_pair[0]) + ", " + str(coordinate_pair[1]) + "]"
            print(coordinate_pair, substring)
#            homopolymer_match = coordinate_pair + " " + substring
            homopolymer_match = match_coords + " " + substring + "\n"
            output_string += homopolymer_match         
    return output_string

def getoriginalGCsubstring(promoter_sequence, pro_length): # Function name
    #total_checked = 0
    #total_matches_found = 0
    gc_output_string = ""
    int_window_size = 50
    gc_match = {}
    final_ranges = {}
    final_gc_content = {}
    total_submatches_found = [] # creating empty list to capture final coordinates
    for y in range(0, pro_length - int_window_size +1): # 0, pro_length (1951), minus 50 bases is 1901
        #print("subsequence starting at %d" %(y))
        promoter_substring = promoter_sequence[y:y+int_window_size:1] #Python substring search = string[begin: end: step], counting up.
#    print(promoter_substring)
        subsequence_gccontent = getGCcontent(promoter_substring)
        if (subsequence_gccontent >=70 or subsequence_gccontent <=30):
            extend_seq = True
            while(extend_seq):
                int_window_size = int_window_size +1
                promoter_substring = promoter_sequence[y:y+int_window_size:1]
                subsequence_gccontent = getGCcontent(promoter_substring)
                if (subsequence_gccontent >=70 or subsequence_gccontent <=30):
                   extend_seq = True
                else:
                   extend_seq = False
        gc_match[y]= subsequence_gccontent
        print(y, subsequence_gccontent,promoter_substring)
        total_submatches_found.append(y)
        range_length = y + int_window_size
        gc_str_match = str(y) + " (" + str(range_length) + ") " + str(subsequence_gccontent) +  " " + str(promoter_substring) + "\n"
        int_window_size = 50
        
    return gc_output_string       

def getGCsubstring(promoter_sequence, pro_length, lowt, hight): # Function name
    #total_checked = 0
    #total_matches_found = 0
    gc_output_string = ""
    int_window_size = 50
    gc_match = {}
    final_ranges = {}
    final_gc_content = {}
    total_submatches_found = [] # creating empty list to capture final coordinates
    xxx = range(0, pro_length - int_window_size +1)
#    print(xxx)
    xxx_iter = iter(xxx)
#    print(xxx_iter)
#    for xxx_iter:
#    for xx in range(0, pro_length - init_window_size): # 0, pro_length (1951), minus 50 bases is 1901
#    for xx in xxx: # 0, pro_length (1951), minus 50 bases is 1901
    while(True):
        try:
            y = next(xxx_iter)
        except StopIteration: 
            break
        else:
            y_start = y
#            next(xxx_iter)
        #print("subsequence starting at %d" %(y))
            promoter_substring = promoter_sequence[y_start:y_start+int_window_size:1] #Python substring search = string[begin: end: step], counting up.
#    print(promoter_substring)
            subsequence_gccontent = getGCcontent(promoter_substring)
            if (subsequence_gccontent >=float(hight) or subsequence_gccontent <=float(lowt)):
                extend_seq = True
                while(extend_seq):
#                    next(xxx_iter)
                    int_window_size = int_window_size +1
                    promoter_substring_last = promoter_substring
                    subsequence_gccontent_last = subsequence_gccontent
                    promoter_substring = promoter_sequence[y_start:y_start+int_window_size:1]
                    subsequence_gccontent = getGCcontent(promoter_substring)
                    if (subsequence_gccontent >=float(hight) or subsequence_gccontent <=float(lowt)):
                        extend_seq = True
                    else:
                       extend_seq = False
                for z in range(0, int_window_size):
                    try:
                        y = next(xxx_iter)
                        print(y)
                    except StopIteration:
                        break

                gc_match[y_start]= subsequence_gccontent_last
                print(y_start, subsequence_gccontent_last,promoter_substring_last)
                total_submatches_found.append(y_start)
                range_length = y_start + int_window_size - 1
                gc_str_match = str(y_start+1) + " (" + str(range_length) + ") " + str(subsequence_gccontent_last) +  " " + str(promoter_substring_last) + "\n"
                gc_output_string = gc_output_string + gc_str_match
#            for z in range(0, 49):
#                next(xxx_iter)
            int_window_size = 50
    return gc_output_string       
            #print("subsequence starting at %d" %(y))
        #print(promoter_substring)
#    range_start = 0
#    range_length = 49
#    last_key = -6
#    range_avg_gc = 0
#    range_tot_gc = 0
#    range_data_points = 1
#    for z in gc_match:
#        key = z
##        range_data_points = range_data_points + 1
##        range_tot_gc = range_tot_gc + gc_match[key]
#        print(key, range_start, range_length, gc_match[key])
#        if key -5 >= last_key: # complete an existing range, and start a new one
#            range_avg_gc = range_tot_gc / range_data_points
#            final_ranges[range_start+1] = range_length + range_start + 1
#            final_gc_content[range_start+1] = range_avg_gc
#            print(range_data_points, range_tot_gc)
#            range_end = range_start + range_length
#            gc_str_seq = promoter_sequence[int(range_start):int(range_end)]
#            gc_str_match = str(range_start) + " (" + str(range_length) + ") " + str(range_avg_gc) + " " + str(gc_str_seq) + "\n"
#            gc_output_string += gc_str_match
#            range_start=key
#            range_length = 49
#            last_key = key
#            range_tot_gc = gc_match[key]
#            range_data_points = 1
#        else: # extend an existing range 
#            range_data_points = range_data_points + 1
#            range_tot_gc = range_tot_gc + gc_match[key]
#            range_length = range_length +(key - last_key)
#            last_key = (last_key + (key - last_key))

#    final_ranges[range_start+1] = range_length + range_start + 1  # complete the last range
#    range_avg_gc = range_tot_gc / range_data_points
#    #final_ranges[range_start] = range_length
#    final_gc_content[range_start+1] = range_avg_gc
#    range_end = range_start + range_length
#    gc_str_seq = promoter_sequence[int(range_start):int(range_end)]
#    gc_str_match = str(range_start) + " (" + str(range_length) + ") " + str(range_avg_gc) +  " " + str(gc_str_seq) + "\n"
#    gc_output_string += gc_str_match
    
#    return gc_output_string
            
#def getGCsubstring(promoter_sequence, pro_length): # Function name
#    #total_checked = 0
#    #total_matches_found = 0
#    gc_output_string = ""
#    gc_match = {}
#    final_ranges = {}
#    final_gc_content = {}
#    total_submatches_found = [] # creating empty list to capture final coordinates
#    for y in range(0, pro_length - 51): # 0, pro_length (1951), minus 50 bases is 1901
#        #print("subsequence starting at %d" %(y))
#        promoter_substring = promoter_sequence[y:y+50:1] #Python substring search = string[begin: end: step], counting up.
##    print(promoter_substring)
#        subsequence_gccontent = getGCcontent(promoter_substring)
#        if (subsequence_gccontent >=70 or subsequence_gccontent <=30):
#           gc_match[y]= subsequence_gccontent
#           print(y, subsequence_gccontent,promoter_substring)
#        total_submatches_found.append(y)
#            #print("subsequence starting at %d" %(y))
#        #print(promoter_substring)
#    range_start = 0
#    range_length = 49
#    last_key = -6
#    range_avg_gc = 0
#    range_tot_gc = 0
#    range_data_points = 1
#    for z in gc_match:
#        key = z
##        range_data_points = range_data_points + 1
##        range_tot_gc = range_tot_gc + gc_match[key]
#        print(key, range_start, range_length, gc_match[key])
#        if key -5 >= last_key: # complete an existing range, and start a new one
#            range_avg_gc = range_tot_gc / range_data_points
#            final_ranges[range_start+1] = range_length + range_start + 1
#            final_gc_content[range_start+1] = range_avg_gc
#            print(range_data_points, range_tot_gc)
#            range_end = range_start + range_length
#            gc_str_seq = promoter_sequence[int(range_start):int(range_end)]
#            gc_str_match = str(range_start) + " (" + str(range_length) + ") " + str(range_avg_gc) + " " + str(gc_str_seq) + "\n"
#            gc_output_string += gc_str_match
#            range_start=key
#            range_length = 49
#            last_key = key
#            range_tot_gc = gc_match[key]
#            range_data_points = 1
#        else: # extend an existing range 
#            range_data_points = range_data_points + 1
#            range_tot_gc = range_tot_gc + gc_match[key]
#            range_length = range_length +(key - last_key)
#            last_key = (last_key + (key - last_key))
#
#    final_ranges[range_start+1] = range_length + range_start + 1  # complete the last range
#    range_avg_gc = range_tot_gc / range_data_points
#    #final_ranges[range_start] = range_length
#    final_gc_content[range_start+1] = range_avg_gc
#    range_end = range_start + range_length
#    gc_str_seq = promoter_sequence[int(range_start):int(range_end)]
#    gc_str_match = str(range_start) + " (" + str(range_length) + ") " + str(range_avg_gc) +  " " + str(gc_str_seq) + "\n"
#    gc_output_string += gc_str_match
#    
#    return gc_output_string
    
    
## this function finds repeated sequences of a given size; reports the first match, and any additional matches
## --- only checks that size ; does not try to expand the length of the match    
def getseqrepeats(promoter_sequence, pro_length): # Function name
    #total_checked = 0
    #total_matches_found = 0
    seq = promoter_sequence
    total_submatches_found = [] # creating empty list to capture final coordinates
    sequenceschecked = []
    for y in range(0, pro_length - 13): # 0, pro_length (1951), minus 50 bases is 1901
        #print("subsequence starting at %d" %(y))
        promoter_substring = promoter_sequence[y:y+12:1]
        if promoter_substring not in sequenceschecked:       
            y = ([[m.start(), m.end()] for m in regex.finditer(promoter_substring, seq, overlapped=False)])
            matches_found = len(y)
            if matches_found > 1:
                total_submatches_found.append(y)
        #print(y)
                for list in y: #list = coordinates
                    coordinates = seq[list[0]:list[1]] #generating the sequence
                    print(list, coordinates)
                sequenceschecked.append(promoter_substring)
    return total_submatches_found

## this function finds repeated sequences starting with a given size, and expanding the size-window 1bp at a time
def getseqrepeats12(promoter_sequence, pro_length): # Function name
    #total_checked = 0
    #total_matches_found = 0
    seq_output_string = ""
    seq = promoter_sequence
    init_window_parameter = 13
    expanded = 0
    total_submatches_found1 = [] # creating empty list to capture final coordinates
    sequenceschecked1 = []
    xxx = range(0, pro_length - init_window_parameter)
#    print(xxx)
    xxx_iter = iter(xxx)
#    print(xxx_iter)
#    for xxx_iter:
#    for xx in range(0, pro_length - init_window_size): # 0, pro_length (1951), minus 50 bases is 1901
#    for xx in xxx: # 0, pro_length (1951), minus 50 bases is 1901
    while(True):
        try:
            xx = next(xxx_iter)
        except StopIteration: 
            break
        else:
#        print("xx = ", xx)
            init_window_size = init_window_parameter
        #print("subsequence starting at %d" %(y))
            promoter_substring = promoter_sequence[xx:xx+init_window_size-1:1]
            print(xx)
            if promoter_substring not in sequenceschecked1:       
                y = ([[m.start(), m.end()] for m in regex.finditer(promoter_substring, seq, overlapped=False)])
                matches_found = len(y)
                if matches_found > 1:
                    total_submatches_found1.append(y)
                    sequenceschecked1.append(promoter_substring)
                    expanded = 0
        # start expanding window size to find max size that still matches somewhere else      
                    while (matches_found > 1):
                        expanded = expanded + 1
                        init_window_size = init_window_size + 1    
                        promoter_substring = promoter_sequence[xx:xx+init_window_size-1:1]
                        if promoter_substring not in sequenceschecked1:       
                            y = ([[m.start(), m.end()] for m in regex.finditer(promoter_substring, seq, overlapped=False)])
                            matches_found = len(y)
#                            print(xx, init_window_size, expanded)
                            if (matches_found > 1):
                                sequenceschecked1.append(promoter_substring)
# add right side min window size substring, so we don't hit this anywhere else in the full seq
                                rt_promoter_substring = promoter_substring[-init_window_parameter]
                                sequenceschecked1.append(rt_promoter_substring)
                                
                                total_submatches_found1.append(y)
                
#                    for list in y: #list = coordinates
#                        foundsequencerepeat = seq[list[0]:list[1]-1] #generating the sequence
#                        print(xx, expanded, list, foundsequencerepeat)
                    for mylist in total_submatches_found1[-1]: #list = coordinates
                        foundsequencerepeat = seq[mylist[0]:mylist[1]] #generating the sequence
                        print(xx, expanded, mylist, foundsequencerepeat)
                        seqoutputsubstring = str(xx) + " " + str(expanded) + " " + str(mylist) + " " + str(foundsequencerepeat) + "\n"
                        seq_output_string += seqoutputsubstring
                else:
                    sequenceschecked1.append(promoter_substring)
#        print("expanded = ", expanded)
#            expanded = expanded + init_window_parameter
            expanded = expanded+1
            while(expanded > 1):
                next(xxx_iter)
                expanded = expanded - 1
    return total_submatches_found1, sequenceschecked1, seq_output_string

    

## this function finds repeated sequences starting with a given size, and expanding the size-window 1bp at a time
def getseqrepeats1(promoter_sequence, pro_length): # Function name
    #total_checked = 0
    #total_matches_found = 0
    seq_output_string = ""
    seq = promoter_sequence   # the sequence to search against (this is anything to the right of the iter position)
    init_window_parameter = 13
    expanded = 0
    total_submatches_found1 = [] # creating empty list to capture final coordinates
    sequenceschecked1 = []
    xxx = range(0, pro_length - init_window_parameter)
#    print(xxx)
    xxx_iter = iter(xxx)
#    print(xxx_iter)
#    for xxx_iter:
#    for xx in range(0, pro_length - init_window_size): # 0, pro_length (1951), minus 50 bases is 1901
#    for xx in xxx: # 0, pro_length (1951), minus 50 bases is 1901
    while(True):
        try:
            xx = next(xxx_iter)
        except StopIteration: 
            break
        else:
#        print("xx = ", xx)
            init_window_size = init_window_parameter
        #print("subsequence starting at %d" %(y))
            promoter_substring = promoter_sequence[xx:xx+init_window_size-1:1]
            seq = promoter_sequence[xx:pro_length:1]   # only want to find repeats to the right of our current position ; to avoid re-finding an already-found repeat
            
            if promoter_substring not in sequenceschecked1:       
                y = ([[m.start(), m.end()] for m in regex.finditer(promoter_substring, seq)])
                matches_found = len(y)
                if matches_found > 1:
                    total_submatches_found1.append(y)
                    sequenceschecked1.append(promoter_substring)
                    expanded = 0
        # start expanding window size to find max size that still matches somewhere else      
                    while (matches_found > 1):
                        expanded = expanded + 1
                        init_window_size = init_window_size + 1    
                        promoter_substring = promoter_sequence[xx:xx+init_window_size-1:1]
                        if promoter_substring not in sequenceschecked1:       
                            y = ([[m.start(), m.end()] for m in regex.finditer(promoter_substring, seq)])
                            matches_found = len(y)
#                            print(xx, init_window_size, expanded)
                            if (matches_found > 1):
                                sequenceschecked1.append(promoter_substring)
# add right side min window size substring, so we don't hit this anywhere else in the full seq
                                rt_promoter_substring = promoter_substring[-init_window_parameter]
                                sequenceschecked1.append(rt_promoter_substring)
                                
                                total_submatches_found1.append(y)
                
#                    for list in y: #list = coordinates
#                        foundsequencerepeat = seq[list[0]:list[1]-1] #generating the sequence
#                        print(xx, expanded, list, foundsequencerepeat)
                    for mylist in total_submatches_found1[-1]: #list = coordinates
                        foundsequencerepeat = seq[mylist[0]:mylist[1]] #generating the sequence
#                        print(xx, expanded, mylist, foundsequencerepeat)
                        realcoord1 = mylist[0]+xx
                        realcoord2 = mylist[1]+xx
                        realcoords = "[" + str(realcoord1) + ", " + str(realcoord2) + "]"
                        print(xx, expanded, realcoords, foundsequencerepeat)
#                        seqoutputsubstring = str(xx) + " " + str(expanded) + " " + str(mylist) + " " + str(foundsequencerepeat) + "\n"
                        seqoutputsubstring = str(xx) + " " + str(expanded) + " " + realcoords + " " + str(foundsequencerepeat) + "\n"
                        seq_output_string += seqoutputsubstring
                else:
                    sequenceschecked1.append(promoter_substring)
#        print("expanded = ", expanded)
#            expanded = expanded + init_window_parameter
#            expanded = expanded+1
            while(expanded > 1):
#                print("inital value = " + str(xx))
#                print("iterating " + str(expanded))
                next(xxx_iter)
                expanded = expanded - 1
    return total_submatches_found1, sequenceschecked1, seq_output_string







### this promoter has sequence matches: 
    # 0-25 (25 T's)  ------    89-115   (26 T's)
    # 25-46          ------    56-77   (21 bp's)
    # 182-195        ------    375-388  (12 bp's)
    # 298-311        ------    345-358   (13 bp's)
    # 467-479        ------    660-672   (12 bp's)


    
#        
#        
#        
##            print(sum())
##        if (i>=70 and i<=30):
##            print("In list")
##        else:
##            print("Not in list")



if __name__ == "__main__":
    application.debug = True
    application.run()
