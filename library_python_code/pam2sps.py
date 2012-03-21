#!/usr/bin/env python

def pam2sps(tree_file, conversion, outfile=None):
    import re
    reg_ex = re.compile('(?<=:)[0-9.]+')
    convert_pam_to_sps = lambda a: str(0.01*float(a.group()))
    convert_sps_to_pam = lambda b: str(100*float(b.group()))

    input_string = open(tree_file).read()
    if conversion == 'pam2sps':
        output_string = reg_ex.sub(convert_pam_to_sps,input_string)
    elif conversion == 'sps2pam':
        output_string = reg_ex.sub(convert_sps_to_pam,input_string)
    else: output_string = reg_ex.sub(convert_pam_to_sps,input_string)

    if outfile:
        writer = open(outfile,'w')
        writer.write(output_string)
        writer.close()
    else: print output_string
