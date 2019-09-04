#!usr/bin/python
#__author__ = "Sujan Mamidi"

######################################################################
## This python program takes a vcf file (phased and bz2 compressed) ##
## and estimates 
# a) Ref allele frequency
# b) Var allele frequency
# c) Minor allele frequency
# d) Expected and observed heterozygosity and then the FIS
# e) Pi for each SNP
# f) number of samples without missing data

## Usage: python SNP_stats_Phased_vcf.py my.vcf.bz2

######################################################################
from __future__ import division
from sys import argv
import bz2

#===================================
def real_main():
    # Defining the program options
    usage = "usage: %prog [SNPfile]"

    snpfile = argv[1]
    outputfile = open((snpfile + ".snpStats"), 'w')

   # Checking the user input
    if len(argv) < 2:
        print("Incorrect number of arguments \n")
        print('%s' % usage)

    else:
        outputfile.write("Chrom pos Ref Var Reffreq Varfreq Pi ObsHet ExpHet Inbreeding_Coeff MAF samples\n")
        source_file = bz2.BZ2File(snpfile,"r")
	for line in source_file:
            if line.startswith("Chrom"):
                continue
            if line.startswith("chrom"):
                continue
            if line.startswith("#"):
                continue
            myline = line.strip().split('\t')
            ncol = len(myline)
            nsamples = ncol - 9
            chrom =  myline[0]
            pos = int(myline[1])
            ref = myline[3]
            var = myline[4]
            countRef = countVar = countmissing = countHet = 0
            for i in range (9,ncol,1):
                if myline[i] == '0|0':
                    countRef += 1
                elif myline[i] == '1|1':
                    countVar += 1
                elif myline[i] == '.|.' :
                    countmissing += 1
                else:
                    countHet += 1

            effsamples = nsamples - countmissing

            try:
                reffreq = round(float(((2 * countRef) + countHet) / (2 * effsamples)),4)
            except ZeroDivisionError:
                reffreq= 0

            try:
                varfreq = round(float(((2 * countVar) + countHet) / (2 * effsamples)),4)
            except ZeroDivisionError:
                varfreq = 0

            try:
                pi = round(float((2 * reffreq * varfreq * 2 * nsamples) / ((2 * nsamples) - 1)),4)
            except ZeroDivisionError:
                pi = 0

            try:
                obsHet = round(float(countHet / effsamples),4)
            except ZeroDivisionError:
                obsHet = 0

            try:
                expHet = round(1 - float(((reffreq * reffreq) + (varfreq * varfreq))),4)
            except ZeroDivisionError:
                expHet = 0

            try:
                FIS = round(float((expHet - obsHet) / expHet),4)
            except ZeroDivisionError:
                FIS = 0

            if reffreq >= varfreq:
                minorallelefreq = varfreq
            else:
                minorallelefreq = reffreq

            outputfile.write(chrom)
            outputfile.write(' ')
            outputfile.write(str(pos))
            outputfile.write(' ')
            outputfile.write(ref)
            outputfile.write(' ')
            outputfile.write(var)
            outputfile.write(' ')
            outputfile.write(str(reffreq))
            outputfile.write(' ')
            outputfile.write(str(varfreq))
            outputfile.write(' ')
            outputfile.write(str(pi))
            outputfile.write(' ')
            outputfile.write(str(obsHet))
            outputfile.write(' ')
            outputfile.write(str(expHet))
            outputfile.write(' ')
            outputfile.write(str(FIS))
            outputfile.write(' ')
            outputfile.write(str(minorallelefreq))
            outputfile.write(' ')
            outputfile.write(str(effsamples))
            outputfile.write('\n')

        outputfile.close()

    return

# ==============================================================
if __name__ == '__main__':
    real_main()




