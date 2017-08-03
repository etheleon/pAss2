#!/usr/bin/env python
#investigate how diverse are the sequences from position 0 to N in the MDR region

#import matplotlib.pyplot as plt
from __future__ import division
import feather
import argparse
import re
import logging
import os

from Bio import AlignIO, Align, SeqIO

import pandas as pd
import numpy as np

import math

from IPython.core.debugger import Tracer

logging.basicConfig(level=logging.INFO)

class mdrSeqAnalysis:
    def __init__(self, ko, root):
        self.ko = ko
        self.mdrFile = "{}/pAss11/{}.fna".format(root, ko)
        self.msaFile= "{}/pAss03/{}.msa".format(root, ko)
        filesExists = os.path.isfile(self.mdrFile) & os.path.isfile(self.mdrFile)
        if filesExists:
            self.mdr      = {}
            self.msa      = {}

            firstRecord   = next(SeqIO.parse(self.mdrFile, "fasta"))
            matches       = re.search("msaStart:(\d+) msaEND:(\d+)", firstRecord.description)
            self.msaStart = int(matches.group(1))
            self.msaEnd   = int(matches.group(2))
            with open(self.mdrFile) as f:
                for record in SeqIO.parse(f, "fasta"):
                    self.mdr[record.id] = str(record.seq)
            logging.info("Total {} mdrs observed".format(len(self.mdr)))
        else:
            raise IOError('{} Files {} and {} does not exist:'.format(self.ko, self.mdrFile, self.msaFile))

    def msaLoc(self):
        contigIDs =[]
        start = []
        end = []
        for record in SeqIO.parse(self.msaFile, "fasta"):
            sequence = str(record.seq)
            #print(record.id)
            headtails   = re.compile("^(-*)\S+[A-Za-z](-*)$")
            regex       = headtails.search(sequence)
            contigIDs.append(record.id)
            start.append(len(regex.group(1)))
            end.append(len(sequence) - len(regex.group(2)))
        return pd.DataFrame({
        "contigID": contigIDs,
        "start": start,
        "end":end,
        "msaS" : self.msaStart,
        "msaE" : self.msaEnd
        })

    def spanDF(self):
        '''
        returns if the contig spans at least 80% (default) of the MDR
        '''
        df = pd.concat([self._isSpanning(r.id, r.seq) for r in SeqIO.parse(self.msaFile, "fasta") if r.id in self.mdr], ignore_index=True).reset_index()
        logging.info("MDR: Total: {} Spanning: {} Not spanning: {}".format(str(len(df)),str(sum(df['spanning'])), str(sum(np.invert(df['spanning'])))))
        df['ko'] = self.ko
        df['contigID']  = df.apply(lambda row : "{}:{}".format(row['ko'], row['id']), axis = 1)
        return df

    def _isSpanning(self, id, seq):
        '''
        returns dataFrame single row id of
        '''
        seq = str(seq)
        #Tracer()()
        mdrRegion = seq[self.msaStart:self.msaEnd]
        b4 = seq[:self.msaStart-1].upper()
        after = seq[self.msaEnd:].upper() #omg major bug....
        frontEmpty = re.search("[ATCG]", b4) is None and mdrRegion[0] == '-'
        backEmpty = re.search("[ATCG]", after) is None and mdrRegion[-1] == '-'
        notSpanning =  frontEmpty or backEmpty

        # if it occupies (threshold eg. 80%) of the MDR
        perc = len(re.sub("-", "", seq[self.msaStart:self.msaEnd])) / (self.msaEnd - self.msaStart)

        def cover80(mdr, frontEmpty, backEmpty, contigID):
            try:
                tailer = frontEmpty and not backEmpty
                header = backEmpty and not frontEmpty
                matches = re.search("^(-*)[ATCGatcg]\S+(-*)$", mdr)
                front = matches.group(1)
                back = matches.group(2)

                if tailer:
                    spanningPerc = len(front) / len(mdr)
                elif header:
                    spanningPerc = len(back) / len(mdr)
                else:
                    spanningPerc = 0
                logging.info("Non-spanning contig: {} spans: {}%".format(contigID, round(spanningPerc*100, 2)))
                return True if spanningPerc > 0.8 else False
            except Exception as err:
                logging.info("This contig: {} error".format(contigID))

        #contig00278
        #problematic = []
        #with open(sa.msaFile, "r") as f:
            #for record in SeqIO.parse(f, "fasta"):
                #if record.id == 'contig00278':
                    #logging.info(record.id)
                    #logging.info(str(record.seq))
                    #problematic.append(str(record.seq))

        return pd.DataFrame({
            "id" : id,
            "spanning": not notSpanning,
            "basePerc": perc,
            "cover80": cover80(mdrRegion, frontEmpty, backEmpty, id) if notSpanning else True,
            "flankingB4": len(re.sub("-", "", b4)),
            "flankingAfter": len(re.sub("-", "", after))
        }, index=[1])

    def shannon(self):
        '''
        Shannon Entropy represents a lower limit of lossless data compression
        H=-\sum_{i=1}^{M} P_i\,log_2\,P_i
        '''
        alignItem = AlignIO.read(self.msaFile, "fasta")
        mdrSeqs = []
        for seq in alignItem:
            if seq.id in self.mdr:
                mdrSeqs.append(seq)
        mdrAlignment = Align.MultipleSeqAlignment(records=mdrSeqs)
        #Tracer()()
        return self._shannon_entropy_list_msa(mdrAlignment)


    def _shannon_entropy_list_msa(self,alignment_file):
        def isGAP(list_input):
            unique_base = set(list_input)
            if '-' in unique_base and len(unique_base) == 1:
                return 1
            else:
                return 0

        def gapPerc(list_input):
            freqDictionary = {i:list_input.count(i) for i in set(list_input)}
            if '-' in freqDictionary:
                return freqDictionary['-']/len(list_input)
            else:
                return 0

        gap_list =[]
        gap_perc =[]
        shannon_entropy_list = []
        for col_no in range(len(list(alignment_file[0]))):
            list_input = list(alignment_file[:, col_no])
            shannon_entropy_list.append(self._shannon_entropy(list_input))
            gap_list.append(isGAP(list_input))
            gap_perc.append(gapPerc(list_input))
        xvalues = [i for i in range(len(shannon_entropy_list))]
        entropyDF = pd.DataFrame({
            "loc": xvalues,
            "shannon Score": shannon_entropy_list,
            "gapPerc": gap_perc,
            "ko": self.ko,
            "allGaps" : gap_list,
        }, index=[xvalues])
        entropyDF['start'] = self.msaStart
        entropyDF['end'] = self.msaEnd
        return entropyDF


    def _shannon_entropy(self,list_input):
        #self.msaFile
        unique_base = set(list_input)                           # Get only the unique bases in a column
        #unique_base = unique_base.discard("-")
        M   =  len(list_input)
        entropy_list = []
        # Number of residues in column
        for base in unique_base:
            n_i = list_input.count(base)                        # Number of residues of type i
            P_i = n_i/float(M)                                  # n_i(Number of residues of type i) / M(Number of residues in column)
            entropy_i = P_i*(math.log(P_i,2))
            entropy_list.append(entropy_i)
        sh_entropy = -(sum(entropy_list))
        #print sh_entropy
        return sh_entropy

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="""
            Class carrying out various MDR sequence analyses
    """)

    parser.add_argument('root', metavar = 'root_dir', help="path to the root directory", default="/w/simulation_fr_the_beginning/reAssemble/everybodyelse/out")
    parser.add_argument('spanfile', metavar = 'span_output', help="Filename of MDR contig span info")
    parser.add_argument('entropyfile', metavar = 'entropy_output', help="Filename of MDR contig shannon info")
    parser.add_argument('--scg', action="store_true", help="To run this on singlecopy genes")
    parser.add_argument('--ko', default="K00927", help="To run this on singlecopy genes")

    args = parser.parse_args()

    spanning = []
    shannon = []
    scg =["K00927","K02316","K02520","K02838","K02867","K02874","K02884","K02887","K02906","K02931","K02935","K02948","K02965","K02982","K02996","K03664","K01937","K02357","K02600","K02863","K02871","K02878","K02886","K02899","K02926","K02933","K02946","K02952","K02967","K02988","K03043"]
    if (args.scg):
        for ko in scg:
            seqAnalysis = mdrSeqAnalysis(ko, args['root'])
            spanningDF = seqAnalysis.spanDF()
            values = seqAnalysis.shannon()
            spanning.append(spanningDF)
            shannon.append(values)
    else:
        seqAnalysis = mdrSeqAnalysis(args.ko, args.root)
        spanningDF = seqAnalysis.spanDF()
        values = seqAnalysis.shannon()
        spanning.append(spanningDF)
        shannon.append(values)

    spanDF = pd.concat(spanning, ignore_index=True)


    feather.write_dataframe(spanDF, 'flanking.raw.feather')
#   library(tidyverse)
#   library(feather)
#   df = read_feather("flanking.raw.feather")
#   summaryDF = df %>% group_by(ko) %>% summarise(mean_after=mean(flankingAfter), sd_after = sd(flankingAfter), mean_b4=mean(flankingB4), sd_b4 = sd(flankingB4)) %>% mutate(afterRibbonUpper = mean_after + sd_after, afterRibbonlower = mean_after - sd_after, b4RibbonUpper = mean_b4 + sd_b4, b4RibbonLower = mean_b4 - sd_b4) %>% select(-sd_after, -sd_b4)
#   plotDF = rbind(summaryDF %>% select(ko, contains('after')) %>% setNames(c("ko", "mean", "upper", "lower")) %>% mutate(loc="after"), summaryDF %>% select(ko, contains('b4')) %>% setNames(c("ko", "mean", "upper", "lower")) %>% mutate(loc = "b4"))
#   { ggplot(plotDF, aes(x=ko, y=mean, ymin=lower, ymax=upper, group=loc, fill = as.factor(loc), color=as.factor(loc))) + geom_line() + geom_ribbon(alpha=0.3) +  ylab("Length") + facet_wrap(~loc, ncol=1) + scale_fill_discrete("Location") + scale_color_discrete("Location") + guides(alpha=FALSE) + theme(axis.text.x=element_text(angle=90)) } %>% ggsave(., file="ribbon.pdf", w=10, h=15)

    spanDF[['id', 'spanning', 'ko', 'basePerc', 'contigID']].to_csv(args.spanfile, index=False, quoting=1)

    shannonDF = pd.concat(shannon, ignore_index=True).reset_index()
    shannonDF.drop('index', axis=1, inplace=True)
    shannonDF.to_csv(args.entropyfile, index=False)

