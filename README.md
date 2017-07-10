# Introduction

Continuation analysis tools from (pAss)[https://github.com/etheleon/pAss], a gene centric assembly tool but written in python

## Shannon Entropy

We calculate the shannon entropy for each base position

```python
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
```

Shannon Entropy applied on the max diversity of a given KO
![shannonplot](https://github.com/etheleon/pAss2/blob/master/shannonPlot.png)

Shannon Entropy applied on the max diversity of 31 Single Copy Genes
![shannonPlotSCG](https://github.com/etheleon/pAss2/blob/master/shannonSCG.png)

## Spanning Analysis

We identify contigs which span the whole length of the max diversity region; spanning is one fo the criteria.

Following plot shows the length of the contigs’ tailing ends before and after the MDR for 31 Single Copy Genes
![ribbon](https://github.com/etheleon/pAss2/blob/master/ribbon.png)
